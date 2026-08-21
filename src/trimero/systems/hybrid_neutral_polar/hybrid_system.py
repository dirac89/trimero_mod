"""
Sistema híbrido tetra-atómico: Rb*(n, l≥3) + Rb(5s) NEUTRO en θ=π + RbCs
POLAR en θ=0. Primera molécula híbrida del proyecto (ronda 2026-08-21).

GEOMETRÍA Y CONVENIOS
---------------------
* Eje de cuantización Z a lo largo de R⃗₂: el RbCs está en +z, a distancia
  R2 del núcleo Rydberg (carga +e en el origen), con su dipolo permanente
  d⃗ apuntando según Z (θ_d = 0). Es EXACTAMENTE el convenio de
  `ChargeDipoleHamiltonian`, así que H_mol se reutiliza sin tocar nada.
* El Rb neutro (perturbador de Fermi) está en R⃗₁ = −R1·ẑ, es decir θ=π,
  SIN compañero en θ=0: la configuración lineal NO es simétrica y el
  `parity_factor` del trímero NO aplica. La transformación del elemento
  individual de θ=0 a θ=π es

      ⟨l₁m|V(θ=π)|l₂m⟩ = (−1)^{l₁+l₂} · ⟨l₁m|V(θ=0)|l₂m⟩,

  derivada con el operador paridad (𝒫∇𝒫†=−∇: los dos signos del término p
  se cancelan) y verificada POR FUERZA BRUTA contra ψ y ∇ψ numéricos en
  (0,0,−R) en tests/systems/hybrid_neutral_polar/test_transformacion_theta_pi.py.
* V_Fermi actúa sólo sobre el electrón: es identidad sobre el rotor
  (δ_{N N'} δ_{M_N M_N'}). H_A es diagonal; H_mol conserva M_J = m_l+M_N por
  simetría axial alrededor de Z; V_Fermi(θ=π) también (el perturbador está
  SOBRE el eje). La conservación de M_J no se asume: se verifica en
  test_conservacion_mj.py con elementos de matriz entre bloques distintos.

PARÁMETROS DEL RbCs (confirmados a mano esta ronda)
---------------------------------------------------
    d = 1.225 D   → d_au = 1.225 × DEBYE_TO_EA0 = 0.481952126075 e·a₀
                    (constante CODATA del módulo polar; con el Debye
                    redondeado 3.33564e-30 C·m sale 0.481951942617:
                    diferencia 4·10⁻⁷ relativa)
    B = 490.17 MHz → B_au = 490.17e6 / HZ_PER_HARTREE = 7.449750e−08 E_h
                    (verificado por dos vías independientes: idéntico)

DOMINIO DE VALIDEZ
------------------
Con n*=35 la tabla A_s/A_p se lee tal cual (remapeo identidad) para
111 ≤ R1 ≤ 2448 a₀, y k² > 0 exige R1 < 2450 a₀. Los R1 representativos de
esta ronda (600, 900, 1100 a₀) están dentro. R2 no tiene restricción de
tabla (H_mol no usa tablas de dispersión), sólo sentido físico R2 ≫ n².

ESTADO DE LA RONDA
------------------
Sólo construcción del módulo y tests de límites (puntos 2-4). Las curvas
BOP del híbrido y ⟨cos θ_d⟩ son Fase 3, ronda siguiente.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Optional

import numpy as np

from trimero.basis.quantum import CoupledBasis, QuantumBasisBlock
from trimero.basis.radial import RadialBasis
from trimero.systems.rb_atom import Atom
from trimero.systems.hybrid_neutral_polar.parity import phase_pi
from trimero.systems.rb_krb_polar.charge_dipole import (
    DEBYE_TO_EA0,
    HZ_PER_HARTREE,
    ChargeDipoleHamiltonian,
    RydbergElectronField,
    State,
    rydberg_diagonal,
)
from trimero.systems.rb_krb_polar.rb_defects import neighbor_levels
from trimero.systems.rb_neutral_perturber.fermi_krb import (
    FermiPseudopotential,
    ScatteringLengths,
    hydrogenic_R,
    hydrogenic_dR,
)
from trimero.systems.rb_neutral_perturber.linear_trimer import NeighborRadials

__all__ = ["GHZ_PER_HARTREE", "D_RBCS_DEBYE", "B_RBCS_MHZ",
           "D_RBCS_AU", "B_RBCS_AU", "HybridNeutralPolar"]

GHZ_PER_HARTREE = HZ_PER_HARTREE / 1.0e9

# Parámetros moleculares del RbCs (estado vibracional de partida).
D_RBCS_DEBYE = 1.225          # momento dipolar permanente, Debye
B_RBCS_MHZ = 490.17           # constante rotacional, MHz

# Conversiones verificadas esta ronda (ver cabecera):
D_RBCS_AU = D_RBCS_DEBYE * DEBYE_TO_EA0          # 0.481952126075 e·a₀
B_RBCS_AU = B_RBCS_MHZ * 1.0e6 / HZ_PER_HARTREE  # 7.449749956417e-08 E_h


@dataclass
class HybridNeutralPolar:
    """
    H(R1, R2) del sistema híbrido Rb*–(Rb neutro en θ=π)–(RbCs en θ=0).

    Construir cuesta unos segundos (base radial); reutilizar la instancia
    para barrer. `hamiltonian` devuelve la matriz del bloque M_J en hartree;
    `eigvals` sus autovalores.

    Args:
        n_manifold: n del manifold cuasi-degenerado (l ≥ l_min). Los vecinos
            individuales son (n+3)s, (n+2)p, (n+1)d — l = 0, 1, 2.
        N_max: rotor rígido del RbCs incluido hasta N = N_max.
        d_debye / B_mhz: parámetros del RbCs; las conversiones a u.a. están
            precalculadas arriba y verificadas.
        include_electron_field: incluye −d⃗·F_elec (campo del electrón
            Rydberg sobre el dipolo). False reproduce H_mol de la ronda del
            ion a secas.
        s_wave / p_wave: canales del pseudopotencial de Fermi.
        radial_source: radiales inyectadas a V_Fermi para l ≤ 2:
            "tabulated" (Coulomb exacta tabulada, como el trímero validado)
            o "hydrogenic" (hidrogenoide entera, camino del sistema polar).
        p_interpolation: interpolación de 1/A_p, como en los otros sistemas.
    """

    n_manifold: int = 35
    N_max: int = 6
    l_min: int = 3
    d_debye: float = D_RBCS_DEBYE
    B_mhz: float = B_RBCS_MHZ
    include_electron_field: bool = True
    s_wave: bool = True
    p_wave: bool = True
    radial_source: str = "tabulated"
    p_interpolation: str = "inverse"

    basis: CoupledBasis = field(init=False)
    radial_basis: RadialBasis = field(init=False)
    neighbors: Optional[NeighborRadials] = field(init=False)
    scattering: ScatteringLengths = field(init=False)
    pseudo: FermiPseudopotential = field(init=False)
    hmol: ChargeDipoleHamiltonian = field(init=False)

    def __post_init__(self):
        if self.radial_source not in ("tabulated", "hydrogenic"):
            raise ValueError(
                f"radial_source debe ser 'tabulated' o 'hydrogenic', "
                f"no {self.radial_source!r}")
        self.l_max = self.n_manifold - 1
        self.neighbor_levels = neighbor_levels(self.n_manifold)
        self.basis = CoupledBasis(
            N_max=self.N_max, manifold_l_min=self.l_min,
            manifold_l_max=self.l_max, neighbor_l=(0, 1, 2))
        self.radial_basis = RadialBasis(
            n_manifold=self.n_manifold, l_min=self.l_min,
            l_max=self.l_max, neighbor_l=(0, 1, 2))
        data_dir = (Path(__file__).resolve().parents[4]
                    / "data" / "Wavefunction")
        self.neighbors = (NeighborRadials(data_dir)
                          if self.radial_source == "tabulated" else None)
        # Convenio del paper: un único k(R), el del manifold; con
        # uniform_n_star == n_star_table el remapeo es la identidad.
        self.scattering = ScatteringLengths(
            data_dir=data_dir, n_star_table=float(self.n_manifold),
            p_interpolation=self.p_interpolation)
        self.pseudo = FermiPseudopotential(
            self.radial_basis, self.scattering,
            s_wave=self.s_wave, p_wave=self.p_wave,
            n_manifold=self.n_manifold, uniform_n_star=float(self.n_manifold),
            radial_fn=self._radial_fn, dradial_fn=self._dradial_fn)
        self.hmol = ChargeDipoleHamiltonian(
            B=B_RBCS_AU, d=D_RBCS_AU,
            electron_field=(RydbergElectronField(self.radial_basis)
                            if self.include_electron_field else None))
        self._blocks: Dict[tuple, np.ndarray] = {}

    @property
    def d_au(self) -> float:
        """Momento dipolar del RbCs en e·a₀."""
        return self.d_debye * DEBYE_TO_EA0

    @property
    def B_au(self) -> float:
        """Constante rotacional del RbCs en E_h."""
        return self.B_mhz * 1.0e6 / HZ_PER_HARTREE

    @property
    def E_manifold(self) -> float:
        """
        Cero de energía: manifold degenerado (los defectos δ_l = 0 para
        l ≥ 3 hacen E_Rb(n,l) = −0.5/n² exacto; verificado esta ronda).
        """
        return Atom(self.n_manifold, 3).E_Rb()

    # ------------------------------------------------------------ radiales
    def _radial_fn(self, l: int, R: float) -> float:
        if self.neighbors is not None and self.neighbors.has(l):
            return self.neighbors.R(l, R)
        return float(hydrogenic_R(self.radial_basis.n_of_l(l), l, R))

    def _dradial_fn(self, l: int, R: float) -> float:
        if self.neighbors is not None and self.neighbors.has(l):
            return self.neighbors.dR(l, R)
        return float(hydrogenic_dR(self.radial_basis.n_of_l(l), l, R))

    # -------------------------------------------------------------- bloque
    def block(self, M_J: int = 0) -> QuantumBasisBlock:
        return self.basis.get_block(M_J)

    def rydberg_diagonal(self, M_J: int = 0) -> np.ndarray:
        key = ("E", M_J)
        if key not in self._blocks:
            self._blocks[key] = rydberg_diagonal(
                self.block(M_J), n_manifold=self.n_manifold)
        return self._blocks[key]

    # --------------------------------------------------- V_Fermi en θ=π
    def fermi_pi_element_unchecked(self, si: State, sj: State,
                                   R1: float) -> float:
        """
        ⟨si|V_Fermi(θ=π)|sj⟩ SIN comprobar M_J ni ortogonalidad del rotor:
        devuelve sólo el factor electrónico fase × elemento(θ=0).

        Existe para VERIFICAR numéricamente que el elemento completo se anula
        entre M_J distintos (test_conservacion_mj.py), igual que hace
        `_matrix_element_unchecked` en el lado polar. No usar en producción.
        """
        l_i, m_i, _, _ = si
        l_j, m_j, _, _ = sj
        if m_i != m_j:
            return 0.0
        return phase_pi(l_i, l_j) * self.pseudo.electron_element(
            l_i, m_i, l_j, m_j, R1)

    def fermi_pi_element(self, si: State, sj: State, R1: float) -> float:
        """
        ⟨si|V_Fermi(θ=π)|sj⟩ completo: factor electrónico (con la fase de
        paridad) × ortogonalidad del rotor δ_{N N'} δ_{M_N M_N'}.

        Raises:
            ValueError: si los dos estados tienen M_J distinto. El perturbador
                está sobre el eje Z, así que V_Fermi conserva m_l y M_J; un
                elemento entre bloques distintos indica error de indexación.
        """
        if si[1] + si[3] != sj[1] + sj[3]:
            raise ValueError(
                f"V_Fermi conserva M_J: elemento pedido entre M_J="
                f"{si[1] + si[3]} (estado {si}) y M_J={sj[1] + sj[3]} "
                f"(estado {sj})")
        if si[2] != sj[2] or si[3] != sj[3]:
            return 0.0
        return self.fermi_pi_element_unchecked(si, sj, R1)

    def fermi_pi_matrix(self, blk: QuantumBasisBlock, R1: float) -> np.ndarray:
        """
        Matriz densa de V_Fermi(θ=π) para un bloque M_J.

        Sólo visita pares con el mismo (m_l, N, M_N): el resto es cero
        exacto (reglas de selección Δm=0, |m|≤1 más identidad rotacional).
        Cada elemento lleva la fase (−1)^{l+l'} verificada por fuerza bruta
        en test_transformacion_theta_pi.py.
        """
        states = blk.states
        dim = len(states)
        groups: Dict[tuple, list] = {}
        for i, st in enumerate(states):
            groups.setdefault((st[1], st[2], st[3]), []).append(i)
        V = np.zeros((dim, dim), dtype=np.float64)
        for idx in groups.values():
            for i in idx:
                li, mi, _, _ = states[i]
                for j in idx:
                    lj, mj, _, _ = states[j]
                    V[i, j] = phase_pi(li, lj) * self.pseudo.electron_element(
                        li, mi, lj, mj, R1)
        return V

    # ---------------------------------------------------------- Hamiltoniano
    def hamiltonian(self, R1: float, R2: float, M_J: int = 0,
                    fermi: bool = True, mol: bool = True) -> np.ndarray:
        """
        H(R1, R2) del bloque M_J en hartree.

        `fermi=False` apaga V_Fermi EXACTAMENTE (matriz nula, límite polar
        puro); `mol=False` apaga H_mol (límite neutro puro). Como en
        BOPSystem, los interruptores existen para MEDIR límites, nunca para
        disimular dominios inválidos.
        """
        blk = self.block(M_J)
        H = np.diag(self.rydberg_diagonal(M_J))
        if mol:
            H = H + self.hmol.build(blk, R2)
        if fermi:
            H = H + self.fermi_pi_matrix(blk, R1)
        return H

    def eigvals(self, R1: float, R2: float, M_J: int = 0,
                fermi: bool = True, mol: bool = True) -> np.ndarray:
        """Autovalores en hartree (sin desplazar)."""
        return np.linalg.eigvalsh(self.hamiltonian(R1, R2, M_J, fermi, mol))

    def spectrum_ghz(self, R1: float, R2: float, M_J: int = 0,
                     fermi: bool = True, mol: bool = True) -> np.ndarray:
        """Autovalores en GHz respecto al manifold libre."""
        return ((self.eigvals(R1, R2, M_J, fermi, mol) - self.E_manifold)
                * GHZ_PER_HARTREE)
