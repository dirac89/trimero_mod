"""
Montaje del sistema Rb*-KRb para curvas BOP, parametrizado por el manifold.

Hasta ahora el manifold n=24 y su vecino 27s estaban repetidos a mano en
`run_bop_curve.py`, en `CoupledBasis` (l ≤ 23), en `RadialBasis` y en
`n_star_of_l`. Este módulo reúne ese montaje en un sitio: se pide un n y sale
todo lo demás.

COMPOSICIÓN DE LA BASE ELECTRÓNICA
----------------------------------
Aguilera-Fernández, Sadeghpour, Schmelcher & González-Férez, J. Phys.: Conf.
Ser. 635, 012023 (2015), arXiv:1507.07972, la define textualmente como

    «the (n, l≥3) degenerate manifold, and the energetically neighboring
     levels (n+1)d, (n+2)p, and (n+3)s»

es decir el manifold MÁS TRES niveles individuales, no uno. Para n=24: 25d, 26p
y 27s. Las rondas anteriores de esta sesión usaron sólo manifold + (n+3)s;
`neighbors=(0,)` reproduce esa base incompleta para poder medir la diferencia.

    H(R) = H_a + H_mol + V_Fermi
         = diag(E_ryd) + [B·N² - d·F_ryd] + V_Fermi

Lo que depende del manifold, y de dónde sale:

  * l máximo del manifold ................ n_manifold - 1   (CoupledBasis)
  * qué niveles son los vecinos .......... rb_defects.neighbor_levels(n)
  * n de las funciones radiales .......... n_manifold, y round(n*) para cada
                                           vecino (RadialBasis; aproximación
                                           documentada en basis/radial.py)
  * energías diagonales de H_a ........... rydberg_diagonal(n_manifold)
  * n* del remapeo k(R) del Fermi ........ FermiPseudopotential(n_manifold)

El cero de energía es siempre E(manifold) + KRb(N=0) = -1/(2 n_manifold²), que
para l ≥ 3 es exacto: `Atom.E_Rb()` no aplica defecto cuántico a l ≥ 3.

⚠️ El límite superior del dominio en R NO es «donde se acaba la física»: es
donde deja de estar definido el REMAPEO SEMICLÁSICO k(R) con el que leemos las
tablas de longitudes de dispersión (k² < 0 pasado el punto de retorno clásico
del par más ligado). H_mol no tiene esa restricción — usa la función de onda
hidrogenoide, que decae pero no se anula ahí. Es una limitación de cómo leemos
las tablas, no del problema.
"""

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from trimero.basis.quantum import CoupledBasis, QuantumBasisBlock
from trimero.basis.radial import RadialBasis
from trimero.systems.rb_krb_polar.charge_dipole import (
    B_KRB_GHZ,
    HZ_PER_HARTREE,
    ChargeDipoleHamiltonian,
    RydbergElectronField,
    rydberg_diagonal,
)
from trimero.systems.rb_neutral_perturber.fermi_krb import (
    FermiPseudopotential,
    ScatteringLengths,
)
from trimero.systems.rb_krb_polar.rb_defects import n_star_nl, neighbor_levels

__all__ = ["GHZ_PER_HARTREE", "BOPSystem"]

GHZ_PER_HARTREE = HZ_PER_HARTREE / 1.0e9


@dataclass
class BOPSystem:
    """
    Sistema completo para un manifold dado. Construir cuesta unos segundos
    (tablas radiales); reutilizar la instancia para todo un barrido en R.

    Args:
        n_manifold: n del manifold cuasi-degenerado (l ≥ l_min).
        n_s: n del estado ns vecino (l=0). Por defecto n_manifold + 3, que es
            la relación del paper (n=24 ↔ 27s, n=25 ↔ 28s).
        delta0_ns: δ₀ de la serie s. `DELTA0_NS_PAPER` para comparar con
            González-Férez 2015; None delega en `Atom.E_Rb()`.
    """

    n_manifold: int = 24
    n_s: Optional[int] = None
    delta0_ns: Optional[float] = None
    N_max: int = 6
    l_min: int = 3
    p_interpolation: str = "inverse"
    # Vecinos individuales incluidos: (0,1,2) = (n+3)s, (n+2)p, (n+1)d, que es
    # la base del paper. (0,) reproduce la base INCOMPLETA de las rondas
    # anteriores, sólo para medir cuánto cambian las cosas al completarla.
    neighbors: tuple = (0, 1, 2)

    basis: CoupledBasis = field(init=False)
    radial: RadialBasis = field(init=False)
    scattering: ScatteringLengths = field(init=False)
    hmol: ChargeDipoleHamiltonian = field(init=False)
    fermi: FermiPseudopotential = field(init=False)

    def __post_init__(self):
        self.levels = neighbor_levels(self.n_manifold)   # {0: n+3, 1: n+2, 2: n+1}
        if self.n_s is not None:
            self.levels[0] = self.n_s
        self.n_s = self.levels[0]
        self.neighbor_l = tuple(sorted(l for l in self.levels
                                       if l in self.neighbors))
        self.levels = {l: self.levels[l] for l in self.neighbor_l}
        self.l_max = self.n_manifold - 1
        self.basis = CoupledBasis(N_max=self.N_max, manifold_l_min=self.l_min,
                                  manifold_l_max=self.l_max,
                                  neighbor_l=self.neighbor_l)
        self.radial = RadialBasis(n_manifold=self.n_manifold, l_min=self.l_min,
                                  l_max=self.l_max, neighbor_l=self.neighbor_l)
        self.n_s_eff = self.radial.n_of_l(0) if 0 in self.neighbor_l else None
        self.scattering = ScatteringLengths(p_interpolation=self.p_interpolation)
        self.hmol = ChargeDipoleHamiltonian(
            electron_field=RydbergElectronField(self.radial))
        self.fermi = FermiPseudopotential(
            self.radial, self.scattering, delta0_ns=self.delta0_ns,
            n_manifold=self.n_manifold, n_s=self.n_s)
        self._blocks = {}

    # -- energías de referencia ----------------------------------------
    def n_star_of_l(self, l: int) -> float:
        """n* exacto del nivel asociado a este l."""
        return self.fermi.n_star_of_l(l)

    def n_star_s(self) -> float:
        """n* del estado ns, con el δ₀ configurado."""
        return n_star_nl(self.levels[0], 0, self.delta0_ns)

    @property
    def E_manifold(self) -> float:
        """Cero de energía: manifold + KRb(N=0), en E_h. Exacto para l ≥ 3."""
        return -0.5 / self.n_manifold**2

    @property
    def E_ns(self) -> float:
        """Energía del estado ns vecino, en E_h."""
        return -0.5 / self.n_star_s() ** 2

    def delta_E_ghz(self, l: int) -> float:
        """ΔE(nivel l) = E(nivel) - E(manifold), en GHz. Umbral rotacional N=0."""
        if l not in self.levels:
            return 0.0
        return (-0.5 / self.n_star_of_l(l) ** 2 - self.E_manifold) * GHZ_PER_HARTREE

    def delta_E_ns_ghz(self) -> float:
        """ΔE(ns) = E(ns) - E(manifold), en GHz. Umbral rotacional N=0."""
        return (self.E_ns - self.E_manifold) * GHZ_PER_HARTREE

    def thresholds_level_ghz(self, l: int, N_values, B_ghz: float = B_KRB_GHZ):
        """Umbrales del nivel l + KRb(N): ΔE(l) + B·N(N+1), en GHz."""
        d = self.delta_E_ghz(l)
        return {int(N): d + B_ghz * N * (N + 1) for N in N_values}

    def thresholds_ns_ghz(self, N_values, B_ghz: float = B_KRB_GHZ):
        """Umbrales ns + KRb(N): ΔE(ns) + B·N(N+1), en GHz."""
        d = self.delta_E_ns_ghz()
        return {int(N): d + B_ghz * N * (N + 1) for N in N_values}

    def thresholds_manifold_ghz(self, N_values, B_ghz: float = B_KRB_GHZ):
        """Umbrales manifold + KRb(N): B·N(N+1), en GHz (el cero es N=0)."""
        return {int(N): B_ghz * N * (N + 1) for N in N_values}

    # -- bloque y Hamiltoniano -----------------------------------------
    def block(self, M_J: int = 0) -> QuantumBasisBlock:
        return self.basis.get_block(M_J)

    def rydberg_diagonal(self, M_J: int = 0) -> np.ndarray:
        key = ("E", M_J)
        if key not in self._blocks:
            self._blocks[key] = rydberg_diagonal(
                self.block(M_J), n_manifold=self.n_manifold, n_s=self.n_s,
                delta0_ns=self.delta0_ns)
        return self._blocks[key]

    def is_manifold(self, M_J: int = 0) -> np.ndarray:
        """Máscara de estados con carácter de manifold (l ≥ l_min)."""
        key = ("m", M_J)
        if key not in self._blocks:
            self._blocks[key] = np.array(
                [st[0] >= self.l_min for st in self.block(M_J).states])
        return self._blocks[key]

    def hamiltonian(self, R: float, M_J: int = 0,
                    fermi: bool = True) -> np.ndarray:
        """
        H(R) del bloque M_J.

        `fermi=False` pone V_Fermi ≡ 0 EXACTAMENTE (matriz nula, no un valor
        devuelto por error ni una extrapolación de las tablas). Sirve para
        MEDIR cuánto aporta el pseudopotencial, y para poder trazar H_a+H_mol
        más allá del dominio del remapeo k(R). Ojo: hacerlo no es gratis, hay
        que justificar en cada caso que V_Fermi sea despreciable ahí — ver
        docs/analysis_extension_dominio_fig1.md, donde se comprueba que en el
        borde NO lo es.
        """
        blk = self.block(M_J)
        H = np.diag(self.rydberg_diagonal(M_J)) + self.hmol.build(blk, R)
        if fermi:
            H = H + self.fermi.build(blk, R)
        return H

    def eigvals(self, R: float, M_J: int = 0, fermi: bool = True) -> np.ndarray:
        """Sólo autovalores: ~2x más rápido, suficiente para curvas adiabáticas."""
        return np.linalg.eigvalsh(self.hamiltonian(R, M_J, fermi))

    def solve(self, R: float, M_J: int = 0, fermi: bool = True):
        return np.linalg.eigh(self.hamiltonian(R, M_J, fermi))

    def character_curve(self, R: float, M_J: int = 0, fermi: bool = True,
                        weight: float = 0.5):
        """
        (E − E_manifold [GHz], k, peso) de la curva adiabática más baja con
        CARÁCTER de manifold en este R.

        No se usa el índice fijo identificado en un extremo: con la base
        completa esa curva cambia de carácter al cruzar los estados de
        (n+1)d y (n+2)p. Ver docs/analysis_base_correcta_3_vecinos.md §5.1.
        """
        w, V = self.solve(R, M_J, fermi)
        mask = self.is_manifold(M_J)
        for k in range(len(w)):
            wm = float(np.sum(V[:, k][mask] ** 2))
            if wm > weight:
                return (w[k] - self.E_manifold) * GHZ_PER_HARTREE, k, wm
        return float("nan"), -1, 0.0

    # -- utilidades de dominio y resonancia -----------------------------
    def domain_bounds(self):
        return self.fermi.domain_bounds()

    def p_resonance_window(self, margin_factor: float = 2.0):
        """
        Ventana de exclusión de la resonancia de forma p EN LA COORDENADA R DE
        ESTE MANIFOLD. La resonancia es la misma en energía (ε ≈ 24.7 meV) pero
        cae en un R distinto para cada n, porque ε = E_manifold + 1/R.
        """
        return self.scattering.p_resonance_window(
            self.fermi.n_star_manifold(), margin_factor=margin_factor)

    def lowest_manifold_index(self, R: float, M_J: int = 0,
                              weight: float = 0.5) -> int:
        """
        Índice k de la curva adiabática más baja con carácter de manifold.

        Por la regla de no cruce dentro de un bloque M_J, este k es constante en
        todo el dominio: basta identificarlo una vez (típicamente en R_max, donde
        la asignación a umbrales es inequívoca).
        """
        w, V = self.solve(R, M_J)
        mask = self.is_manifold(M_J)
        for k in range(len(w)):
            if float(np.sum(V[:, k][mask] ** 2)) > weight:
                return k
        raise RuntimeError("ningún autoestado con peso de manifold suficiente")
