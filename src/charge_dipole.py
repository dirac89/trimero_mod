"""
Hamiltoniano de acoplamiento carga-dipolo para el trímero Rb*-KRb.

    H_mol = B·N² - d·F_ryd(R,r)                                  (Ec. 3)
    F_ryd(R,r) = e·R/R³ + e·(r-R)/|r-R|³                         (Ec. 4)

Base científica: González-Férez, Sadeghpour & Schmelcher,
                 New J. Phys. 17, 013021 (2015).

ALCANCE DE ESTE MÓDULO (ronda actual)
-------------------------------------
Implementado:
  * B·N²                       — término rotacional, diagonal.
  * -d·(e·R/R³)                — campo del ion Rb⁺, primer término de Ec. 4.

  * -d·(e·(r-R)/|r-R|³)        — campo del electrón Rydberg, vía la expansión
                                 multipolar completa (Ec. A.6-A.10). Opcional:
                                 sólo se incluye si se pasa `electron_field`.

Geometría y convenio de signos
------------------------------
El ion Rb⁺ está en el origen, la molécula KRb en R⃗, y el eje de cuantización
Z se toma a lo largo de R⃗ (el mismo eje que define M_J = m_l + M_N). Con e=1
en unidades atómicas, el campo del ion en la posición de la molécula es

    F_ion = R⃗/R³ = (1/R²) Ẑ,

de modo que la energía del dipolo permanente d⃗ de KRb en ese campo es

    -d⃗·F_ion = -(d/R²)·cos(θ_d),

con θ_d el ángulo polar del eje molecular respecto de Ẑ. Por tanto este
término es diagonal en (l, m_l) y en M_N, y conecta N con N±1.

Todas las magnitudes están en unidades atómicas (E_h, a₀, e·a₀).
"""

from typing import List, Sequence, Tuple

import numpy as np

from angular_algebra import gaunt
from atom import Atom
from quantum_basis import QuantumBasisBlock
from rydberg_radial import RadialBasis

__all__ = [
    "HZ_PER_HARTREE",
    "DEBYE_TO_EA0",
    "B_KRB_GHZ",
    "D_KRB_DEBYE",
    "B_KRB_AU",
    "D_KRB_AU",
    "ChargeDipoleHamiltonian",
    "RydbergElectronField",
    "rydberg_diagonal",
]

# --- constantes -------------------------------------------------------
# CODATA 2018 (NIST): 1 E_h = 6.579683920502e15 Hz
HZ_PER_HARTREE = 6.579683920502e15
# Factor estándar de conversión de momento dipolar
DEBYE_TO_EA0 = 0.393430307

# B(KRb): Ni et al., Phys. Chem. Chem. Phys. 11, 9626 (2009), citado en el paper
B_KRB_GHZ = 1.114
# d(KRb): Ni et al., Science 322, 231 (2008), citado en el paper
D_KRB_DEBYE = 0.566

B_KRB_AU = B_KRB_GHZ * 1.0e9 / HZ_PER_HARTREE   # ≈ 1.693e-7 E_h
D_KRB_AU = D_KRB_DEBYE * DEBYE_TO_EA0           # ≈ 0.2227 e·a₀

# Estado = (l, m_l, N, M_N)
State = Tuple[int, int, int, int]


class RydbergElectronField:
    """
    Componentes esféricas del campo del ELECTRÓN Rydberg en la posición de la
    molécula,  F⃗_elec = e·(r⃗-R⃗)/|r⃗-R⃗|³, sobre la base electrónica {|l m_l⟩}.

    Derivación (eje Z ∥ R⃗, el mismo convenio de la ronda del ion):

        F⃗_elec/e = ∇_R (1/|r⃗-R⃗|)
                 = ẑ ∂f/∂R  -  ρ̂ (1/R) ∂f/∂γ,     f = Σ_k g_k(r,R) P_k(cos γ)

    con γ el ángulo entre r⃗ y R⃗ (aquí el ángulo polar del electrón) y ρ̂ el
    unitario radial cilíndrico. De ahí, usando dP_k(cosγ)/dγ = -sinγ P_k'(cosγ)
    y P_k^1 = -sqrt(1-x²) P_k' (Condon-Shortley):

        F_z   = Σ_k (∂g_k/∂R) P_k(cos θ)
        F_ρ   = -(1/R) Σ_k g_k P_k^1(cos θ)

    Pasando a componentes esféricas F_{±1} = ∓(F_x ± i F_y)/√2 y usando
    P_k(cosθ) = sqrt(4π/(2k+1)) Y_{k0},
    P_k^1 e^{±iφ} = ± sqrt(4π k(k+1)/(2k+1)) Y_{k,±1}:

        ⟨l1 m1|F_0   |l2 m2⟩ = Σ_k Z^k sqrt(4π/(2k+1))          ·C(l1m1|k, 0|l2m2)
        ⟨l1 m1|F_{±1}|l2 m2⟩ = (1/(√2 R)) Σ_k G^k
                                 sqrt(4π k(k+1)/(2k+1))          ·C(l1m1|k,±1|l2m2)

    con C el coeficiente de Gaunt ∫Y*_{l1m1} Y_{kq} Y_{l2m2} dΩ.

    La suma en k NO es una truncación: el Gaunt anula todo k fuera de
    |l1-l2| ≤ k ≤ l1+l2 con l1+l2+k par, así que cada elemento de matriz es
    una suma FINITA y EXACTA.
    """

    def __init__(self, radial: "RadialBasis" = None):
        self.radial = radial if radial is not None else RadialBasis()
        self._cache = {}

    def k_values(self, l1: int, l2: int, q: int, k_max: int = None):
        """Multipolos que contribuyen a ⟨l1 m1|F_q|l2 m2⟩ (Gaunt no nulo)."""
        k_hi = l1 + l2 if k_max is None else min(l1 + l2, k_max)
        k_lo = max(abs(l1 - l2), abs(q))
        return tuple(k for k in range(k_lo, k_hi + 1) if (l1 + l2 + k) % 2 == 0)

    def field_element(
        self, l1: int, m1: int, l2: int, m2: int, mu: int, R: float, k_max: int = None
    ) -> float:
        """
        ⟨l1 m1| F_mu |l2 m2⟩ del campo del electrón, en unidades atómicas.

        `k_max` trunca la expansión multipolar; sólo para los tests del límite
        dipolar (k_max=1) y del monopolo (k_max=0). En producción se deja None,
        que es exacto.
        """
        if m1 != m2 + mu:
            return 0.0
        key = (l1, m1, l2, m2, mu, R, k_max)
        hit = self._cache.get(key)
        if hit is not None:
            return hit

        ks = self.k_values(l1, l2, mu, k_max)
        if not ks:
            self._cache[key] = 0.0
            return 0.0

        ang = np.array([gaunt(l1, m1, k, mu, l2, m2) for k in ks])
        karr = np.array(ks, dtype=float)
        if mu == 0:
            rad = self.radial.dg_integrals(l1, l2, ks, R)
            val = float(np.sum(rad * np.sqrt(4.0 * np.pi / (2.0 * karr + 1.0)) * ang))
        else:
            rad = self.radial.g_integrals(l1, l2, ks, R)
            coef = np.sqrt(4.0 * np.pi * karr * (karr + 1.0) / (2.0 * karr + 1.0))
            val = float(np.sum(rad * coef * ang) / (np.sqrt(2.0) * R))

        self._cache[key] = val
        return val


class ChargeDipoleHamiltonian:
    """
    H_mol = B·N² - d·F_ion(R)  sobre la base acoplada {(l, m_l, N, M_N)}
    bloqueada por M_J = m_l + M_N.

    El término del campo del electrón Rydberg NO está incluido todavía.
    """

    def __init__(
        self,
        B: float = B_KRB_AU,
        d: float = D_KRB_AU,
        electron_field: "RydbergElectronField" = None,
    ):
        """
        Args:
            B: constante rotacional de KRb en E_h.
            d: momento dipolar permanente de KRb en e·a₀.
            electron_field: si se pasa, se incluye el término del campo del
                electrón Rydberg. Si es None (por defecto) el Hamiltoniano es
                exactamente el verificado en la ronda del ion — la ausencia de
                este argumento no cambia ningún resultado anterior.
        """
        self.B = B
        self.d = d
        self.electron_field = electron_field

    # -- bloques elementales -------------------------------------------
    @staticmethod
    def cos_theta_element(N1: int, M_N1: int, N2: int, M_N2: int) -> float:
        """
        ⟨N1 M_N1| cos(θ_d) |N2 M_N2⟩ para el rotor rígido.

        cos θ = sqrt(4π/3)·Y_10, de donde ΔM_N = 0 y ΔN = ±1:

            ⟨N-1 M|cosθ|N M⟩ = ⟨N M|cosθ|N-1 M⟩
                             = sqrt( (N² - M²) / ((2N-1)(2N+1)) )

        Se escribe con N_max = max(N1,N2) para que las dos ramas sean la
        MISMA expresión: la simetría es exacta bit a bit, no aproximada.
        """
        if M_N1 != M_N2:
            return 0.0
        if abs(N1 - N2) != 1:
            return 0.0
        n = max(N1, N2)
        m = M_N1
        return np.sqrt((n * n - m * m) / ((2.0 * n - 1.0) * (2.0 * n + 1.0)))

    def rotational_element(self, N1: int, M_N1: int, N2: int, M_N2: int) -> float:
        """
        ⟨N1 M_N1| B·N² |N2 M_N2⟩ = δ_{N1 N2} δ_{M_N1 M_N2} · B·N(N+1).

        N² es diagonal en la base |N M_N⟩ con autovalor N(N+1) (en ℏ²=1).
        """
        if N1 != N2 or M_N1 != M_N2:
            return 0.0
        return self.B * N1 * (N1 + 1)

    def ion_field_element(
        self, state_i: State, state_j: State, R: float
    ) -> float:
        """
        ⟨i| -d⃗·F_ion(R) |j⟩ con F_ion = R⃗/R³ = (1/R²)Ẑ.

        Diagonal en (l, m_l) porque el operador no actúa sobre el electrón
        Rydberg; ΔM_N = 0 y ΔN = ±1 por cos(θ_d).
        """
        l_i, m_i, N_i, MN_i = state_i
        l_j, m_j, N_j, MN_j = state_j
        if l_i != l_j or m_i != m_j:
            return 0.0
        ang = self.cos_theta_element(N_i, MN_i, N_j, MN_j)
        if ang == 0.0:
            return 0.0
        return -self.d * ang / (R * R)

    # -- elemento de matriz completo -----------------------------------
    def electron_field_element(
        self, state_i: State, state_j: State, R: float
    ) -> float:
        """
        ⟨i| -d⃗·F_elec(R) |j⟩ con F_elec = e·(r⃗-R⃗)/|r⃗-R⃗|³.

        Producto escalar en componentes esféricas: -d⃗·F⃗ = -Σ_μ (-1)^μ d_μ F_{-μ},
        con d_μ = d·sqrt(4π/3)·Y_{1μ}(Ω_d) actuando sobre el rotor. De ahí:

            ⟨i|-d⃗·F⃗|j⟩ = -d·sqrt(4π/3) Σ_μ (-1)^μ
                            ⟨N' M'|Y_{1μ}|N M⟩ · ⟨l' m'|F_{-μ}|l m⟩

        Reglas de selección resultantes: ΔM_N = +μ y Δm_l = -μ, luego ΔM_J = 0
        pero M_N NO se conserva por separado. ΔN = ±1 porque d⃗ es un operador
        de rango 1 sobre el rotor. Δl es libre (limitado por el Gaunt).

        Devuelve 0.0 si no hay `electron_field` configurado.
        """
        if self.electron_field is None:
            return 0.0
        l_i, m_i, N_i, MN_i = state_i
        l_j, m_j, N_j, MN_j = state_j
        total = 0.0
        for mu in (-1, 0, 1):
            if MN_i != MN_j + mu:
                continue
            rot = gaunt(N_i, MN_i, 1, mu, N_j, MN_j)
            if rot == 0.0:
                continue
            fe = self.electron_field.field_element(l_i, m_i, l_j, m_j, -mu, R)
            if fe == 0.0:
                continue
            total += (-1.0 if mu % 2 else 1.0) * rot * fe
        return -self.d * np.sqrt(4.0 * np.pi / 3.0) * total

    def _matrix_element_unchecked(
        self, state_i: State, state_j: State, R: float
    ) -> float:
        """
        ⟨i| H_mol |j⟩ SIN comprobar la regla de selección en M_J.

        Existe para poder VERIFICAR numéricamente que el acoplamiento se
        anula entre M_J distintos, en lugar de asumirlo por construcción.
        No usar en producción: usar `matrix_element`.
        """
        l_i, m_i, N_i, MN_i = state_i
        l_j, m_j, N_j, MN_j = state_j
        # B·N² actúa sólo sobre el rotor, pero el elemento en la base
        # producto exige además ortogonalidad de la parte electrónica:
        # ⟨l m_l| l' m_l'⟩ = δ_{l l'} δ_{m_l m_l'}.
        if l_i == l_j and m_i == m_j:
            rot = self.rotational_element(N_i, MN_i, N_j, MN_j)
        else:
            rot = 0.0
        ion = self.ion_field_element(state_i, state_j, R)
        ele = self.electron_field_element(state_i, state_j, R)
        return rot + ion + ele

    def matrix_element(self, state_i: State, state_j: State, R: float) -> float:
        """
        ⟨i| H_mol |j⟩ con i, j = (l, m_l, N, M_N).

        Raises:
            ValueError: si los dos estados tienen M_J = m_l + M_N distinto.
                El acoplamiento carga-dipolo conserva M_J (simetría axial en
                torno a R⃗), así que un elemento entre M_J distintos indica un
                error de indexación, no un cero físico legítimo.
        """
        MJ_i = state_i[1] + state_i[3]
        MJ_j = state_j[1] + state_j[3]
        if MJ_i != MJ_j:
            raise ValueError(
                f"H_mol conserva M_J: elemento pedido entre M_J={MJ_i} "
                f"(estado {state_i}) y M_J={MJ_j} (estado {state_j})"
            )
        return self._matrix_element_unchecked(state_i, state_j, R)

    # -- construcción de la matriz del bloque --------------------------
    def rotational_diagonal(self, block: QuantumBasisBlock) -> np.ndarray:
        """Vector de energías rotacionales B·N(N+1) en el orden del bloque."""
        return np.array(
            [self.B * N * (N + 1) for (_, _, N, _) in block.states],
            dtype=np.float64,
        )

    def build_reference(self, block: QuantumBasisBlock, R: float) -> np.ndarray:
        """
        Matriz densa recorriendo los dim² pares. O(dim²) llamadas: sólo para
        bloques pequeños y para comprobar `build` en los tests.

        Llama a `matrix_element` de forma independiente para (i,j) y (j,i): la
        simetría NO se impone a posteriori, se comprueba en los tests.
        """
        states: Sequence[State] = block.states
        dim = len(states)
        H = np.zeros((dim, dim), dtype=np.float64)
        for i, si in enumerate(states):
            for j, sj in enumerate(states):
                H[i, j] = self.matrix_element(si, sj, R)
        return H

    def build(self, block: QuantumBasisBlock, R: float) -> np.ndarray:
        """
        Matriz densa de H_mol para un bloque M_J.

        Sólo visita los pares que las reglas de selección permiten: el término
        rotacional es diagonal, y el acoplamiento dipolar es de rango 1 sobre el
        rotor, luego ΔN = ±1 con ΔM_N = μ ∈ {-1,0,1} y Δm_l = -μ. El resto de
        pares es cero exactamente y no se calcula.

        Cada elemento se calcula de forma independiente en (i,j) y (j,i); la
        equivalencia con `build_reference` y la simetría se comprueban en los
        tests, no se asumen.
        """
        states: Sequence[State] = block.states
        index = {st: i for i, st in enumerate(states)}
        l_values = (
            self.electron_field.radial.l_values
            if self.electron_field is not None
            else None
        )
        H = np.zeros((len(states), len(states)), dtype=np.float64)
        for i, si in enumerate(states):
            H[i, i] = self.matrix_element(si, si, R)
            l_i, m_i, N_i, MN_i = si
            partner_l = l_values if l_values is not None else (l_i,)
            for mu in (-1, 0, 1):
                m_j = m_i + mu
                MN_j = MN_i - mu
                for N_j in (N_i - 1, N_i + 1):
                    if N_j < 0 or abs(MN_j) > N_j:
                        continue
                    for l_j in partner_l:
                        if abs(m_j) > l_j:
                            continue
                        j = index.get((l_j, m_j, N_j, MN_j))
                        if j is None:
                            continue
                        H[i, j] = self.matrix_element(si, states[j], R)
        return H


# --- H_a: energías Rydberg (diagonal) ---------------------------------
def rydberg_diagonal(
    block: QuantumBasisBlock, n_manifold: int = 24, n_s: int = 27
) -> np.ndarray:
    """
    Diagonal de H_a en el orden del bloque: energías del Rb* sin perturbar.

    Manifold cuasi-degenerado n=24, l≥3, más el estado vecino 27s (l=0),
    tal como se usa en González-Férez 2015. Las energías salen de
    `atom.Atom.E_Rb()` (defectos cuánticos), no se duplican aquí.

    NOTA: esto es el H_a *libre*. El pseudopotencial de Fermi Rydberg-KRb
    (fermi_potentials.py) no está incluido en esta ronda.
    """
    cache = {}
    out = np.empty(len(block.states), dtype=np.float64)
    for k, (l, _, _, _) in enumerate(block.states):
        if l not in cache:
            n = n_s if l == 0 else n_manifold
            cache[l] = Atom(n, l).E_Rb()
        out[k] = cache[l]
    return out
