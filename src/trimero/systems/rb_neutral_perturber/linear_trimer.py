"""
Trímero Rydberg ultralargo LINEAL SIMÉTRICO en campo eléctrico DC.

Aguilera-Fernández, Schmelcher & González-Férez,
J. Phys. B 49, 124002 (2016) — arXiv:1601.05049.

    H = H₀ + F·r + V(r⃗, R⃗₁) + V(r⃗, R⃗₂)                        (Ec. 2-4)

    V(r⃗, R⃗ᵢ) = 2π A_s[k(Rᵢ)] δ³(r⃗-R⃗ᵢ)
              + 6π A_p³[k(Rᵢ)] ∇⃖ δ³(r⃗-R⃗ᵢ) ∇⃗                    (Ec. 1)

    k²(R)/2 = 1/R − 1/2n²,   F = F Ẑ

Geometría **simétrica lineal** (Fig. 1(a) del paper): los dos átomos neutros
sobre el eje Z del LFF, a lados opuestos del core y a la misma distancia:

    R₁ = R₂ = R,   θ₁ = 0,   θ₂ = π

Base (§II del paper): manifold degenerado Rb(n=35, l≥3) más los tres niveles
vecinos 38s, 37p y 36d. El defecto cuántico del 35f se desprecia.
El cero de energía es E(n=35, l≥3) = −1/(2·35²).

SIMETRÍA DEL SISTEMA — m_l, NO M_J
----------------------------------
Aquí no hay KRb ni rotor: no existen (N, M_N) y por tanto M_J = m_l + M_N no
es el buen número cuántico de las rondas de Rb*-KRb. El eje Z es a la vez el
eje de los dos perturbadores y la dirección del campo, así que **los tres
términos del Hamiltoniano conservan m_l**:

  * V_Fermi con el perturbador en θ=0 o π sólo conecta m₁=m₂ (los armónicos
    esféricos con |m|≥2 y sus gradientes se anulan sobre el eje),
  * F·r = F·z tiene Δm=0, Δl=±1.

Luego H es bloque-diagonal en m_l **incluso con campo**, y la nomenclatura
molecular del paper es literalmente m_l: Σ ≡ m_l=0, Π ≡ |m_l|=1. A partir de
|m_l| ≥ 2 el pseudopotencial es idénticamente cero y sólo queda el manifold
desplazado por el campo.

Se reutiliza `CoupledBasis` con `N_max=0`, que degenera a estados (l, m_l, 0, 0)
y a bloques etiquetados por M_J = m_l. No es un apaño: es la afirmación de que
este sistema es el límite «sin rotor» de la base compartida.

LOS DOS PERTURBADORES SE REDUCEN A UN FACTOR DE PARIDAD
------------------------------------------------------
Como R⃗₂ = −R⃗₁ y ψ_{lm}(−r⃗) = (−1)^l ψ_{lm}(r⃗), también
(∇ψ)(−r⃗) = (−1)^{l+1} (∇ψ)(r⃗). Los dos términos del pseudopotencial son
bilineales en ψ o en ∇ψ, luego

    ⟨l₁m|V(R⃗₂)|l₂m⟩ = (−1)^{l₁+l₂} ⟨l₁m|V(R⃗₁)|l₂m⟩

para AMBAS ondas (los dos signos extra del gradiente se cancelan). Y como
R₁ = R₂, las longitudes de dispersión son las mismas. Por tanto

    V_total = [1 + (−1)^{l₁+l₂}] · V(θ=0)

que vale 2·V(θ=0) si l₁+l₂ es par y CERO si es impar. Ésa es exactamente la
separación gerade/ungerade que el paper cita de su Ref. [15]: en ausencia de
campo el bloque m_l se parte en l pares y l impares, sin acoplamiento entre
ellos. El término F·r, con Δl=±1, es el único que los conecta: por eso el
campo «acopla los estados gerade y ungerade» (§III.A).

APROXIMACIONES EXPLÍCITAS
-------------------------
1. **Funciones radiales de los vecinos.** 38s, 37p y 36d tienen n* no entero.
   Se usan las funciones de Coulomb TABULADAS (`rvsR38s/37p/36d.dat`),
   interpoladas con spline cúbico; la derivada sale del mismo spline. La
   alternativa hidrogenoide de n entero (`radial_source="hydrogenic"`) acierta
   la amplitud al ~1 % pero se desfasa: ver `docs/analysis_trimero_lineal_campo_dc.md`.

2. **Convenio de signo de los vecinos.** El signo relativo entre las radiales
   tabuladas y los elementos dipolares ⟨38s|r|37p⟩, ⟨37p|r|36d⟩, ⟨36d|r|35f⟩
   de `exp_val_r.txt` se ASUME consistente (mismo generador original). Para
   l ≥ 3 sí está verificado: `exp_val_r.txt` reproduce exactamente, signo
   incluido, ∫u_{35,l} r u_{35,l+1} dr con el convenio de `hydrogenic_R`.

3. **Rango en R.** Las tablas A_s/A_p llegan a 2448 a₀ ≈ 2n² = 2450 a₀, el
   punto de retorno clásico de n=35. Más allá k² < 0 y no hay lectura posible
   sin extrapolar al límite k→0. Las figuras del paper llegan a 2600 a₀; ese
   tramo queda fuera.

BUGS ENCONTRADOS EN LOS DATOS DEL LEGADO (no se propagan aquí)
-------------------------------------------------------------
* `rvsDR38s.dat` vale 10⁻⁶ veces la derivada real de `rvsR38s.dat` (las de
  37p y 36d sí son correctas al 0.1 %). Aquí no se usa: la derivada del 38s
  sale del spline.
* `rvsR38s.dat`/`rvsDR38s.dat` tienen 780 filas (malla uniforme 111..2448 con
  paso 3) mientras `rvsAS/rvsAP/rvsR37p/rvsR36d` tienen 776 (les faltan 4
  filas sobre la resonancia p, hueco 750→765 a₀). `trimer.py` indexa las tres
  por número de fila, así que a partir de la fila 214 lee el 38s con un
  desfase de 12 a₀. Aquí se interpola por valor de R, no por índice.
"""

from pathlib import Path
from typing import Dict, Sequence

import numpy as np
from scipy.interpolate import CubicSpline

from trimero.basis.quantum import CoupledBasis
from trimero.basis.radial import RadialBasis
from trimero.systems.rb_atom import Atom
from trimero.systems.rb_krb_polar.rb_defects import neighbor_levels
from trimero.systems.rb_neutral_perturber.fermi_krb import (
    FermiPseudopotential,
    ScatteringLengths,
    hydrogenic_R,
    hydrogenic_dR,
)

__all__ = [
    "HARTREE_TO_GHZ",
    "AU_FIELD_IN_V_PER_M",
    "field_au",
    "NeighborRadials",
    "SymmetricLinearTrimer",
]

# CODATA 2018: 1 E_h = 6.579 683 920 502e15 Hz.
# OJO: `trimer.py` usa 6.579683920729e9 y lo llama «EhtoGHz». Está 10³ alto
# (es Hz/10⁶, no GHz). Bug del legado congelado por los golden files.
HARTREE_TO_GHZ = 6.579683920502e6

# CODATA 2018: 1 u.a. de campo eléctrico = 5.142 206 747 63e11 V/m.
AU_FIELD_IN_V_PER_M = 5.14220674763e11


def field_au(F_V_per_m: float) -> float:
    """Intensidad de campo en unidades atómicas a partir de V/m."""
    return F_V_per_m / AU_FIELD_IN_V_PER_M


class NeighborRadials:
    """
    R_{nl}(r) y dR_{nl}/dr de los niveles vecinos, de las tablas del legado.

    `rvsR38s.dat`, `rvsR37p.dat`, `rvsR36d.dat` son las funciones de Coulomb
    con n* no entero, tabuladas en r. Se interpolan con spline cúbico y la
    derivada se toma del propio spline (error O(h⁴) con h=3 a₀, ~10⁻⁵
    relativo), en vez de leer `rvsDR*.dat`: la del 38s está mal por 10⁻⁶.
    """

    FILES = {0: "rvsR38s", 1: "rvsR37p", 2: "rvsR36d"}

    def __init__(self, data_dir=None):
        if data_dir is None:
            data_dir = Path(__file__).resolve().parents[4] / "data" / "Wavefunction"
        self.data_dir = Path(data_dir)
        self._spline: Dict[int, CubicSpline] = {}
        for l, stem in self.FILES.items():
            path = self.data_dir / f"{stem}.dat"
            if not path.exists():
                raise FileNotFoundError(f"falta la función de onda tabulada {path}")
            tab = np.loadtxt(path)
            self._spline[l] = CubicSpline(tab[:, 0], tab[:, 1])
        self.r_min = float(min(s.x[0] for s in self._spline.values()))
        self.r_max = float(max(s.x[-1] for s in self._spline.values()))

    def has(self, l: int) -> bool:
        return l in self._spline

    def R(self, l: int, r: float) -> float:
        return float(self._spline[l](r))

    def dR(self, l: int, r: float) -> float:
        return float(self._spline[l](r, 1))


class SymmetricLinearTrimer:
    """
    H(R, F) del trímero lineal simétrico Rb(5s)–Rb*(n,l≥3)–Rb(5s).

    Un objeto por sistema; `hamiltonian(R, F_au, m_l)` devuelve la matriz del
    bloque m_l en hartree, y `energies` sus autovalores en GHz respecto al
    manifold libre de campo.
    """

    def __init__(
        self,
        n_manifold: int = 35,
        s_wave: bool = True,
        p_wave: bool = True,
        radial_source: str = "tabulated",
        n_perturbers: int = 2,
        neighbor_sign: Dict[int, float] = None,
        data_dir=None,
        scattering: ScatteringLengths = None,
    ):
        """
        `radial_source`: "tabulated" (por defecto) usa las funciones de Coulomb
            tabuladas para l ≤ 2; "hydrogenic" usa la hidrogenoide de n entero
            de `RadialBasis`, que es lo que hace el camino de Rb*-KRb.
        `n_perturbers`: 2 = trímero lineal simétrico (θ=0 y θ=π), que es el
            sistema del paper. 1 = dímero con el átomo en θ=0, que es la curva
            de referencia punteada de las Figs. 3, 5 y 7. Nada más cambia:
            misma base, mismo campo, mismo cero de energía.
        `neighbor_sign`: {l: ±1} para invertir el convenio de signo de una
            radial tabulada. Existe para MEDIR la sensibilidad a la
            aproximación 2 de la cabecera, no para producción.
        """
        if n_perturbers not in (1, 2):
            raise ValueError(f"n_perturbers debe ser 1 o 2, no {n_perturbers}")
        if radial_source not in ("tabulated", "hydrogenic"):
            raise ValueError(f"radial_source debe ser 'tabulated' o 'hydrogenic', "
                             f"no {radial_source!r}")
        self.n_manifold = n_manifold
        self.radial_source = radial_source
        self.n_perturbers = n_perturbers
        self.neighbor_sign = dict(neighbor_sign or {})
        if data_dir is None:
            data_dir = Path(__file__).resolve().parents[4] / "data" / "Wavefunction"
        self.data_dir = Path(data_dir)

        self.l_max = n_manifold - 1
        self.neighbor_levels = neighbor_levels(n_manifold)   # {0: n+3, 1: n+2, 2: n+1}

        # Base electrónica pura: el rotor de KRb se apaga con N_max=0, y
        # entonces M_J = m_l exactamente. Ver cabecera del módulo.
        self.basis = CoupledBasis(
            N_max=0, manifold_l_min=3, manifold_l_max=self.l_max, neighbor_l=(0, 1, 2)
        )
        self.radial_basis = RadialBasis(
            n_manifold=n_manifold, l_min=3, l_max=self.l_max, neighbor_l=(0, 1, 2)
        )

        self.neighbors = (NeighborRadials(self.data_dir)
                          if radial_source == "tabulated" else None)

        # Un único k(R) por punto, el del manifold: es el convenio del paper,
        # A_s[k(R)] y A_p[k(R)] con k del número cuántico principal Rydberg.
        # Con n*=35 == n* de la tabla, el remapeo es la identidad.
        self.scattering = scattering if scattering is not None else ScatteringLengths(
            data_dir=self.data_dir, n_star_table=float(n_manifold)
        )
        self.pseudo = FermiPseudopotential(
            self.radial_basis, self.scattering,
            s_wave=s_wave, p_wave=p_wave, n_manifold=n_manifold,
            uniform_n_star=float(n_manifold),
            radial_fn=self._radial, dradial_fn=self._dradial,
        )

        # Energías de un electrón (hartree) por l. l ≥ 3: manifold degenerado.
        self.E_manifold = Atom(n_manifold, 3).E_Rb()
        self._E_of_l = {
            l: Atom(self.neighbor_levels[l], l).E_Rb() if l in self.neighbor_levels
            else self.E_manifold
            for l in range(0, self.l_max + 1)
        }

        # ⟨l|r|l+1⟩, índice l = 0..n-2.
        #
        # Las TRES PRIMERAS son las cruzadas entre vecinos —⟨38s|r|37p⟩,
        # ⟨37p|r|36d⟩, ⟨36d|r|35f⟩— y sólo existen tabuladas: sus n* no son
        # enteros y no hay forma cerrada. Se leen de `exp_val_r.txt`.
        #
        # De l ≥ 3 en adelante son intra-manifold y tienen forma cerrada
        # EXACTA, ⟨n,l|r|n,l+1⟩ = −(3/2)n√(n²−(l+1)²) con el convenio de signo
        # de `hydrogenic_R`. El fichero las reproduce sólo a 3e-4 relativo
        # (residuo de la integración numérica original, ver test T8b), así que
        # se usa la forma cerrada. La diferencia es ~4e-3 GHz a 500 V/m.
        dip = np.loadtxt(self.data_dir / "exp_val_r.txt")
        if dip.shape != (n_manifold - 1,):
            raise ValueError(
                f"exp_val_r.txt tiene {dip.shape} valores; se esperaban "
                f"{n_manifold - 1} (⟨l|r|l+1⟩ para l=0..{n_manifold - 2})"
            )
        self.dipole_table = dip
        self.dipole = dip.copy()
        n = float(n_manifold)
        for l in range(3, n_manifold - 1):
            self.dipole[l] = -1.5 * n * np.sqrt(n**2 - (l + 1) ** 2)
        self._atom = Atom(n_manifold, 0)

    # ------------------------------------------------------------ radiales
    def _neighbor_scale(self, l: int) -> float:
        return float(self.neighbor_sign.get(l, 1.0))

    def _radial(self, l: int, R: float) -> float:
        if self.neighbors is not None and self.neighbors.has(l):
            return self._neighbor_scale(l) * self.neighbors.R(l, R)
        return float(hydrogenic_R(self.radial_basis.n_of_l(l), l, R))

    def _dradial(self, l: int, R: float) -> float:
        if self.neighbors is not None and self.neighbors.has(l):
            return self._neighbor_scale(l) * self.neighbors.dR(l, R)
        return float(hydrogenic_dR(self.radial_basis.n_of_l(l), l, R))

    # -------------------------------------------------------------- bloque
    def l_values(self, m_l: int) -> Sequence[int]:
        """
        Los l del bloque m_l, en orden creciente. Es la base del bloque.

        `CoupledBasis` los enumera manifold-primero (3..n-1, 0, 1, 2); aquí se
        ordenan por l porque la matriz de campo es tridiagonal en l y así se
        lee. El orden es interno y consistente en `diagonal`,
        `pseudopotential_matrix` y `field_matrix`, que llaman todas aquí.
        """
        return sorted(l for (l, m, _N, _MN) in self.basis.get_block(m_l).states
                      if m == m_l)

    # ---------------------------------------------------------- términos H
    def parity_factor(self, l1: int, l2: int) -> float:
        """
        Peso geométrico del par (l₁,l₂): 1 + (−1)^{l₁+l₂} con los dos
        perturbadores —2 si l₁+l₂ es par, 0 si impar, ver cabecera— y 1 con
        uno solo, donde no hay nada que cancelar.
        """
        if self.n_perturbers == 1:
            return 1.0
        return 1.0 + (-1.0) ** (l1 + l2)

    def pseudopotential_matrix(self, R: float, m_l: int) -> np.ndarray:
        """V(R⃗₁) + V(R⃗₂) en el bloque m_l, hartree."""
        ls = self.l_values(m_l)
        V = np.zeros((len(ls), len(ls)))
        for i, l1 in enumerate(ls):
            for j, l2 in enumerate(ls):
                f = self.parity_factor(l1, l2)
                if f == 0.0:
                    continue
                V[i, j] = f * self.pseudo.electron_element(l1, m_l, l2, m_l, R)
        return V

    def field_matrix(self, F_au: float, m_l: int) -> np.ndarray:
        """F·z en el bloque m_l, hartree. Δl=±1, Δm=0."""
        ls = self.l_values(m_l)
        H = np.zeros((len(ls), len(ls)))
        if F_au == 0.0:
            return H
        for i, l1 in enumerate(ls):
            for j, l2 in enumerate(ls):
                if abs(l1 - l2) != 1:
                    continue
                H[i, j] = self._atom.Vfield(
                    l1, l2, m_l, m_l, self.dipole[min(l1, l2)], F_au
                )
        return H

    def diagonal(self, m_l: int) -> np.ndarray:
        """E_{nl} de cada estado del bloque, hartree."""
        return np.array([self._E_of_l[l] for l in self.l_values(m_l)])

    def hamiltonian(self, R: float, F_au: float = 0.0, m_l: int = 0) -> np.ndarray:
        """H del bloque m_l a distancia R y campo F (hartree, simétrica)."""
        H = np.diag(self.diagonal(m_l))
        H += self.pseudopotential_matrix(R, m_l)
        H += self.field_matrix(F_au, m_l)
        return H

    # ------------------------------------------------------------ espectro
    def energies(self, R: float, F_au: float = 0.0, m_l: int = 0) -> np.ndarray:
        """
        Autovalores del bloque m_l en GHz, relativos al manifold sin campo.

        Orden ascendente, como devuelve `eigvalsh`.
        """
        w = np.linalg.eigvalsh(self.hamiltonian(R, F_au, m_l))
        return (w - self.E_manifold) * HARTREE_TO_GHZ

    def spectrum(self, R: float, F_au: float = 0.0, m_l: int = 0):
        """
        (E, W): energías en GHz y peso de manifold de cada autovector.

        W = Σ_{l≥3} |c_l|². Sirve para separar los estados que evolucionan del
        manifold degenerado —los que dibuja el paper— de los que evolucionan de
        38s, 37p y 36d, que caen 20, 102 y 54 GHz por debajo y dominarían
        cualquier ordenación por energía.
        """
        w, v = np.linalg.eigh(self.hamiltonian(R, F_au, m_l))
        ls = np.array(self.l_values(m_l))
        W = (np.abs(v[ls >= 3, :]) ** 2).sum(axis=0)
        return (w - self.E_manifold) * HARTREE_TO_GHZ, W

    def curves(self, R_values, F_au: float = 0.0, m_l: int = 0) -> np.ndarray:
        """Matriz (len(R), dim) de energías en GHz. Filas = puntos de R."""
        return np.array([self.energies(float(R), F_au, m_l) for R in R_values])

    def spectra(self, R_values, F_au: float = 0.0, m_l: int = 0):
        """(E, W) apilados sobre R: dos matrices (len(R), dim)."""
        out = [self.spectrum(float(R), F_au, m_l) for R in R_values]
        return np.array([e for e, _ in out]), np.array([w for _, w in out])

    def table_R(self, R_min: float = None, R_max: float = None) -> np.ndarray:
        """
        La malla nativa de `rvsAS.dat`/`rvsAP.dat`, opcionalmente recortada.

        Evaluar EN los nodos de la tabla elimina toda interpolación de A_s y
        A_p: es la misma elección que hace el camino legado (`R = As[row,0]`).
        """
        R = self.scattering.R_table
        if R_min is not None:
            R = R[R >= R_min]
        if R_max is not None:
            R = R[R <= R_max]
        return R
