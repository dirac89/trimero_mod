"""
Pseudopotencial de Fermi V(r) para Rb*-KRb sobre la base acoplada CoupledBasis.

González-Férez, Sadeghpour & Schmelcher, NJP 17, 013021 (2015).

    V(r⃗) = 2π A_s[k(R)] δ³(r⃗-R⃗)  +  6π A_p[k(R)] δ³(r⃗-R⃗) ∇⃖·∇⃗

POR QUÉ UN MÓDULO NUEVO Y NO `fermi.FermiPotentials`
----------------------------------------------------
`FermiPotentials` se escribió para el modelo de Aguilera-Fernández 2016 (Rb* con
DOS perturbadores neutros) y no es reutilizable aquí sin arrastrar problemas:

  1. Usa `mathlib.special.Spherical`, que es `lpmv` CRUDO: le falta la
     normalización sqrt((2l+1)/4π·(l-m)!/(l+m)!) del armónico esférico. Sus
     términos `Vs` y `VpA` van sin normalizar mientras `VpB`/`VpC` (vía `DOlm`
     y `DPhilm`) sí la llevan. Es internamente inconsistente y el factor
     depende de l, así que no puede absorberse fuera (`trimer.py` suma
     `Vsp()` directamente).
  2. Su lógica de ramas `l<=2` está atada a funciones de onda externas
     (`wave1`, `wave2`) del manifold n=35 del código de 2016.
  3. Corregirla cambiaría los resultados de `trimer.py`, que está congelado por
     los golden files y que el usuario pidió no tocar.

Este módulo implementa el pseudopotencial de cero desde la geometría de esta
base, y se valida contra ψ y ∇ψ evaluados numéricamente (test F2).

GEOMETRÍA: PERTURBADOR SOBRE EL EJE
-----------------------------------
El eje de cuantización es Z ∥ R⃗ (mismo convenio de las rondas del ion y del
electrón), luego el perturbador está en θ=0, donde
Y_lm(0,φ) = sqrt((2l+1)/4π)·δ_{m0}. Con ∇ψ_{lm}|_{θ=0} no nulo sólo para
|m| ≤ 1, las formas cerradas son:

    ⟨l₁ 0|V_s|l₂ 0⟩   = (A_s/2)·R_{l₁}(R)·R_{l₂}(R)·√((2l₁+1)(2l₂+1))
    ⟨l₁ 0|V_p|l₂ 0⟩   = (3/2)·A_p·R'_{l₁}(R)·R'_{l₂}(R)·√((2l₁+1)(2l₂+1))
    ⟨l₁ ±1|V_p|l₂ ±1⟩ = (3/4)·A_p·R_{l₁}(R)R_{l₂}(R)/R²
                          ·√(l₁(l₁+1)l₂(l₂+1)(2l₁+1)(2l₂+1))

y cero en cualquier otro caso (Δm_l ≠ 0, o |m_l| ≥ 2, o V_s con m ≠ 0).

V es un operador PURAMENTE ELECTRÓNICO: diagonal en (N, M_N). Conserva M_J por
ausencia de acoplamiento al rotor, no por una regla de selección de momento
angular como en H_mol. Mezcla l fuertemente: ahí está la ligadura ULRM.

APROXIMACIÓN 1 — KRb como perturbador puntual
---------------------------------------------
El pseudopotencial trata a KRb como un dispersor SIN estructura interna,
caracterizado por una única pareja (A_s, A_p): la aproximación de "core
molecular puntual". El paper de 2015 no discute explícitamente su validez para
un dímero polar. Marcado como aproximación explícita, no dado por hecho.

APROXIMACIÓN 2 — remapeo k(R) de las tablas de n=35
---------------------------------------------------
Ver `ScatteringLengths`. Decisión del usuario: opción (i).
"""

from pathlib import Path
from typing import Dict, Sequence, Tuple

import numpy as np
from scipy.special import eval_genlaguerre, gammaln

from trimero.basis.quantum import QuantumBasisBlock
from trimero.basis.radial import RadialBasis
from trimero.systems.rb_krb_polar.rb_defects import n_star_nl, neighbor_levels

__all__ = [
    "N_STAR_TABLE",
    "E_TABLE_HARTREE",
    "n_star_of_l",
    "hydrogenic_R",
    "hydrogenic_dR",
    "ScatteringLengths",
    "FermiPseudopotential",
]

# La tabla original se generó sobre el manifold n=35 del código de 2016.
N_STAR_TABLE = 35.0
E_TABLE_HARTREE = -0.5 / N_STAR_TABLE**2

# Manifold por defecto: el de González-Férez 2015 (n=24, con 27s/26p/25d).
N_MANIFOLD_DEFAULT = 24

State = Tuple[int, int, int, int]


def n_star_of_l(l: int, delta0_ns: float = None,
                n_manifold: int = N_MANIFOLD_DEFAULT,
                n_s: int = None) -> float:
    """
    n* efectivo del estado electrónico asociado a cada l de la base.

    l = 0, 1, 2 son los tres vecinos individuales —(n+3)s, (n+2)p, (n+1)d— y se
    usa su n* con defecto cuántico EXACTO (23.8679 para 27s, 23.3512 para 26p,
    23.6538 para 25d), no el entero. Ojo: `RadialBasis.n_of_l` sí devuelve un
    entero, porque ahí la función radial se aproxima por una hidrogenoide de n
    entero (limitación ya documentada en `basis/radial.py`). Son dos cosas
    distintas a propósito: el remapeo k(R) usa la energía precisa, la forma
    radial usa la aproximada.

    Para l ≥ 3 el defecto cuántico es despreciable: n* = n_manifold exacto.
    `n_s` fuerza el nivel de l=0 (compatibilidad); si es None sale de
    `rb_defects.neighbor_levels`.
    """
    levels = neighbor_levels(n_manifold)
    if n_s is not None:
        levels[0] = n_s
    if l not in levels:
        return float(n_manifold)
    return n_star_nl(levels[l], l, delta0_ns)


# ----------------------------------------------------------- radiales
def hydrogenic_R(n: int, l: int, r):
    """R_{nl}(r) hidrogenoide (Z=1), evaluada en espacio logarítmico."""
    r = np.asarray(r, dtype=float)
    rho = 2.0 * r / n
    log_norm = 0.5 * (
        3.0 * np.log(2.0 / n) + gammaln(n - l) - np.log(2.0 * n) - gammaln(n + l + 1)
    )
    out = np.exp(log_norm + l * np.log(rho) - rho / 2.0) * eval_genlaguerre(
        n - l - 1, 2 * l + 1, rho
    )
    return float(out) if out.ndim == 0 else out


def hydrogenic_dR(n: int, l: int, r):
    """
    dR_{nl}/dr analítica.

        dR/dρ = N e^{-ρ/2} ρ^{l-1} [ (l - ρ/2)·L^{2l+1}_{n-l-1} - ρ·L^{2l+2}_{n-l-2} ]
        dR/dr = (2/n)·dR/dρ            (ρ = 2r/n,  dL^α_k/dρ = -L^{α+1}_{k-1})
    """
    r = np.asarray(r, dtype=float)
    rho = 2.0 * r / n
    log_norm = 0.5 * (
        3.0 * np.log(2.0 / n) + gammaln(n - l) - np.log(2.0 * n) - gammaln(n + l + 1)
    )
    lag1 = eval_genlaguerre(n - l - 1, 2 * l + 1, rho)
    # Para l = n-1 el segundo Laguerre tendría grado -1: no existe, vale 0.
    lag2 = (
        eval_genlaguerre(n - l - 2, 2 * l + 2, rho)
        if n - l - 2 >= 0
        else np.zeros_like(rho)
    )
    pref = np.exp(log_norm - rho / 2.0 + (l - 1) * np.log(rho))
    out = (2.0 / n) * pref * ((l - rho / 2.0) * lag1 - rho * lag2)
    return float(out) if out.ndim == 0 else out


# ----------------------------------------------------------- tablas A_s, A_p
class ScatteringLengths:
    """
    Longitudes de dispersión A_s(k), A_p(k) leídas de `rvsAS.dat`/`rvsAP.dat`
    con REMAPEO SEMICLÁSICO en k.

    ⚠️ APROXIMACIÓN EXPLÍCITA, NO EQUIVALENCIA EXACTA ⚠️

    Las tablas están tabuladas en R para el manifold n=35 del código de 2016
    (su R máximo, 2448 a₀, coincide con el punto de retorno clásico externo
    2n² = 2450 a₀ de n=35). La longitud de dispersión depende del momento del
    electrón incidente, no de R, y la relación entre ambos,

        k²(R)/2 = E_n + 1/R,     E_n = -1/(2 n*²)

    mete n* dentro de la tabla. Para usarlas con n*=24 / n*=23.869 se invierte:

        R  --(n* actual)-->  k  --(n*=35)-->  R' = 1/(k²/2 - E_35)

    y se interpola A_s(R'), A_p(R') en la tabla original.

    Esto ASUME que la columna 2 es A(k) calculada sólo con la energía
    semiclásica, sin dependencia adicional en l ni en la estructura del
    perturbador. La procedencia de las tablas sigue sin verificarse
    (`docs/analysis_procedencia_rvsAS_rvsAP.md`): el script generador nunca
    apareció. Es una inferencia razonada, no un hecho comprobado.

    PARA SUSTITUIRLO por la opción (ii) (tablas regeneradas para n=24/27),
    basta con reemplazar `remap_R` por la identidad y apuntar `data_dir` a las
    tablas nuevas; nada más en el módulo depende del remapeo.

    `enabled=False` fuerza A_s = A_p = 0 (usado por el test del límite de
    acoplamiento nulo).

    INTERPOLACIÓN DE A_p A TRAVÉS DE LA RESONANCIA
    ----------------------------------------------
    `A_p(k) = -tan δ_p(k)/k³` DIVERGE donde δ_p = π/2 (resonancia de forma).
    La expansión de rango efectivo de onda p dice que

        k³ cot δ_p(k) = -1/a_p + ½ r_p k² + O(k⁴)

    es analítica en k², y como `1/A_p = -k³ cot δ_p`, resulta que **1/A_p ES la
    función de rango efectivo**: suave a través de la resonancia, donde tiene un
    CERO simple en vez de un polo. Con forma Breit-Wigner,
    `tan δ_p = (Γ/2)/(E_r-E)`:

        A_p   = -(Γ/2)/[k³(E_r-E)]     -> polo simple en E=E_r
        1/A_p = -k³(E_r-E)/(Γ/2)       -> cero simple, lineal en E

    Por eso `p_interpolation="inverse"` (por defecto) interpola 1/A_p en la
    energía cinética local del electrón ε = k²/2, y luego invierte. Comprobado
    empíricamente sobre 20 nodos que abarcan la resonancia: ajuste lineal en ε
    da R² = 0.988 para 1/A_p frente a R² = 0.294 para A_p.

    `p_interpolation="linear"` recupera el comportamiento anterior (interpolar
    A_p directamente en R'), para poder comparar.

    A_s NO se invierte nunca: no tiene polo (max|A_s| = 16 a₀) pero SÍ cruza
    cero (mínimo de Ramsauer en R'=429 a₀), donde 1/A_s sería singular.
    """

    def __init__(self, data_dir=None, enabled: bool = True,
                 n_star_table: float = N_STAR_TABLE,
                 p_interpolation: str = "inverse"):
        if p_interpolation not in ("inverse", "linear"):
            raise ValueError(f"p_interpolation debe ser 'inverse' o 'linear', "
                             f"no {p_interpolation!r}")
        self.enabled = enabled
        self.p_interpolation = p_interpolation
        self.n_star_table = n_star_table
        self.E_table = -0.5 / n_star_table**2
        if data_dir is None:
            data_dir = Path(__file__).resolve().parents[4] / "data" / "Wavefunction"
        self.data_dir = Path(data_dir)
        As = np.loadtxt(self.data_dir / "rvsAS.dat")
        Ap = np.loadtxt(self.data_dir / "rvsAP.dat")
        if not np.array_equal(As[:, 0], Ap[:, 0]):
            raise ValueError("rvsAS.dat y rvsAP.dat no comparten la malla en R")
        self.R_table = As[:, 0]
        self.A_s_table = As[:, 1]
        self.A_p_table = Ap[:, 1]
        if np.any(self.A_p_table == 0.0):
            raise ValueError(
                "rvsAP.dat contiene A_p = 0 exactamente: 1/A_p sería singular. "
                "Usa p_interpolation='linear' o corrige la tabla."
            )
        # Energía cinética local del electrón en cada nodo: eps = k²/2.
        # Es la variable en la que 1/A_p es analítica (rango efectivo).
        self.eps_table = self.E_table + 1.0 / self.R_table
        order = np.argsort(self.eps_table)
        self._eps_asc = self.eps_table[order]
        self._inv_ap_asc = 1.0 / self.A_p_table[order]

    # -- relación semiclásica -----------------------------------------
    @staticmethod
    def energy_of_n_star(n_star: float) -> float:
        return -0.5 / n_star**2

    def k_of_R(self, R: float, n_star: float) -> float:
        """k(R) del manifold actual. NaN si R es clásicamente prohibido."""
        val = 2.0 * (self.energy_of_n_star(n_star) + 1.0 / R)
        return float(np.sqrt(val)) if val > 0.0 else float("nan")

    def remap_R_from_energy(self, R: float, E: float) -> float:
        """
        R' en la tabla n=35 con el mismo k. Inf si no hay solución.

        Si la energía pedida ES la de la tabla (caso n*=35, el manifold para el
        que se generó), el remapeo es la IDENTIDAD algebraica R' = 1/(1/R) = R.
        Se cortocircuita para devolver R exactamente: sin el atajo, la ida y
        vuelta en coma flotante mueve el extremo R=2448 fuera de la tabla por
        1 ulp y `scattering_from_energy` lo rechazaría. No es un cambio de
        física, es evitar reconstruir un número que ya se tiene.
        """
        if E == self.E_table:
            return float(R)
        denom = E + 1.0 / R - self.E_table
        return float(1.0 / denom) if denom > 0.0 else float("inf")

    def remap_R(self, R: float, n_star: float) -> float:
        return self.remap_R_from_energy(R, self.energy_of_n_star(n_star))

    def in_domain(self, R: float, n_star: float) -> bool:
        """True si R es clásicamente permitido Y su R' cae dentro de la tabla."""
        if not np.isfinite(self.k_of_R(R, n_star)):
            return False
        Rp = self.remap_R(R, n_star)
        return bool(self.R_table.min() <= Rp <= self.R_table.max())

    # -- lectura -------------------------------------------------------
    def scattering_from_energy(self, R: float, E: float) -> Tuple[float, float]:
        """(A_s, A_p) para una energía electrónica efectiva E. Nunca extrapola."""
        if not self.enabled:
            return 0.0, 0.0
        if E + 1.0 / R <= 0.0:
            raise ValueError(
                f"R={R} es clásicamente prohibido para E={E:.6e} E_h "
                f"(punto de retorno externo en R={-1.0/E:.1f} a0): k² < 0, "
                "no hay longitud de dispersión definida"
            )
        Rp = self.remap_R_from_energy(R, E)
        lo, hi = self.R_table.min(), self.R_table.max()
        if not (lo <= Rp <= hi):
            raise ValueError(
                f"el remapeo de R={R} da R'={Rp:.2f}, fuera de la tabla "
                f"[{lo:.0f}, {hi:.0f}] a0. No se extrapola."
            )
        A_s = float(np.interp(Rp, self.R_table, self.A_s_table))
        if self.p_interpolation == "linear":
            A_p = float(np.interp(Rp, self.R_table, self.A_p_table))
        else:
            # 1/A_p es suave a través del polo; se interpola y se invierte.
            inv = float(np.interp(E + 1.0 / R, self._eps_asc, self._inv_ap_asc))
            A_p = float("inf") if inv == 0.0 else 1.0 / inv
        return A_s, A_p

    def scattering(self, R: float, n_star: float) -> Tuple[float, float]:
        return self.scattering_from_energy(R, self.energy_of_n_star(n_star))

    def scattering_pair(self, R: float, n_star_1: float, n_star_2: float):
        """
        (A_s, A_p) para un elemento de matriz entre dos estados de energías
        distintas. Se usa la energía MEDIA: es simétrica en i↔j, así que la
        matriz sale exactamente simétrica, y en la diagonal se reduce al valor
        del propio estado. Es una elección, no un resultado del paper; sólo
        importa entre el 27s y el manifold, que difieren 1 % en n*.
        """
        E = 0.5 * (self.energy_of_n_star(n_star_1) + self.energy_of_n_star(n_star_2))
        return self.scattering_from_energy(R, E)

    # -- resonancia de forma p: localización y ventana de exclusión ------
    def p_resonance_window(self, n_star: float,
                           margin_factor: float = 2.0) -> Dict[str, float]:
        """
        Ventana de R a excluir de cualquier comparación cuantitativa, por caer
        sobre la resonancia de forma p, donde el pseudopotencial de rango cero
        diverge (limitación física conocida del modelo, no implementamos la
        corrección de rango efectivo de Omont 1977).

        Criterio, todo calculado de la tabla, nada a mano:

        1. **Centro** = polo de A_p, definido como el CERO del interpolante
           lineal de `1/A_p` en ε=k²/2 entre los dos nodos que cambian de signo.
           Es literalmente el polo que ve el modo `p_interpolation="inverse"`,
           no un punto medio del hueco de malla.
        2. **Anchura** = FWHM de |A_p|: extensión en R de los nodos con
           |A_p| >= max|A_p|/2. Ojo: max|A_p| es finito sólo porque la malla no
           cae sobre el polo, así que esta FWHM es una escala operativa, no la
           anchura de la resonancia en sentido estricto (Γ).
        3. **Ventana** = R_polo ± margin_factor × FWHM.

        Devuelve las tres cosas para poder documentarlas, no sólo los límites.
        """
        E = self.energy_of_n_star(n_star)
        eps, inv = self._eps_asc, self._inv_ap_asc
        A_p = 1.0 / inv

        def R_of_eps(e):
            if e <= E:
                raise ValueError(f"ε={e:.6e} <= E({n_star})={E:.6e}: sin R real")
            return 1.0 / (e - E)

        cross = np.flatnonzero(np.diff(np.sign(A_p)) != 0)
        if len(cross) != 1:
            raise ValueError(f"se esperaba UN cambio de signo de A_p, hay {len(cross)}")
        i = int(cross[0])
        # cero del interpolante lineal de 1/A_p en ε
        eps_pole = eps[i] - inv[i] * (eps[i + 1] - eps[i]) / (inv[i + 1] - inv[i])

        mask = np.abs(A_p) >= 0.5 * np.abs(A_p).max()
        R_hi_fwhm = R_of_eps(float(eps[mask].min()))   # ε menor -> R mayor
        R_lo_fwhm = R_of_eps(float(eps[mask].max()))
        fwhm = R_hi_fwhm - R_lo_fwhm

        R_pole = R_of_eps(eps_pole)
        return {
            "eps_pole": float(eps_pole),
            "R_pole": float(R_pole),
            "fwhm": float(fwhm),
            "R_fwhm_lo": float(R_lo_fwhm),
            "R_fwhm_hi": float(R_hi_fwhm),
            "margin_factor": float(margin_factor),
            "R_lo": float(R_pole - margin_factor * fwhm),
            "R_hi": float(R_pole + margin_factor * fwhm),
        }


# ----------------------------------------------------------- pseudopotencial
class FermiPseudopotential:
    """V(r) sobre CoupledBasis. Ver la cabecera del módulo para las fórmulas."""

    def __init__(self, radial: RadialBasis, scattering: ScatteringLengths,
                 s_wave: bool = True, p_wave: bool = True,
                 delta0_ns: float = None,
                 n_manifold: int = N_MANIFOLD_DEFAULT,
                 n_s: int = None,
                 uniform_n_star: float = None,
                 radial_fn=None, dradial_fn=None):
        """
        `uniform_n_star`: si se da, TODOS los elementos de matriz usan ese n*
            para evaluar k(R), en vez de la energía media del par (l1,l2).
            Es el convenio del paper: A_s[k(R)], A_p[k(R)] con un único k por
            R, el del número cuántico principal del estado Rydberg. Con
            `uniform_n_star == scattering.n_star_table` el remapeo es la
            identidad y la tabla se lee tal cual fue generada.

        `radial_fn(l, R)` / `dradial_fn(l, R)`: reemplazan R_{n(l),l}(R) y su
            derivada por funciones externas. Sirven para inyectar funciones de
            Coulomb tabuladas para los niveles vecinos, cuyo n* no es entero y
            que la hidrogenoide de n entero sólo aproxima en amplitud (la fase
            se desvía). None -> hidrogenoide, comportamiento de siempre.
        """
        self.radial = radial
        self.scattering = scattering
        self.s_wave = s_wave
        self.p_wave = p_wave
        self.uniform_n_star = uniform_n_star
        self.radial_fn = radial_fn
        self.dradial_fn = dradial_fn
        # None -> n* del ns tal como está en rb_atom.py. Ver systems/rb_krb_polar/rb_defects.py
        self.delta0_ns = delta0_ns
        # Manifold sobre el que se evalúa k(R). Por defecto n=24 + 27s.
        self.n_manifold = n_manifold
        self.n_s = n_s
        self._cache: Dict[tuple, float] = {}

    def n_star_of_l(self, l: int) -> float:
        """n* del estado de este l, con el manifold configurado."""
        return n_star_of_l(l, self.delta0_ns, self.n_manifold, self.n_s)

    def n_star_manifold(self) -> float:
        """n* del manifold (l ≥ 3), que es el que fija la coordenada R."""
        return float(self.n_manifold)

    def large_gaps(self, factor: float = 3.0):
        """
        Intervalos (R'_lo, R'_hi) donde la malla de la tabla tiene un paso
        anómalamente grande.

        Importa porque `np.interp` los cruza en línea recta sin avisar. En
        `rvsAP.dat` hay un hueco de 15 a0 (frente a 3 a0 de mediana) justo
        sobre el polo de la resonancia de forma p: interpolar ahí une dos
        puntos a lados opuestos de una divergencia, y el valor resultante no
        es una lectura de la tabla.
        """
        R = self.scattering.R_table
        steps = np.diff(R)
        thr = factor * float(np.median(steps))
        return [(float(R[i]), float(R[i + 1]))
                for i in np.flatnonzero(steps > thr)]

    def bridges_gap(self, R: float, factor: float = 3.0) -> bool:
        """True si el R' remapeado de este R cae dentro de un hueco de malla."""
        ns = [self.n_star_of_l(l) for l in self.radial.l_values]
        gaps = self.large_gaps(factor)
        # Basta con que UN par (l1,l2) caiga en el hueco: ese elemento de
        # matriz ya lleva un A_p interpolado a través de la divergencia.
        for a in ns:
            for b in ns:
                E = 0.5 * (-0.5 / a**2 + -0.5 / b**2)
                Rp = self.scattering.remap_R_from_energy(R, E)
                if any(lo < Rp < hi for lo, hi in gaps):
                    return True
        return False

    def domain_bounds(self):
        """
        Ventana [R_min, R_max] en la que TODOS los elementos de matriz tienen
        longitud de dispersión definida.

        No basta con mirar el manifold: cada par (l1,l2) usa la energía media
        de sus dos estados, y el par más ligado (27s-27s) tiene el punto de
        retorno clásico más corto de todos. El límite lo fija el par más
        restrictivo, no uno elegido a mano.
        """
        ns = [self.n_star_of_l(l) for l in self.radial.l_values]
        E_t = self.scattering.E_table
        lo_t, hi_t = self.scattering.R_table.min(), self.scattering.R_table.max()
        R_min, R_max = -np.inf, np.inf
        for a in ns:
            for b in ns:
                E = 0.5 * (-0.5 / a**2 + -0.5 / b**2)
                # R' = 1/(E + 1/R - E_t) debe caer en [lo_t, hi_t]
                R_max = min(R_max, 1.0 / (1.0 / hi_t + E_t - E))
                R_min = max(R_min, 1.0 / (1.0 / lo_t + E_t - E))
        # Margen de 1e-9 relativo (~1e-6 a0, físicamente irrelevante): el
        # extremo exacto da R' = 2448.0000000001 al reconstruirlo en coma
        # flotante y `scattering_from_energy` lo rechazaría por el borde.
        return float(R_min * (1.0 + 1e-9)), float(R_max * (1.0 - 1e-9))

    def R_of(self, l: int, R: float) -> float:
        """R_{n(l),l}(R) en el punto del perturbador."""
        if self.radial_fn is not None:
            return float(self.radial_fn(l, R))
        return float(hydrogenic_R(self.radial.n_of_l(l), l, R))

    def dR_of(self, l: int, R: float) -> float:
        """dR_{n(l),l}/dr evaluada en el punto del perturbador."""
        if self.dradial_fn is not None:
            return float(self.dradial_fn(l, R))
        return float(hydrogenic_dR(self.radial.n_of_l(l), l, R))

    def electron_element(self, l1: int, m1: int, l2: int, m2: int, R: float) -> float:
        """⟨l₁ m₁| V |l₂ m₂⟩, parte electrónica."""
        if m1 != m2 or abs(m1) > 1:
            return 0.0
        key = (l1, m1, l2, m2, R)
        hit = self._cache.get(key)
        if hit is not None:
            return hit

        if self.uniform_n_star is not None:
            A_s, A_p = self.scattering.scattering(R, self.uniform_n_star)
        else:
            A_s, A_p = self.scattering.scattering_pair(
                R, self.n_star_of_l(l1), self.n_star_of_l(l2)
            )
        norm = np.sqrt((2.0 * l1 + 1.0) * (2.0 * l2 + 1.0))
        val = 0.0
        if m1 == 0:
            if self.s_wave:
                val += 0.5 * A_s * self.R_of(l1, R) * self.R_of(l2, R) * norm
            if self.p_wave:
                val += 1.5 * A_p * self.dR_of(l1, R) * self.dR_of(l2, R) * norm
        else:
            if self.p_wave:
                val += (
                    0.75 * A_p
                    * self.R_of(l1, R) * self.R_of(l2, R) / (R * R)
                    * np.sqrt(l1 * (l1 + 1.0) * l2 * (l2 + 1.0)) * norm
                )
        val = float(val)
        self._cache[key] = val
        return val

    def _matrix_element_unchecked(self, si: State, sj: State, R: float) -> float:
        """Sin comprobar M_J: existe para VERIFICAR el cero, no para producción."""
        l1, m1, N1, MN1 = si
        l2, m2, N2, MN2 = sj
        if N1 != N2 or MN1 != MN2:      # V es puramente electrónico
            return 0.0
        return self.electron_element(l1, m1, l2, m2, R)

    def matrix_element(self, si: State, sj: State, R: float) -> float:
        MJ_i, MJ_j = si[1] + si[3], sj[1] + sj[3]
        if MJ_i != MJ_j:
            raise ValueError(
                f"V_Fermi conserva M_J: elemento pedido entre M_J={MJ_i} "
                f"({si}) y M_J={MJ_j} ({sj})"
            )
        return self._matrix_element_unchecked(si, sj, R)

    def build_reference(self, block: QuantumBasisBlock, R: float) -> np.ndarray:
        """Barrido dim² completo. Para bloques pequeños y para chequear `build`."""
        states: Sequence[State] = block.states
        H = np.zeros((len(states), len(states)), dtype=np.float64)
        for i, si in enumerate(states):
            for j, sj in enumerate(states):
                H[i, j] = self.matrix_element(si, sj, R)
        return H

    def build(self, block: QuantumBasisBlock, R: float) -> np.ndarray:
        """
        Sólo visita los pares permitidos: mismo (m_l, N, M_N) con |m_l| ≤ 1 y
        l libre. El resto es cero exactamente.
        """
        states: Sequence[State] = block.states
        by_key: Dict[tuple, list] = {}
        for i, (l, m, N, MN) in enumerate(states):
            if abs(m) <= 1:
                by_key.setdefault((m, N, MN), []).append(i)
        H = np.zeros((len(states), len(states)), dtype=np.float64)
        for idxs in by_key.values():
            for i in idxs:
                for j in idxs:
                    H[i, j] = self.matrix_element(states[i], states[j], R)
        return H
