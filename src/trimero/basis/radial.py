"""
Base radial del electrón Rydberg e integrales radiales de la expansión
multipolar del campo e·(r-R)/|r-R|³.

González-Férez, Sadeghpour & Schmelcher, NJP 17, 013021 (2015), Ec. A.6-A.10.

Con  g_k(r,R) = r_<^k / r_>^{k+1}  (r_< = min(r,R), r_> = max(r,R)), la
expansión de Legendre de 1/|r⃗-R⃗| da lugar a dos familias de integrales
radiales sobre u_a(r) = r·R_{n_a l_a}(r):

    G^k_{ab}(R) = ∫ u_a u_b  g_k(r,R)      dr
    Z^k_{ab}(R) = ∫ u_a u_b  ∂g_k/∂R (r,R) dr

Ambas se escriben en términos de dos piezas acotadas (ningún factor crece sin
control, porque (r/R)^k ≤ 1 dentro y (R/r)^k ≤ 1 fuera):

    Gin_k = (1/R) ∫_{r<R} u_a u_b (r/R)^k dr
    Gout_k =      ∫_{r>R} u_a u_b (R/r)^k dr/r

    G^k  = Gin_k + Gout_k
    Z^k  = [ -(k+1)·Gin_k + k·Gout_k ] / R

APROXIMACIÓN ABIERTA (documentada): los niveles vecinos individuales —(n+3)s,
(n+2)p, (n+1)d— se representan con hidrogenoides de número cuántico principal
efectivo ENTERO, el más próximo a su n*, en lugar de funciones de Coulomb con n*
no entero. Para n=24: 27s → n_eff=24 (n*=23.869), 26p → n_eff=23 (n*=23.351),
25d → n_eff=24 (n*=23.654). El error en la extensión radial es del orden del
1-3 %. La ENERGÍA de esos niveles sí usa el defecto cuántico exacto
(`atom.Atom.E_Rb()`); sólo se aproxima la forma radial.
`scipy.special.hyperu` resulta inestable (NaN) para l=0 con n* no entero, así
que la función de Coulomb exacta queda pendiente.
"""

from typing import Dict, List, Sequence, Tuple

import numpy as np
from scipy.special import eval_genlaguerre, gammaln

from trimero.systems.rb_krb_polar.rb_defects import n_star_nl, neighbor_levels

__all__ = ["RadialBasis"]


class RadialBasis:
    """Funciones radiales u_{nl}(r) en malla e integrales G^k, Z^k cacheadas."""

    def __init__(
        self,
        n_manifold: int = 24,
        l_min: int = 3,
        l_max: int = None,
        neighbor_l=(0, 1, 2),
        neighbor_n_eff: Dict[int, int] = None,
        r_max: float = 3500.0,
        n_points: int = 12000,
    ):
        """
        `neighbor_l` son los niveles individuales de la base: l=0 → (n+3)s,
        l=1 → (n+2)p, l=2 → (n+1)d. Cada uno se representa por la hidrogenoide
        del ENTERO más próximo a su n* (misma aproximación que se venía usando
        para el 27s: n_eff = 24, no 27). `neighbor_n_eff` permite fijar esos
        enteros a mano; si es None se calculan de los defectos cuánticos.
        """
        self.n_manifold = n_manifold
        self.l_min = l_min
        self.l_max = n_manifold - 1 if l_max is None else l_max
        self.neighbor_l = tuple(sorted(neighbor_l))
        if neighbor_n_eff is None:
            levels = neighbor_levels(n_manifold)
            neighbor_n_eff = {l: int(round(n_star_nl(levels[l], l)))
                              for l in self.neighbor_l}
        self.neighbor_n_eff = dict(neighbor_n_eff)
        self.r_max = r_max
        self.n_points = n_points

        # Malla uniforme en √r: resuelve las oscilaciones cerca del núcleo sin
        # desperdiciar puntos en la cola exponencial.
        x = np.linspace(0.0, np.sqrt(r_max), n_points)
        self.r = x * x

        self.l_values: List[int] = list(self.neighbor_l) + list(
            range(l_min, self.l_max + 1))
        self._u_cache: Dict[int, np.ndarray] = {}
        self._parts_cache: Dict[Tuple[int, int, float], Tuple[np.ndarray, np.ndarray, tuple]] = {}

    @property
    def n_s_eff(self) -> int:
        """n entero de la radial del ns. Se conserva por compatibilidad."""
        return self.neighbor_n_eff[0]

    # -- funciones de onda ---------------------------------------------
    def n_of_l(self, l: int) -> int:
        """Número cuántico principal (efectivo) asociado a cada l de la base."""
        return self.neighbor_n_eff.get(l, self.n_manifold)

    def u(self, l: int) -> np.ndarray:
        """u_{n(l),l}(r) = r·R_{n(l),l}(r), normalizada: ∫u² dr = 1."""
        if l not in self._u_cache:
            self._u_cache[l] = self._hydrogenic_u(self.n_of_l(l), l, self.r)
        return self._u_cache[l]

    @staticmethod
    def _hydrogenic_u(n: int, l: int, r: np.ndarray) -> np.ndarray:
        """
        u_{nl}(r) hidrogenoide (Z=1), evaluada en espacio logarítmico para
        evitar el desbordamiento de rho^l con l=23 y rho~250.
        """
        rho = 2.0 * r / n
        log_norm = 0.5 * (
            3.0 * np.log(2.0 / n) + gammaln(n - l) - np.log(2.0 * n) - gammaln(n + l + 1)
        )
        out = np.zeros_like(r)
        pos = r > 0.0
        out[pos] = (
            r[pos]
            * np.exp(log_norm + l * np.log(rho[pos]) - rho[pos] / 2.0)
            * eval_genlaguerre(n - l - 1, 2 * l + 1, rho[pos])
        )
        return out

    # -- integrales radiales -------------------------------------------
    def _parts(self, l1: int, l2: int, R: float, k_values: Sequence[int]):
        """(Gin_k, Gout_k) para los k pedidos. Cacheado por (l1,l2,R,ks)."""
        key = (min(l1, l2), max(l1, l2), float(R), tuple(k_values))
        hit = self._parts_cache.get(key)
        if hit is not None:
            return hit

        r = self.r
        f = self.u(l1) * self.u(l2)
        ks = np.asarray(k_values, dtype=float)

        idx = int(np.searchsorted(r, R))
        f_R = float(np.interp(R, r, f))

        # Se inserta R como nodo exacto para que el corte r≶R no introduzca
        # error de media celda en el punto donde g_k tiene su vértice.
        r_in = np.append(r[:idx], R)
        f_in = np.append(f[:idx], f_R)
        ratio_in = r_in / R                        # ≤ 1
        Gin = np.trapezoid(f_in * ratio_in[None, :] ** ks[:, None], r_in, axis=1) / R

        if idx < len(r):
            r_out = np.concatenate(([R], r[idx:]))
            f_out = np.concatenate(([f_R], f[idx:]))
            ratio_out = R / r_out                  # ≤ 1
            Gout = np.trapezoid(
                f_out * ratio_out[None, :] ** ks[:, None] / r_out, r_out, axis=1
            )
        else:
            Gout = np.zeros_like(Gin)

        self._parts_cache[key] = (Gin, Gout, tuple(k_values))
        return self._parts_cache[key]

    def g_integrals(self, l1: int, l2: int, k_values: Sequence[int], R: float) -> np.ndarray:
        """G^k_{l1 l2}(R) = ∫ u1 u2 g_k dr, vector sobre k_values."""
        Gin, Gout, _ = self._parts(l1, l2, R, k_values)
        return Gin + Gout

    def dg_integrals(self, l1: int, l2: int, k_values: Sequence[int], R: float) -> np.ndarray:
        """Z^k_{l1 l2}(R) = ∫ u1 u2 (∂g_k/∂R) dr, vector sobre k_values."""
        Gin, Gout, _ = self._parts(l1, l2, R, k_values)
        ks = np.asarray(k_values, dtype=float)
        return (-(ks + 1.0) * Gin + ks * Gout) / R
