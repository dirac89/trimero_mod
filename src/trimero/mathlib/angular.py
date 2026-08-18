"""
Álgebra de momento angular: símbolos 3j de Wigner y coeficientes de Gaunt.

Necesarios para la expansión en armónicos esféricos del campo del electrón
Rydberg (Ec. A.6-A.10 de González-Férez, Sadeghpour & Schmelcher, NJP 17,
013021, 2015) y para el elemento de matriz del rotor.

Se implementan aquí (y no se toma de una librería) porque el proyecto no tiene
sympy entre sus dependencias. Ambas funciones se validan en
`test_rydberg_field.py` contra valores analíticos cerrados y contra la regla de
suma de ortogonalidad.

Convenio: armónicos esféricos complejos con fase de Condon-Shortley.
"""

from functools import lru_cache

import numpy as np
from scipy.special import gammaln

__all__ = ["wigner_3j", "gaunt"]


@lru_cache(maxsize=None)
def wigner_3j(j1: int, j2: int, j3: int, m1: int, m2: int, m3: int) -> float:
    """
    Símbolo 3j de Wigner  (j1 j2 j3 ; m1 m2 m3)  por la fórmula de Racah.

    La suma alternante se evalúa en punto flotante con log-factoriales
    (`gammaln`); la precisión resultante se comprueba en los tests con la
    regla de suma de ortogonalidad, que es sensible a cancelaciones.
    """
    if m1 + m2 + m3 != 0:
        return 0.0
    if abs(m1) > j1 or abs(m2) > j2 or abs(m3) > j3:
        return 0.0
    if j3 < abs(j1 - j2) or j3 > j1 + j2:
        return 0.0
    if (j1 + j2 + j3) % 2 == 1 and m1 == 0 and m2 == 0 and m3 == 0:
        return 0.0

    t_min = max(0, j2 - j3 - m1, j1 - j3 + m2)
    t_max = min(j1 + j2 - j3, j1 - m1, j2 + m2)
    if t_max < t_min:
        return 0.0

    # prefactor: sqrt( Δ(j1j2j3) · Π (j±m)! )
    log_pref = 0.5 * (
        gammaln(j1 + j2 - j3 + 1)
        + gammaln(j1 - j2 + j3 + 1)
        + gammaln(-j1 + j2 + j3 + 1)
        - gammaln(j1 + j2 + j3 + 2)
        + gammaln(j1 + m1 + 1)
        + gammaln(j1 - m1 + 1)
        + gammaln(j2 + m2 + 1)
        + gammaln(j2 - m2 + 1)
        + gammaln(j3 + m3 + 1)
        + gammaln(j3 - m3 + 1)
    )

    t = np.arange(t_min, t_max + 1)
    log_den = (
        gammaln(t + 1)
        + gammaln(j3 - j2 + t + m1 + 1)
        + gammaln(j3 - j1 + t - m2 + 1)
        + gammaln(j1 + j2 - j3 - t + 1)
        + gammaln(j1 - t - m1 + 1)
        + gammaln(j2 - t + m2 + 1)
    )
    # Se factoriza el término dominante antes de sumar para limitar la
    # cancelación catastrófica de la serie alternante.
    shift = log_den.min()
    terms = np.where(t % 2 == 0, 1.0, -1.0) * np.exp(shift - log_den)
    s = terms.sum()
    if s == 0.0:
        return 0.0

    sign = -1.0 if (j1 - j2 - m3) % 2 else 1.0
    return sign * np.sign(s) * np.exp(log_pref - shift + np.log(abs(s)))


@lru_cache(maxsize=None)
def gaunt(l1: int, m1: int, l2: int, m2: int, l3: int, m3: int) -> float:
    """
    Coeficiente de Gaunt  ∫ Y*_{l1 m1} Y_{l2 m2} Y_{l3 m3} dΩ.

        = (-1)^{m1} sqrt((2l1+1)(2l2+1)(2l3+1)/4π)
          · (l1 l2 l3; 0 0 0) · (l1 l2 l3; -m1 m2 m3)

    Nótese el complejo conjugado en el PRIMER armónico (convenio de elemento
    de matriz ⟨l1 m1| O |l3 m3⟩).
    """
    if m2 + m3 != m1:
        return 0.0
    w0 = wigner_3j(l1, l2, l3, 0, 0, 0)
    if w0 == 0.0:
        return 0.0
    wm = wigner_3j(l1, l2, l3, -m1, m2, m3)
    if wm == 0.0:
        return 0.0
    pref = np.sqrt((2 * l1 + 1) * (2 * l2 + 1) * (2 * l3 + 1) / (4.0 * np.pi))
    sign = -1.0 if m1 % 2 else 1.0
    return sign * pref * w0 * wm
