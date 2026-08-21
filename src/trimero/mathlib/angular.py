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

    Orden canónico de argumentos: la fórmula de Racah no es simétrica en
    punto flotante bajo permutación de columnas (cada orden produce una
    serie alternante distinta y, por tanto, un redondeo distinto). Para
    garantizar que cualquier permutación devuelva EXACTAMENTE el mismo
    valor (bit a bit) se reordenan los argumentos según una regla única,
    definida sobre los pares (j, m) y no sobre posiciones:

      1. Orientación global de los m: se comparan lexicográficamente las
         secuencias ordenadas de los pares (j, m) con m tal cual y con m
         volteado; la mayor define la orientación canónica. Voltear
         cuesta la fase exacta (-1)^{j1+j2+j3} de la simetría m->-m.
      2. Orden ascendente de los pares (j, m); el par mayor queda en la
         posición 3. Colocar el mayor j al final minimiza j1+j2-j3, es
         decir, acorta la serie alternante y reduce la cancelación
         catastrófica.
      3. La fase de permutación de columnas se cuenta como el número de
         inversiones ESTRICTAS de la secuencia de pares ya orientada (un
         intercambio de pares idénticos no cuenta): impar -> fase
         (-1)^{j1+j2+j3}, aplicada con aritmética entera (factor exacto).

    Con ello la forma canónica y las fases dependen sólo del conjunto de
    pares con su orientación global, no del orden de entrada: cualquier
    permutación de columnas o volteo simultáneo de los tres m comparte la
    MISMA evaluación de Racah bit a bit, y las relaciones entre las 12
    formas equivalentes se satisfacen exactamente. Las simetrías de Regge
    quedan fuera a propósito.
    """
    if m1 + m2 + m3 != 0:
        return 0.0
    if abs(m1) > j1 or abs(m2) > j2 or abs(m3) > j3:
        return 0.0
    if j3 < abs(j1 - j2) or j3 > j1 + j2:
        return 0.0
    if (j1 + j2 + j3) % 2 == 1 and m1 == 0 and m2 == 0 and m3 == 0:
        return 0.0

    pares = ((j1, m1), (j2, m2), (j3, m3))
    seq_p = tuple(sorted(pares))
    pares_m = tuple((j, -m) for (j, m) in pares)
    seq_m = tuple(sorted(pares_m))
    fase_extra = 1.0
    if seq_m > seq_p:
        usar, seq = pares_m, seq_m
        fase_extra = (-1.0) ** (j1 + j2 + j3)
    else:
        usar, seq = pares, seq_p
    (J1, M1), (J2, M2), (J3, M3) = seq
    impar = (
        sum(1 for a in range(3) for b in range(a + 1, 3) if usar[a] > usar[b]) % 2
        == 1
    )

    t_min = max(0, J2 - J3 - M1, J1 - J3 + M2)
    t_max = min(J1 + J2 - J3, J1 - M1, J2 + M2)
    if t_max < t_min:
        return 0.0

    # prefactor: sqrt( Δ(j1j2j3) · Π (j±m)! )
    log_pref = 0.5 * (
        gammaln(J1 + J2 - J3 + 1)
        + gammaln(J1 - J2 + J3 + 1)
        + gammaln(-J1 + J2 + J3 + 1)
        - gammaln(J1 + J2 + J3 + 2)
        + gammaln(J1 + M1 + 1)
        + gammaln(J1 - M1 + 1)
        + gammaln(J2 + M2 + 1)
        + gammaln(J2 - M2 + 1)
        + gammaln(J3 + M3 + 1)
        + gammaln(J3 - M3 + 1)
    )

    t = np.arange(t_min, t_max + 1)
    log_den = (
        gammaln(t + 1)
        + gammaln(J3 - J2 + t + M1 + 1)
        + gammaln(J3 - J1 + t - M2 + 1)
        + gammaln(J1 + J2 - J3 - t + 1)
        + gammaln(J1 - t - M1 + 1)
        + gammaln(J2 - t + M2 + 1)
    )
    # Se factoriza el término dominante antes de sumar para limitar la
    # cancelación catastrófica de la serie alternante.
    shift = log_den.min()
    terms = np.where(t % 2 == 0, 1.0, -1.0) * np.exp(shift - log_den)
    s = terms.sum()
    if s == 0.0:
        return 0.0

    sign = -1.0 if (J1 - J2 - M3) % 2 else 1.0
    val = sign * np.sign(s) * np.exp(log_pref - shift + np.log(abs(s)))
    if impar:
        val *= (-1.0) ** (j1 + j2 + j3)
    return val * fase_extra


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
