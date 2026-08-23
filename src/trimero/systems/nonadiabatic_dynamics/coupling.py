"""
Acoplamiento no adiabático de derivada ⟨Ψᵢ(R)|d/dR|Ψⱼ(R)⟩, genérico.

Módulo de la Fase 1 de `docs/PLAN_nonadiabatic_dynamics.md`. No conoce ningún
sistema físico concreto: recibe una función `solve(R) -> (w, V)` (autovalores
ordenados y autovectores como columnas, mismo orden de base en todo R) o una
malla de autovectores ya calculada, y devuelve el acoplamiento de derivada por
diferencias finitas centradas.

EL PASO MÁS DELICADO: continuidad de fase
------------------------------------------
`numpy.linalg.eigh` no garantiza continuidad de signo entre diagonalizaciones
vecinas — el signo de cada autovector es arbitrario en cada llamada
independiente. Si no se fija antes de derivar, un simple cambio de signo entre
R y R+h se leería como un acoplamiento de derivada gigante y espurio, en vez
del acoplamiento físico real. La fijación se hace de forma SECUENCIAL: en cada
paso, el signo de cada autovector se elige maximizando el solapamiento con el
mismo autovector (misma columna) en el paso anterior (`fix_eigenvector_signs`).
Esto asume no-cruce entre autovalores vecinos del subespacio elegido dentro de
un paso de malla — si dos autovalores se cruzan exactamente entre R y R+h, la
correspondencia por solapamiento columna-a-columna deja de ser válida; con un
paso h suficientemente fino (el que ya usa `trace_curve` para lo mismo) no
ocurre en la práctica.

FÓRMULA Y SU ORDEN DE ERROR
----------------------------
Con Ψⱼ(R) de fase continua, dΨⱼ/dR se aproxima por diferencias centradas:

    dΨⱼ/dR|_{Rₖ} ≈ [Ψⱼ(Rₖ₊₁) − Ψⱼ(Rₖ₋₁)] / (2h)     (orden h²)

y A_ij(Rₖ) = Ψᵢ(Rₖ) · dΨⱼ/dR|_{Rₖ}. La antisimetría A_ij = −A_ji y la
anulación de la diagonal A_ii = 0 son identidades EXACTAS del operador
d/dR sobre una base ortonormal para todo R (se derivan de d/dR⟨Ψᵢ|Ψⱼ⟩=0);
esta fórmula discreta las respeta sólo hasta O(h²), el mismo orden que la
diferencia centrada — por eso los tests de antisimetría y diagonal nula usan
una tolerancia atada a h, no cero exacto, y hay un test de convergencia
explícito que verifica que el residuo cae ~4x al reducir h a la mitad. Ver
docs/analysis_fase1_acoplamiento_derivada.md para la verificación numérica.
"""
from typing import Callable, Optional, Tuple

import numpy as np

__all__ = ["fix_eigenvector_signs", "eigenbasis_along_R", "derivative_coupling"]


def fix_eigenvector_signs(V: np.ndarray, V_ref: np.ndarray) -> np.ndarray:
    """
    Copia de V con el signo de cada columna invertido si su solapamiento con
    la columna correspondiente de V_ref es negativo.

    Args:
        V: (dim, n_states), autovectores a fijar (columnas).
        V_ref: (dim, n_states), referencia de fase (paso anterior de R).

    Returns:
        (dim, n_states) con signo por columna elegido para que
        sum(V[:, j] * V_ref[:, j]) >= 0 para todo j.
    """
    overlap = np.sum(V * V_ref, axis=0)
    signs = np.where(overlap < 0.0, -1.0, 1.0)
    return V * signs


def eigenbasis_along_R(
    R: np.ndarray, solve: Callable[[float], Tuple[np.ndarray, np.ndarray]]
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Diagonaliza `solve(R_k)` en cada punto de la malla y fija la continuidad
    de fase secuencialmente (ver docstring del módulo).

    Args:
        R: malla de R, monótona (creciente o decreciente).
        solve: R -> (w, V); w autovalores, V autovectores en columnas. Mismo
            orden y dimensión de base en todo R.

    Returns:
        R (tal cual se pasó), W (n_R, n_states), V (n_R, dim, n_states) con
        fase continua a lo largo de R.
    """
    R = np.asarray(R, dtype=float)
    if len(R) == 0:
        raise ValueError("R está vacío")

    w0, V0 = solve(float(R[0]))
    w0 = np.asarray(w0, dtype=float)
    V0 = np.asarray(V0, dtype=float)
    dim, n_states = V0.shape

    W = np.empty((len(R), n_states))
    V = np.empty((len(R), dim, n_states))
    W[0], V[0] = w0, V0

    for k in range(1, len(R)):
        wk, Vk = solve(float(R[k]))
        V[k] = fix_eigenvector_signs(np.asarray(Vk, dtype=float), V[k - 1])
        W[k] = np.asarray(wk, dtype=float)

    return R, W, V


def derivative_coupling(
    R: np.ndarray,
    V: Optional[np.ndarray] = None,
    solve: Optional[Callable[[float], Tuple[np.ndarray, np.ndarray]]] = None,
    states: Optional[np.ndarray] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    ⟨Ψᵢ(R)|d/dR|Ψⱼ(R)⟩ por diferencias finitas centradas, en los puntos
    interiores de una malla UNIFORME de R (se necesitan Rₖ₋₁ y Rₖ₊₁).

    Args:
        R: malla de R, uniforme (se verifica explícitamente), len(R) >= 3.
        V: (n_R, dim, n_states) autovectores ya calculados en cada R de la
            malla (columnas = autoestados). Si se da, la continuidad de fase
            se fija aquí igualmente (idempotente si ya venía fijada) — no se
            asume que quien la generó lo hiciera. Mutuamente excluyente con
            `solve`.
        solve: R -> (w, V); se usa junto con `eigenbasis_along_R` para generar
            la malla de autovectores internamente. Mutuamente excluyente
            con `V`.
        states: índices (en la base de `V`/`solve`) de los estados a incluir.
            Por defecto, todos. Restringir reduce el coste si sólo interesa
            un subespacio pequeño (p.ej. dos estados cerca de un cruce).

    Returns:
        R_mid: R[1:-1].
        A: (n_R-2, n_states, n_states), A[k, i, j] =
            ⟨Ψᵢ(R_mid[k])|d/dR|Ψⱼ(R_mid[k])⟩.
    """
    R = np.asarray(R, dtype=float)
    if len(R) < 3:
        raise ValueError(
            "se necesitan al menos 3 puntos de R para diferencias centradas")

    steps = np.diff(R)
    if not np.allclose(steps, steps[0], rtol=1e-9, atol=1e-12):
        raise ValueError(
            "derivative_coupling requiere una malla de R uniforme "
            f"(pasos observados: {steps})")
    h = float(steps[0])

    if V is None and solve is None:
        raise ValueError("hay que dar V precalculado o una función solve")
    if V is not None and solve is not None:
        raise ValueError("V y solve son mutuamente excluyentes")

    if V is None:
        _, _, V = eigenbasis_along_R(R, solve)
    else:
        V = np.asarray(V, dtype=float).copy()
        for k in range(1, len(R)):
            V[k] = fix_eigenvector_signs(V[k], V[k - 1])

    if states is not None:
        V = V[:, :, states]

    n_states = V.shape[2]
    n_mid = len(R) - 2
    A = np.empty((n_mid, n_states, n_states))
    for k in range(n_mid):
        dpsi = (V[k + 2] - V[k]) / (2.0 * h)   # dΨ/dR en R[k+1], orden h²
        A[k] = V[k + 1].T @ dpsi

    return R[1:-1], A
