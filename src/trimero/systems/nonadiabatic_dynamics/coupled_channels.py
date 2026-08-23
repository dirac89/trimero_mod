"""
Ecuación radial de canales acoplados en representación adiabática (Ec. 8 de
Mellado-Alcedo, Guttridge, Cornish, Sadeghpour & González-Férez, PRA 110,
013314, 2024), Fase 2 de `docs/PLAN_nonadiabatic_dynamics.md`.

Módulo GENÉRICO (no conoce ningún sistema físico concreto), como
`coupling.py` de la Fase 1. Unidades atómicas: energía en Hartree, R en a₀,
masa en m_e, ħ=1 (parámetro por si se quiere usar con otras unidades).

ECUACIÓN QUE SE RESUELVE
-------------------------
La función de onda molecular total se expande en la base electrónica
adiabática {Ψᵢ(r;R)} (autovectores de H_el(R), ya calculados en la Fase 1):

    Ψ(r,R) = Σᵢ χᵢ(R) Ψᵢ(r;R)

Sustituyendo en la ecuación de Schrödinger completa y proyectando sobre
⟨Ψᵢ(r;R)|, el operador cinético nuclear T = -ħ²/(2μ) d²/dR² actúa sobre el
producto χⱼ(R)Ψⱼ(r;R) y genera, además del término diagonal esperado, los
acoplamientos no adiabáticos de la Fase 1:

    d²/dR² [χⱼΨⱼ] = χⱼ''Ψⱼ + 2χⱼ'Ψⱼ' + χⱼΨⱼ''

Proyectando sobre Ψᵢ y usando A_ij = ⟨Ψᵢ|Ψⱼ'⟩ (Fase 1) y B_ij = ⟨Ψᵢ|Ψⱼ''⟩
(nuevo en este módulo, `second_derivative_coupling`):

    -ħ²/(2μ) Σⱼ [δ_ij χⱼ'' + 2A_ij χⱼ' + B_ij χⱼ] + Vᵢ(R) χᵢ = E χᵢ

Nótese que B_ii ≠ 0 EN GENERAL, aunque A_ii=0 (Fase 1, test b): derivando dos
veces la normalización ⟨Ψᵢ|Ψᵢ⟩=1 sale ⟨Ψᵢ'|Ψᵢ'⟩ + ⟨Ψᵢ|Ψᵢ''⟩ = 0, es decir
B_ii = -Σₖ|A_ik|² ≤ 0: es la corrección diagonal de Born-Oppenheimer (DBOC),
un término real que NO se descarta aquí (regla del proyecto: "no se fuerza
ninguna limitación conocida a desaparecer"). `second_derivative_coupling` la
calcula igual que cualquier otro elemento de B, sin tratamiento especial.

ELECCIÓN DE DISCRETIZACIÓN: DIFERENCIAS FINITAS, NO DVR
----------------------------------------------------------
Se decide ANTES de implementar, por tres razones:

1. Consistencia de malla. A_ij y B_ij de este módulo (y de la Fase 1) sólo
   existen como diferencias finitas centradas sobre una malla de R UNIFORME:
   no son funciones continuas que se puedan proyectar sobre una base DVR sin
   interpolar (e introducir un error adicional no caracterizado). Con
   diferencias finitas, el operador cinético nuclear T, el operador de
   primera derivada D1 y los acoplamientos A_ij, B_ij viven en la MISMA malla,
   sin pasos de interpolación.
2. Lo que realmente se conoce de Vᵢ(R) son valores en R discretos, sacados de
   `sysm.solve(R)` para cada R de la malla (una diagonalización electrónica
   por punto) — no una función analítica que un DVR pudiera explotar mejor
   para converger espectralmente. La ventaja característica de un DVR
   (convergencia espectral con potenciales suaves conocidos en todo punto)
   no aplica aquí: la información de entrada ya es discreta.
3. Simplicidad y verificabilidad para esta primera implementación (N=2
   canales): la matriz de energía cinética de 3 puntos es el operador mejor
   caracterizado en métodos numéricos, permite un test de límite directo
   contra un oscilador armónico ANALÍTICO (§test 2a) y un test de
   convergencia en el paso de malla exactamente igual al de la Fase 1.

Si un test de convergencia futuro (Fase 3, diagrama de estabilización) muestra
que diferencias finitas necesita una malla prohibitivamente fina para la
precisión que se busca, se reconsiderará DVR — no antes, y con la misma
disciplina: se documenta, no se cambia en silencio.

SIMETRIZACIÓN
--------------
El operador de acoplamiento -A_ij d/dR - ½A_ij' es Hermítico en el continuo
(se verifica por integración por partes: ∫f(-A_ij g' - ½A_ij'g)dR =
∫g(-A_ji f' - ½A_ji'f)dR usando A_ji=-A_ij), pero la discretización directa
diag(A_ij) @ D1 no lo respeta exactamente — el mismo error O(h²) de
diferencias finitas centradas ya caracterizado en la Fase 1 (antisimetría de
A_ij, test 1a). La matriz completa se simetriza explícitamente,
H = (H + Hᵀ)/2, antes de diagonalizar: es la forma estándar de tratar este
artefacto de discretización sin ocultarlo (queda documentado aquí, no en un
ajuste de tolerancia posterior).
"""
from typing import Callable, Optional, Tuple

import numpy as np
from scipy.linalg import eigh

from trimero.systems.nonadiabatic_dynamics.coupling import (
    derivative_coupling,
    eigenbasis_along_R,
    fix_eigenvector_signs,
)

__all__ = [
    "second_derivative_coupling",
    "radial_kinetic_matrix",
    "first_derivative_matrix",
    "build_coupled_hamiltonian",
    "adiabatic_coupled_hamiltonian",
    "solve_coupled_channels",
]


# --------------------------------------------------------- acoplamiento B_ij
def second_derivative_coupling(
    R: np.ndarray,
    V: Optional[np.ndarray] = None,
    solve: Optional[Callable[[float], Tuple[np.ndarray, np.ndarray]]] = None,
    states: Optional[np.ndarray] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    ⟨Ψᵢ(R)|d²/dR²|Ψⱼ(R)⟩ por diferencias finitas centradas (3 puntos), en los
    puntos interiores de una malla uniforme. Réplica exacta de la interfaz de
    `coupling.derivative_coupling`, para el mismo tipo de entrada (V
    precalculado o `solve`), y el mismo tratamiento de continuidad de fase
    (`fix_eigenvector_signs`, secuencial).

    Returns:
        R_mid: R[1:-1].
        B: (n_R-2, n_states, n_states), B[k,i,j] = ⟨Ψᵢ(R_mid[k])|d²/dR²|Ψⱼ(R_mid[k])⟩.
    """
    R = np.asarray(R, dtype=float)
    if len(R) < 3:
        raise ValueError(
            "se necesitan al menos 3 puntos de R para diferencias centradas")

    steps = np.diff(R)
    if not np.allclose(steps, steps[0], rtol=1e-9, atol=1e-12):
        raise ValueError(
            "second_derivative_coupling requiere una malla de R uniforme "
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
    B = np.empty((n_mid, n_states, n_states))
    for k in range(n_mid):
        d2psi = (V[k + 2] - 2.0 * V[k + 1] + V[k]) / (h * h)
        B[k] = V[k + 1].T @ d2psi

    return R[1:-1], B


# ------------------------------------------------- operadores en la malla R
def radial_kinetic_matrix(n: int, h: float, mu: float, hbar: float = 1.0) -> np.ndarray:
    """T = -ħ²/(2μ) d²/dR², diferencias centradas de 3 puntos, condiciones de
    contorno de Dirichlet (χ=0 fuera de la malla, "partícula en una caja")."""
    T = np.zeros((n, n))
    diag = hbar**2 / (mu * h * h)
    off = -hbar**2 / (2.0 * mu * h * h)
    np.fill_diagonal(T, diag)
    idx = np.arange(n - 1)
    T[idx, idx + 1] = off
    T[idx + 1, idx] = off
    return T


def first_derivative_matrix(n: int, h: float) -> np.ndarray:
    """d/dR, diferencias centradas, antisimétrica exacta, Dirichlet en los
    bordes (el vecino fuera de la malla se toma como χ=0, así que en el
    primer/último punto el estencil sólo usa el vecino interior)."""
    D1 = np.zeros((n, n))
    idx = np.arange(n - 1)
    D1[idx, idx + 1] = 1.0 / (2.0 * h)
    D1[idx + 1, idx] = -1.0 / (2.0 * h)
    return D1


def build_coupled_hamiltonian(
    R_grid: np.ndarray,
    V: np.ndarray,
    A: Optional[np.ndarray] = None,
    B: Optional[np.ndarray] = None,
    mu: float = 1.0,
    hbar: float = 1.0,
) -> np.ndarray:
    """
    Matriz de canales acoplados de la ecuación del docstring del módulo,
    ensamblada en bloques (n_ch·n_R, n_ch·n_R), simetrizada.

    Args:
        R_grid: malla UNIFORME de R (los puntos donde viven V, A, B), len n_R.
        V: (n_R, n_ch), potenciales adiabáticos Vᵢ(R) (diagonal).
        A: (n_R, n_ch, n_ch) o None (= sin acoplamiento de primera derivada).
        B: (n_R, n_ch, n_ch) o None (= sin acoplamiento de segunda derivada,
            incluida la corrección diagonal de Born-Oppenheimer).
        mu: masa reducida nuclear, en m_e (unidades atómicas).
        hbar: por defecto 1.0 (unidades atómicas).

    Returns:
        H, (n_ch·n_R, n_ch·n_R), Hermítica (real simétrica).
    """
    R_grid = np.asarray(R_grid, dtype=float)
    n_R = len(R_grid)
    V = np.asarray(V, dtype=float)
    n_ch = V.shape[1]
    steps = np.diff(R_grid)
    if not np.allclose(steps, steps[0], rtol=1e-9, atol=1e-12):
        raise ValueError("build_coupled_hamiltonian requiere una malla uniforme")
    h = float(steps[0])

    T = radial_kinetic_matrix(n_R, h, mu, hbar)
    D1 = first_derivative_matrix(n_R, h) if A is not None else None

    H = np.zeros((n_ch * n_R, n_ch * n_R))
    for i in range(n_ch):
        for j in range(n_ch):
            block = np.zeros((n_R, n_R))
            if i == j:
                block += T
                block += np.diag(V[:, i])
            if B is not None:
                block += -(hbar**2 / (2.0 * mu)) * np.diag(B[:, i, j])
            if i != j and A is not None:
                block += -(hbar**2 / mu) * (np.diag(A[:, i, j]) @ D1)
            H[i * n_R:(i + 1) * n_R, j * n_R:(j + 1) * n_R] = block

    return 0.5 * (H + H.T)


# ------------------------------------------------- pipeline adiabático (Fase 1 -> 2)
def adiabatic_coupled_hamiltonian(
    R: np.ndarray,
    solve: Callable[[float], Tuple[np.ndarray, np.ndarray]],
    states: Optional[np.ndarray] = None,
    mu: float = 1.0,
    hbar: float = 1.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Conecta la Fase 1 (`solve(R)->(w,V)`) con `build_coupled_hamiltonian`:
    calcula Vᵢ(R) (autovalores), A_ij (`coupling.derivative_coupling`) y B_ij
    (`second_derivative_coupling`) sobre la misma malla de fase continua, y
    ensambla H.

    Returns:
        R_mid, H.
    """
    R_out, W, V = eigenbasis_along_R(R, solve)
    if states is not None:
        W = W[:, states]
        V = V[:, :, states]
    R_mid, A = derivative_coupling(R_out, V=V)
    _, B = second_derivative_coupling(R_out, V=V)
    V_mid = W[1:-1]
    H = build_coupled_hamiltonian(R_mid, V_mid, A=A, B=B, mu=mu, hbar=hbar)
    return R_mid, H


def solve_coupled_channels(
    R_grid: np.ndarray,
    V: np.ndarray,
    A: Optional[np.ndarray] = None,
    B: Optional[np.ndarray] = None,
    mu: float = 1.0,
    hbar: float = 1.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Ensambla (`build_coupled_hamiltonian`) y diagonaliza. Devuelve TODOS los
    autovalores/autovectores (el filtrado de qué es un estado ligado real
    frente a un estado de caja es la Fase 3, método de estabilización).

    Returns:
        E: (n_ch·n_R,) autovalores ordenados.
        chi: (n_ch·n_R, n_ch·n_R) autovectores (columnas); `chi[:, k]` se
            reshape a (n_ch, n_R) para leer χᵢ(R) del autoestado k.
    """
    H = build_coupled_hamiltonian(R_grid, V, A=A, B=B, mu=mu, hbar=hbar)
    E, chi = eigh(H)
    return E, chi
