"""
Método de estabilización (Hazi & Taylor, Phys. Rev. A 1, 1109, 1970), Fase 3
de `docs/PLAN_nonadiabatic_dynamics.md`.

Módulo GENÉRICO (no conoce ningún sistema físico concreto ni la ecuación de
canales acoplados de la Fase 2 directamente): recibe un
`hamiltonian_builder(L) -> (R_grid, H)` que ensambla un Hamiltoniano en una
caja de tamaño L (típicamente `coupled_channels.build_coupled_hamiltonian`
con una malla que se extiende con L), y encuentra qué autovalores son
estados ligados REALES frente a estados de caja.

LA IDEA FÍSICA
---------------
Un estado ligado real vive dentro del pozo y decae exponencialmente fuera de
él: una vez que la caja es más grande que el alcance de esa cola, su energía
prácticamente no cambia al seguir agrandando la caja. Un "estado de caja"
(la discretización del continuo por el confinamiento artificial) sí depende
fuertemente de L — para una caja infinita 1D de longitud L, el n-ésimo nivel
escala como Eₙ(L) ≈ n²π²ħ²/(2μL²), así que dEₙ/dL = -2Eₙ/L. Según L crece,
los estados de caja "barren" hacia abajo en energía y cruzan (en realidad,
por la regla de no-cruce, se acoplan levemente y EVITAN cruzar) la
trayectoria plana del estado ligado real — el "diagrama de estabilización"
clásico: una línea horizontal atravesada por un enjambre de líneas que caen.

CRITERIO AUTOMÁTICO (no visual)
---------------------------------
En vez de identificar la meseta a ojo, se compara |dEₙ/dL| medido en cada
trayectoria contra la escala de un estado de caja EN ESA MISMA ENERGÍA Y
LONGITUD, 2|E|/L (la fórmula de arriba, despejada). Un punto de la
trayectoria se clasifica ESTABLE si

    |dE/dL| < umbral · (2|E|/L)

con `umbral` (por defecto 0.1) fijando cuánto más lenta que un estado de
caja típico debe ser la variación para contar como "plana". Una trayectoria
completa se acepta como estado ligado si es estable en al menos una fracción
`min_fraction` de los L escaneados Y en el L más grande (el punto donde más
se puede confiar en que la cola exponencial ya cabe dentro de la caja).

EMPAREJAMIENTO DE TRAYECTORIAS ENTRE CAJAS VECINAS
-----------------------------------------------------
A diferencia de la Fase 1 (donde el parámetro es R y se sigue la fase de los
AUTOVECTORES), aquí la dimensión del Hamiltoniano cambia con L (la malla
crece), así que no hay autovectores de dimensión constante que comparar.
En su lugar se sigue la trayectoria por PROXIMIDAD EN ENERGÍA entre cajas de
tamaño vecino (asignación óptima, algoritmo húngaro,
`scipy.optimize.linear_sum_assignment`): válido porque, para un Hamiltoniano
real que depende de un único parámetro continuo, autovalores de la misma
simetría no se cruzan exactamente (regla de no-cruce) — con un paso de L
suficientemente fino, la trayectoria más cercana en energía es la física
correcta. Es el mismo principio de no-cruce que ya usa
`trimero.simulation.bop_tracking.trace_curve` para seguir curvas en R, sólo
que aquí el emparejamiento es por energía en vez de por solapamiento de
autovector (porque la dimensión cambia).
"""
from typing import Callable, List, Tuple

import numpy as np
from scipy.optimize import linear_sum_assignment

__all__ = [
    "stabilization_scan",
    "track_trajectories",
    "classify_stability",
    "stable_state_summary",
]


def stabilization_scan(
    L_values: np.ndarray,
    hamiltonian_builder: Callable[[float], Tuple[np.ndarray, np.ndarray]],
    n_keep: int = 20,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Diagonaliza `hamiltonian_builder(L)` para cada L y guarda los `n_keep`
    autovalores más bajos.

    Returns:
        L_values, E: (n_L, n_keep), autovalores ordenados en cada fila (SIN
        emparejar todavía entre L vecinos — eso es `track_trajectories`).
    """
    L_values = np.asarray(L_values, dtype=float)
    E = np.empty((len(L_values), n_keep))
    for i, L in enumerate(L_values):
        _, H = hamiltonian_builder(float(L))
        w = np.linalg.eigvalsh(H)
        if len(w) < n_keep:
            raise ValueError(
                f"n_keep={n_keep} pero el Hamiltoniano en L={L} sólo tiene "
                f"{len(w)} autovalores")
        E[i] = w[:n_keep]
    return L_values, E


def track_trajectories(E: np.ndarray) -> np.ndarray:
    """
    Reordena cada fila de E (autovalores por caja, sin correspondencia entre
    filas) en trayectorias continuas por L, emparejando cajas vecinas por
    proximidad en energía (asignación óptima, algoritmo húngaro). Ver
    docstring del módulo para la justificación (regla de no-cruce).

    Returns:
        trajectories: misma forma que E, columna k = una trayectoria continua.
    """
    n_L, n_keep = E.shape
    trajectories = np.empty_like(E)
    trajectories[0] = np.sort(E[0])
    for i in range(1, n_L):
        cost = np.abs(trajectories[i - 1][:, None] - E[i][None, :])
        row, col = linear_sum_assignment(cost)
        trajectories[i, row] = E[i, col]
    return trajectories


def classify_stability(
    L_values: np.ndarray, trajectories: np.ndarray, threshold_ratio: float = 0.1
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    dE/dL por diferencias centradas en los puntos interiores de L_values
    (malla UNIFORME, se verifica). Clasifica ESTABLE cada punto donde
    |dE/dL| < threshold_ratio · (2|E|/L) (ver docstring del módulo).

    Returns:
        L_mid: L_values[1:-1].
        stable: (n_L-2, n_keep) bool.
        slope: (n_L-2, n_keep), dE/dL medido (para inspección/diagnóstico).
    """
    L_values = np.asarray(L_values, dtype=float)
    if len(L_values) < 3:
        raise ValueError("se necesitan al menos 3 valores de L")
    steps = np.diff(L_values)
    if not np.allclose(steps, steps[0], rtol=1e-9, atol=1e-12):
        raise ValueError(
            "classify_stability requiere una malla de L uniforme "
            f"(pasos observados: {steps})")
    h = float(steps[0])

    slope = (trajectories[2:] - trajectories[:-2]) / (2.0 * h)
    E_mid = trajectories[1:-1]
    L_mid = L_values[1:-1]
    box_scale = 2.0 * np.abs(E_mid) / L_mid[:, None]
    stable = np.abs(slope) < threshold_ratio * box_scale
    return L_mid, stable, slope


def stable_state_summary(
    L_mid: np.ndarray,
    trajectories_mid: np.ndarray,
    stable: np.ndarray,
    min_fraction: float = 0.5,
) -> List[dict]:
    """
    Resumen por trayectoria: se acepta como estado ligado real si es estable
    en al menos `min_fraction` de los L escaneados Y en el L más grande (el
    punto donde más se puede confiar en que la caja ya contiene la cola
    exponencial del estado real).

    Returns:
        lista de dicts {index, fraction_stable, is_bound, E_mean}, E_mean es
        el promedio de la trayectoria SÓLO en los tramos estables (nan si
        nunca es estable).
    """
    n_L, n_keep = stable.shape
    summaries = []
    for k in range(n_keep):
        frac = float(stable[:, k].mean())
        is_bound = bool(frac >= min_fraction and stable[-1, k])
        if stable[:, k].any():
            E_mean = float(trajectories_mid[stable[:, k], k].mean())
        else:
            E_mean = float("nan")
        summaries.append({
            "index": k, "fraction_stable": frac, "is_bound": is_bound,
            "E_mean": E_mean,
        })
    return summaries
