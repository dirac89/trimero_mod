"""
Seguimiento de una curva de potencial Born-Oppenheimer (BOP) a través de R.

El problema que resuelve: elegir el estado en cada R por un criterio local
("el más bajo con carácter de manifold") hace que la etiqueta salte entre
estados en regiones de casi-degeneración, produciendo curvas discontinuas que
no son comparables con una figura publicada.

Aquí se sigue UNA curva por SOLAPAMIENTO máximo del autovector con el paso
anterior, |⟨v_R | v_{R-h}⟩|², con refinamiento adaptativo del paso cuando el
solapamiento cae por debajo de un umbral.

NOTA FÍSICA: dentro de un mismo bloque M_J las curvas adiabáticas NO se cruzan
(regla de no cruce), así que el índice del autovalor ordenado ya es de por sí
una etiqueta adiabática global. La excepción son los cruces genuinos entre las
dos especies de reflexión σ_v que conviven en M_J=0, donde sí puede haber
degeneración exacta. Por eso `trace_curve` devuelve también el índice ordenado
en cada paso: si se mantiene constante, ambos criterios coinciden y la curva es
la adiabática; si cambia, el seguimiento por solapamiento ha pasado por un
cruce real y hay que mirarlo.
"""

from typing import Callable, Dict, List

import numpy as np

__all__ = ["trace_curve"]


def trace_curve(
    solve: Callable[[float], tuple],
    R_start: float,
    R_end: float,
    step: float,
    select_initial: Callable[[np.ndarray, np.ndarray], int],
    overlap_threshold: float = 0.7,
    max_refine: int = 6,
    n_context: int = 0,
    progress: Callable[[str], None] = None,
) -> Dict:
    """
    Args:
        solve: R -> (w, V) autovalores ordenados y autovectores (columnas).
        step: paso base en R.
        select_initial: (w, V) -> índice del estado del que se parte.
        overlap_threshold: si |⟨v|v_prev⟩|² cae por debajo, se parte el paso
            por la mitad y se reintenta, en vez de aceptar el salto.
        max_refine: número máximo de bisecciones consecutivas.
        n_context: cuántos autovalores más bajos guardar en cada R para dibujar
            el fondo de curvas vecinas (0 = ninguno).

    Returns:
        dict con R, E, index, overlap, weight_slot, context, refinements.
    """
    cache: Dict[float, tuple] = {}

    def solve_cached(R):
        key = round(R, 9)
        if key not in cache:
            cache[key] = solve(R)
        return cache[key]

    w, V = solve_cached(R_start)
    k0 = select_initial(w, V)
    v_prev = V[:, k0]

    R_list = [R_start]
    E_list = [float(w[k0])]
    idx_list = [int(k0)]
    ov_list = [1.0]
    ctx_list = [np.array(w[:n_context]) if n_context else None]
    refinements: List[Dict] = []

    R_cur = R_start
    while R_cur < R_end - 1e-9:
        h = step
        depth = 0
        while True:
            target = min(R_cur + h, R_end)
            w, V = solve_cached(target)
            ov = (v_prev @ V) ** 2
            idx = int(np.argmax(ov))
            best = float(ov[idx])
            if best >= overlap_threshold or depth >= max_refine:
                break
            depth += 1
            refinements.append(
                {"R_from": R_cur, "R_try": target, "overlap": best,
                 "depth": depth, "new_step": h / 2.0}
            )
            if progress:
                progress(f"    refinando en R={R_cur:.2f}->{target:.2f}: "
                         f"solapamiento {best:.3f}, paso {h:.3f}->{h/2:.3f}")
            h /= 2.0

        v_prev = V[:, idx]
        R_cur = target
        R_list.append(target)
        E_list.append(float(w[idx]))
        idx_list.append(idx)
        ov_list.append(best)
        ctx_list.append(np.array(w[:n_context]) if n_context else None)

    return {
        "R": np.array(R_list),
        "E": np.array(E_list),
        "index": np.array(idx_list),
        "overlap": np.array(ov_list),
        "context": np.array(ctx_list) if n_context else None,
        "refinements": refinements,
        "n_solves": len(cache),
    }
