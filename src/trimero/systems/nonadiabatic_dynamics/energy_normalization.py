"""
Normalización en energía de la función de onda del continuo, Fase 4 de
`docs/PLAN_nonadiabatic_dynamics.md`.

⚠️ RECONSTRUCCIÓN PROPIA, NO RÉPLICA DE LA REFERENCIA ORIGINAL. El método
que se necesita está descrito en González-Férez, Weidemüller & Schmelcher,
Phys. Rev. A 76, 023402 (2007) — la referencia que cita el paper de 2024
para esto exactamente ("A numerical way to obtain these energy normalized
wave functions is described in Ref. [2007]"). Ese texto NO se consiguió
(seis vías de acceso agotadas y documentadas en
`docs/analysis_fase4_busqueda_referencia.md`: DOI/APS 403, arXiv sin
preprint, página del grupo de Heidelberg de esa época caída, repositorio de
la UGR sin este artículo, Semantic Scholar limitado por rate-limiting).

Lo que sigue es una reconstrucción DESDE PRIMEROS PRINCIPIOS de una técnica
estándar de teoría de colisiones (normalización de estados del continuo
cuantizados en caja, ver p.ej. Sakurai & Napolitano, *Modern Quantum
Mechanics*, o Landau & Lifshitz, *Quantum Mechanics*, §21 sobre
normalización del continuo), verificada aquí contra un caso analítico
cerrado (`test_energy_normalization.py`). NO se ha contrastado contra el
método real del paper de 2007 — pendiente si se consigue acceso en el
futuro (ver `docs/analysis_fase4_normalizacion_y_decaimiento.md`).

DERIVACIÓN
-----------
Los autoestados de una caja finita (Dirichlet, `build_coupled_hamiltonian`
de la Fase 2) forman un conjunto ORTONORMAL DISCRETO: ⟨χₙ|χₘ⟩ = δₙₘ, y
satisfacen la resolución de la identidad

    Σₙ |χₙ⟩⟨χₙ| = 1̂                                              (1)

Contar estados: el número de autoestados con energía en [E, E+dE] es
dn = ρ(E) dE, con ρ(E) = dn/dE la densidad de estados — LOCAL, en esa
misma caja de tamaño L (no la deriva en L de la Fase 3; aquí L está fija y
se cuenta a lo largo de n, el índice del autoestado). Agrupando estados
vecinos en bins de energía dE (válido cuando dE es grande frente al
espaciado local ΔEₙ pero pequeño frente a la escala en la que varía la
física), la suma discreta de (1) se reescribe como integral:

    Σₙ |χₙ⟩⟨χₙ|  =  ∫ dE ρ(E) |χ_{n(E)}⟩⟨χ_{n(E)}|                (2)

Se quiere que esto coincida con la resolución de la identidad del continuo,
construida con estados normalizados en ENERGÍA (⟨χ_E|χ_E'⟩ = δ(E-E')):

    ∫ dE |χ_E⟩⟨χ_E| = 1̂                                          (3)

Comparando (2) y (3) término a término en dE:

    |χ_E⟩⟨χ_E| = ρ(E) |χ_{n(E)}⟩⟨χ_{n(E)}|
    χ_E = √ρ(E) · χ_{n(E)} = χ_{n(E)} / √(ΔEₙ)                     (4)

con ΔEₙ = 1/ρ(E) el espaciado LOCAL entre autovalores de caja vecinos a esa
energía, EN LA MISMA CAJA (estimado por diferencia central: ΔEₙ ≈
(E_{n+1}-E_{n-1})/2). Esta es la Ec. (4) que implementa
`energy_normalize_box_states`.

CONVENCIÓN DE NORMALIZACIÓN discreta ↔ continua
--------------------------------------------------
`numpy.linalg.eigh`/`scipy.linalg.eigh` devuelven autovectores v con
Σᵢvᵢ²=1 (norma discreta, sin peso de malla). La función de onda
CONTINUAMENTE normalizada (∫|χ|²dR=1) en los puntos de la malla es
χₙ(R_k) = v_k/√h (h = paso de malla): Σₖ(v_k/√h)² h = Σv_k² = 1. Combinando
con (4): χ_E(R_k) = v_k / √h / √(ΔEₙ).

CONVENCIÓN DE LA PARED DURA: R empieza en h, no en 0
-------------------------------------------------------
`radial_kinetic_matrix`/`build_coupled_hamiltonian` (Fase 2) imponen
Dirichlet en el punto INMEDIATAMENTE FUERA de la malla (χ=0 en el vecino
ficticio, no en el primer punto de la malla). Si la malla se arma con
`R = linspace(0, L, n)` (incluyendo R=0 como punto real), la pared física
queda en R=-h, no en R=0, y comparar contra una fórmula analítica escrita
para una pared en R=0 introduce un desfase de fase de orden h que NO
converge como el resto del código (se detectó exactamente este efecto en
esta misma ronda: con esa malla el error frente a la solución exacta caía
sólo como O(h), no O(h²), y desaparecía por completo al mover la malla a
`R = h, 2h, ..., nh` — pared física exactamente en R=0). Toda función de
este módulo asume esta segunda convención; se documenta aquí porque costó
una iteración real detectarlo (ver
docs/analysis_fase4_normalizacion_y_decaimiento.md).
"""
from typing import Optional

import numpy as np

__all__ = ["local_level_spacing", "energy_normalize_box_states"]


def local_level_spacing(E: np.ndarray) -> np.ndarray:
    """
    ΔEₙ ≈ (E_{n+1}-E_{n-1})/2, diferencia central en el ÍNDICE n de
    autovalores de UNA sola caja (a L fija) — no confundir con dE/dL de
    `stabilization.classify_stability`, que deriva en el tamaño de caja, no
    en el índice de nivel a caja fija.

    Args:
        E: autovalores ORDENADOS de una única diagonalización.

    Returns:
        mismo tamaño que E; extremos con diferencia hacia el vecino único.
    """
    E = np.asarray(E, dtype=float)
    if len(E) < 2:
        raise ValueError("se necesitan al menos 2 autovalores")
    dE = np.empty_like(E)
    dE[0] = E[1] - E[0]
    dE[-1] = E[-1] - E[-2]
    if len(E) > 2:
        dE[1:-1] = (E[2:] - E[:-2]) / 2.0
    return dE


def energy_normalize_box_states(
    R_grid: np.ndarray, chi_box: np.ndarray, E: np.ndarray, index: Optional[int] = None
) -> np.ndarray:
    """
    Convierte autovectores de caja (L²-normalizados, Σv²=1) a normalización
    en energía, Ec. (4) del docstring del módulo.

    Args:
        R_grid: malla UNIFORME de R (se verifica), empezando en R=h (ver
            docstring del módulo — la pared física debe caer en R=0).
        chi_box: (n_R, n_states), autovectores en columnas (convención
            estándar de eigh, Σv²=1 por columna).
        E: (n_states,), autovalores de ESA MISMA caja, mismo orden que las
            columnas de chi_box.
        index: si se da, devuelve sólo esa columna ya normalizada en
            energía (1D). Si no, devuelve todas (2D, misma forma que
            chi_box).

    Returns:
        χ_E en los puntos de R_grid (1D si `index` se dio, si no 2D).
    """
    R_grid = np.asarray(R_grid, dtype=float)
    steps = np.diff(R_grid)
    if not np.allclose(steps, steps[0], rtol=1e-9, atol=1e-12):
        raise ValueError(
            "energy_normalize_box_states requiere una malla de R uniforme "
            f"(pasos observados: {steps})")
    h = float(steps[0])

    chi_box = np.asarray(chi_box, dtype=float)
    dE = local_level_spacing(np.asarray(E, dtype=float))

    chi_R = chi_box / np.sqrt(h)
    chi_E = chi_R / np.sqrt(dE)[None, :]

    if index is not None:
        return chi_E[:, index]
    return chi_E
