"""
Tasas de decaimiento no adiabático (predisociación), Fase 4 de
`docs/PLAN_nonadiabatic_dynamics.md`.

Réplica de la Ec. 15 de Mellado-Alcedo, Guttridge, Cornish, Sadeghpour &
González-Férez, PRA 110, 013314 (2024) — ecuación que SÍ se pudo leer
completa (arxiv.org/html/2401.09618):

    Γᵢ = (2π/ħ) |⟨χᵢ^d|A_du|χⱼ^u⟩|²                                (15)

donde χᵢ^d es un estado LIGADO del canal "d" (discreto, L²-normalizado a 1)
y χⱼ^u es un estado del CONTINUO del canal "u", normalizado en ENERGÍA
(`energy_normalization.py` — reconstrucción propia, ver advertencia en ese
módulo). A_du es el operador de acoplamiento de la Ec. 8 (Fase 2,
`coupled_channels.py`):

    A_du = -ħ²/μ · A_du(R) d/dR - ħ²/(2μ) · B_du(R)

Se trata Γᵢ como REGLA DE ORO DE FERMI de orden más bajo: χᵢ^d y χⱼ^u se
resuelven cada uno en su propio canal DESACOPLADO (orden cero, A=B=0), y
A_du hace de perturbación entre ellos — es exactamente la estructura que
implica citar la regla de oro (density-of-states × |elemento de matriz|²) y
la que usa el paper de 2024 para "ilustrar el acoplamiento inducido en los
cruces evitados" (texto del paper, Sec. IV).

CONSISTENCIA DE REPRESENTACIÓN (por qué no hace falta pesar por h)
----------------------------------------------------------------------
`build_coupled_hamiltonian` actúa sobre autovectores v en la convención
Σv²=1 (no en la función de onda continuamente normalizada χ=v/√h). Esa
misma matriz, restringida al bloque fuera de la diagonal (i≠j), ES el
operador A_du discretizado, en la MISMA convención. Para el elemento de
matriz físico ∫χᵢ^d(R)[A_du χⱼ^u](R)dR, expresando todo en la convención
v (factor √h en cada función de onda) los dos factores 1/√h se cancelan
contra el peso h de la cuadratura de Riemann, y queda

    ⟨χᵢ^d|A_du|χⱼ^u⟩ = v_d · (Op_du · v_u) / √(ΔEⱼ)

con v_d, v_u los autovectores PUROS (Σv²=1) de los dos canales
desacoplados, Op_du el bloque fuera de la diagonal de
`build_coupled_hamiltonian` para un sistema de 2 canales, y ΔEⱼ el
espaciado local de `energy_normalization.local_level_spacing` del estado
continuo — la única normalización "extra" que hace falta aplicar aquí es
justo el factor 1/√(ΔEⱼ) de la Ec. (4) de `energy_normalization.py`; el
resto ya lo captura la propia matriz.
"""
from typing import Optional, Tuple

import numpy as np

from trimero.systems.nonadiabatic_dynamics.coupled_channels import build_coupled_hamiltonian
from trimero.systems.nonadiabatic_dynamics.energy_normalization import local_level_spacing

__all__ = ["coupling_operator_block", "coupling_matrix_element", "decay_rate",
           "nonadiabatic_decay_rate"]


def coupling_operator_block(
    R_grid: np.ndarray, A_du: np.ndarray, B_du: np.ndarray, mu: float = 1.0, hbar: float = 1.0
) -> np.ndarray:
    """
    Bloque fuera de la diagonal (canal d -> canal u) de
    `build_coupled_hamiltonian` para un sistema de 2 canales, aislado. Se
    obtiene ensamblando un H de 2 canales con V=0 (el bloque fuera de la
    diagonal no depende de V) y extrayendo el bloque (0,1).
    """
    n_R = len(R_grid)
    V = np.zeros((n_R, 2))
    A = np.zeros((n_R, 2, 2))
    A[:, 0, 1], A[:, 1, 0] = A_du, -A_du
    B = np.zeros((n_R, 2, 2))
    B[:, 0, 1], B[:, 1, 0] = B_du, B_du
    H = build_coupled_hamiltonian(R_grid, V, A=A, B=B, mu=mu, hbar=hbar)
    return H[:n_R, n_R:]


def coupling_matrix_element(
    v_d: np.ndarray, v_u: np.ndarray, Op_du: np.ndarray, delta_E_u: float
) -> float:
    """⟨χᵢ^d|A_du|χⱼ^u⟩, ver docstring del módulo. v_d, v_u: autovectores
    PUROS (Σv²=1) de los canales d, u DESACOPLADOS. delta_E_u: espaciado
    local del estado continuo (`energy_normalization.local_level_spacing`)."""
    raw = v_d @ (Op_du @ v_u)
    return raw / np.sqrt(delta_E_u)


def decay_rate(matrix_element: float, hbar: float = 1.0) -> float:
    """Γ = (2π/ħ)|elemento de matriz|², Ec. 15."""
    return (2.0 * np.pi / hbar) * abs(matrix_element) ** 2


def nonadiabatic_decay_rate(
    R_grid: np.ndarray,
    V_d: np.ndarray,
    V_u: np.ndarray,
    A_du: np.ndarray,
    B_du: np.ndarray,
    bound_index_d: int,
    mu: float = 1.0,
    hbar: float = 1.0,
) -> Tuple[float, dict]:
    """
    Ensambla y diagonaliza los DOS canales DESACOPLADOS por separado
    (orden cero de la regla de oro), localiza el estado del continuo de u
    más cercano en energía al estado ligado `bound_index_d` de d, y calcula
    Γ vía la Ec. 15.

    Returns:
        Γᵢ, y un dict de diagnóstico {E_d, E_u, continuum_index,
        delta_E_u, matrix_element} para inspección/tests.
    """
    H_d = build_coupled_hamiltonian(R_grid, V_d[:, None], mu=mu, hbar=hbar)
    H_u = build_coupled_hamiltonian(R_grid, V_u[:, None], mu=mu, hbar=hbar)
    E_d, chi_d = np.linalg.eigh(H_d)
    E_u, chi_u = np.linalg.eigh(H_u)

    E_bound = E_d[bound_index_d]
    if E_bound <= E_u[0]:
        raise ValueError(
            f"E_bound={E_bound:.6e} está por debajo del umbral del canal u "
            f"({E_u[0]:.6e}): no hay predisociación posible, revisa V_u.")
    continuum_index = int(np.argmin(np.abs(E_u - E_bound)))

    Op_du = coupling_operator_block(R_grid, A_du, B_du, mu=mu, hbar=hbar)
    dE_u = local_level_spacing(E_u)[continuum_index]

    v_d = chi_d[:, bound_index_d]
    v_u = chi_u[:, continuum_index]
    me = coupling_matrix_element(v_d, v_u, Op_du, dE_u)
    gamma = decay_rate(me, hbar=hbar)

    return gamma, {
        "E_d": float(E_bound), "E_u": float(E_u[continuum_index]),
        "continuum_index": continuum_index, "delta_E_u": float(dE_u),
        "matrix_element": float(me),
    }
