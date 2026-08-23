"""
Tests de la normalización en energía del continuo
(`trimero.systems.nonadiabatic_dynamics.energy_normalization`), Fase 4 de
`docs/PLAN_nonadiabatic_dynamics.md`.

⚠️ Este módulo es una RECONSTRUCCIÓN PROPIA (ver docstring de
`energy_normalization.py`), no una réplica del método de González-Férez,
Weidemüller & Schmelcher, PRA 76, 023402 (2007) — ese texto no se consiguió
(`docs/analysis_fase4_busqueda_referencia.md`). Se verifica aquí contra un
caso analítico CERRADO, independiente del propio método: la partícula libre
1D radial (pared dura en R=0), cuya solución normalizada en energía es
conocida en forma cerrada,

    u_E(R) = √(2μ/(πħ²k)) · sin(kR),    k = √(2μE)

(normalización estándar de la onda parcial libre a δ(E-E'), p.ej. Sakurai &
Napolitano §A.5, o Landau & Lifshitz §21 y su normalización de la onda
esférica libre l=0).

NOTA DE ESTA RONDA (dejar registrado — costó una iteración real, ver
docs/analysis_fase4_normalizacion_y_decaimiento.md §2): la primera versión
de este test usaba `R = np.linspace(0, L, n)` (incluyendo R=0 como punto de
malla) y comparaba directamente contra u_E(R) — el error caía sólo como
O(h), no O(h²) como en el resto del proyecto. La causa: con esa malla la
pared de Dirichlet real de `build_coupled_hamiltonian` (que impone χ=0 en el
vecino FUERA del array, no en el primer punto del array) cae en R=-h, no en
R=0, así que la fase de la solución numérica está desplazada un paso de
malla respecto a u_E(R). Corregido armando la malla como
`R = h, 2h, ..., n·h` (pared física exactamente en R=0, sin ningún punto en
R=0 dentro del array) — la misma convención que ya usan
`radial_kinetic_matrix`/`build_coupled_hamiltonian`, aplicada aquí de forma
explícita y consistente.
"""
import numpy as np
import pytest
from scipy.linalg import eigh_tridiagonal

from trimero.systems.nonadiabatic_dynamics.coupled_channels import build_coupled_hamiltonian
from trimero.systems.nonadiabatic_dynamics.energy_normalization import (
    energy_normalize_box_states,
    local_level_spacing,
)


def free_particle_box(mu, L, h):
    """Caja radial libre (V=0), pared dura en R=0 (malla R=h,2h,...,nh) y en
    R=L+h. Usa `eigh_tridiagonal` (exacto para esta matriz, mucho más rápido
    que la ruta densa de `build_coupled_hamiltonian` para las cajas grandes
    que hacen falta aquí para que R<<L)."""
    n = int(round(L / h))
    R = np.arange(1, n + 1) * h
    diag = np.full(n, 1.0 / (mu * h * h))
    off = np.full(n - 1, -1.0 / (2.0 * mu * h * h))
    E, chi = eigh_tridiagonal(diag, off)
    return R, E, chi


def analytic_energy_normalized_free_particle(R, mu, k):
    """u_E(R) = sqrt(2 mu/(pi k)) sin(kR), hbar=1."""
    return np.sqrt(2.0 * mu / (np.pi * k)) * np.sin(k * R)


# ---------------------------------------------------- local_level_spacing
def test_local_level_spacing_matches_infinite_box_formula():
    """Para la caja libre infinita, ΔEₙ analítico = π k/(μL) (derivado de
    Eₙ=n²π²/(2μL²)). Referencia independiente del propio código: fórmula
    cerrada, no la propia diferencia central aplicada a otra cosa."""
    mu, L, h = 1.0, 2000.0, 0.2
    R, E, chi = free_particle_box(mu, L, h)
    idx = np.argmin(np.abs(E - 0.01))
    dE = local_level_spacing(E)
    k = np.sqrt(2.0 * mu * E[idx])
    dE_exact = np.pi * k / (mu * L)
    np.testing.assert_allclose(dE[idx], dE_exact, rtol=2e-3)


def test_local_level_spacing_rejects_too_short_array():
    with pytest.raises(ValueError):
        local_level_spacing(np.array([1.0]))


# --------------------------------- verificación central: contra u_E(R) exacto
def test_energy_normalization_matches_analytic_free_particle():
    """
    Núcleo de la Fase 4 (normalización): en una caja libre grande, la
    reconstrucción χ_E de `energy_normalize_box_states` debe coincidir con
    u_E(R) = √(2μ/(πk))sin(kR) en la región R≪L, con el error cayendo como
    O(h²) — igual que todos los demás operadores de diferencias finitas del
    proyecto (Fases 1-3).
    """
    mu, L, E0 = 1.0, 600.0, 0.01
    errors = []
    for h in (0.2, 0.1, 0.05):
        R, E, chi = free_particle_box(mu, L, h)
        idx = np.argmin(np.abs(E - E0))
        k = np.sqrt(2.0 * mu * E[idx])

        chi_E = energy_normalize_box_states(R, chi, E, index=idx)
        u_E = analytic_energy_normalized_free_particle(R, mu, k)

        # fase global arbitraria del autovector: fijar signo en un punto
        i0 = np.argmax(R > 5.0 / k)
        if chi_E[i0] * u_E[i0] < 0:
            chi_E = -chi_E

        mask = (R > 3.0) & (R < 150.0)  # R << L=600
        err = np.max(np.abs(chi_E[mask] - u_E[mask]))
        errors.append(err)

    print(f"\n  errores vs h=(0.2,0.1,0.05): {errors}")
    assert errors[0] > errors[1] > errors[2]
    # convergencia O(h^2) limpia: medido errors[0]/errors[1]=4.00,
    # errors[1]/errors[2]=4.00 (ver docs/analysis_fase4_normalizacion_y_decaimiento.md)
    assert errors[0] / errors[1] > 3.5
    assert errors[1] / errors[2] > 3.5
    # margen sobre el residuo medido con h=0.05 (8.40e-5), no un valor
    # adivinado antes de ejecutar
    assert errors[2] < 2e-4


# ------------------------------------------------------- validación de malla
def test_energy_normalize_rejects_nonuniform_grid():
    R = np.array([0.1, 0.2, 0.35])
    chi = np.zeros((3, 2))
    E = np.array([0.1, 0.2])
    with pytest.raises(ValueError):
        energy_normalize_box_states(R, chi, E)
