"""
Tests de los factores de Franck-Condon (`franck_condon.py`), Fase 5 de
`docs/PLAN_nonadiabatic_dynamics.md`.

Referencia analítica CERRADA: el solapamiento entre dos gaussianas
normalizadas es un resultado estándar de libro de texto,

    ⟨φ_a|φ_b⟩ = √(2σₐσᵦ/(σₐ²+σᵦ²)) · exp[-(Rₐ-Rᵦ)²/(2(σₐ²+σᵦ²))]

Se usa el estado fundamental de un oscilador armónico (ya verificado contra
la energía analítica en la Fase 2) como estado FINAL χᵢᵏ — es exactamente
gaussiano, con anchura σ_HO=1/√(μω) — y se solapa contra el paquete de
ondas gaussiano `gaussian_wavepacket` (el modelo de estado inicial de esta
fase) centrado en un punto distinto. Dos caminos numéricos independientes
del solapamiento (cuadratura discreta de `franck_condon_factor` vs. fórmula
cerrada de dos gaussianas) deben coincidir, convergiendo O(h²) con el paso
de malla — igual que el resto del proyecto.
"""
import numpy as np
import pytest

from trimero.systems.nonadiabatic_dynamics.coupled_channels import build_coupled_hamiltonian
from trimero.systems.nonadiabatic_dynamics.franck_condon import (
    franck_condon_factor,
    gaussian_wavepacket,
)

MU, OMEGA = 1.0, 0.05
R_CENTER = 100.0  # centro del oscilador armónico, lejos de la pared R=0
SIGMA_HO = 1.0 / np.sqrt(MU * OMEGA)


def two_gaussian_overlap(Ra, sa, Rb, sb):
    return (np.sqrt(2 * sa * sb / (sa**2 + sb**2))
            * np.exp(-(Ra - Rb) ** 2 / (2 * (sa**2 + sb**2))))


def ho_ground_state(L, h):
    n = int(round(L / h))
    R = np.arange(1, n + 1) * h
    V = 0.5 * MU * OMEGA**2 * (R - R_CENTER) ** 2
    H = build_coupled_hamiltonian(R, V[:, None], mu=MU)
    E, chi = np.linalg.eigh(H)
    v0 = chi[:, 0]
    i0 = np.argmin(np.abs(R - R_CENTER))
    if v0[i0] < 0:
        v0 = -v0
    return R, v0, E[0]


# --------------------------------------------------------- gaussian_wavepacket
def test_gaussian_wavepacket_is_normalized():
    R = np.arange(1, 4000) * 0.05
    u = gaussian_wavepacket(R, R0=110.0, sigma=6.0)
    h = R[1] - R[0]
    np.testing.assert_allclose(np.sum(u**2) * h, 1.0, rtol=1e-8)


def test_gaussian_wavepacket_rejects_nonuniform_grid():
    R = np.array([1.0, 2.0, 3.5])
    with pytest.raises(ValueError):
        gaussian_wavepacket(R, R0=2.0, sigma=1.0)


# ----------------------------------------------- solapamiento vs. cerrado
def test_franck_condon_matches_two_gaussian_overlap_closed_form():
    R0, sigma_scat = 110.0, 6.0
    F_exact = two_gaussian_overlap(R_CENTER, SIGMA_HO, R0, sigma_scat)

    errors = []
    for h in (0.2, 0.1, 0.05):
        R, v0, E0 = ho_ground_state(250.0, h)
        u_scat = gaussian_wavepacket(R, R0, sigma_scat)
        F_num = franck_condon_factor(R, v0, u_scat)
        errors.append(abs(F_num - F_exact))

    print(f"\n  F_exacto = {F_exact:.8e}")
    print(f"  errores vs h=(0.2,0.1,0.05): {errors}")
    assert errors[0] > errors[1] > errors[2]
    assert errors[0] / errors[1] > 3.5
    assert errors[1] / errors[2] > 3.5
    assert errors[2] < 1e-4


def test_franck_condon_coupling_weight_scales_result():
    """Un peso electrónico C(R)=c constante debe escalar F por c exactamente
    (linealidad de la integral) — no depende del solapamiento en sí."""
    R, v0, _ = ho_ground_state(250.0, 0.1)
    u_scat = gaussian_wavepacket(R, R0=110.0, sigma=6.0)
    F_unweighted = franck_condon_factor(R, v0, u_scat)
    c = 3.7
    F_weighted = franck_condon_factor(R, v0, u_scat, coupling_weight=np.full(len(R), c))
    np.testing.assert_allclose(F_weighted, c * F_unweighted, rtol=1e-12)


def test_franck_condon_rejects_nonuniform_grid():
    R = np.array([1.0, 2.0, 3.5])
    v = np.zeros(3)
    u = np.zeros(3)
    with pytest.raises(ValueError):
        franck_condon_factor(R, v, u)


# -------------------------------------------- decae a cero lejos del solapamiento
def test_franck_condon_vanishes_for_far_separated_states():
    """Si el paquete inicial está muchas anchuras lejos del estado final, F
    debe ser despreciable (cola gaussiana) — sanity check de orden de
    magnitud, no de precisión."""
    R, v0, _ = ho_ground_state(400.0, 0.1)
    u_scat_near = gaussian_wavepacket(R, R0=110.0, sigma=6.0)
    u_scat_far = gaussian_wavepacket(R, R0=300.0, sigma=6.0)  # ~45 sigma_HO lejos
    F_near = abs(franck_condon_factor(R, v0, u_scat_near))
    F_far = abs(franck_condon_factor(R, v0, u_scat_far))
    print(f"\n  F_near={F_near:.4e}  F_far={F_far:.4e}")
    assert F_far < 1e-10 * F_near
