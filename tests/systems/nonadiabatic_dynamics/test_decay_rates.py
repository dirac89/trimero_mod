"""
Tests de la tasa de decaimiento no adiabático (`decay_rates.py`), Fase 4 de
`docs/PLAN_nonadiabatic_dynamics.md`.

Modelo de juguete tipo Breit-Wigner, con referencia analítica INDEPENDIENTE
(cuadratura de funciones de onda cerradas, no el propio método de
diferencias finitas):

  - Canal "d" (discreto): pozo cuadrado RADIAL (pared dura en R=0, atractivo
    -V0 para R<a) — la misma familia que la Fase 3, restringida a R≥0 (como
    toda R física del proyecto). Único estado ligado con solución analítica
    CERRADA: ψ_d(R) = A sin(kR) dentro, A sin(ka)e^{-κ(R-a)} fuera, con A
    normalizado analíticamente (sin cuadratura). Su energía coincide con la
    rama impar de la ecuación trascendente ya usada y verificada en la Fase
    3 (`test_stabilization.py`): E_d = -0.05482888535439445 (V0=0.1, a=8.0,
    μ=1.0).
  - Canal "u" (continuo): partícula libre desplazada, V_u(R) = V_shift =
    -0.15 (constante, por debajo de E_d, para que el estado ligado de d
    quede embebido en el continuo de u — condición necesaria para que haya
    predisociación). Solución cerrada u_E(R) = √(2μ/(πk))sin(kR).
  - Acoplamiento: A_du(R) = A0 constante, B_du=0 (se deja B fuera para que
    la integral de solapamiento analítica sea tratable en forma cerrada).

Con esto, Γ_analítico = (2π)|-(1/μ)·A0·∫ψ_d(R)·dχ_u/dR dR|² se calcula por
cuadratura (`scipy.integrate.quad`) de las funciones cerradas — un camino
totalmente distinto de `nonadiabatic_decay_rate` (que ensambla y diagonaliza
Hamiltonianos de diferencias finitas). Concuerdan dentro de ~1%.

Advertencia sobre el orden de convergencia (documentada, no ocultada): a
diferencia de las Fases 1-3, Γ NO muestra una convergencia O(h²) limpia y
monótona al refinar h a L fija. La razón: además del error de discretización
espacial, hay un segundo efecto — el estado del continuo más cercano en
energía a E_d rara vez cae EXACTAMENTE en E_d (la caja sólo aproxima un
continuo con niveles espaciados ΔE_u), así que Γ se evalúa en un E_u
ligeramente distinto de E_d en cada malla, y Γ(E) varía con la energía. Este
"ruido de emparejamiento en energía" no decrece monótonamente con h a L
fija. Por eso el test verifica acuerdo cuantitativo (dentro de un margen
generoso) y la tendencia exacta Γ∝A0² (una identidad algebraica de la
fórmula, no sujeta a este efecto), en vez de exigir un orden de convergencia
formal — tal como pidió el usuario ("confirma orden de magnitud y
tendencia", no una prueba de convergencia).
"""
import numpy as np
import pytest
from scipy.integrate import quad

from trimero.systems.nonadiabatic_dynamics.coupled_channels import build_coupled_hamiltonian
from trimero.systems.nonadiabatic_dynamics.decay_rates import (
    coupling_matrix_element,
    coupling_operator_block,
    decay_rate,
    nonadiabatic_decay_rate,
)

MU, V0, A_WELL = 1.0, 0.1, 8.0
V_SHIFT = -0.15
E_D_EXACT = -0.05482888535439445  # rama impar del pozo cuadrado, Fase 3


def analytic_gamma(E_d, A0, mu=MU, well_depth=V0, well_width=A_WELL, threshold=V_SHIFT):
    """Γ por cuadratura de las funciones de onda cerradas (ver docstring)."""
    k_d = np.sqrt(2.0 * mu * (E_d + well_depth))
    kappa_d = np.sqrt(-2.0 * mu * E_d)
    norm2 = (well_width / 2 - np.sin(2 * k_d * well_width) / (4 * k_d)
             + np.sin(k_d * well_width) ** 2 / (2 * kappa_d))
    amp_d = 1.0 / np.sqrt(norm2)

    def psi_d(R):
        R = np.asarray(R, dtype=float)
        return np.where(R < well_width, amp_d * np.sin(k_d * R),
                        amp_d * np.sin(k_d * well_width) * np.exp(-kappa_d * (R - well_width)))

    k_u = np.sqrt(2.0 * mu * (E_d - threshold))

    def dchi_u_dR(R):
        return np.sqrt(2.0 * mu / (np.pi * k_u)) * k_u * np.cos(k_u * R)

    I1, _ = quad(lambda R: psi_d(R) * dchi_u_dR(R), 0.0, well_width, limit=200)
    I2, _ = quad(lambda R: psi_d(R) * dchi_u_dR(R), well_width, well_width + 300.0, limit=400)
    matrix_el = -(1.0 / mu) * A0 * (I1 + I2)
    return 2.0 * np.pi * abs(matrix_el) ** 2


def toy_grid(L, h):
    n = int(round(L / h))
    R = np.arange(1, n + 1) * h
    V_d = np.where(R < A_WELL, -V0, 0.0)
    V_u = np.full(n, V_SHIFT)
    return R, V_d, V_u


def find_bound_index(R, V_d):
    H_d = build_coupled_hamiltonian(R, V_d[:, None], mu=MU)
    E_d = np.linalg.eigvalsh(H_d)
    idx = int(np.argmin(np.abs(E_d - E_D_EXACT)))
    assert E_d[idx] < 0.0
    return idx, E_d[idx]


# --------------------------------------------------------- sanidad del toy
def test_toy_well_reproduces_known_bound_energy():
    """El pozo radial (pared dura en R=0) tiene el mismo único estado ligado
    ya verificado en la Fase 3 (rama impar): confirma que el modelo de
    juguete está bien planteado antes de calcular nada de Γ."""
    R, V_d, V_u = toy_grid(150.0, 0.05)
    idx, E_d = find_bound_index(R, V_d)
    np.testing.assert_allclose(E_d, E_D_EXACT, atol=3e-4)


# ---------------------------------------- núcleo: comparación con analítico
def test_decay_rate_matches_analytic_quadrature_reference():
    A0 = 0.01
    R, V_d, V_u = toy_grid(150.0, 0.1)
    bound_idx, E_d_numeric = find_bound_index(R, V_d)
    A_du = np.full(len(R), A0)
    B_du = np.zeros(len(R))

    gamma, diag = nonadiabatic_decay_rate(R, V_d, V_u, A_du, B_du, bound_idx, mu=MU)
    gamma_exact = analytic_gamma(E_D_EXACT, A0)

    print(f"\n  E_d numérico={diag['E_d']:.8e}  (exacto {E_D_EXACT:.8e})")
    print(f"  E_u (continuo más cercano)={diag['E_u']:.8e}  ΔE_u={diag['delta_E_u']:.4e}")
    print(f"  Gamma numérico  = {gamma:.6e}")
    print(f"  Gamma analítico = {gamma_exact:.6e}  (cuadratura, E_d exacto)")
    print(f"  diferencia relativa = {abs(gamma - gamma_exact) / gamma_exact:.4%}")

    # acuerdo observado ~0.1-1.3% en varias mallas; 2% deja margen sin ser
    # tan laxo que un error de factor/signo real pase desapercibido.
    assert abs(gamma - gamma_exact) / gamma_exact < 0.02


# --------------------------------------------------- tendencia: Gamma ~ A0^2
def test_decay_rate_scales_quadratically_with_coupling_strength():
    """Identidad algebraica de la Ec. 15 (M lineal en A0 con B_du=0 y A_du
    constante ⇒ Γ=|M|² cuadrático en A0), no sujeta al ruido de
    emparejamiento en energía del test anterior: debe cumplirse con muy
    poco margen."""
    R, V_d, V_u = toy_grid(150.0, 0.1)
    bound_idx, _ = find_bound_index(R, V_d)
    ratios = []
    for A0 in (0.005, 0.01, 0.02):
        A_du = np.full(len(R), A0)
        B_du = np.zeros(len(R))
        gamma, _ = nonadiabatic_decay_rate(R, V_d, V_u, A_du, B_du, bound_idx, mu=MU)
        ratios.append(gamma / A0**2)
    print(f"\n  Gamma/A0^2 para A0=(0.005,0.01,0.02): {ratios}")
    np.testing.assert_allclose(ratios, ratios[0], rtol=1e-9)


# --------------------------------------------------- sanidad: umbral, forma
def test_nonadiabatic_decay_rate_rejects_bound_state_below_continuum_threshold():
    """Si V_u no está por debajo de E_d, no hay predisociación posible:
    debe fallar explícitamente, no devolver un número sin sentido."""
    R, V_d, _ = toy_grid(150.0, 0.1)
    bound_idx, _ = find_bound_index(R, V_d)
    V_u_too_high = np.full(len(R), 1.0)  # muy por encima de E_d < 0
    A_du = np.zeros(len(R))
    B_du = np.zeros(len(R))
    with pytest.raises(ValueError):
        nonadiabatic_decay_rate(R, V_d, V_u_too_high, A_du, B_du, bound_idx, mu=MU)


def test_coupling_operator_block_is_finite_and_correct_shape():
    R, V_d, V_u = toy_grid(50.0, 0.2)
    A_du = np.full(len(R), 0.01)
    B_du = np.zeros(len(R))
    Op = coupling_operator_block(R, A_du, B_du, mu=MU)
    assert Op.shape == (len(R), len(R))
    assert np.all(np.isfinite(Op))
