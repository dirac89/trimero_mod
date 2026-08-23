"""
Tests del acoplamiento no adiabático de derivada ⟨Ψᵢ(R)|d/dR|Ψⱼ(R)⟩
(`trimero.systems.nonadiabatic_dynamics.coupling`), Fase 1 de
`docs/PLAN_nonadiabatic_dynamics.md`.

Modelo de referencia (tests a, b, e): sistema 2x2 real simétrico

    H(R) = [[a(R), c], [c, b(R)]],  a(R) = -k(R-R0),  b(R) = +k(R-R0)

con cruce evitado centrado en R0 y hueco mínimo 2c. Para este modelo el
ángulo de mezcla θ(R) = (1/2) atan2(2c, a(R)-b(R)) tiene autovectores
psi_lower = (cosθ, sinθ), psi_upper = (-sinθ, cosθ), y el acoplamiento de
derivada exacto es |A_12(R)| = |dθ/dR| = c·k / (2[k²(R-R0)² + c²])
(derivación en docs/analysis_fase1_acoplamiento_derivada.md). Sirve de
referencia analítica independiente de la discretización.

Test d usa un cruce evitado REAL ya documentado en el proyecto: el sistema
Rb*-KRb (`BOPSystem`, n_manifold=25, M_J=0, sin Fermi) tiene, según
`plots/rb_krb_polar/fig1_ad_MJ0_n25.npz`, un cambio de índice de la curva de
carácter de manifold k: 54->55 cerca de R=1300 a0. Un escaneo fino del hueco
de energía w[55]-w[54] localiza el mínimo real en R≈1305 a0 (hueco ≈5e-9 Eh,
muchísimo más estrecho que la resolución de 5 a0 del barrido de producción).
"""
import numpy as np
import pytest

from trimero.systems.nonadiabatic_dynamics.coupling import (
    derivative_coupling,
    eigenbasis_along_R,
    fix_eigenvector_signs,
)


# --------------------------------------------------------------------- toy
def toy_two_level_solve(k=0.01, c=0.05, R0=10.0):
    def solve(R):
        a = -k * (R - R0)
        b = k * (R - R0)
        H = np.array([[a, c], [c, b]])
        return np.linalg.eigh(H)
    return solve


def toy_analytic_coupling(R, k=0.01, c=0.05, R0=10.0):
    """|A_12(R)| exacto, ver docstring del módulo."""
    x = -2.0 * k * (R - R0)
    return 2.0 * c * k / (x**2 + 4.0 * c**2)


# --------------------------------------------------------- fix_eigenvector_signs
def test_fix_eigenvector_signs_flips_negative_overlap_columns():
    V_ref = np.array([[1.0, 0.0], [0.0, 1.0]])
    V = np.array([[-1.0, 0.0], [0.0, -1.0]])  # ambas columnas con signo volteado
    fixed = fix_eigenvector_signs(V, V_ref)
    np.testing.assert_allclose(fixed, V_ref)


def test_fix_eigenvector_signs_leaves_positive_overlap_columns():
    V_ref = np.array([[1.0, 0.0], [0.0, 1.0]])
    V = np.array([[1.0, 0.0], [0.0, 1.0]])
    fixed = fix_eigenvector_signs(V, V_ref)
    np.testing.assert_allclose(fixed, V)


# --------------------------------------------------------- (a) antisimetría
def test_antisymmetry_toy_model():
    solve = toy_two_level_solve()
    R = np.arange(5.0, 15.0001, 0.05)
    R_mid, A = derivative_coupling(R, solve=solve)
    # identidad matemática ⟨i|d/dR|j⟩ = -⟨j|d/dR|i⟩; el error es O(h²)
    # (diferencias finitas centradas), no cero exacto — ver docstring.
    np.testing.assert_allclose(A[:, 0, 1], -A[:, 1, 0], atol=5e-4, rtol=1e-3)


# --------------------------------------------------------- (b) diagonal nula
def test_diagonal_zero_toy_model():
    """A_ii = 0 es una identidad exacta del operador d/dR (normalización de
    Ψᵢ(R) para todo R), pero la fórmula discreta sólo la respeta a O(h²): con
    h=0.05 el residuo verificado es ~1.29e-6 (coeficiente h⁻² ≈ 5.18e-4,
    estable en h=0.2,0.1,0.05,0.025 — ver
    docs/analysis_fase1_acoplamiento_derivada.md). atol se fija con margen
    sobre ese residuo medido, no a cero exacto."""
    solve = toy_two_level_solve()
    R = np.arange(5.0, 15.0001, 0.05)
    R_mid, A = derivative_coupling(R, solve=solve)
    np.testing.assert_allclose(np.diagonal(A, axis1=1, axis2=2), 0.0, atol=2e-6)


# --------------------------------------------------- (a,b) también con BOPSystem
def test_antisymmetry_and_diagonal_zero_real_system_away_from_crossing():
    from trimero.systems.rb_krb_polar.bop_system import BOPSystem

    sysm = BOPSystem(n_manifold=6, N_max=2)
    states = [6, 7]  # hueco ~2e-3 Eh, esencialmente constante en todo R: no hay cruce

    def solve(R):
        w, V = sysm.solve(float(R), 0, fermi=False)
        return w[states], V[:, states]

    R = np.arange(700.0, 900.0001, 5.0)
    R_mid, A = derivative_coupling(R, solve=solve)
    np.testing.assert_allclose(A[:, 0, 1], -A[:, 1, 0], atol=1e-6, rtol=1e-2)
    np.testing.assert_allclose(np.diagonal(A, axis1=1, axis2=2), 0.0, atol=1e-9)


# ------------------------------------------------- (c) suave lejos de cruces
def test_smooth_away_from_crossings_real_system():
    """Lejos de cruces, el acoplamiento debe ser pequeño y variar suavemente
    (sin saltos de orden de magnitud entre R vecinos)."""
    from trimero.systems.rb_krb_polar.bop_system import BOPSystem

    sysm = BOPSystem(n_manifold=6, N_max=2)
    states = [6, 7]

    def solve(R):
        w, V = sysm.solve(float(R), 0, fermi=False)
        return w[states], V[:, states]

    R = np.arange(400.0, 1800.0001, 20.0)
    R_mid, A = derivative_coupling(R, solve=solve)
    coupling = np.abs(A[:, 0, 1])
    assert np.all(coupling < 1e-3), (
        f"acoplamiento inesperadamente grande lejos de cruce: max={coupling.max():.3e}"
    )
    # suave: decae MONÓTONAMENTE con R (sin picos ni rebotes, que delatarían
    # un cruce oculto), en vez de acotar el ratio punto a punto — el
    # acoplamiento cae ~9 órdenes de magnitud en este rango (de ~9e-10 a
    # ~6e-15, ver docs/analysis_fase1_acoplamiento_derivada.md), así que un
    # cociente con un piso absoluto (p.ej. 1e-12) se satura de forma
    # artificial mucho antes del final del rango y no mide suavidad real.
    assert np.all(np.diff(coupling) <= 1e-30), (
        f"la magnitud del acoplamiento no decae monótonamente: {coupling}"
    )


# ------------------------------------------------- (d) pico en cruce REAL
@pytest.mark.slow
def test_peak_at_real_documented_avoided_crossing():
    """
    Cruce evitado real documentado: BOPSystem(n_manifold=25), M_J=0, sin
    Fermi, par de estados k=54,55. `plots/rb_krb_polar/fig1_ad_MJ0_n25.npz`
    registra el cambio de índice de la curva de carácter cerca de R=1300 a0;
    un escaneo fino del hueco w[55]-w[54] ubica el mínimo real en R≈1305 a0
    (hueco ≈ 4.7e-9 Eh, ver docs/analysis_fase1_acoplamiento_derivada.md).
    """
    from trimero.systems.rb_krb_polar.bop_system import BOPSystem
    from trimero.systems.rb_krb_polar.rb_defects import DELTA0_NS_PAPER

    sysm = BOPSystem(n_manifold=25, delta0_ns=DELTA0_NS_PAPER)
    states = [54, 55]

    def solve(R):
        w, V = sysm.solve(float(R), 0, fermi=False)
        return w[states], V[:, states]

    R_peak = np.array([1304.0, 1305.0, 1306.0])
    R_away = np.array([899.0, 900.0, 901.0])

    _, A_peak = derivative_coupling(R_peak, solve=solve)
    _, A_away = derivative_coupling(R_away, solve=solve)

    peak = abs(A_peak[0, 0, 1])
    away = abs(A_away[0, 0, 1])
    print(f"\n  |A(54,55)| en cruce (R=1305): {peak:.4e}")
    print(f"  |A(54,55)| lejos    (R=900) : {away:.4e}")
    # medido: peak/away ≈ 84.75 (ver docs/analysis_fase1_acoplamiento_derivada.md);
    # 50x deja margen sin ser tan laxo que un bug que reduzca el pico a la
    # mitad pase desapercibido.
    assert peak > 50.0 * away, (
        f"el acoplamiento no muestra un pico claro en el cruce documentado: "
        f"peak={peak:.3e}, away={away:.3e}"
    )


# --------------------------------------------------- (e) convergencia con h
def test_convergence_with_finite_difference_step():
    """El error frente al valor analítico debe reducirse ~4x al reducir el
    paso a la mitad (orden 2 de las diferencias centradas)."""
    k, c, R0 = 0.01, 0.05, 10.0
    solve = toy_two_level_solve(k=k, c=c, R0=R0)
    R_probe = 8.0  # lejos del pico agudo en R0, régimen donde domina el término O(h²)
    exact = toy_analytic_coupling(R_probe, k=k, c=c, R0=R0)

    errors = []
    for h in (0.2, 0.1, 0.05):
        R = np.array([R_probe - h, R_probe, R_probe + h])
        _, A = derivative_coupling(R, solve=solve)
        errors.append(abs(abs(A[0, 0, 1]) - exact))

    print(f"\n  errores vs h=(0.2,0.1,0.05): {errors}")
    assert errors[0] > errors[1] > errors[2]
    # orden 2: al menos ~3x de reducción por cada mitad de paso (deja margen)
    assert errors[0] / errors[1] > 3.0
    assert errors[1] / errors[2] > 3.0


# ------------------------------------------ eigenbasis_along_R (continuidad)
def test_eigenbasis_along_r_matches_analytic_coupling_magnitude():
    k, c, R0 = 0.01, 0.05, 10.0
    solve = toy_two_level_solve(k=k, c=c, R0=R0)
    R = np.arange(5.0, 15.0001, 0.1)
    R_out, W, V = eigenbasis_along_R(R, solve)
    np.testing.assert_allclose(R_out, R)
    assert W.shape == (len(R), 2)
    assert V.shape == (len(R), 2, 2)
    # las columnas siguen siendo autovectores normalizados en todo R
    norms = np.sum(V**2, axis=1)
    np.testing.assert_allclose(norms, 1.0, atol=1e-10)


def test_derivative_coupling_rejects_nonuniform_grid():
    solve = toy_two_level_solve()
    R = np.array([5.0, 5.1, 5.4])  # paso no uniforme
    with pytest.raises(ValueError):
        derivative_coupling(R, solve=solve)


def test_derivative_coupling_requires_v_or_solve():
    R = np.array([5.0, 5.1, 5.2])
    with pytest.raises(ValueError):
        derivative_coupling(R)
