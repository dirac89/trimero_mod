"""
Tests de la ecuación radial de canales acoplados
(`trimero.systems.nonadiabatic_dynamics.coupled_channels`), Fase 2 de
`docs/PLAN_nonadiabatic_dynamics.md`.

Cada test usa una referencia INDEPENDIENTE del propio método (nunca el
método comparado contra sí mismo con el acoplamiento a cero):

  - test 2a: oscilador armónico, energías E_n=(n+½)ħω ANALÍTICAS.
  - test 2a: dos canales con acoplamiento A=B=0 deben reproducir Los DOS
    solves de un solo canal, resueltos independientemente.
  - test 2b: doble pozo de dos canales armónicos con acoplamiento diabático
    CONSTANTE Δ. Se resuelve dos veces por caminos completamente distintos:
    (i) representación DIABÁTICA — matriz de bloques con acoplamiento local
    Δ·I, sin ninguna maquinaria de la Fase 1/2 (T + V + acoplamiento
    multiplicativo, el problema mejor caracterizado que hay); (ii)
    representación ADIABÁTICA — diagonalizando el Hamiltoniano diabático en
    cada R para obtener V_±(R) y las autofunciones, y construyendo A_ij, B_ij
    con la maquinaria real de este módulo. Las energías DEBEN coincidir: es
    un principio físico exacto (independencia de representación), no una
    aproximación — así que sirve de referencia fuerte sin depender de ningún
    resultado "ya conocido" de la literatura.
  - test de aplicación a curvas reales: cruce evitado real de BOPSystem
    (Fase 1d, R≈1305 a₀, estados 54/55). No se afirma que las energías
    resultantes sean niveles vibracionales físicos reales del sistema (la
    ventana de R usada, 50 a₀, es mucho menor que el dominio real de la
    molécula — ver docs/analysis_fase2_canales_acoplados.md): sólo se
    verifica que el acoplamiento tiene un efecto real y medible sobre el
    espectro frente a la solución desacoplada, consistente con el pico de
    A_ij ya encontrado en la Fase 1.
"""
import numpy as np
import pytest

from trimero.systems.nonadiabatic_dynamics.coupling import (
    derivative_coupling,
    eigenbasis_along_R,
)
from trimero.systems.nonadiabatic_dynamics.coupled_channels import (
    adiabatic_coupled_hamiltonian,
    build_coupled_hamiltonian,
    radial_kinetic_matrix,
    second_derivative_coupling,
    solve_coupled_channels,
)

GHZ_PER_HARTREE = 6.579683920502e6


# --------------------------------------------------------- 2a: 1 canal vs HO
def test_single_channel_matches_analytic_harmonic_oscillator():
    mu, omega = 1.0, 0.05
    L, n = 40.0, 400
    R = np.linspace(-L, L, n)
    V = 0.5 * mu * omega**2 * R**2

    E, _ = solve_coupled_channels(R, V[:, None], mu=mu)
    exact = np.array([(k + 0.5) * omega for k in range(6)])
    np.testing.assert_allclose(E[:6], exact, rtol=2e-3)


def test_single_channel_convergence_with_grid_spacing():
    """El error frente al oscilador armónico exacto debe caer ~4x al doblar
    la resolución de la malla (orden 2 de la diferencia centrada, igual que
    en la Fase 1)."""
    mu, omega = 1.0, 0.05
    L = 40.0
    exact0 = 0.5 * omega
    errors = []
    for n in (100, 200, 400):
        R = np.linspace(-L, L, n)
        V = 0.5 * mu * omega**2 * R**2
        E, _ = solve_coupled_channels(R, V[:, None], mu=mu)
        errors.append(abs(E[0] - exact0))
    print(f"\n  errores (n=100,200,400): {errors}")
    assert errors[0] > errors[1] > errors[2]
    assert errors[0] / errors[1] > 3.5
    assert errors[1] / errors[2] > 3.5


# ------------------------------------------------ 2a: acoplamiento a cero
def test_zero_coupling_recovers_two_independent_single_channel_solves():
    mu, omega = 1.0, 0.05
    L, n = 40.0, 400
    R = np.linspace(-L, L, n)
    V1 = 0.5 * mu * omega**2 * R**2
    V2 = 0.5 * mu * omega**2 * (R - 5.0) ** 2  # centro distinto, canal independiente
    V = np.stack([V1, V2], axis=1)

    E_coupled, _ = solve_coupled_channels(R, V, A=None, B=None, mu=mu)

    E1, _ = solve_coupled_channels(R, V1[:, None], mu=mu)
    E2, _ = solve_coupled_channels(R, V2[:, None], mu=mu)
    expected = np.sort(np.concatenate([E1[:8], E2[:8]]))[:8]

    np.testing.assert_allclose(np.sort(E_coupled[:8]), expected, rtol=1e-10)


# --------------------------------------------------- 2b: doble pozo, dos vías
def test_double_well_diabatic_and_adiabatic_representations_agree():
    mu, omega = 50.0, 0.02
    R1, R2 = -8.0, 8.0
    Delta = 0.01
    L, n = 30.0, 300
    R = np.linspace(-L, L, n)
    h = R[1] - R[0]

    V1 = 0.5 * mu * omega**2 * (R - R1) ** 2
    V2 = 0.5 * mu * omega**2 * (R - R2) ** 2

    # (i) referencia diabática: acoplamiento local constante, sin ninguna
    # maquinaria de la Fase 1/2, sólo T + V + Delta*I por bloques.
    T = radial_kinetic_matrix(n, h, mu)
    H_dia = np.zeros((2 * n, 2 * n))
    H_dia[:n, :n] = T + np.diag(V1)
    H_dia[n:, n:] = T + np.diag(V2)
    H_dia[:n, n:] = Delta * np.eye(n)
    H_dia[n:, :n] = Delta * np.eye(n)
    E_dia, _ = np.linalg.eigh(H_dia)

    # (ii) representación adiabática: se diagonaliza el mismo H_dia(R) en cada
    # R con la maquinaria de la Fase 1/2 (A_ij, B_ij por diferencias finitas).
    def solve(Rk):
        Hk = np.array([[0.5 * mu * omega**2 * (Rk - R1) ** 2, Delta],
                       [Delta, 0.5 * mu * omega**2 * (Rk - R2) ** 2]])
        return np.linalg.eigh(Hk)

    R_mid, H_adia = adiabatic_coupled_hamiltonian(R, solve, mu=mu)
    E_adia, _ = np.linalg.eigh(H_adia)

    print(f"\n  diabática  lowest 6: {E_dia[:6]}")
    print(f"  adiabática lowest 6: {E_adia[:6]}")
    print(f"  diferencia máxima  : {np.max(np.abs(E_adia[:6] - E_dia[:6])):.3e}")
    # dos vías de cálculo INDEPENDIENTES del mismo problema físico; deben
    # coincidir por un principio exacto (independencia de representación),
    # no por construcción del test.
    np.testing.assert_allclose(E_adia[:6], E_dia[:6], atol=1e-7)


# ---------------------------------------------------------- sanidad general
def test_build_coupled_hamiltonian_is_symmetric_and_finite():
    mu = 1.0
    n = 50
    R = np.linspace(-10.0, 10.0, n)
    V = np.stack([0.01 * R**2, 0.01 * (R - 1.0) ** 2], axis=1)
    rng_A = 1e-3 * np.sin(0.1 * R)[:, None, None] * np.array([[0.0, 1.0], [-1.0, 0.0]])
    B = np.zeros((n, 2, 2))
    H = build_coupled_hamiltonian(R, V, A=rng_A, B=B, mu=mu)
    assert np.all(np.isfinite(H))
    np.testing.assert_allclose(H, H.T, atol=1e-12)


# ------------------------------------------- aplicación a un cruce evitado real
@pytest.mark.slow
def test_real_avoided_crossing_coupling_shifts_the_spectrum():
    """
    Cruce evitado real de la Fase 1 (BOPSystem n_manifold=25, M_J=0, sin
    Fermi, estados 54/55, mínimo de hueco en R≈1305 a₀). Ventana de R
    limitada (50 a₀) por coste computacional (51 diagonalizaciones de
    1113x1113, ~65 s) — ver docs/analysis_fase2_canales_acoplados.md para la
    advertencia explícita sobre qué SÍ y qué NO se puede concluir de esta
    ventana. No se afirma que las energías resultantes sean niveles
    vibracionales físicos reales (para eso, Fase 3 y un dominio de R
    realista); sólo se verifica que la maquinaria produce una matriz sana
    (Hermítica, finita) y que el acoplamiento tiene un efecto real y no
    trivial sobre el espectro, tal como predice el pico de A_ij de la Fase 1.
    """
    from trimero.systems.rb_krb_polar.bop_system import BOPSystem
    from trimero.systems.rb_krb_polar.rb_defects import DELTA0_NS_PAPER

    sysm = BOPSystem(n_manifold=25, delta0_ns=DELTA0_NS_PAPER)
    states = [54, 55]

    def solve(R):
        w, V = sysm.solve(float(R), 0, fermi=False)
        return w[states], V[:, states]

    mu = 113522.0  # ~ reduced mass Rb-RbCs en m_e; ver docs/analysis_fase2_canales_acoplados.md
    R = np.arange(1280.0, 1330.0001, 1.0)

    R_out, W, V = eigenbasis_along_R(R, solve)
    R_mid, A = derivative_coupling(R_out, V=V)
    _, B = second_derivative_coupling(R_out, V=V)
    V_mid = W[1:-1]

    H_coupled = build_coupled_hamiltonian(R_mid, V_mid, A=A, B=B, mu=mu)
    H_decoupled = build_coupled_hamiltonian(R_mid, V_mid, A=None, B=None, mu=mu)
    assert np.all(np.isfinite(H_coupled))
    np.testing.assert_allclose(H_coupled, H_coupled.T, atol=1e-12)

    E_c, _ = np.linalg.eigh(H_coupled)
    E_d, _ = np.linalg.eigh(H_decoupled)

    shift_ghz = np.abs(E_c[:6] - E_d[:6]) * GHZ_PER_HARTREE
    print(f"\n  acoplado   lowest 6 (GHz, rel. a min desacoplado): "
          f"{(E_c[:6] - E_d[0]) * GHZ_PER_HARTREE}")
    print(f"  desacoplado lowest 6 (GHz, rel. a min desacoplado): "
          f"{(E_d[:6] - E_d[0]) * GHZ_PER_HARTREE}")
    print(f"  |desplazamiento| por nivel (GHz): {shift_ghz}")
    # medido: max ~0.38 GHz (ver docs/analysis_fase2_canales_acoplados.md);
    # se exige un efecto claramente no trivial (>50 MHz) y físicamente
    # acotado (<10 GHz, la escala de la propia curva en esta ventana).
    assert np.max(shift_ghz) > 0.05
    assert np.max(shift_ghz) < 10.0
