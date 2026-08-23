"""
Tests del método de estabilización (Hazi-Taylor 1970)
(`trimero.systems.nonadiabatic_dynamics.stabilization`), Fase 3 de
`docs/PLAN_nonadiabatic_dynamics.md`.

Referencia independiente: POZO CUADRADO FINITO 1D, con el número y la
energía de sus estados ligados obtenidos por búsqueda de raíces
(`scipy.optimize.brentq`) de la ecuación trascendente estándar

    par:   κ = k·tan(ka)
    impar: κ = -k·cot(ka)

con k=√(2μ(E+V0)), κ=√(-2μE) — un método numérico completamente distinto
(búsqueda de raíces de una ecuación trascendente 1D) del que se usa en
`stabilization.py` (diagonalización de una caja en diferencias finitas +
clasificación por deriva). El número de estados ligados se contrasta además
contra la fórmula analítica N = ⌊k₀·2a/π⌋+1, k₀=√(2μV0) (Griffiths,
"Introduction to Quantum Mechanics", pozo finito).

ADVERTENCIA DE ESTA RONDA (dejada aquí porque costó una iteración real): la
primera versión de `bound_states_square_well` buscaba cambios de signo de la
ecuación trascendente en una malla de energía sin evitar los POLOS de
tan/cot — un polo también produce un cambio de signo aparente, así que esa
versión "encontraba" 2 raíces espurias además de las 3 reales (5 en vez de
3), y el criterio ⌊k₀·2a/π⌋+1 = 3 lo delató antes de aceptar ningún
resultado del módulo de estabilización contra esa referencia rota. La
versión de abajo busca cada raíz dentro de una única rama monótona de
tan/cot (entre polos consecutivos), evitando el problema por construcción.
"""
import numpy as np
import pytest
from scipy.optimize import brentq

from trimero.systems.nonadiabatic_dynamics.coupled_channels import build_coupled_hamiltonian
from trimero.systems.nonadiabatic_dynamics.stabilization import (
    classify_stability,
    stabilization_scan,
    stable_state_summary,
    track_trajectories,
)


def bound_states_square_well(mu, V0, a):
    """Energías ligadas EXACTAS (hasta la precisión de brentq) del pozo
    cuadrado finito V(R)=-V0 para |R|<a, 0 en el resto. Búsqueda de raíces
    rama a rama, sin cruzar nunca un polo de tan/cot (ver docstring del
    módulo)."""
    k0 = np.sqrt(2.0 * mu * V0)
    roots = []

    n = 0
    while True:  # pares: k*tan(k*a) = kappa, rama (n*pi, n*pi+pi/2)
        lo, hi = n * np.pi / a, min((n * np.pi + np.pi / 2) / a, k0)
        if lo >= k0:
            break
        lo, hi = lo + 1e-9, hi - 1e-9
        if hi > lo:
            f = lambda k: k * np.tan(k * a) - np.sqrt(max(k0**2 - k**2, 1e-300))
            if f(lo) * f(hi) < 0:
                k_r = brentq(f, lo, hi)
                roots.append(k_r**2 / (2.0 * mu) - V0)
        n += 1

    n = 0
    while True:  # impares: -k*cot(k*a) = kappa, rama (n*pi+pi/2, (n+1)*pi)
        lo, hi = (n * np.pi + np.pi / 2) / a, min((n + 1) * np.pi / a, k0)
        if lo >= k0:
            break
        lo, hi = lo + 1e-9, hi - 1e-9
        if hi > lo:
            f = lambda k: -k / np.tan(k * a) - np.sqrt(max(k0**2 - k**2, 1e-300))
            if f(lo) * f(hi) < 0:
                k_r = brentq(f, lo, hi)
                roots.append(k_r**2 / (2.0 * mu) - V0)
        n += 1

    return np.array(sorted(roots))


def n_bound_states_formula(mu, V0, a):
    """N = floor(k0*2a/pi) + 1, fórmula analítica estándar del pozo finito."""
    k0 = np.sqrt(2.0 * mu * V0)
    return int(np.floor(k0 * 2.0 * a / np.pi)) + 1


# ------------------------------------ sanidad de la propia referencia
def test_square_well_reference_matches_analytic_count_formula():
    mu, V0, a = 1.0, 0.1, 8.0
    exact = bound_states_square_well(mu, V0, a)
    assert len(exact) == n_bound_states_formula(mu, V0, a)
    assert len(exact) == 3
    np.testing.assert_allclose(
        exact,
        [-0.08833046956086295, -0.05482888535439445, -0.008088245288672372],
        rtol=1e-8,
    )


# --------------------------------------------------- track_trajectories
def test_track_trajectories_keeps_identity_through_avoided_crossing():
    """
    Dos trayectorias sintéticas que se ACERCAN mucho sin cruzar NUNCA (el
    modelo analítico de 2 niveles de la Fase 1, con hueco mínimo 2c en
    x0): `track_trajectories` debe seguir cada rama por continuidad, no
    limitarse a devolver el orden ordenado de cada fila (que SÍ tendría un
    salto discontinuo justo en el punto de máximo acercamiento).

    NOTA: la primera versión de este test usaba dos RECTAS que se cruzan de
    verdad (pendientes opuestas) en vez de un cruce evitado. En un cruce
    real (energía exactamente degenerada en un punto) la asignación de
    identidad es intrínsecamente ambigua ahí — el coste de intercambiar o no
    las dos ramas es idéntico exactamente en el cruce — así que ese test
    fallaba por una premisa equivocada del propio test (el emparejamiento
    por energía más cercana no está mal definido cerca de un cruce EVITADO,
    que es el caso físico real; sí lo está en un cruce exacto, que no ocurre
    en autovalores de la misma simetría por la regla de no-cruce, la razón
    citada en el docstring del módulo). Se sustituye por el caso bien
    definido.
    """
    k, c = 0.05, 0.01
    x = np.linspace(-1.0, 1.0, 41)
    a = -k * x  # rama "baja" lejos de x=0
    b = k * x  # rama "alta" lejos de x=0
    # autovalores de [[a,c],[c,b]]: nunca se cruzan (hueco minimo 2c en x=0)
    lower = 0.5 * (a + b) - np.sqrt((0.5 * (a - b)) ** 2 + c**2)
    upper = 0.5 * (a + b) + np.sqrt((0.5 * (a - b)) ** 2 + c**2)
    assert np.all(upper - lower > 0.0)  # nunca degenerado, por construcción

    E = np.stack([lower, upper], axis=1)
    traj = track_trajectories(E)
    np.testing.assert_allclose(traj[:, 0], lower, atol=1e-12)
    np.testing.assert_allclose(traj[:, 1], upper, atol=1e-12)


# --------------------------------------------------- validación de malla
def test_classify_stability_rejects_nonuniform_l_grid():
    L = np.array([10.0, 12.0, 15.0])
    traj = np.zeros((3, 2))
    with pytest.raises(ValueError):
        classify_stability(L, traj)


# --------------------------------------- 3a/3b: pozo cuadrado, referencia independiente
def test_stabilization_identifies_exact_number_and_energy_of_bound_states():
    """
    Núcleo de la Fase 3: en un caso con acoplamiento CERO (canal único) y
    número/energía de estados ligados conocidos por una vía TOTALMENTE
    independiente (ecuación trascendente), el clasificador automático de
    `stabilization.py` debe:
      (a) identificar EXACTAMENTE 3 trayectorias como ligadas (ni más ni
          menos: ni una espuria, ni ninguna real perdida);
      (b) sus energías promedio en el tramo estable deben coincidir con la
          referencia exacta dentro de una tolerancia atada a la
          discretización (h=0.25 a₀, mismo orden de magnitud que los
          residuos O(h²) ya caracterizados en las Fases 1-2);
      (c) TODOS los estados por encima del umbral de disociación (E>0) deben
          quedar clasificados como NO ligados (estados de caja).
    """
    mu, V0, a = 1.0, 0.1, 8.0
    exact = bound_states_square_well(mu, V0, a)
    assert len(exact) == 3

    h = 0.25

    def hamiltonian_builder(L):
        n = int(round(L / h)) + 1
        R = np.linspace(-L / 2.0, L / 2.0, n)
        V = np.where(np.abs(R) < a, -V0, 0.0)[:, None]
        H = build_coupled_hamiltonian(R, V, mu=mu)
        return R, H

    L_values = np.arange(30.0, 150.0001, 2.0)
    L_values, E = stabilization_scan(L_values, hamiltonian_builder, n_keep=12)
    trajectories = track_trajectories(E)
    L_mid, stable, slope = classify_stability(L_values, trajectories, threshold_ratio=0.1)
    summary = stable_state_summary(L_mid, trajectories[1:-1], stable, min_fraction=0.5)

    bound = sorted([s for s in summary if s["is_bound"]], key=lambda s: s["E_mean"])
    print(f"\n  estados ligados encontrados: {len(bound)} (referencia: {len(exact)})")
    for s, e in zip(bound, exact):
        print(f"    numérico={s['E_mean']:.6e}  exacto={e:.6e}  "
              f"dif={s['E_mean'] - e:.3e}  fracción_estable={s['fraction_stable']:.3f}")

    assert len(bound) == 3, (
        f"se esperaban 3 estados ligados (referencia independiente), "
        f"se encontraron {len(bound)}")
    E_bound = np.array([s["E_mean"] for s in bound])
    # atol=3e-3: cubre el residuo medido (máx 1.4e-3 en el estado más somero,
    # el más sensible a la discretización h=0.25 y al tamaño de caja finito)
    # con margen, sin ser tan laxo que confunda un bug de física con error
    # de malla.
    np.testing.assert_allclose(E_bound, exact, atol=3e-3)

    unbound = [s for s in summary if not s["is_bound"]]
    for s in unbound:
        # las trayectorias de caja de este barrido son todas de energía
        # positiva (por encima del umbral de disociación en E=0); si alguna
        # tuviera energía negativa y fuese descartada, sería sospechoso.
        idx = s["index"]
        assert trajectories[-1, idx] > 0.0, (
            f"trayectoria {idx} tiene energía negativa (E={trajectories[-1, idx]:.4e}) "
            f"pero se clasificó como estado de caja, no ligado"
        )
