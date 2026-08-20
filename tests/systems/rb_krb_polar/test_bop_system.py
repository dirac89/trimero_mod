"""
Tests del montaje parametrizado por manifold (`systems/rb_krb_polar/bop_system.py`).

Lo que hay que proteger aquí son dos cosas distintas:

  1. Que la base es la del paper —manifold (n,l≥3) + (n+1)d + (n+2)p + (n+3)s—
     y que se ve exactamente en qué se diferencia de la base incompleta que se
     usó antes de leer arXiv:1507.07972 (test S1).
  2. Que todo lo que depende del manifold sale de la física y no de números
     escritos a mano: qué nivel fija el dominio, umbrales rotacionales y el
     desplazamiento en R del polo de la resonancia p (tests S2, S3, S4).

Ejecutar:  poetry run pytest tests/simulation/test_bop_system.py -s
"""

import numpy as np

from trimero.systems.rb_krb_polar.charge_dipole import B_KRB_GHZ
from trimero.systems.rb_krb_polar.bop_system import GHZ_PER_HARTREE, BOPSystem
from trimero.systems.rb_krb_polar.rb_defects import DELTA0_NS_PAPER as D0, n_star_ns


# ================================================================ TEST S1
def test_s1_n24_base_correcta():
    """
    n=24 con la base del paper: manifold (24, l≥3) + 25d + 26p + 27s.

    Los números de referencia de las rondas anteriores (dim 1016, dominio hasta
    1138.92 a0) correspondían a la base INCOMPLETA, sólo manifold + 27s. Aquí se
    comprueban los dos, y que la diferencia es exactamente la esperada.
    """
    s = BOPSystem(n_manifold=24, delta0_ns=D0)
    old = BOPSystem(n_manifold=24, delta0_ns=D0, neighbors=(0,))
    print(f"\n   niveles vecinos: {s.levels}   (l -> n)")
    print(f"   dim(M_J=0): {len(old.block(0))} (incompleta) -> {len(s.block(0))}")
    print(f"   n enteros de las radiales: {s.radial.neighbor_n_eff}")
    lo, hi = s.domain_bounds()
    lo0, hi0 = old.domain_bounds()
    print(f"   dominio: [{lo0:.4f}, {hi0:.4f}] -> [{lo:.4f}, {hi:.4f}] a0")
    print(f"   ΔE(27s) = {s.delta_E_ns_ghz():.5f} GHz")

    assert s.levels == {0: 27, 1: 26, 2: 25}
    assert len(s.block(0)) == 1064 and len(old.block(0)) == 1016
    assert s.basis.total_dimension() == 576 * 49
    assert s.l_max == 23
    assert s.radial.n_of_l(0) == 24 and s.radial.n_of_l(1) == 23
    assert s.radial.n_of_l(2) == 24 and s.radial.n_of_l(5) == 24
    # el umbral N=0 del 27s no cambia: es energía atómica, no depende de la base
    assert abs(s.delta_E_ns_ghz() - (-63.40)) < 0.01
    assert np.isclose(s.delta_E_ns_ghz(), old.delta_E_ns_ghz())
    assert np.isclose(s.E_manifold, -0.5 / 24**2)
    # la ventana de la resonancia p tampoco: depende del manifold, no de los
    # vecinos (docs/analysis_ventana_exclusion_resonancia.md)
    w = s.p_resonance_window()
    assert abs(w["R_pole"] - 562.771) < 0.01
    assert abs(w["fwhm"] - 13.187) < 0.01


# ================================================================ TEST S2
def test_s2_dominio_lo_fija_el_np():
    """
    Completar la base ACORTA el dominio del remapeo k(R): el par más ligado deja
    de ser (n+3)s y pasa a ser (n+2)p, que tiene n* menor.

    Ojo con la lectura: ese tope es el límite del REMAPEO SEMICLÁSICO con el que
    leemos las tablas de dispersión, no el punto donde se acaba la física.
    """
    for n, l_binding in ((24, 1), (25, 1)):
        s = BOPSystem(n_manifold=n, delta0_ns=D0)
        lo, hi = s.domain_bounds()
        turns = {l: 2.0 * s.n_star_of_l(l) ** 2 for l in (0, 1, 2, 3)}
        print(f"\n   n={n}: dominio [{lo:.4f}, {hi:.4f}] a0")
        print("     puntos de retorno clásicos 2 n*^2 por nivel: "
              + ", ".join(f"l={l}: {t:.1f}" for l, t in turns.items()))

        # el nivel más ligado de la base es (n+2)p
        assert min(turns, key=turns.get) == l_binding
        # el tope del dominio queda justo por debajo de su retorno clásico
        assert hi < turns[l_binding]
        assert turns[l_binding] - hi < 1.0
        # y coincide con lo que impone el borde de la tabla, sin ajustar nada
        E = -0.5 / s.n_star_of_l(l_binding) ** 2
        R_max_teo = 1.0 / (1.0 / s.scattering.R_table.max()
                           + s.scattering.E_table - E)
        assert abs(hi - R_max_teo) < 1e-5

    # n=25 concreto: la base incompleta llegaba a 1236.3, la correcta no
    s25 = BOPSystem(n_manifold=25, delta0_ns=D0)
    old25 = BOPSystem(n_manifold=25, delta0_ns=D0, neighbors=(0,))
    print(f"\n   n=25: R_max {old25.domain_bounds()[1]:.2f} (incompleta) -> "
          f"{s25.domain_bounds()[1]:.2f} a0 (correcta)")
    assert s25.domain_bounds()[1] < old25.domain_bounds()[1]
    assert 1180.0 < s25.domain_bounds()[1] < 1190.0
    assert np.isclose(old25.n_star_s(), n_star_ns(28, D0))


# ================================================================ TEST S3
def test_s3_n25_umbrales():
    """Umbrales 28s + KRb(N): ΔE(28s) + B·N(N+1), con B = 1.114 GHz."""
    s = BOPSystem(n_manifold=25, delta0_ns=D0)
    dE = s.delta_E_ns_ghz()
    thr = s.thresholds_ns_ghz(range(0, 7))
    print(f"\n   ΔE(28s) = {dE:.4f} GHz")
    for N in range(7):
        print(f"     N={N}  B·N(N+1) = {B_KRB_GHZ*N*(N+1):8.4f}  "
              f"umbral = {thr[N]:9.4f} GHz")

    # ΔE es la diferencia de energías, sin rotación: umbral N=0
    assert np.isclose(thr[0], dE)
    assert np.isclose(dE, (-0.5 / s.n_star_s() ** 2 + 0.5 / 25**2) * GHZ_PER_HARTREE)
    # los dos que se piden explícitamente
    assert np.isclose(thr[5], dE + 30.0 * B_KRB_GHZ)
    assert np.isclose(thr[6], dE + 42.0 * B_KRB_GHZ)
    # el 28s está MENOS por debajo del manifold que el 27s del suyo
    s24 = BOPSystem(n_manifold=24, delta0_ns=D0)
    print(f"   ΔE(27s) = {s24.delta_E_ns_ghz():.4f} GHz  ->  |ΔE| decrece con n")
    assert dE < 0.0 and dE > s24.delta_E_ns_ghz()


# ================================================================ TEST S4
def test_s4_resonancia_misma_energia_distinto_R():
    """
    La resonancia de forma p es una propiedad de e⁻+Rb(5S): misma ENERGÍA para
    los dos manifolds, distinto R, porque eps = E_manifold + 1/R.
    """
    s24 = BOPSystem(n_manifold=24, delta0_ns=D0)
    s25 = BOPSystem(n_manifold=25, delta0_ns=D0)
    w24, w25 = s24.p_resonance_window(), s25.p_resonance_window()
    print(f"\n   eps_polo  n=24: {w24['eps_pole']:.10f}   n=25: {w25['eps_pole']:.10f}")
    print(f"   R_polo    n=24: {w24['R_pole']:.3f}      n=25: {w25['R_pole']:.3f}")
    print(f"   FWHM      n=24: {w24['fwhm']:.3f}       n=25: {w25['fwhm']:.3f}")

    # misma energía, exactamente: no depende del manifold
    assert np.isclose(w24["eps_pole"], w25["eps_pole"], rtol=0, atol=0)
    # distinto R, y con la relación exacta R = 1/(eps - E_manifold)
    assert np.isclose(w25["R_pole"], 1.0 / (w25["eps_pole"] - s25.E_manifold))
    assert w25["R_pole"] > w24["R_pole"]
    # a n mayor, el manifold está menos ligado -> el mismo eps cae a R mayor
    # y la estructura se ensancha en R
    assert w25["fwhm"] > w24["fwhm"]
    # la ventana sigue siendo simétrica alrededor del polo
    assert np.isclose(w25["R_lo"], w25["R_pole"] - 2.0 * w25["fwhm"])
    assert np.isclose(w25["R_hi"], w25["R_pole"] + 2.0 * w25["fwhm"])


# ================================================================ TEST S5
def test_s5_fermi_off_es_cero_exacto():
    """
    `fermi=False` pone V_Fermi ≡ 0 EXACTAMENTE, y el hueco que deja NO es
    pequeño en el borde del dominio.

    Lo segundo es el resultado que descartó extender la Fig. 1 más allá del
    remapeo: ver docs/analysis_extension_dominio_fig1.md §1.
    """
    s = BOPSystem(n_manifold=25, delta0_ns=D0)
    blk = s.block(0)
    R = 900.0
    H_full = s.hamiltonian(R)
    H_nof = s.hamiltonian(R, fermi=False)
    V = s.fermi.build(blk, R)

    # H con Fermi es EXACTAMENTE H sin Fermi más la matriz del pseudopotencial
    # (bit a bit: es la misma suma que hace `hamiltonian`, no una aproximación).
    # Ojo: (H+V)-H no es exacto en coma flotante, así que se compara al revés.
    assert np.array_equal(H_full, H_nof + V)
    assert np.linalg.norm(V) > 0.0

    # y en el borde del dominio el pseudopotencial NO es despreciable
    hi = s.domain_bounds()[1]
    e_full = s.character_curve(hi - 0.5)[0]
    e_nof = s.character_curve(hi - 0.5, fermi=False)[0]
    jump = e_full - e_nof
    print(f"\n   R_max = {hi:.2f} a0   E con Fermi = {e_full:.4f} GHz   "
          f"sin Fermi = {e_nof:.4f} GHz   salto = {jump:+.4f} GHz "
          f"({100*abs(jump/e_full):.1f} %)")
    assert abs(jump) > 5.0                      # GHz: muy lejos de despreciable
    assert abs(jump) > 0.25 * abs(e_full)       # >25 % de la ligadura
