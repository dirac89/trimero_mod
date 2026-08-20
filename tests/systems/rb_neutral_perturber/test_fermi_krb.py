"""
Tests del pseudopotencial de Fermi V(r) para Rb*-KRb sobre CoupledBasis.

González-Férez, Sadeghpour & Schmelcher, NJP 17, 013021 (2015).

Geometría: ion Rb⁺ en el origen, KRb en R⃗, eje de cuantización Z ∥ R⃗ (mismo
convenio de las rondas anteriores). El perturbador está POR TANTO sobre el eje
(θ=0), lo que restringe fuertemente el acoplamiento.

La referencia independiente de esta suite es `brute_force_fermi_element`:
evalúa ψ y ∇ψ numéricamente en el punto cartesiano (0,0,R) por diferencias
finitas de la función de onda completa, sin usar ninguna de las formas cerradas
del módulo de producción.

Ejecutar:  poetry run pytest tests/systems/rb_neutral_perturber/test_fermi_krb.py -s
"""

import numpy as np
import pytest
from scipy.special import gammaln, lpmv

from trimero.basis.quantum import CoupledBasis
from trimero.basis.radial import RadialBasis
from trimero.systems.rb_krb_polar.charge_dipole import (
    ChargeDipoleHamiltonian,
    RydbergElectronField,
    rydberg_diagonal,
)
from trimero.systems.rb_neutral_perturber.fermi_krb import (
    E_TABLE_HARTREE,
    N_STAR_TABLE,
    FermiPseudopotential,
    ScatteringLengths,
    hydrogenic_R,
    hydrogenic_dR,
    n_star_of_l,
)

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")


# ---------------------------------------------------------------- referencia
def _Ylm(l, m, theta, phi):
    """Y_lm complejo con fase de Condon-Shortley, construido desde lpmv."""
    am = abs(m)
    c = np.sqrt(
        (2 * l + 1) / (4 * np.pi) * np.exp(gammaln(l - am + 1) - gammaln(l + am + 1))
    )
    val = c * lpmv(am, l, np.cos(theta)) * np.exp(1j * am * phi)
    if m < 0:
        val = (-1) ** am * np.conj(val)
    return val


def _psi(n, l, m, x, y, z):
    r = np.sqrt(x * x + y * y + z * z)
    return hydrogenic_R(n, l, r) * _Ylm(l, m, np.arccos(z / r), np.arctan2(y, x))


def _grad_psi(n, l, m, R, h):
    """∇ψ en (0,0,R) por diferencias finitas centradas en cartesianas."""
    g = []
    for axis in range(3):
        p = [0.0, 0.0, R]
        q = [0.0, 0.0, R]
        p[axis] += h
        q[axis] -= h
        g.append((_psi(n, l, m, *p) - _psi(n, l, m, *q)) / (2.0 * h))
    return np.array(g)


def brute_force_fermi_element(n_i, l_i, m_i, n_j, l_j, m_j, R, A_s, A_p, h=0.05):
    """
    ⟨i| 2π A_s δ³(r-R) + 6π A_p δ³(r-R) ∇·∇ |j⟩ evaluado numéricamente.

    No usa ninguna forma cerrada del módulo de producción: sólo ψ y ∇ψ.
    """
    psi_i = _psi(n_i, l_i, m_i, 0.0, 0.0, R)
    psi_j = _psi(n_j, l_j, m_j, 0.0, 0.0, R)
    vs = 2.0 * np.pi * A_s * np.conj(psi_i) * psi_j
    gi = _grad_psi(n_i, l_i, m_i, R, h)
    gj = _grad_psi(n_j, l_j, m_j, R, h)
    vp = 6.0 * np.pi * A_p * np.dot(np.conj(gi), gj)
    return complex(vs + vp)


# ---------------------------------------------------------------- fixtures
BASIS = CoupledBasis(N_max=6, manifold_l_min=3)
RADIAL = RadialBasis()
SCAT = ScatteringLengths()


# ================================================================ TEST F0
def test_f0_radial_derivative():
    """dR_nl/dr analítica contra diferencias finitas centradas."""
    print("\n  dR_nl/dr analítica vs diferencias finitas centradas (h=1e-4):")
    print("    n   l      r    |   analítica      |  dif. finitas    |   rel")
    print("   ---------------- |------------------|------------------|--------")
    worst = 0.0
    h = 1e-4
    for (n, l) in [(24, 0), (24, 3), (24, 12), (24, 23)]:
        for r in (120.0, 600.0, 900.0, 1400.0):
            ana = hydrogenic_dR(n, l, r)
            num = (hydrogenic_R(n, l, r + h) - hydrogenic_R(n, l, r - h)) / (2 * h)
            scale = max(abs(ana), abs(num))
            rel = abs(ana - num) / scale if scale > 0 else 0.0
            worst = max(worst, rel)
            print(f"   {n:3d} {l:3d} {r:7.1f} | {ana:+.10e} | {num:+.10e} | {rel:.1e}")
    print(f"   peor diferencia relativa = {worst:.2e}")
    assert worst < 1e-6, worst


# ================================================================ TEST F1
def test_f1_remap_domain():
    """Remapeo k(R): dominio válido y detección de casos fuera de tabla."""
    print(f"\n  Tabla original: n*={N_STAR_TABLE}, E={E_TABLE_HARTREE:.6e} E_h, "
          f"R ∈ [{SCAT.R_table.min():.0f}, {SCAT.R_table.max():.0f}] a0")
    print(f"  n* usado: l≥3 -> {n_star_of_l(3):.6f} ;  27s (l=0) -> {n_star_of_l(0):.6f}")
    print("\n      R    | n*      |   k(R)    |    R'     | en tabla | motivo")
    print("   --------|---------|-----------|-----------|----------|--------")
    for R in (300.0, 600.0, 900.0, 1100.0, 1152.0, 1500.0, 3000.0):
        for l in (3, 0):
            ns = n_star_of_l(l)
            k = SCAT.k_of_R(R, ns)
            ok = SCAT.in_domain(R, ns)
            Rp = SCAT.remap_R(R, ns)
            why = "" if ok else ("clásicamente prohibido" if not np.isfinite(k)
                                 else "R' fuera de [111,2448]")
            print(f"   {R:8.1f}| {ns:8.5f}| {k:9.6f} | {Rp:9.2f} | {str(ok):^8s} | {why}")
        print()
    # el rango que pediste (600-3000) NO está entero en dominio: se reporta
    assert SCAT.in_domain(600.0, n_star_of_l(3))
    assert SCAT.in_domain(900.0, n_star_of_l(3))
    assert not SCAT.in_domain(1500.0, n_star_of_l(3))
    with pytest.raises(ValueError):
        SCAT.scattering(1500.0, n_star_of_l(3))
    print("   scattering(R=1500) lanza ValueError en vez de extrapolar ✓")

    A_s, A_p = SCAT.scattering(900.0, n_star_of_l(3))
    print(f"\n   A_s(R=900, n*=24) = {A_s:+.6f} a0      "
          f"A_p(R=900, n*=24) = {A_p:+.4f} a0³")


# ================================================================ TEST F2
def test_f2_expansion_vs_brute_force():
    """Formas cerradas del módulo contra ψ y ∇ψ numéricos en (0,0,R)."""
    fermi = FermiPseudopotential(RADIAL, SCAT)
    R = 900.0
    A_s, A_p = SCAT.scattering(R, n_star_of_l(3))
    casos = [(3, 0, 3, 0), (5, 0, 5, 0), (5, 0, 8, 0), (23, 0, 22, 0),
             (5, 1, 5, 1), (6, 1, 9, 1), (5, -1, 7, -1), (5, 2, 5, 2),
             (5, 0, 5, 1), (7, 3, 7, 3)]
    print(f"\n  R = {R:g} a0,  A_s = {A_s:.6f},  A_p = {A_p:.4f}")
    print("   l1 m1 l2 m2 |   cerrada        |  fuerza bruta    |   rel")
    print("   ------------|------------------|------------------|--------")
    vals = []
    for (l1, m1, l2, m2) in casos:
        n1, n2 = RADIAL.n_of_l(l1), RADIAL.n_of_l(l2)
        closed = fermi.electron_element(l1, m1, l2, m2, R)
        bf = brute_force_fermi_element(n1, l1, m1, n2, l2, m2, R, A_s, A_p)
        assert abs(bf.imag) < 1e-12 * max(abs(bf), 1e-30), f"parte imaginaria {bf}"
        vals.append((l1, m1, l2, m2, closed, bf.real))
    # Los casos |m_l| >= 2 dan CERO EXACTO en la forma cerrada; la diferencia
    # finita deja ruido de redondeo (~1e-49). Un cociente relativo frente a ese
    # ruido no significa nada, así que la escala es la del mayor elemento real.
    floor = 1e-12 * max(abs(c) for *_, c, _ in vals)
    worst = 0.0
    for (l1, m1, l2, m2, closed, bf) in vals:
        rel = abs(closed - bf) / max(abs(closed), abs(bf), floor)
        worst = max(worst, rel)
        nota = "  <- cero exacto (|m|>=2 o Δm!=0)" if abs(closed) < floor else ""
        print(f"   {l1:2d}{m1:3d}{l2:3d}{m2:3d} | {closed:+.10e} | {bf:+.10e} | {rel:.1e}{nota}")
    print(f"   escala de corte (cero numérico) = {floor:.2e}")
    print(f"   peor diferencia relativa = {worst:.2e}")
    assert worst < 1e-5, worst

    # El residuo ~1e-7 es de la diferencia finita, no de la forma cerrada:
    # se COMPRUEBA reduciendo el paso h.
    print("\n   convergencia en el paso h (elemento (5,0)<-(8,0)):")
    closed = fermi.electron_element(5, 0, 8, 0, R)
    prev = None
    for h in (0.4, 0.2, 0.1, 0.05):
        bf = brute_force_fermi_element(24, 5, 0, 24, 8, 0, R, A_s, A_p, h=h).real
        rel = abs(bf / closed - 1.0)
        print(f"     h = {h:5.2f}:  bruta = {bf:+.12e}   rel = {rel:.2e}")
        if prev is not None:
            assert rel < prev, "la diferencia finita no converge a la forma cerrada"
        prev = rel
    print("   converge monótonamente hacia la forma cerrada ✓")


# ================================================================ TEST F3c
def test_f3c_mj_conservation():
    """V(r) conserva M_J; se verifica, no se asume."""
    fermi = FermiPseudopotential(RADIAL, SCAT)
    R = 900.0
    print("\n  V_Fermi es un operador puramente ELECTRÓNICO: diagonal en (N, M_N).")
    print("  Con el perturbador sobre el eje (θ=0) sólo sobreviven |m_l| ≤ 1 y")
    print("  Δm_l = 0, luego ΔM_J = 0 — pero por AUSENCIA de acoplamiento al")
    print("  rotor, no por una regla de selección de momento angular como en H_mol.")

    pares, worst = 0, 0.0
    for MJa, MJb in [(0, 1), (0, -2), (2, 5)]:
        for xa in BASIS.get_block(MJa).states[:40]:
            for yb in BASIS.get_block(MJb).states[:40]:
                v = fermi._matrix_element_unchecked(xa, yb, R)
                worst = max(worst, abs(v))
                pares += 1
    print(f"\n   {pares} pares con M_J distinto SIN comprobación: max|elem| = {worst:.1e}")
    assert worst == 0.0

    with pytest.raises(ValueError):
        fermi.matrix_element((5, 3, 3, -3), (5, 3, 3, -2), R)
    print("   matrix_element entre M_J distintos lanza ValueError ✓")

    # diagonal en N y M_N, pero NO en l
    si = (5, 0, 3, 0)
    print(f"\n   ΔN≠0:   ⟨{si}|V|{(5,0,4,0)}⟩ = "
          f"{fermi.matrix_element(si,(5,0,4,0),R):.3e}  (debe ser 0)")
    assert fermi.matrix_element(si, (5, 0, 4, 0), R) == 0.0
    v_l = fermi.matrix_element(si, (8, 0, 3, 0), R)
    print(f"   Δl≠0:   ⟨{si}|V|{(8,0,3,0)}⟩ = {v_l:.6e}  (debe ser ≠ 0: mezcla l)")
    assert v_l != 0.0
    v_m2 = fermi.matrix_element((5, 2, 3, -2), (7, 2, 3, -2), R)
    print(f"   |m_l|=2:⟨(5,2,3,-2)|V|(7,2,3,-2)⟩ = {v_m2:.3e}  (debe ser 0: θ=0)")
    assert v_m2 == 0.0


# ================================================================ TEST F3b
def test_f3b_hermiticity_full():
    """Hermiticidad de H_a + H_mol + V(r)."""
    efield = RydbergElectronField(RADIAL)
    h_mol = ChargeDipoleHamiltonian(electron_field=efield)
    fermi = FermiPseudopotential(RADIAL, SCAT)
    R = 900.0

    # V sólo tiene soporte en |m_l| <= 1. Un bloque de |M_J| grande no contiene
    # NINGÚN estado así, y barrerlo daría max|V| = 0: el test sería vacuo. Se
    # barre justo el subespacio donde V vive.
    blk0 = BASIS.get_block(0)
    support = [st for st in blk0.states if abs(st[1]) <= 1]
    worst = scale = 0.0
    for x in support:
        for y in support:
            a = fermi.matrix_element(x, y, R)
            b = fermi.matrix_element(y, x, R)
            worst = max(worst, abs(a - b))
            scale = max(scale, abs(a))
    print(f"\n  Barrido completo del soporte de V (|m_l|<=1) en M_J=0: "
          f"{len(support)} estados, {len(support)**2} pares")
    print(f"    max|V_ij - V_ji| = {worst:.3e}  (max|V_ij| = {scale:.3e})")
    assert scale > 0.0, "el barrido de simetría es vacuo: V es idénticamente nulo ahí"
    # No es cero exacto: al intercambiar i<->j el producto
    # l1(l1+1)·l2(l2+1)·(2l1+1)(2l2+1) se agrupa en otro orden y la
    # multiplicación en coma flotante no es asociativa. ~4e-16 relativo.
    print(f"    relativo = {worst/scale:.2e}")
    assert worst <= 1e-15 * scale, (worst, scale)

    blk0 = BASIS.get_block(0)
    H = (np.diag(rydberg_diagonal(blk0)) + h_mol.build(blk0, R)
         + fermi.build(blk0, R))
    asym = np.linalg.norm(H - H.T)
    print(f"\n  H = H_a + H_mol + V_Fermi, bloque M_J=0 {H.shape}:")
    print(f"    ||H||_F = {np.linalg.norm(H):.9e}")
    print(f"    ||H - H.T||_F = {asym:.3e}  (relativo {asym/np.linalg.norm(H):.1e})")
    assert asym <= 1e-14 * np.linalg.norm(H)

    # build vs barrido dim² completo, en un bloque CON soporte de V
    blk7 = BASIS.get_block(7)
    Hb = fermi.build(blk7, R)
    Hr = fermi.build_reference(blk7, R)
    print(f"\n  build vs build_reference, bloque M_J=7 (dim {len(blk7)}):")
    print(f"    max|V| = {np.abs(Hr).max():.3e}  (no vacuo)")
    print(f"    allclose(rtol=1e-12) = {np.allclose(Hb, Hr, rtol=1e-12, atol=0.0)}, "
          f"max|dif| = {np.abs(Hb-Hr).max():.1e}")
    assert np.abs(Hr).max() > 0.0, "comparación vacua: V es nulo en este bloque"
    assert np.allclose(Hb, Hr, rtol=1e-12, atol=0.0)


# ================================================================ TEST F3a
def test_f3a_zero_scattering_limit():
    """A_s = A_p = 0  =>  el espectro vuelve a H_a + H_mol ya verificado."""
    efield = RydbergElectronField(RADIAL)
    h_mol = ChargeDipoleHamiltonian(electron_field=efield)
    fermi_off = FermiPseudopotential(RADIAL, ScatteringLengths(enabled=False))
    blk = BASIS.get_block(0)
    R = 900.0

    V = fermi_off.build(blk, R)
    print(f"\n  Con ScatteringLengths(enabled=False): ||V||_F = {np.linalg.norm(V):.3e}")
    assert np.linalg.norm(V) == 0.0

    H_ref = np.diag(rydberg_diagonal(blk)) + h_mol.build(blk, R)
    H_off = H_ref + V
    ev_ref = np.sort(np.linalg.eigvalsh(H_ref))
    ev_off = np.sort(np.linalg.eigvalsh(H_off))
    print(f"  max|E(A=0) - E(H_a+H_mol)| = {np.abs(ev_off-ev_ref).max():.3e} E_h")
    print(f"  allclose(rtol=1e-12) = {np.allclose(ev_off, ev_ref, rtol=1e-12, atol=0.0)}")
    assert np.allclose(ev_off, ev_ref, rtol=1e-12, atol=0.0)

    # control positivo: con el remapeo ACTIVO el espectro SÍ cambia
    fermi_on = FermiPseudopotential(RADIAL, SCAT)
    ev_on = np.sort(np.linalg.eigvalsh(H_ref + fermi_on.build(blk, R)))
    d = np.abs(ev_on - ev_ref).max()
    print(f"  control positivo (remapeo activo): max|ΔE| = {d:.6e} E_h "
          f"({d*6.579683920502e15/1e9:.3f} GHz)")
    assert d > 1e-10


# ================================================================ TEST F3d
def test_f3d_magnitude_vs_charge_dipole():
    """Magnitud de V_Fermi frente a H_mol. Se reporta tal cual salga."""
    efield = RydbergElectronField(RADIAL)
    h_mol = ChargeDipoleHamiltonian(electron_field=efield)
    h_rot = ChargeDipoleHamiltonian()
    fermi = FermiPseudopotential(RADIAL, SCAT)
    blk = BASIS.get_block(0)
    E_ryd = rydberg_diagonal(blk)
    E_rot = h_rot.rotational_diagonal(blk)
    ref = np.sort(E_ryd + E_rot)

    print("\n      R    |  ||V_cd||_F   | ||V_Fermi||_F |  cociente  | "
          "ΔE_cd [GHz] | ΔE_Fermi [GHz]")
    print("   --------|---------------|---------------|------------|"
          "-------------|---------------")
    GHZ = 6.579683920502e15 / 1e9
    # R=1100 se cayó al completar la base: el par 26p-26p (n*=23.3512) tiene
    # su punto de retorno clásico en 1090.6 a0, más corto que el del 27s. Con
    # la base incompleta ese par no existía y el dominio llegaba a 1138.9.
    for R in (400.0, 500.0, 600.0, 800.0, 900.0, 1080.0):
        V_cd = h_mol.build(blk, R) - np.diag(E_rot)
        V_f = fermi.build(blk, R)
        H_cd = np.diag(E_ryd) + h_mol.build(blk, R)
        H_f = np.diag(E_ryd) + np.diag(E_rot) + V_f
        d_cd = np.abs(np.sort(np.linalg.eigvalsh(H_cd)) - ref).max()
        d_f = np.abs(np.sort(np.linalg.eigvalsh(H_f)) - ref).max()
        print(f"   {R:8.0f}| {np.linalg.norm(V_cd):.7e} | {np.linalg.norm(V_f):.7e} | "
              f"{np.linalg.norm(V_f)/np.linalg.norm(V_cd):10.2f} | "
              f"{d_cd*GHZ:11.4f} | {d_f*GHZ:13.4f}")


# ================================================================ TEST F4
def test_f4_p_resonance_window():
    """
    La ventana de exclusión de la resonancia de forma p: que el centro sea de
    verdad el polo, que la anchura sea la FWHM y que el margen se aplique.

    El polo se define como el CERO del interpolante de 1/A_p, que es justo el
    que usa `p_interpolation="inverse"`; el test lo comprueba evaluando A_p a
    ambos lados: debe cambiar de signo y ser enorme.
    """
    ns = n_star_of_l(3)
    w = SCAT.p_resonance_window(ns, margin_factor=2.0)
    print(f"\n   polo:  eps = {w['eps_pole']:.8f} E_h = "
          f"{w['eps_pole']*27211.386245988:.3f} meV  ->  R = {w['R_pole']:.3f} a0")
    print(f"   FWHM de |A_p|: R in [{w['R_fwhm_lo']:.3f}, {w['R_fwhm_hi']:.3f}] a0"
          f"  ->  {w['fwhm']:.3f} a0")
    print(f"   ventana (±2×FWHM): [{w['R_lo']:.2f}, {w['R_hi']:.2f}] a0")

    # el centro cae dentro de la FWHM
    assert w["R_fwhm_lo"] < w["R_pole"] < w["R_fwhm_hi"]
    # el margen se aplica tal cual, simétrico
    assert np.isclose(w["R_lo"], w["R_pole"] - 2.0 * w["fwhm"])
    assert np.isclose(w["R_hi"], w["R_pole"] + 2.0 * w["fwhm"])
    assert np.isclose(w["margin_factor"], 2.0)
    # escala en el sitio: valores de regresión de docs/analysis_resonancia_onda_p.md
    assert abs(w["R_pole"] - 562.77) < 0.5, w["R_pole"]
    assert abs(w["fwhm"] - 13.2) < 0.5, w["fwhm"]

    # a ambos lados del polo, A_p debe cambiar de signo y dispararse
    A_p_lo = SCAT.scattering(w["R_pole"] - 0.5, ns)[1]
    A_p_hi = SCAT.scattering(w["R_pole"] + 0.5, ns)[1]
    print(f"   A_p(R_polo-0.5) = {A_p_lo:+.4e} a0^3")
    print(f"   A_p(R_polo+0.5) = {A_p_hi:+.4e} a0^3")
    assert A_p_lo * A_p_hi < 0.0
    assert min(abs(A_p_lo), abs(A_p_hi)) > 1.0e6

    # fuera de la ventana A_p ya es "normal": órdenes de magnitud por debajo
    for R in (w["R_lo"] - 5.0, w["R_hi"] + 5.0):
        A_p = SCAT.scattering(R, ns)[1]
        print(f"   A_p(R={R:.1f}) = {A_p:+.4e} a0^3")
        assert abs(A_p) < 1.0e5

    # el margen escala como se anuncia
    w4 = SCAT.p_resonance_window(ns, margin_factor=4.0)
    assert np.isclose(w4["R_hi"] - w4["R_lo"], 2.0 * (w["R_hi"] - w["R_lo"]))
