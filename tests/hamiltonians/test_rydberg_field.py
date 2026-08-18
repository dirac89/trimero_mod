#!/usr/bin/env python3
"""
Tests analíticos del término del campo del ELECTRÓN Rydberg en F_ryd:

    -d·(e·(r-R)/|r-R|³)                                  (2º término de Ec. 4)

González-Férez, Sadeghpour & Schmelcher, NJP 17, 013021 (2015).

El término del ion Rb⁺ ya fue verificado en la ronda anterior
(test_charge_dipole.py) y aquí se reutiliza SIN redefinir su convenio.

La pieza central de esta suite es TEST 6: una implementación de REFERENCIA por
cuadratura 2D directa sobre (r, cosθ) del campo crudo (r-R)/|r-R|³, sin ninguna
expansión en armónicos esféricos, sin símbolos 3j y sin integrales radiales
separadas. Si la expansión multipolar de producción coincide con ella, el
álgebra de las Ec. A.6-A.10 es correcta.

Ejecutar:  poetry run python src/test_rydberg_field.py
"""

import time

import numpy as np
from scipy.special import gammaln, lpmv, roots_legendre

from angular_algebra import gaunt, wigner_3j
from charge_dipole import (
    B_KRB_AU,
    D_KRB_AU,
    ChargeDipoleHamiltonian,
    RydbergElectronField,
    rydberg_diagonal,
)
from quantum_basis import CoupledBasis
from rydberg_radial import RadialBasis



# ======================================================================
# REFERENCIA INDEPENDIENTE: cuadratura 2D del campo crudo
# ======================================================================
def _norm_legendre(l, m, x):
    """N_lm(x) tal que Y_lm(θ,φ) = N_lm(cosθ)·e^{imφ}  (fase Condon-Shortley)."""
    am = abs(m)
    c = np.sqrt((2 * l + 1) / (4 * np.pi) * np.exp(gammaln(l - am + 1) - gammaln(l + am + 1)))
    val = c * lpmv(am, l, x)
    if m < 0:
        val = val * (-1) ** am
    return val


def brute_force_field_element(radial, l1, m1, l2, m2, mu, R, n_x=800):
    """
    ⟨l1 m1| F_μ |l2 m2⟩ del campo del electrón, por cuadratura 2D directa.

        F_z = (r cosθ - R)/D³ ,  F_ρ = r sinθ/D³ ,  D = |r⃗-R⃗|
        F_{+1} = -(F_x + iF_y)/√2 = -F_ρ e^{iφ}/√2
        F_{-1} = +(F_x - iF_y)/√2 = +F_ρ e^{-iφ}/√2

    La integral en φ es analítica y da 2π·δ_{m1, m2+μ}. Lo que queda se
    integra numéricamente: Gauss-Legendre en x=cosθ, trapecio en r sobre la
    malla de `radial`. NO usa expansión multipolar ni símbolos 3j.
    """
    if m1 != m2 + mu:
        return 0.0
    r = radial.r
    f_r = radial.u(l1) * radial.u(l2)
    x, w = roots_legendre(n_x)
    D3 = (r[:, None] ** 2 + R**2 - 2.0 * R * np.outer(r, x)) ** 1.5
    ang = _norm_legendre(l1, m1, x) * _norm_legendre(l2, m2, x)
    if mu == 0:
        kern = (np.outer(r, x) - R) / D3
        pref = 1.0
    else:
        kern = np.outer(r, np.sqrt(1.0 - x * x)) / D3
        pref = -1.0 / np.sqrt(2.0) if mu == 1 else 1.0 / np.sqrt(2.0)
    inner = kern @ (ang * w)
    return pref * 2.0 * np.pi * np.trapezoid(f_r * inner, r)


# Objetos compartidos
RADIAL = RadialBasis()
EFIELD = RydbergElectronField(RADIAL)
BASIS = CoupledBasis(N_max=6, manifold_l_min=3)


# ======================================================================
# TEST 0a — álgebra angular contra valores analíticos conocidos
# ======================================================================
def test_0a_angular_algebra():
    casos = [
        ((0, 0, 0, 0, 0, 0), 1.0, "(000;000)"),
        ((1, 1, 1, 0, 0, 0), 0.0, "(111;000) suma impar"),
        ((1, 1, 0, 0, 0, 0), -1.0 / np.sqrt(3.0), "(110;000)"),
        ((1, 1, 2, 0, 0, 0), np.sqrt(2.0 / 15.0), "(112;000)"),
        ((1, 1, 0, 1, -1, 0), 1.0 / np.sqrt(3.0), "(110;1-10)"),
        ((2, 1, 1, 0, 0, 0), np.sqrt(2.0 / 15.0), "(211;000) cíclico de (112)"),
        ((2, 2, 2, 0, 0, 0), -np.sqrt(2.0 / 35.0), "(222;000)"),
    ]
    print("   símbolo 3j        |    calculado    |     esperado    |   diff")
    print("  -------------------|-----------------|-----------------|---------")
    for args, want, label in casos:
        got = wigner_3j(*args)
        print(f"   {label:18s}| {got:15.12f} | {want:15.12f} | {abs(got - want):.1e}")
        assert abs(got - want) < 1e-13, f"{label}: {got} != {want}"

    # fórmula cerrada de Racah para (j1 j2 j3; 0 0 0) con J=j1+j2+j3 par:
    #   (-1)^g sqrt((2g-2j1)!(2g-2j2)!(2g-2j3)!/(2g+1)!) · g!/((g-j1)!(g-j2)!(g-j3)!)
    from math import factorial
    print("\n   contra la fórmula cerrada de (j1 j2 j3; 0 0 0):")
    worst = 0.0
    n_cmp = 0
    for j1 in range(0, 12):
        for j2 in range(0, 12):
            for j3 in range(abs(j1 - j2), min(j1 + j2, 20) + 1):
                if (j1 + j2 + j3) % 2:
                    continue
                g = (j1 + j2 + j3) // 2
                want = ((-1) ** g) * np.sqrt(
                    factorial(2 * g - 2 * j1) * factorial(2 * g - 2 * j2)
                    * factorial(2 * g - 2 * j3) / factorial(2 * g + 1)
                ) * factorial(g) / (factorial(g - j1) * factorial(g - j2) * factorial(g - j3))
                got = wigner_3j(j1, j2, j3, 0, 0, 0)
                worst = max(worst, abs(got - want))
                n_cmp += 1
    print(f"     {n_cmp} tripletes comparados, max|diferencia| = {worst:.2e}")
    assert worst < 1e-13, worst

    # fórmula cerrada (l 1 l+1; 000) = (-1)^{l+1} sqrt((l+1)/((2l+1)(2l+3)))
    print("\n   (l 1 l+1; 000) contra fórmula cerrada:")
    for l in range(0, 12):
        got = wigner_3j(l, 1, l + 1, 0, 0, 0)
        want = ((-1.0) ** (l + 1)) * np.sqrt(
            (l + 1.0) / ((2.0 * l + 1.0) * (2.0 * l + 3.0))
        )
        assert abs(got - want) < 1e-13, f"l={l}: {got} != {want}"
    print(f"     l = 0..11 ✓  (ej. l=7: {wigner_3j(7,1,8,0,0,0):.12f})")

    # ortogonalidad: con m3 FIJO,  Σ_{m1} (2j3+1)·(j1 j2 j3; m1, -m1-m3, m3)² = 1
    print("\n   regla de suma Σ_{m1} (2j3+1)·3j² = 1  (m3 fijo):")
    for (j1, j2, j3, m3) in [(2, 3, 4, 0), (5, 4, 3, 2), (1, 1, 2, -1),
                             (6, 6, 6, 3), (23, 22, 45, 0), (23, 23, 46, -5)]:
        acc = 0.0
        for m1 in range(-j1, j1 + 1):
            acc += (2 * j3 + 1) * wigner_3j(j1, j2, j3, m1, -m1 - m3, m3) ** 2
        print(f"     ({j1},{j2},{j3}) m3={m3:3d}: {acc:.14f}   diff={abs(acc - 1.0):.1e}")
        assert abs(acc - 1.0) < 1e-11, (j1, j2, j3, m3, acc)

    # Gaunt: ∫Y*_lm Y_00 Y_lm dΩ = 1/sqrt(4π)
    print("\n   Gaunt ∫Y*_lm Y_00 Y_lm dΩ = 1/√(4π) = %.12f:" % (1 / np.sqrt(4 * np.pi)))
    for l, m in [(0, 0), (3, 2), (7, -5), (23, 23)]:
        got = gaunt(l, m, 0, 0, l, m)
        assert abs(got - 1 / np.sqrt(4 * np.pi)) < 1e-13, (l, m, got)
    print("     (l,m) = (0,0),(3,2),(7,-5),(23,23) ✓")


# ======================================================================
# TEST 0b — funciones radiales: normalización y ⟨r⟩ exacto
# ======================================================================
def test_0b_radial_basis():
    print(f"  malla: {RADIAL.n_points} puntos, r ∈ [0, {RADIAL.r_max:g}] a0 (uniforme en √r)")
    print("\n    l  |  n  |    <u|u>       |    <r> num.   |   <r> exacto  |   rel")
    print("   ----|-----|----------------|---------------|---------------|--------")
    for l in [0, 3, 5, 12, 23]:
        n = RADIAL.n_of_l(l)
        u = RADIAL.u(l)
        nrm = np.trapezoid(u * u, RADIAL.r)
        r_num = np.trapezoid(u * u * RADIAL.r, RADIAL.r)
        r_exact = (3.0 * n * n - l * (l + 1)) / 2.0
        rel = abs(r_num / r_exact - 1.0)
        print(f"   {l:4d}| {n:4d}| {nrm:.12f} | {r_num:13.6f} | {r_exact:13.6f} | {rel:.1e}")
        assert abs(nrm - 1.0) < 1e-10, f"l={l} no normalizada: {nrm}"
        assert rel < 1e-9, f"l={l}: <r> = {r_num} != {r_exact}"

    tail = abs(RADIAL.u(0)[-1]), abs(RADIAL.u(23)[-1])
    print(f"\n   |u| en r_max: l=0 -> {tail[0]:.2e},  l=23 -> {tail[1]:.2e}  (cola despreciable)")
    assert max(tail) < 1e-9


# ======================================================================
# TEST 0c — cierre con la ronda anterior: el elemento de cos(θ_d) obtenido
# vía Gaunt debe COINCIDIR con cos_theta_element ya verificado.
# ======================================================================
def test_0c_rotor_matches_round1():
    h = ChargeDipoleHamiltonian()
    worst = 0.0
    n = 0
    for N1 in range(0, 7):
        for M1 in range(-N1, N1 + 1):
            for N2 in range(0, 7):
                for M2 in range(-N2, N2 + 1):
                    via_gaunt = np.sqrt(4.0 * np.pi / 3.0) * gaunt(N1, M1, 1, 0, N2, M2)
                    ronda1 = h.cos_theta_element(N1, M1, N2, M2)
                    worst = max(worst, abs(via_gaunt - ronda1))
                    n += 1
    print(f"  ⟨N' M'|cosθ|N M⟩ vía Gaunt/3j  vs  cos_theta_element (ronda 1)")
    print(f"    {n} pares comparados,  max|diferencia| = {worst:.3e}")
    assert worst < 1e-14, worst
    print("    el convenio de la ronda 1 se reproduce desde el álgebra angular general ✓")


# ======================================================================
# TEST 6 — expansión multipolar vs CUADRATURA 2D BRUTA (validación clave)
# ======================================================================
def test_6_expansion_vs_brute_force():
    casos = [
        (5, 2, 5, 2, 0), (6, 2, 5, 2, 0), (23, 0, 22, 0, 0),
        (4, 1, 7, 1, 0), (6, 3, 5, 2, 1), (6, 1, 5, 2, -1),
        (12, -5, 11, -4, -1), (3, 0, 3, 0, 0), (0, 0, 1, 0, 0),
    ]
    for R in (2500.0, 1500.0, 900.0):
        print(f"\n  --- R = {R:g} a0 " + "-" * 46)
        print("   l1 m1  l2 m2  μ |   expansión      |   fuerza bruta   |   rel")
        print("   ---------------|------------------|------------------|--------")
        worst = 0.0
        for (l1, m1, l2, m2, mu) in casos:
            if l1 not in RADIAL.l_values or l2 not in RADIAL.l_values:
                continue
            exp_ = EFIELD.field_element(l1, m1, l2, m2, mu, R)
            bf = brute_force_field_element(RADIAL, l1, m1, l2, m2, mu, R)
            scale = max(abs(exp_), abs(bf))
            rel = abs(exp_ - bf) / scale if scale > 0 else 0.0
            worst = max(worst, rel)
            print(f"   {l1:2d}{m1:3d} {l2:3d}{m2:3d} {mu:2d} | {exp_:+.10e} | {bf:+.10e} | {rel:.1e}")
        tol = 1e-6 if R > 2000 else 1e-3
        print(f"   peor diferencia relativa = {worst:.2e}   (tolerancia {tol:.0e})")
        assert worst < tol, f"expansión != fuerza bruta a R={R}: {worst}"

    # A R=900 la densidad electrónica solapa con r≈R, donde el campo crudo
    # 1/D³ casi diverge y la cuadratura de fuerza bruta pierde precisión. Se
    # COMPRUEBA (no se afirma) que el residuo es de la cuadratura: al refinar
    # los nodos, la fuerza bruta converge HACIA el valor de la expansión.
    print("\n  --- convergencia de la fuerza bruta a R=900 (el residuo es suyo) ---")
    l1, m1, l2, m2, mu = 3, 0, 3, 0, 0
    exp_ = EFIELD.field_element(l1, m1, l2, m2, mu, 900.0)
    print(f"   expansión (analítica en el corte r≶R) = {exp_:.12e}")
    prev = None
    for n_x in (400, 800, 1600, 3200):
        bf = brute_force_field_element(RADIAL, l1, m1, l2, m2, mu, 900.0, n_x=n_x)
        rel = abs(bf / exp_ - 1.0)
        print(f"     n_x = {n_x:5d}:  bruta = {bf:.12e}   rel = {rel:.2e}")
        if prev is not None:
            assert rel < prev, "la fuerza bruta no converge hacia la expansión"
        prev = rel
    print("   converge monótonamente hacia la expansión ✓")


# ======================================================================
# TEST 3 — conservación de M_J con ΔM_N ≠ 0
# ======================================================================
def test_3_mj_conserved_mn_not():
    h = ChargeDipoleHamiltonian(electron_field=EFIELD)
    R = 1500.0

    # ΔM_N = +1, Δm_l = -1  =>  ΔM_J = 0.  Debe ser NO NULO.
    si = (5, 2, 3, -2)     # M_J = 0
    sj = (6, 1, 2, -1)     # M_J = 0;  Δm_l = -1, ΔM_N = +1  =>  ΔM_J = 0
    v = h.matrix_element(si, sj, R)
    print(f"  ΔM_N ≠ 0 con ΔM_J = 0:")
    print(f"    i = {si}  (m_l={si[1]}, M_N={si[3]}, M_J={si[1]+si[3]})")
    print(f"    j = {sj}  (m_l={sj[1]}, M_N={sj[3]}, M_J={sj[1]+sj[3]})")
    print(f"    Δm_l = {sj[1]-si[1]:+d}, ΔM_N = {sj[3]-si[3]:+d}, ΔM_J = 0")
    print(f"    ⟨i|H_mol|j⟩ = {v:+.9e}  (debe ser ≠ 0)")
    assert v != 0.0, "el término del electrón no acopla ΔM_N≠0: expansión mal"

    # el término del ion SOLO no puede producir esto
    h_ion = ChargeDipoleHamiltonian()
    v_ion = h_ion.matrix_element(si, sj, R)
    print(f"    sólo ion: {v_ion:+.9e}  (debe ser 0: el ion no mezcla M_N) ✓")
    assert v_ion == 0.0

    # ΔN debe seguir siendo ±1 (d⃗ es rango 1 sobre el rotor)
    sk = (6, 1, 5, -1)   # M_J=0, ΔN = +2
    print(f"\n    ΔN=±2: j={sk} -> ⟨i|H|j⟩ = {h.matrix_element(si, sk, R):+.3e} (debe ser 0)")
    assert h.matrix_element(si, sk, R) == 0.0

    # M_J distinto: sigue lanzando, y el cálculo crudo sigue dando cero exacto
    bad = (6, 1, 2, 0)   # M_J = 1
    try:
        h.matrix_element(si, bad, R)
    except ValueError as exc:
        print(f"\n    M_J distinto -> ValueError ✓")
    else:
        raise AssertionError("no detectó M_J distinto")

    pares, worst = 0, 0.0
    for MJa, MJb in [(0, 1), (0, -2), (2, 5)]:
        for xa in BASIS.get_block(MJa).states[:40]:
            for yb in BASIS.get_block(MJb).states[:40]:
                v = h._matrix_element_unchecked(xa, yb, R)
                worst = max(worst, abs(v))
                pares += 1
    print(f"    {pares} pares con M_J distinto SIN comprobación: max|elem| = {worst:.1e} ✓")
    assert worst == 0.0


# ======================================================================
# TEST 2 — hermiticidad con AMBOS términos
# ======================================================================
def test_2_hermiticity_both_terms():
    h = ChargeDipoleHamiltonian(electron_field=EFIELD)
    R = 1500.0

    small_MJ = max(
        (mj for mj in BASIS.block_M_J_values() if len(BASIS.get_block(mj)) <= 200),
        key=lambda mj: len(BASIS.get_block(mj)),
    )
    blk = BASIS.get_block(small_MJ)
    print(f"  Barrido completo: bloque M_J={small_MJ}, dim={len(blk)} ({len(blk)**2} pares)")
    worst = 0.0
    scale = 0.0
    for x in blk.states:
        for y in blk.states:
            a = h.matrix_element(x, y, R)
            b = h.matrix_element(y, x, R)
            worst = max(worst, abs(a - b))
            scale = max(scale, abs(a))
    print(f"    max|H_ij - H_ji| = {worst:.3e}   (max|H_ij| = {scale:.3e})")
    assert worst <= 1e-16 * max(scale, 1e-30), f"asimetría {worst} frente a escala {scale}"

    # `build` sólo visita los pares permitidos por las reglas de selección:
    # se comprueba que no pierde nada frente al barrido dim² completo.
    Hb = h.build(blk, R)
    Hr = h.build_reference(blk, R)
    print(f"    build vs build_reference (dim {len(blk)}): "
          f"np.array_equal = {np.array_equal(Hb, Hr)}, "
          f"max|dif| = {np.abs(Hb - Hr).max():.1e}")
    assert np.array_equal(Hb, Hr), "build pierde elementos frente al barrido completo"

    blk0 = BASIS.get_block(0)
    t0 = time.perf_counter()
    H = h.build(blk0, R)
    print(f"\n  build(M_J=0, R={R:g}) -> {H.shape} en {time.perf_counter()-t0:.2f} s")
    asym = np.linalg.norm(H - H.T)
    print(f"    ||H||_F = {np.linalg.norm(H):.9e}")
    print(f"    ||H - H.T||_F = {asym:.3e}   (relativo {asym/np.linalg.norm(H):.1e})")
    assert asym <= 1e-14 * np.linalg.norm(H), asym

    off = np.abs(H - np.diag(np.diag(H))).max()
    print(f"    max|fuera de diagonal| = {off:.6e}  (no vacuo)")
    assert off > 0.0
    ev = np.linalg.eigvalsh(H)
    print(f"    autovalores reales: min={ev.min():.9e} max={ev.max():.9e}")


# ======================================================================
# TEST 4 — límite dipolar: expansión completa vs truncada a 1er orden en r/R
# ======================================================================
def test_4_first_order_dipole_limit():
    print("  Campo del electrón a 1er orden en r/R (r < R en todo el soporte):")
    print("    F_elec ≈ -ẑ/R²  -  (2cosγ ẑ - sinγ ρ̂)·r/R³")
    print("    = monopolo (carga -e en el origen) + dipolo.")
    print("    Corresponde EXACTAMENTE a truncar la expansión multipolar en k ≤ 1.")
    print("\n   elemento          |     R     |   completa     |   k<=1        |   rel")
    print("  -------------------|-----------|----------------|---------------|--------")
    casos = [(5, 2, 5, 2, 0), (6, 2, 5, 2, 0), (6, 3, 5, 2, 1)]
    for (l1, m1, l2, m2, mu) in casos:
        prev = None
        for R in (1500.0, 3000.0, 6000.0, 12000.0):
            full = EFIELD.field_element(l1, m1, l2, m2, mu, R)
            appr = EFIELD.field_element(l1, m1, l2, m2, mu, R, k_max=1)
            rel = abs(full - appr) / abs(full) if full != 0 else 0.0
            print(f"   ({l1},{m1})<-({l2},{m2}) μ={mu:2d} | {R:9.0f} | {full:+.7e} | {appr:+.7e} | {rel:.2e}")
            if prev is not None:
                assert rel < prev, (
                    f"la diferencia relativa no decrece con R: {rel} >= {prev}"
                )
            prev = rel
        print("  " + "-" * 70)
    print("  La diferencia relativa decrece monótonamente con R en los 3 casos ✓")


# ======================================================================
# TEST 5 — consistencia con el término del ion: neutralidad del átomo Rydberg
# ======================================================================
def test_5_monopole_cancels_ion():
    print("  El término k=0 del campo del electrón vale exactamente -1/R² cuando")
    print("  toda la densidad electrónica está dentro de R, y CANCELA el campo del")
    print("  ion (+1/R²): el átomo Rydberg es neutro. Esto liga este término con el")
    print("  de la ronda anterior sin necesidad de un test artificial.\n")
    print("      R    |  F_ion=1/R²  | ⟨lm|F_0^elec|lm⟩ |   suma        | |suma|/F_ion | ratio·R²")
    print("   --------|--------------|------------------|---------------|-------------|----------")
    l, m = 5, 2
    ratios, scaled = [], []
    for R in (3000.0, 6000.0, 12000.0, 24000.0):
        f_ion = 1.0 / R**2
        f_el = EFIELD.field_element(l, m, l, m, 0, R)
        tot = f_ion + f_el
        ratio = abs(tot) / f_ion
        ratios.append(ratio)
        scaled.append(ratio * R**2)
        print(f"   {R:8.0f}| {f_ion:.6e} | {f_el:+.9e} | {tot:+.6e} | {ratio:11.4e} | {ratio*R**2:.5e}")
    assert all(ratios[i + 1] < ratios[i] for i in range(len(ratios) - 1)), ratios
    print("\n   El residuo NO es cero: el monopolo se cancela exactamente, pero queda")
    print("   el CUADRUPOLO del electrón, que cae como 1/R⁴. Por eso el cociente")
    print("   residuo/ion cae como 1/R², y ratio·R² es constante:")
    drift = abs(scaled[-1] / scaled[-2] - 1.0)
    print(f"     ratio·R² : {scaled[-2]:.6e} -> {scaled[-1]:.6e}   deriva = {drift:.2e}")
    assert drift < 0.02, f"el residuo no escala como 1/R⁴: deriva {drift}"

    # sólo el término k=0 debe valer exactamente -1/R²
    R = 6000.0
    k0 = EFIELD.field_element(l, m, l, m, 0, R, k_max=0)
    print(f"\n   sólo k=0 a R={R:g}: {k0:.12e}   vs  -1/R² = {-1/R**2:.12e}")
    rel_k0 = abs(k0 + 1.0 / R**2) * R**2
    print(f"     diferencia relativa = {rel_k0:.2e}")
    assert rel_k0 < 1e-9, (k0, -1 / R**2)

    # regresión exacta: sin electron_field, la matriz es BIT A BIT la de la ronda 1
    blk = BASIS.get_block(7)
    H_ion_now = ChargeDipoleHamiltonian().build(blk, 1500.0)
    H_ref = np.zeros_like(H_ion_now)
    h_ref = ChargeDipoleHamiltonian()
    for i, si in enumerate(blk.states):
        for j, sj in enumerate(blk.states):
            H_ref[i, j] = (
                h_ref.rotational_element(si[2], si[3], sj[2], sj[3])
                if (si[0] == sj[0] and si[1] == sj[1]) else 0.0
            ) + h_ref.ion_field_element(si, sj, 1500.0)
    print(f"\n   regresión ion-solo (bloque M_J=7, dim {len(blk)}):")
    print(f"     np.array_equal(H_nuevo, H_ronda1) = {np.array_equal(H_ion_now, H_ref)}")
    assert np.array_equal(H_ion_now, H_ref)


# ======================================================================
# TEST 1 — límite R→∞ con AMBOS términos
# ======================================================================
def test_1_large_R_both_terms():
    h = ChargeDipoleHamiltonian(electron_field=EFIELD)
    blk = BASIS.get_block(0)
    E_ryd = rydberg_diagonal(blk)
    ref = np.sort(E_ryd + h.rotational_diagonal(blk))
    print(f"  Bloque M_J=0, dim = {len(blk)}")

    R_big = 20000.0
    H = np.diag(E_ryd) + h.build(blk, R_big)
    off = np.abs(H - np.diag(np.diag(H))).max()
    print(f"\n   R = {R_big:g} a0:  max|fuera de diagonal| = {off:.6e} E_h  (no vacuo)")
    assert off > 0.0

    ev = np.sort(np.linalg.eigvalsh(H))
    dmax = np.abs(ev - ref).max()
    print(f"   max|autovalor - (E_Ryd + B·N(N+1))| = {dmax:.6e} E_h   tol = 1e-10")
    for k in range(5):
        print(f"     {k}: {ev[k]:.15e}  vs {ref[k]:.15e}   Δ={ev[k]-ref[k]:+.3e}")
    assert dmax < 1e-10, dmax

    R_mid = 1500.0
    Hm = np.diag(E_ryd) + h.build(blk, R_mid)
    dmid = np.abs(np.sort(np.linalg.eigvalsh(Hm)) - ref).max()
    print(f"\n   Control positivo R={R_mid:g}: max|E - desacoplado| = {dmid:.6e} E_h "
          f"({dmid*6.579683920502e15/1e6:.3f} MHz)")
    assert dmid > 1e-10


# ======================================================================
