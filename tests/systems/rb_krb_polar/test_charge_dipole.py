#!/usr/bin/env python3
"""
Tests analíticos del Hamiltoniano de acoplamiento carga-dipolo

    H_mol = B·N² - d·F_ryd(R,r),    F_ryd = e·R/R³ + e·(r-R)/|r-R|³

Base científica: González-Férez, Sadeghpour & Schmelcher,
                 New J. Phys. 17, 013021 (2015), Ec. 3-4.

ALCANCE DE ESTA RONDA: sólo el término B·N² y el término del campo del ion
Rb⁺ (e·R/R³). El término del campo del electrón Rydberg (Ec. A.6-A.10) NO
está implementado todavía y estos tests NO lo cubren.

Estos tests son condiciones NECESARIAS que cualquier implementación correcta
debe cumplir. Se escribieron y se vieron fallar ANTES de la implementación.

Ejecutar:  poetry run python src/test_charge_dipole.py    (desde raíz)
       o:  poetry run python test_charge_dipole.py         (desde src/)
"""

import time

import numpy as np

from trimero.basis.quantum import CoupledBasis
from trimero.systems.rb_krb_polar.charge_dipole import (
    B_KRB_AU,
    D_KRB_AU,
    ChargeDipoleHamiltonian,
    rydberg_diagonal,
)


# ----------------------------------------------------------------------


# Base compartida por todos los tests (construirla es barato, ~0.1 s)
BASIS = CoupledBasis(N_max=6, manifold_l_min=3)


# ======================================================================
# TEST 3 (del enunciado) — se ejecuta PRIMERO porque es la verificación
# "a mano" que debe pasar ANTES de integrarse en la matriz completa.
#   ⟨N,M_N| B·N² |N,M_N⟩ = B·N(N+1)
# ======================================================================
def test_3_rotational_diagonal_element():
    h = ChargeDipoleHamiltonian()

    print(f"  B(KRb) = {B_KRB_AU:.9e} E_h   ({B_KRB_AU * 6.579683920502e15 / 1e9:.4f} GHz)")
    print(f"  d(KRb) = {D_KRB_AU:.9e} e·a0")
    print()

    # --- caso de mano, par (N,M_N) arbitrario ---------------------------
    N, M_N = 4, -2
    got = h.rotational_element(N, M_N, N, M_N)
    want = B_KRB_AU * N * (N + 1)
    print(f"  Caso de mano: (N,M_N) = ({N},{M_N})")
    print(f"    ⟨N M_N|B·N²|N M_N⟩ calculado = {got:.15e}")
    print(f"    B·N(N+1) = B·{N}·{N + 1} = {want:.15e}")
    print(f"    diferencia absoluta          = {abs(got - want):.3e}")
    assert got == want, f"{got} != {want}"

    # --- todos los (N,M_N) de la base -----------------------------------
    print("\n    N | M_N |   B·N(N+1) [E_h]   | ok")
    print("   ---|-----|--------------------|----")
    for N in range(0, BASIS.N_max + 1):
        for M_N in range(-N, N + 1):
            got = h.rotational_element(N, M_N, N, M_N)
            want = B_KRB_AU * N * (N + 1)
            assert got == want, f"(N={N},M_N={M_N}): {got} != {want}"
        print(f"   {N:3d} | all | {B_KRB_AU * N * (N + 1):.12e} | ✓")

    # --- N² es diagonal: nada fuera de la diagonal ----------------------
    n_off = 0
    for N1 in range(0, BASIS.N_max + 1):
        for M1 in range(-N1, N1 + 1):
            for N2 in range(0, BASIS.N_max + 1):
                for M2 in range(-N2, N2 + 1):
                    if (N1, M1) == (N2, M2):
                        continue
                    v = h.rotational_element(N1, M1, N2, M2)
                    assert v == 0.0, f"⟨{N1},{M1}|B N²|{N2},{M2}⟩ = {v} != 0"
                    n_off += 1
    print(f"\n  B·N² estrictamente diagonal: {n_off} elementos fuera de diagonal, todos == 0.0 ✓")


# ======================================================================
# TEST 4 (del enunciado) — regla de selección M_J.
# El acoplamiento carga-dipolo NO mezcla M_J distintos.
#   (a) matrix_element() DETECTA el error y lanza ValueError.
#   (b) el cálculo crudo (sin la comprobación) da EXACTAMENTE 0.0.
# ======================================================================
def test_4_mj_selection_rule():
    h = ChargeDipoleHamiltonian()
    R = 1500.0

    # --- (a) el código detecta el error y lanza -------------------------
    si = (5, 3, 3, -3)   # M_J = 3 + (-3) = 0
    sj = (5, 3, 3, -2)   # M_J = 3 + (-2) = 1
    print(f"  (a) estado i = {si}  ->  M_J = {si[1] + si[3]}")
    print(f"      estado j = {sj}  ->  M_J = {sj[1] + sj[3]}")
    try:
        h.matrix_element(si, sj, R)
    except ValueError as exc:
        print(f"      matrix_element() lanzó ValueError: {exc} ✓")
    else:
        raise AssertionError("matrix_element() NO detectó el cruce de M_J")

    # simétrico
    try:
        h.matrix_element(sj, si, R)
    except ValueError:
        print("      matrix_element(j,i) también lanza ✓")
    else:
        raise AssertionError("matrix_element(j,i) NO detectó el cruce de M_J")

    # el mismo par pero con M_J igual NO debe lanzar
    ok = h.matrix_element(si, si, R)
    print(f"      matrix_element(i,i) (M_J igual) no lanza, vale {ok:.9e} ✓")

    # --- (b) verificación explícita de que el valor crudo es cero -------
    # NO se asume "cero por construcción": se calcula sin la comprobación.
    pairs = 0
    worst = 0.0
    for MJ_a, MJ_b in [(0, 1), (0, 2), (0, -1), (3, -4), (0, 7)]:
        sa = BASIS.get_block(MJ_a).states[:60]
        sb = BASIS.get_block(MJ_b).states[:60]
        for x in sa:
            for y in sb:
                assert x[1] + x[3] != y[1] + y[3]
                v = h._matrix_element_unchecked(x, y, R)
                worst = max(worst, abs(v))
                assert v == 0.0, f"elemento no nulo entre M_J distintos: {x} {y} -> {v}"
                pairs += 1
    print(f"\n  (b) {pairs} pares (l,m_l,N,M_N) con M_J distinto evaluados SIN la comprobación")
    print(f"      max |⟨i|H_cd|j⟩| = {worst:.1e}  (exigido exactamente 0.0) ✓")

    # control positivo: dentro de un mismo M_J el operador NO es nulo
    blk = BASIS.get_block(0)
    nz = 0
    for x in blk.states[:400]:
        for y in blk.states[:400]:
            if h._matrix_element_unchecked(x, y, R) != 0.0:
                nz += 1
    print(f"      control positivo: {nz} elementos NO nulos dentro del bloque M_J=0 ✓")
    assert nz > 0, "el operador es idénticamente nulo: el test (b) sería vacuo"


# ======================================================================
# TEST 2 (del enunciado) — hermiticidad / simetría real.
# ======================================================================
def test_2_hermiticity():
    h = ChargeDipoleHamiltonian()
    R = 1500.0

    # --- barrido COMPLETO O(dim²) elemento a elemento en un bloque chico
    small_MJ = max(
        mj for mj in BASIS.block_M_J_values() if len(BASIS.get_block(mj)) <= 200
    )
    # ese `max` da un bloque pequeño; buscamos el mayor bloque <= 200 estados
    small_MJ = max(
        (mj for mj in BASIS.block_M_J_values() if len(BASIS.get_block(mj)) <= 200),
        key=lambda mj: len(BASIS.get_block(mj)),
    )
    blk = BASIS.get_block(small_MJ)
    print(f"  Barrido completo elemento a elemento: bloque M_J={small_MJ}, dim={len(blk)}")
    worst = 0.0
    for x in blk.states:
        for y in blk.states:
            a = h.matrix_element(x, y, R)
            b = h.matrix_element(y, x, R)
            worst = max(worst, abs(a - b))
    print(f"    max |⟨i|H|j⟩ - ⟨j|H|i⟩| = {worst:.3e}  sobre {len(blk)**2} pares")
    assert worst == 0.0, f"matrix_element no es simétrico: {worst}"

    # --- matriz completa del bloque M_J=0 -------------------------------
    blk0 = BASIS.get_block(0)
    t0 = time.perf_counter()
    H = h.build(blk0, R)
    print(f"\n  build(M_J=0, R={R}) -> {H.shape}  en {time.perf_counter() - t0:.2f} s")
    assert H.dtype == np.float64, H.dtype
    asym = np.linalg.norm(H - H.T)
    print(f"    ||H||_F        = {np.linalg.norm(H):.9e}")
    print(f"    ||H - H.T||_F  = {asym:.3e}")
    print(f"    np.array_equal(H, H.T) = {np.array_equal(H, H.T)}")
    assert asym == 0.0, f"H no es simétrica: ||H-H.T|| = {asym}"

    # no vacuo: debe haber acoplamiento fuera de la diagonal
    offdiag = H - np.diag(np.diag(H))
    print(f"    max |elemento fuera de diagonal| = {np.abs(offdiag).max():.6e}")
    assert np.abs(offdiag).max() > 0.0, "H es diagonal: el test de simetría sería vacuo"

    # autovalores reales (eigvalsh sobre matriz real simétrica)
    ev = np.linalg.eigvalsh(H)
    print(f"    autovalores de H_mol(R=1500): min={ev.min():.9e}  max={ev.max():.9e}")
    assert np.all(np.isreal(ev))


# ======================================================================
# TEST 1 (del enunciado) — límite R -> ∞.
# El acoplamiento -d·F_ion cae como 1/R²; a R muy grande el espectro debe
# recuperar E_Rydberg(l) + B·N(N+1) SIN mezcla.
# ======================================================================
def test_1_large_R_limit():
    h = ChargeDipoleHamiltonian()
    blk = BASIS.get_block(0)
    print(f"  Bloque M_J=0, dim = {len(blk)}")

    E_ryd = rydberg_diagonal(blk)
    E_rot = h.rotational_diagonal(blk)
    ref = np.sort(E_ryd + E_rot)
    print(f"    E_Rydberg: {len(np.unique(E_ryd))} valores distintos, "
          f"[{E_ryd.min():.9e}, {E_ryd.max():.9e}] E_h")
    print(f"    E_rot    : {len(np.unique(E_rot))} valores distintos, "
          f"[{E_rot.min():.9e}, {E_rot.max():.9e}] E_h")

    R_big = 10000.0
    H = np.diag(E_ryd) + h.build(blk, R_big)

    # el término de acoplamiento debe seguir estando presente (no vacuo)
    offdiag = np.abs(H - np.diag(np.diag(H))).max()
    print(f"\n    R = {R_big:g} a0:  max|off-diag| = {offdiag:.6e} E_h  "
          f"(d/R² = {D_KRB_AU / R_big**2:.3e})")
    assert offdiag > 0.0, "sin acoplamiento a R grande: el test sería vacuo"

    ev = np.sort(np.linalg.eigvalsh(H))
    dmax = np.abs(ev - ref).max()
    print(f"    max |autovalor - (E_Ryd + B N(N+1))| = {dmax:.6e} E_h")
    print(f"    tolerancia                            = 1.000000e-10 E_h")
    print("\n    primeros 6 autovalores vs referencia desacoplada:")
    for k in range(6):
        print(f"      {k}: {ev[k]:.15e}   {ref[k]:.15e}   Δ={ev[k] - ref[k]:+.3e}")
    assert dmax < 1e-10, f"el límite R->∞ no se recupera: Δmax = {dmax}"

    # ------- CONTROL POSITIVO -------------------------------------------
    # A R intermedio el espectro SÍ debe apartarse del desacoplado; si no,
    # el test anterior pasaría también con un acoplamiento idénticamente
    # nulo (bug silencioso).
    R_mid = 1500.0
    Hm = np.diag(E_ryd) + h.build(blk, R_mid)
    evm = np.sort(np.linalg.eigvalsh(Hm))
    dmid = np.abs(evm - ref).max()
    print(f"\n    Control positivo R = {R_mid:g} a0:")
    print(f"      max |autovalor - desacoplado| = {dmid:.6e} E_h  "
          f"({dmid * 6.579683920502e15 / 1e6:.3f} MHz)")
    assert dmid > 1e-10, (
        "a R=1500 el espectro coincide con el desacoplado: el acoplamiento "
        "no está haciendo nada y el test del límite R->∞ es vacuo"
    )


# ======================================================================
