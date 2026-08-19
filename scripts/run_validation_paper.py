#!/usr/bin/env python3
"""
Validación contra González-Férez 2015 con delta0(ns) = 3.13180 (override local).

  2) Convergencia en N: bloque M_J=0 con N_max=6 vs N_max=8.
  3) <cos theta_d> del autoestado más bajo derivado del manifold n=24, vs R.

Ejecutar:  poetry run python scripts/run_validation_paper.py
"""
import time
import numpy as np

from trimero.basis.quantum import CoupledBasis
from trimero.basis.radial import RadialBasis
from trimero.hamiltonians.charge_dipole import (
    HZ_PER_HARTREE, ChargeDipoleHamiltonian, RydbergElectronField, rydberg_diagonal,
)
from trimero.hamiltonians.fermi_krb import (
    FermiPseudopotential, ScatteringLengths, n_star_of_l,
)
from trimero.systems.rb_defects import DELTA0_NS_PAPER

GHZ = HZ_PER_HARTREE / 1.0e9
D0 = DELTA0_NS_PAPER
E_MANIFOLD = -0.5 / 24.0**2

RADIAL = RadialBasis()
SCAT = ScatteringLengths()
EFIELD = RydbergElectronField(RADIAL)
HMOL = ChargeDipoleHamiltonian(electron_field=EFIELD)
FERMI = FermiPseudopotential(RADIAL, SCAT, delta0_ns=D0)


def total_H(block, R):
    return (np.diag(rydberg_diagonal(block, delta0_ns=D0))
            + HMOL.build(block, R) + FERMI.build(block, R))


def cos_theta_matrix(block):
    """<i|cos(theta_d)|j> = delta_ll' delta_mm' <N M|cos|N' M'>. Independiente de R."""
    states = block.states
    index = {st: i for i, st in enumerate(states)}
    C = np.zeros((len(states), len(states)))
    for i, (l, m, N, MN) in enumerate(states):
        for Np in (N - 1, N + 1):
            if Np < 0 or abs(MN) > Np:
                continue
            j = index.get((l, m, Np, MN))
            if j is not None:
                C[i, j] = HMOL.cos_theta_element(N, MN, Np, MN)
    return C


def main():
    print("=" * 88)
    print(f"VALIDACIÓN vs González-Férez 2015   —   delta0(ns) = {D0}")
    print("=" * 88)

    # ---------------- 2) convergencia en N ----------------
    print("\n2) CONVERGENCIA EN N  (bloque M_J=0)")
    print("-" * 88)
    b6 = CoupledBasis(N_max=6, manifold_l_min=3).get_block(0)
    b8 = CoupledBasis(N_max=8, manifold_l_min=3).get_block(0)
    print(f"   dim(N_max=6) = {len(b6)}    dim(N_max=8) = {len(b8)}")
    for R in (600.0, 900.0):
        t0 = time.perf_counter()
        e6 = np.sort(np.linalg.eigvalsh(total_H(b6, R)))
        e8 = np.sort(np.linalg.eigvalsh(total_H(b8, R)))
        print(f"\n   R = {R:g} a0   ({time.perf_counter()-t0:.1f} s)")
        print("     k |     E(N=6) [E_h]     |     E(N=8) [E_h]     |  |ΔE| [MHz] |"
              "  |ΔE|/|E|   | |ΔE|/|E-E_man|")
        print("    ---|----------------------|----------------------|-------------|"
              "-------------|---------------")
        worst_abs = worst_rel = worst_bind = 0.0
        for k in range(10):
            d = abs(e8[k] - e6[k])
            rel = d / abs(e6[k])
            bind = abs(e6[k] - E_MANIFOLD)
            relb = d / bind if bind > 0 else np.nan
            worst_abs = max(worst_abs, d * GHZ * 1e3)
            worst_rel = max(worst_rel, rel)
            worst_bind = max(worst_bind, relb)
            print(f"    {k:3d}| {e6[k]:.15e} | {e8[k]:.15e} | {d*GHZ*1e3:11.5f} | "
                  f"{rel:11.4e} | {relb:13.4e}")
        print(f"     peor sobre los 10 primeros: {worst_abs:.5f} MHz, "
              f"rel(E) = {worst_rel:.3e}, rel(ligadura) = {worst_bind:.3e}")
        print(f"     criterio del paper: rel < 2e-6  ->  "
              f"{'CUMPLE' if worst_rel < 2e-6 else 'NO CUMPLE'} en rel(E)")

    # ---------------- 3) <cos theta_d> ----------------
    print("\n\n3) ORIENTACIÓN <cos theta_d>  —  autoestado más bajo del manifold n=24")
    print("-" * 88)
    C = cos_theta_matrix(b6)
    l_of = np.array([st[0] for st in b6.states])
    is_manifold = l_of >= 3
    print(f"   operador cos(theta_d) sobre M_J=0: ||C||_F = {np.linalg.norm(C):.4f}, "
          f"simétrico = {np.array_equal(C, C.T)}")
    # El dominio del remapeo lo fija el par MÁS ligado presente en la matriz
    # (media 27s+manifold), no el manifold solo: su punto de retorno es más corto.
    E_pair = 0.5 * (-0.5 / n_star_of_l(0, D0) ** 2 + -0.5 / n_star_of_l(3, D0) ** 2)
    R_turn = -1.0 / E_pair
    print(f"   dominio: punto de retorno del par 27s-manifold = {R_turn:.1f} a0")

    print("\n      R    |  E [E_h]              | E-E_man [GHz] | peso manif. | "
          "<cos> (más bajo) | <cos> (adiabático)")
    print("   --------|-----------------------|---------------|-------------|"
          "-----------------|-------------------")
    Rs = np.concatenate([np.arange(110.0, 400.0, 30.0), np.arange(400.0, 1140.0, 40.0)])
    v_prev = None
    rows, skipped = [], []
    for R in Rs:
        try:
            H = total_H(b6, R)
        except ValueError as exc:
            skipped.append((R, str(exc).split(":")[0]))
            continue
        w, V = np.linalg.eigh(H)
        # (a) más bajo con carácter de manifold
        pick = None
        for idx in np.argsort(w):
            if float(np.sum(V[:, idx][is_manifold] ** 2)) > 0.5:
                pick = idx
                break
        v_low = V[:, pick]
        cos_low = float(v_low @ C @ v_low)
        # (b) seguimiento adiabático por solapamiento máximo con el paso anterior
        if v_prev is None:
            idx_ad = pick
        else:
            idx_ad = int(np.argmax(np.abs(v_prev @ V)))
        v_ad = V[:, idx_ad]
        cos_ad = float(v_ad @ C @ v_ad)
        v_prev = v_ad
        wm = float(np.sum(v_low[is_manifold] ** 2))
        rows.append((R, w[pick], wm, cos_low, w[idx_ad], cos_ad))
        print(f"   {R:8.1f}| {w[pick]:+.15e} | {(w[pick]-E_MANIFOLD)*GHZ:13.4f} | "
              f"{wm:11.4f} | {cos_low:+16.6f} | {cos_ad:+18.6f}")
    for R, why in skipped:
        print(f"   {R:8.1f}|  --- fuera de dominio: {why}")

    smooth = [r for r in rows if r[0] <= 560.0]
    m_low = max(rows, key=lambda r: abs(r[3]))
    m_ad = max(rows, key=lambda r: abs(r[5]))
    m_sm = max(smooth, key=lambda r: abs(r[3]))
    print(f"\n   máx |<cos>| rama suave (R<=560, monótona): {abs(m_sm[3]):.6f} en R={m_sm[0]:.0f} a0")
    print(f"   máx |<cos>| criterio 'más bajo':            {abs(m_low[3]):.6f} en R={m_low[0]:.0f} a0")
    print(f"   máx |<cos>| seguimiento adiabático:         {abs(m_ad[5]):.6f} en R={m_ad[0]:.0f} a0")
    cruce = [r for r in smooth if abs(r[3]) < 0.78]
    if cruce:
        print(f"   la rama suave cruza 0.78 entre R={smooth[smooth.index(cruce[0])-1][0]:.0f} "
              f"y R={cruce[0][0]:.0f} a0")
    print(f"   paper (M_J=0): 0.78")
    print("=" * 88)


if __name__ == "__main__":
    main()
