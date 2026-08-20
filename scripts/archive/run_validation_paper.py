#!/usr/bin/env python3
"""
Validación contra González-Férez 2015 con delta0(ns) = 3.13180 (override local).

  2) Convergencia en N: bloque M_J=0 con N_max=6 vs N_max=8.
  3) <cos theta_d> del autoestado más bajo derivado del manifold, vs R.

⚠️ SUPERADO por `scripts/run_basis_correction_check.py`, que hace lo mismo con la
base CORRECTA (manifold + (n+1)d + (n+2)p + (n+3)s) y además compara contra la
base incompleta que se usaba antes. Este script se conserva porque su salida está
citada en docs/analysis_verificacion_tabla_I.md §10-11; ahora usa `BOPSystem`, así
que sus números son ya los de la base completa.

Ejecutar:  poetry run python scripts/run_validation_paper.py
"""
import time
import numpy as np

from trimero.systems.rb_krb_polar.bop_system import GHZ_PER_HARTREE as GHZ, BOPSystem
from trimero.systems.rb_krb_polar.rb_defects import DELTA0_NS_PAPER

D0 = DELTA0_NS_PAPER
N_MANIFOLD = 24                      # el manifold del paper; el resto se deriva
SYS6 = BOPSystem(n_manifold=N_MANIFOLD, delta0_ns=D0, N_max=6)
SYS8 = BOPSystem(n_manifold=N_MANIFOLD, delta0_ns=D0, N_max=8)
E_MANIFOLD = SYS6.E_manifold
HMOL = SYS6.hmol


def total_H(system, block, R):
    return system.hamiltonian(R)


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
    b6, b8 = SYS6.block(0), SYS8.block(0)
    print(f"   dim(N_max=6) = {len(b6)}    dim(N_max=8) = {len(b8)}")
    for R in (600.0, 900.0):
        t0 = time.perf_counter()
        e6 = np.sort(SYS6.eigvals(R))
        e8 = np.sort(SYS8.eigvals(R))
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
    print(f"\n\n3) ORIENTACIÓN <cos theta_d>  —  más bajo del manifold n={N_MANIFOLD}")
    print("-" * 88)
    C = cos_theta_matrix(b6)
    is_manifold = SYS6.is_manifold(0)
    print(f"   operador cos(theta_d) sobre M_J=0: ||C||_F = {np.linalg.norm(C):.4f}, "
          f"simétrico = {np.array_equal(C, C.T)}")
    # El dominio del remapeo lo fija el par MÁS ligado presente en la matriz.
    # Con la base completa ése es el (n+2)p, no el (n+3)s.
    turns = {l: 2.0 * SYS6.n_star_of_l(l) ** 2 for l in sorted(SYS6.levels)}
    print(f"   dominio del remapeo: {SYS6.domain_bounds()}")
    print(f"   puntos de retorno clásicos 2n*² por vecino: "
          + ", ".join(f"l={l}: {t:.1f}" for l, t in turns.items()))

    print("\n      R    |  E [E_h]              | E-E_man [GHz] | peso manif. | "
          "<cos> (más bajo) | <cos> (adiabático)")
    print("   --------|-----------------------|---------------|-------------|"
          "-----------------|-------------------")
    Rs = np.concatenate([np.arange(110.0, 400.0, 30.0), np.arange(400.0, 1140.0, 40.0)])
    v_prev = None
    rows, skipped = [], []
    for R in Rs:
        try:
            H = SYS6.hamiltonian(float(R))
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
