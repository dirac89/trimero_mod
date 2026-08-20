#!/usr/bin/env python3
"""
Qué cambia al completar la base electrónica con los tres vecinos del paper.

Aguilera-Fernández, Sadeghpour, Schmelcher & González-Férez, J. Phys.: Conf.
Ser. 635, 012023 (2015), arXiv:1507.07972, define la base como

    «the (n, l≥3) degenerate manifold, and the energetically neighboring
     levels (n+1)d, (n+2)p, and (n+3)s»

Las rondas anteriores usaron sólo manifold + (n+3)s. Este script mide, para
n=24, qué cambia al añadir (n+1)d y (n+2)p:

  [1] Composición y tamaño de la base, y coste de diagonalizar.
  [2] Tabla I: qué niveles son los vecinos y cuáles son puntos de contraste.
  [3] Convergencia en N (N_max=6 vs N_max=8), base completa vs incompleta.
  [4] Orientación <cos theta_d>, base completa vs incompleta.

Ejecutar:  poetry run python scripts/run_basis_correction_check.py
"""
import argparse
import time

import numpy as np

from trimero.systems.rb_krb_polar.charge_dipole import B_KRB_GHZ
from trimero.systems.rb_krb_polar.bop_system import GHZ_PER_HARTREE as GHZ, BOPSystem
from trimero.systems.rb_krb_polar.rb_defects import (
    DELTA0_NS_PAPER as D0, energy_rb, n_star_nl, neighbor_levels,
)

RULE = "=" * 88
L_NAME = {0: "s", 1: "p", 2: "d", 3: "f"}


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-manifold", type=int, default=24)
    ap.add_argument("--skip-n8", action="store_true",
                    help="salta la convergencia N=8 (es lo más caro)")
    return ap.parse_args()


def cos_theta_matrix(system, block):
    """<i|cos(theta_d)|j> = delta_ll' delta_mm' <N M|cos|N' M'>. No depende de R."""
    states = block.states
    index = {st: i for i, st in enumerate(states)}
    C = np.zeros((len(states), len(states)))
    for i, (l, m, N, MN) in enumerate(states):
        for Np in (N - 1, N + 1):
            if Np < 0 or abs(MN) > Np:
                continue
            j = index.get((l, m, Np, MN))
            if j is not None:
                C[i, j] = system.hmol.cos_theta_element(N, MN, Np, MN)
    return C


def orientation_curve(system, block, C, Rs):
    """(R, E, peso de manifold, <cos>) del más bajo con carácter de manifold."""
    mask = system.is_manifold(0)
    out = []
    for R in Rs:
        try:
            H = system.hamiltonian(float(R))
        except ValueError:
            continue
        w, V = np.linalg.eigh(H)
        pick = None
        for idx in np.argsort(w):
            if float(np.sum(V[:, idx][mask] ** 2)) > 0.5:
                pick = idx
                break
        if pick is None:
            continue
        v = V[:, pick]
        out.append((float(R), float(w[pick]),
                    float(np.sum(v[mask] ** 2)), float(v @ C @ v)))
    return out


def main():
    args = parse_args()
    n = args.n_manifold
    print(RULE)
    print(f"BASE CORREGIDA (manifold + 3 vecinos) FRENTE A LA INCOMPLETA — n = {n}")
    print(RULE)

    full = BOPSystem(n_manifold=n, delta0_ns=D0)
    part = BOPSystem(n_manifold=n, delta0_ns=D0, neighbors=(0,))

    # ------------------------------------------------ [1] composición
    lv = neighbor_levels(n)
    print("\n[1] COMPOSICIÓN Y TAMAÑO DE LA BASE\n")
    print(f"    manifold: n={n}, l = {full.l_min}..{full.l_max}"
          f"   ({sum(2*l+1 for l in range(full.l_min, full.l_max+1))} estados "
          "electrónicos)")
    print("    vecinos individuales (l identifica el nivel):")
    for l in (2, 1, 0):
        ns = n_star_nl(lv[l], l, D0 if l == 0 else None)
        print(f"      l={l}  ->  {lv[l]}{L_NAME[l]}   n* = {ns:.6f}   "
              f"radial hidrogenoide n = {full.radial.n_of_l(l)}   "
              f"({2*l+1} estados m_l)")
    print(f"\n    base electrónica total: "
          f"{sum(2*l+1 for l in range(full.l_min, full.l_max+1))} + 5 + 3 + 1 "
          f"= {full.basis.total_dimension()//(( full.N_max+1)**2)} = {n}²")
    print(f"    dim total (x rotor {(full.N_max+1)**2}): "
          f"{part.basis.total_dimension()} -> {full.basis.total_dimension()}")
    print(f"    dim(bloque M_J=0, N_max={full.N_max}): "
          f"{len(part.block(0))} (incompleta) -> {len(full.block(0))} (correcta)"
          f"   [+{len(full.block(0))-len(part.block(0))}, "
          f"+{100*(len(full.block(0))/len(part.block(0))-1):.1f} %]")

    print("\n    coste de UN punto del barrido (R nuevo cada vez, sin cachés "
          "reutilizables;\n    media de 3 R distintos):")
    for name, s in (("incompleta", part), ("correcta", full)):
        s.hamiltonian(600.0)                       # calienta las tablas radiales
        tb = te = 0.0
        for R in (611.0, 623.0, 637.0):
            t0 = time.perf_counter()
            H = s.hamiltonian(R)
            t1 = time.perf_counter()
            np.linalg.eigvalsh(H)
            t2 = time.perf_counter()
            tb += t1 - t0
            te += t2 - t1
        tb /= 3.0
        te /= 3.0
        print(f"      {name:11s} dim={len(s.block(0)):5d}   "
              f"construir H: {tb:.3f} s   eigvalsh: {te:.3f} s   "
              f"total {tb+te:.3f} s")

    lo_p, hi_p = part.domain_bounds()
    lo_f, hi_f = full.domain_bounds()
    print("\n    dominio del remapeo k(R) — OJO, es el límite de cómo LEEMOS las")
    print("    tablas de dispersión, no donde se acaba la física:")
    print(f"      incompleta: [{lo_p:.2f}, {hi_p:.2f}] a0")
    print(f"      correcta:   [{lo_f:.2f}, {hi_f:.2f}] a0")
    turns = {l: 2.0 * full.n_star_of_l(l) ** 2 for l in (0, 1, 2, 3)}
    print("      puntos de retorno clásicos 2n*²: "
          + ", ".join(f"{lv.get(l, n)}{L_NAME[l]}: {t:.1f}" for l, t in turns.items()))
    l_bind = min(turns, key=turns.get)
    print(f"      el par más ligado pasa a ser {lv[l_bind]}{L_NAME[l_bind]}-"
          f"{lv[l_bind]}{L_NAME[l_bind]}: el dominio se ACORTA "
          f"{hi_p-hi_f:.1f} a0")

    # ------------------------------------------------ [2] Tabla I
    print("\n" + "-" * 88)
    print("[2] TABLA I: QUÉ NIVELES SON VECINOS DE LA BASE Y CUÁLES SON CONTRASTE")
    print("-" * 88)
    E_ref = energy_rb(n, 3, D0)
    tabla = [(n + 1, 3), (n + 4, 0), (n + 2, 2), (n + 3, 1),
             (n, 3), (n + 3, 0), (n + 1, 2), (n + 2, 1)]
    print(f"    referencia de la tabla: E({n}, l=3) = {E_ref:.15e} E_h\n")
    print("      n | l |     E_nl [E_h]      |    n*      | ΔE vs manifold [GHz] | "
          "¿en la base?")
    print("    ----|---|---------------------|------------|----------------------|"
          "-------------")
    en_base = 0
    for nn, ll in sorted(tabla, key=lambda t: -energy_rb(t[0], t[1], D0)):
        E = energy_rb(nn, ll, D0)
        ns = n_star_nl(nn, ll, D0 if ll == 0 else None)
        dE = abs(E - E_ref) * GHZ
        if ll >= 3 and nn == n:
            tag = "SÍ (manifold)"
        elif lv.get(ll) == nn:
            tag = f"SÍ (vecino l={ll})"
            en_base += 1
        else:
            tag = "no (contraste)"
        print(f"    {nn:3d} | {ll:1d} | {E:.15e} | {ns:10.7f} | {dE:20.5f} | {tag}")
    print(f"\n    vecinos de la base presentes en la Tabla I: {en_base} de 3")
    print(f"    -> los tres vecinos son {lv[2]}d, {lv[1]}p, {lv[0]}s, y los tres")
    print("       están en la tabla. Los otros cuatro niveles ({}, {}, {}, {}) son"
          .format(f"{n+1},l=3", f"{n+2},l=2", f"{n+3},l=1", f"{n+4},l=0"))
    print("       puntos de contraste: el siguiente de cada serie, NO entran en la base.")
    print("    (la Tabla I son energías atómicas: su valor no depende de la base;")
    print("     lo que se verifica aquí es que entendimos qué es cada fila)")

    # ------------------------------------------------ [3] convergencia en N
    print("\n" + "-" * 88)
    print("[3] CONVERGENCIA EN N (N_max=6 vs N_max=8)")
    print("-" * 88)
    if args.skip_n8:
        print("    saltado por --skip-n8")
    else:
        for name, neigh in (("incompleta (manifold+ns)", (0,)),
                            ("correcta (manifold+3)", (0, 1, 2))):
            s6 = BOPSystem(n_manifold=n, delta0_ns=D0, neighbors=neigh, N_max=6)
            s8 = BOPSystem(n_manifold=n, delta0_ns=D0, neighbors=neigh, N_max=8)
            print(f"\n    base {name}: dim(N=6) = {len(s6.block(0))}, "
                  f"dim(N=8) = {len(s8.block(0))}")
            print("      criterio       | R [a0] | peor |ΔE| [MHz] | peor |ΔE|/|E| | "
                  "peor |ΔE|/|E-E_man| | ¿rel(E) < 2e-6?")
            print("      ---------------|--------|----------------|---------------|"
                  "---------------------|----------------")
            for R in (600.0, 900.0):
                w6, V6 = s6.solve(R)
                w8, V8 = s8.solve(R)
                k6 = s6.lowest_manifold_index(R)
                k8 = s8.lowest_manifold_index(R)
                for crit, e6, e8 in (
                    ("10 más bajos", w6[:10], w8[:10]),
                    (f"manifold k={k6}", w6[k6:k6 + 1], w8[k8:k8 + 1]),
                ):
                    d = np.abs(e8 - e6)
                    rel = d / np.abs(e6)
                    relb = d / np.abs(e6 - s6.E_manifold)
                    print(f"      {crit:14s} | {R:6.0f} | {d.max()*GHZ*1e3:14.5f} | "
                          f"{rel.max():13.4e} | {relb.max():20.4e} | "
                          f"{'CUMPLE' if rel.max() < 2e-6 else 'NO CUMPLE':>15s}")
            # qué son los 10 más bajos: cambia al completar la base
            w, V = s6.solve(900.0)
            mask = s6.is_manifold(0)
            l_of = np.array([st[0] for st in s6.block(0).states])
            print("      composición de los 10 más bajos en R=900 (l dominante, "
                  "E-E_man [GHz]):")
            comp = []
            for k in range(10):
                c2 = V[:, k] ** 2
                l_dom = int(l_of[np.argmax(c2)])
                comp.append(f"l={l_dom}:{(w[k]-s6.E_manifold)*GHZ:.1f}")
            print("        " + "  ".join(comp))

    # ------------------------------------------------ [4] orientación
    print("\n" + "-" * 88)
    print("[4] ORIENTACIÓN <cos theta_d>  (más bajo con carácter de manifold)")
    print("-" * 88)
    Rs = np.concatenate([np.arange(110.0, 400.0, 30.0), np.arange(400.0, 600.0, 40.0)])
    curves = {}
    for name, s in (("incompleta", part), ("correcta", full)):
        C = cos_theta_matrix(s, s.block(0))
        curves[name] = {r[0]: r for r in orientation_curve(s, s.block(0), C, Rs)}
        print(f"    ||C||_F ({name}) = {np.linalg.norm(C):.4f}")
    print("\n      R [a0] | <cos> incompleta | <cos> correcta |   Δ<cos>   | "
          "E-E_man incompl. | E-E_man correcta [GHz]")
    print("      -------|------------------|----------------|------------|"
          "------------------|----------------------")
    worst = 0.0
    for R in Rs:
        a = curves["incompleta"].get(float(R))
        b = curves["correcta"].get(float(R))
        if a is None or b is None:
            continue
        dc = b[3] - a[3]
        worst = max(worst, abs(dc))
        print(f"      {R:6.1f} | {a[3]:+16.6f} | {b[3]:+14.6f} | {dc:+10.2e} | "
              f"{(a[1]-part.E_manifold)*GHZ:16.4f} | "
              f"{(b[1]-full.E_manifold)*GHZ:20.4f}")
    print(f"\n    peor |Δ<cos>| en la rama suave = {worst:.3e}")
    for name in ("incompleta", "correcta"):
        rows = sorted(curves[name].values())
        cruces = [r for r in rows if abs(r[3]) < 0.78]
        if cruces and rows.index(cruces[0]) > 0:
            prev = rows[rows.index(cruces[0]) - 1]
            print(f"    base {name}: cruza 0.78 entre R = {prev[0]:.0f} y "
                  f"{cruces[0][0]:.0f} a0   (máx en el borde: {rows[0][3]:+.6f} "
                  f"en R = {rows[0][0]:.0f})")
    print(RULE)


if __name__ == "__main__":
    main()
