#!/usr/bin/env python3
"""
Fig. 1(a) de Aguilera-Fernández et al. 2015: curvas BOP M_J=0 para n=24 y n=25.

    poetry run python scripts/run_fig1_curves.py

Base: manifold (n, l>=3) + (n+1)d + (n+2)p + (n+3)s (arXiv:1507.07972).
Curva: la más baja con CARÁCTER de manifold (peso > 50 %) en cada R — no el
índice fijo, que con la base completa cambia de carácter (ver
docs/analysis_base_correcta_3_vecinos.md §5.1).

ORDEN DE LO QUE HACE, que es también el orden en que hay que leerlo:

  [1] COMPROBACIÓN PREVIA, antes de extender nada: mide ||V_Fermi||, la parte
      de H_mol que depende de R, y sobre todo cuánto se mueve la curva de
      carácter si se pone V_Fermi = 0, todo cerca del borde del dominio del
      remapeo k(R). Si el salto en el empalme no es pequeño, la extensión
      "H_a+H_mol solo" NO está justificada y el script lo dice.
  [2] Barrido con H completo dentro del dominio del remapeo.
  [3] Barrido con V_Fermi ≡ 0 más allá del dominio, SÓLO como diagnóstico de
      cuánto se despega: se dibuja marcado como rechazado si [1] falla.
  [4] Figura con las dos curvas en R ∈ [400, 1800] a0, con la ventana de la
      resonancia p y el final del dominio marcados explícitamente.
  [5] Comparación cualitativa n=24 vs n=25 contra las tres tendencias que
      describe el texto del paper.
"""
import argparse
import os
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from trimero.systems.rb_krb_polar.bop_system import GHZ_PER_HARTREE as GHZ, BOPSystem
from trimero.systems.rb_krb_polar.rb_defects import DELTA0_NS_PAPER as D0

RULE = "=" * 96
# Criterio de aceptación de la extensión, fijado ANTES de mirar los números:
# el salto en el empalme al poner V_Fermi = 0 debe ser pequeño frente a la
# escala vertical de la figura que queremos comparar (~decenas de GHz).
JUMP_TOL_GHZ = 1.0
JUMP_TOL_REL = 0.05          # ... y frente a la propia ligadura de la curva
# Umbral para delimitar la cola de la resonancia p: mismo criterio que en
# docs/analysis_ventana_exclusion_resonancia.md §3.
TAIL_GHZ = 120.0


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n", type=int, nargs="+", default=[24, 25])
    ap.add_argument("--rmin", type=float, default=400.0)
    ap.add_argument("--rmax", type=float, default=1800.0)
    ap.add_argument("--step", type=float, default=5.0)
    ap.add_argument("--step-ext", type=float, default=10.0,
                    help="paso en la región sin V_Fermi (más barata y más lisa)")
    ap.add_argument("--ymin", type=float, default=-60.0)
    ap.add_argument("--ymax", type=float, default=5.0)
    ap.add_argument("--reuse", action="store_true")
    ap.add_argument("--out", default="plots/archive/fig1_MJ0_n24_n25.png")
    return ap.parse_args()


# ------------------------------------------------------------- paso 1
def edge_check(sysm, n_points=8, span=150.0):
    """
    ¿Es V_Fermi despreciable en el borde del dominio? Devuelve (ok, filas).

    Tres medidas: norma de V_Fermi, norma de la parte de H_mol que depende de
    R, y el salto en la curva de carácter al poner V_Fermi = 0 — que es el
    error que introduciría la extrapolación justo en el empalme.
    """
    blk = sysm.block(0)
    hi = sysm.domain_bounds()[1]
    E_rot = sysm.hmol.rotational_diagonal(blk)
    # densidad electrónica del manifold, normalizada al máximo del lóbulo externo
    u2 = sysm.radial.u(sysm.l_min) ** 2
    r = sysm.radial.r
    dens_max = float(u2[r > 0.5 * 2 * sysm.n_manifold**2].max())

    rows = []
    for R in np.linspace(hi - span, hi - 0.5, n_points):
        Vf = sysm.fermi.build(blk, float(R))
        Vm = sysm.hmol.build(blk, float(R)) - np.diag(E_rot)
        e_full = sysm.character_curve(float(R))[0]
        e_nof = sysm.character_curve(float(R), fermi=False)[0]
        rows.append({
            "R": float(R),
            "nf": float(np.linalg.norm(Vf)),
            "nm": float(np.linalg.norm(Vm)),
            "e_full": e_full,
            "e_nof": e_nof,
            "jump": e_full - e_nof,
            "dens": float(np.interp(R, r, u2)) / dens_max,
        })
    last = rows[-1]
    ok = (abs(last["jump"]) < JUMP_TOL_GHZ
          and abs(last["jump"]) < JUMP_TOL_REL * abs(last["e_full"]))
    return ok, rows


def print_edge_check(sysm, rows, ok):
    hi = sysm.domain_bounds()[1]
    print(f"\n  n = {sysm.n_manifold}:  borde del dominio del remapeo "
          f"R_max = {hi:.2f} a0   (2n² = {2*sysm.n_manifold**2} a0)")
    print("     R    | ||V_F||_F  | ||H_mol-d||_F | V_F/H_mol | E completo | "
          "E sin Fermi |  SALTO   | |psi(R)|²/max")
    print("   -------|------------|---------------|-----------|------------|"
          "-------------|----------|--------------")
    for x in rows:
        print(f"   {x['R']:7.1f}| {x['nf']:.4e} | {x['nm']:.7e} | "
              f"{x['nf']/x['nm']:9.3f} | {x['e_full']:10.4f} | {x['e_nof']:11.4f} | "
              f"{x['jump']:+8.4f} | {x['dens']:13.4f}")
    last = rows[-1]
    print(f"\n   criterio fijado de antemano: |salto| < {JUMP_TOL_GHZ} GHz  Y  "
          f"< {JUMP_TOL_REL*100:.0f} % de |E|")
    print(f"   medido en el empalme: |salto| = {abs(last['jump']):.4f} GHz = "
          f"{100*abs(last['jump']/last['e_full']):.1f} % de |E| = "
          f"{abs(last['e_full']):.4f} GHz")
    print(f"   -> {'PASA' if ok else 'NO PASA'}: la extensión con V_Fermi = 0 "
          f"{'está' if ok else 'NO está'} justificada")
    return ok


# ------------------------------------------------------------- barridos
def sweep(sysm, R_values, fermi, label):
    E = np.full(len(R_values), np.nan)
    K = np.full(len(R_values), -1, dtype=int)
    W = np.zeros(len(R_values))
    t0 = time.perf_counter()
    for i, R in enumerate(R_values):
        E[i], K[i], W[i] = sysm.character_curve(float(R), fermi=fermi)
        if i % 40 == 0:
            print(f"      {label} R = {R:8.2f} ({i+1}/{len(R_values)})", flush=True)
    print(f"      {label}: {len(R_values)} puntos en "
          f"{time.perf_counter()-t0:.1f} s")
    return E, K, W


def local_minima(R, y):
    return [i for i in range(1, len(y) - 1)
            if np.isfinite(y[i - 1:i + 2]).all()
            and y[i] < y[i - 1] and y[i] < y[i + 1]]


def local_maxima(R, y):
    return [i for i in range(1, len(y) - 1)
            if np.isfinite(y[i - 1:i + 2]).all()
            and y[i] > y[i - 1] and y[i] > y[i + 1]]


def main():
    args = parse_args()
    print(RULE)
    print("Fig. 1(a) — curvas BOP M_J=0 de Rb*-KRb para varios manifolds")
    print(RULE)

    systems = {n: BOPSystem(n_manifold=n, delta0_ns=D0) for n in args.n}

    # ---------------------------------------------------- [1] comprobación
    print("\n[1] COMPROBACIÓN PREVIA: ¿es V_Fermi despreciable en el borde?")
    print("-" * 96)
    print("    Se mide antes de extender nada. El SALTO es lo que decide: es el")
    print("    error que la extrapolación V_Fermi=0 mete justo en el empalme.")
    verdict = {}
    for n, s in systems.items():
        ok, rows = edge_check(s)
        verdict[n] = print_edge_check(s, rows, ok)
    extend_ok = all(verdict.values())
    print("\n" + "-" * 96)
    if extend_ok:
        print("    RESULTADO: los dos manifolds pasan. Se extiende con H_a+H_mol.")
    else:
        print("    RESULTADO: NO pasa. La región R > R_max se calcula igualmente,")
        print("    pero se marca como EXTENSIÓN RECHAZADA y no se usa para ninguna")
        print("    conclusión cuantitativa. Ver el documento de análisis.")
    print("-" * 96)

    # ---------------------------------------------------- [2]/[3] barridos
    data = {}
    for n, s in systems.items():
        npz = f"plots/archive/fig1_MJ0_n{n}.npz"
        if args.reuse and os.path.exists(npz):
            d = np.load(npz)
            data[n] = {k: d[k] for k in d.files}
            print(f"\n[2] n={n}: reutilizando {npz}")
            continue
        hi = s.domain_bounds()[1]
        R_in = np.arange(args.rmin, hi, args.step)
        R_in = np.append(R_in, hi)
        R_ex = np.arange(hi + args.step_ext, args.rmax + 1e-9, args.step_ext)
        print(f"\n[2] n={n}: barrido CON V_Fermi en R ∈ [{args.rmin}, {hi:.2f}] "
              f"({len(R_in)} puntos)")
        E_in, K_in, W_in = sweep(s, R_in, True, f"n={n} Fermi")
        print(f"[3] n={n}: barrido SIN V_Fermi en R ∈ ({hi:.2f}, {args.rmax}] "
              f"({len(R_ex)} puntos)")
        E_ex, K_ex, W_ex = sweep(s, R_ex, False, f"n={n} sin-Fermi")
        data[n] = {"R_in": R_in, "E_in": E_in, "K_in": K_in, "W_in": W_in,
                   "R_ex": R_ex, "E_ex": E_ex, "K_ex": K_ex, "W_ex": W_ex,
                   "R_max": np.array([hi])}
        np.savez(npz, **data[n])
        print(f"    datos en {npz}")

    # ---------------------------------------------------- [5] tendencias
    print("\n" + "-" * 96)
    print("[5] FORMA DE LAS CURVAS Y COMPARACIÓN CON EL TEXTO DEL PAPER")
    print("-" * 96)
    summary = {}
    for n, s in systems.items():
        d = data[n]
        R, y = d["R_in"], d["E_in"]
        win = s.p_resonance_window()
        # Además de la ventana, se quita la BANDA DE COLA butterfly: el tramo
        # justo a la derecha en el que la curva sigue a cientos de GHz por la
        # resonancia. Sin quitarlo, el primer punto tras la ventana se cuela
        # como "pozo" de 600 GHz, que no es un pozo de la curva comparable.
        tail_mask = (R > win["R_hi"]) & (np.abs(y) > TAIL_GHZ)
        R_tail = float(R[tail_mask].max()) if tail_mask.any() else win["R_hi"]
        keep = (R < win["R_lo"]) | (R > R_tail)
        Rk, yk = R[keep], y[keep]
        mins = local_minima(Rk, yk)
        maxs = local_maxima(Rk, yk)
        # profundidad de cada pozo respecto del máximo inmediatamente a su derecha
        wells = []
        for i in mins:
            nxt = [m for m in maxs if m > i]
            depth = (yk[nxt[0]] - yk[i]) if nxt else float("nan")
            wells.append((float(Rk[i]), float(yk[i]), float(depth)))
        summary[n] = {"wells": wells, "win": win, "R_tail": R_tail,
                      "R_max": float(s.domain_bounds()[1])}
        print(f"\n    n = {n}   ventana p excluida: "
              f"[{win['R_lo']:.1f}, {win['R_hi']:.1f}] a0; "
              f"cola butterfly hasta {R_tail:.1f} a0; "
              f"dominio hasta {s.domain_bounds()[1]:.1f} a0")
        print("      #  |  R_pozo [a0] | E-E_man [GHz] | profundidad [GHz] | R/2n²")
        print("      ---|--------------|---------------|-------------------|------")
        for j, (Rw, Ew, dep) in enumerate(wells):
            dtxt = f"{dep:17.4f}" if np.isfinite(dep) else " " * 13 + "---"
            print(f"      {j:3d}| {Rw:12.2f} | {Ew:13.4f} | {dtxt} | "
                  f"{Rw/(2*n*n):.3f}")

    if len(summary) == 2:
        a, b = sorted(summary)          # a = n menor, b = n mayor
        wa, wb = summary[a]["wells"], summary[b]["wells"]
        print("\n    Tres tendencias que describe el texto para el conjunto de curvas:")
        # 1) pozos a mayor R al aumentar n
        if wa and wb:
            print(f"\n    (1) ¿los pozos se desplazan a MAYOR R al aumentar n?")
            print(f"        pozo más externo:  n={a}: R = {wa[-1][0]:.1f} a0 "
                  f"({wa[-1][0]/(2*a*a):.3f}·2n²)   "
                  f"n={b}: R = {wb[-1][0]:.1f} a0 ({wb[-1][0]/(2*b*b):.3f}·2n²)")
            n_pair = min(len(wa), len(wb))
            for j in range(1, n_pair + 1):
                print(f"        pozo -{j} (contando desde fuera): "
                      f"n={a}: {wa[-j][0]:.1f} a0   n={b}: {wb[-j][0]:.1f} a0   "
                      f"Δ = {wb[-j][0]-wa[-j][0]:+.1f} a0")
            print(f"\n    (2) ¿la PROFUNDIDAD de los pozos decrece con n?")
            for j in range(1, n_pair + 1):
                da, db = wa[-j][2], wb[-j][2]
                if not (np.isfinite(da) and np.isfinite(db)):
                    print(f"        pozo -{j}: profundidad NO medible (no hay "
                          "máximo a su derecha dentro del dominio)")
                    continue
                print(f"        pozo -{j}: n={a}: {da:.4f} GHz   n={b}: {db:.4f} GHz"
                      f"   {'DECRECE' if db < da else 'CRECE'}")
            print(f"\n    (3) ¿la amplitud de la oscilación DECRECE con R?")
            for n in (a, b):
                w = summary[n]["wells"]
                deps = [d for _, _, d in w if np.isfinite(d)]
                if len(deps) >= 2:
                    print(f"        n={n}: profundidades de dentro a fuera = "
                          + ", ".join(f"{d:.3f}" for d in deps)
                          + ("  -> decreciente" if all(
                              deps[i] >= deps[i + 1] for i in range(len(deps) - 1))
                             else "  -> NO monótona"))
                else:
                    print(f"        n={n}: sólo {len(deps)} pozo(s) con "
                          "profundidad medible: sin base para juzgar")

    # ---------------------------------------------------- [4] figura
    fig, axes = plt.subplots(2, 1, figsize=(11.0, 9.0), sharex=True)
    colors = {24: "crimson", 25: "tab:blue"}
    for ax, (ylo, yhi), ttl in (
        (axes[0], (args.ymin, args.ymax),
         "Ventana de energía de la Fig. 1(a)"),
        (axes[1], (-160.0, 20.0),
         "Rango ampliado — se ve el fondo del pozo butterfly recortado"),
    ):
        for n, s in systems.items():
            d = data[n]
            win = summary[n]["win"] if n in summary else s.p_resonance_window()
            R, y = d["R_in"], d["E_in"]
            keep = (R < win["R_lo"]) | (R > win["R_hi"])
            c = colors.get(n, None)
            ax.plot(R, np.where(keep, y, np.nan), color=c, lw=2.6, zorder=4,
                    label=f"$n={n}$: $H_a+H_{{mol}}+V_{{Fermi}}$ "
                          f"(dominio del remapeo)")
            ax.axvspan(win["R_lo"], win["R_hi"], color=c, alpha=0.18, zorder=1)
            R_tail = summary.get(n, {}).get("R_tail", win["R_hi"])
            if R_tail > win["R_hi"]:
                ax.axvspan(win["R_hi"], R_tail, color=c, alpha=0.07,
                           hatch="///", ec=c, lw=0.0, zorder=1)
            hi = float(d["R_max"][0])
            ax.axvline(hi, color=c, lw=1.4, ls=":", zorder=3)
            ax.plot(d["R_ex"], d["E_ex"], color=c, lw=1.4, ls=(0, (1, 2)),
                    alpha=0.55, zorder=2,
                    label=f"$n={n}$: extensión $V_{{Fermi}}\\equiv0$ "
                          "— RECHAZADA" if not extend_ok else
                          f"$n={n}$: extensión $H_a+H_{{mol}}$")
            # el salto en el empalme, dibujado
            if np.isfinite(y[-1]) and len(d["E_ex"]):
                ax.plot([hi, hi], [y[-1], d["E_ex"][0]], color=c, lw=1.2,
                        ls="-", alpha=0.9, zorder=5)
                ax.plot([hi], [d["E_ex"][0]], "v", color=c, ms=6, zorder=6)
        ax.axhline(0, color="k", lw=0.8, ls=":", zorder=2)
        ax.set_xlim(args.rmin, args.rmax)
        ax.set_ylim(ylo, yhi)
        ax.grid(alpha=0.25)
        ax.set_ylabel(r"$E - E_{n,l\geq3}$  [GHz]")
        ax.set_title(ttl, fontsize=10)
        ax.legend(loc="lower right", fontsize=8)
    for n, s in systems.items():
        for ax in axes:
            ax.axvline(2.0 * n * n, color=colors.get(n), lw=1.0, ls="-.",
                       alpha=0.45, zorder=1)
    axes[1].set_xlabel(r"$R$  [$a_0$]")
    txt = ("banda llena = ventana excluida de la resonancia p   ·   "
           "banda rayada = cola butterfly (pozos no comparables)\n"
           "punteado vertical = fin del dominio del remapeo $k(R)$   ·   "
           "raya-punto = $2n^2 a_0$   ·   "
           r"línea continua + triángulo = SALTO al poner $V_{Fermi}\equiv0$")
    if not extend_ok:
        txt += ("\n⚠️ la extensión punteada NO está justificada: el salto en el "
                "empalme es de 12-14 GHz (ver docs/analysis_extension_dominio_fig1.md)")
    fig.suptitle(
        "Rb*-KRb, curvas BOP $M_J=0$ — comparación con la Fig. 1(a) de "
        "Aguilera-Fernández et al. 2015\n"
        "base: manifold $(n,l\\geq3)$ + $(n{+}1)d$ + $(n{+}2)p$ + $(n{+}3)s$   ·   "
        "curva de carácter manifold (>50 % de peso)\n" + txt,
        fontsize=9.5)
    fig.tight_layout(rect=[0, 0, 1, 0.86])
    fig.savefig(args.out, dpi=150)
    print(f"\n    PNG guardado en {args.out}")
    print(RULE)


if __name__ == "__main__":
    main()
