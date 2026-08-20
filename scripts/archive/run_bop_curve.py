#!/usr/bin/env python3
"""
Curva BOP de M_J=0 para Rb*-KRb: la más baja que evoluciona del manifold n=24.

H = H_a + H_mol + V_Fermi, delta0_ns = 3.13180 (comparación con
González-Férez, Sadeghpour & Schmelcher, NJP 17, 013021, 2015, Fig. 2).

Ejecutar:  poetry run python scripts/run_bop_curve.py
"""
import argparse
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from trimero.systems.rb_krb_polar.bop_system import GHZ_PER_HARTREE, BOPSystem
from trimero.simulation.bop_tracking import trace_curve
from trimero.systems.rb_krb_polar.rb_defects import DELTA0_NS_PAPER as D0

# El manifold ya no está hardcodeado: se monta en BOPSystem. Este script es el
# caso n=24 + 27s del paper; para n=25 + 28s ver scripts/run_bop_manifold.py.
SYSTEM = BOPSystem(n_manifold=24, n_s=27, delta0_ns=D0)

GHZ = GHZ_PER_HARTREE
E_MAN = SYSTEM.E_manifold         # manifold n=24 + KRb(N=0): cero de energía
N_CONTEXT = 60

BLOCK = SYSTEM.block(0)
FERMI = SYSTEM.fermi
IS_MANIFOLD = SYSTEM.is_manifold(0)


def solve(R):
    return SYSTEM.solve(R)


def select_lowest_manifold(w, V):
    for k in np.argsort(w):
        if float(np.sum(V[:, k][IS_MANIFOLD] ** 2)) > 0.5:
            return int(k)
    raise RuntimeError("ningún autoestado con >50% de peso en el manifold")


def local_minima(R, y):
    """Mínimos locales estrictos en una malla no uniforme."""
    out = []
    for i in range(1, len(y) - 1):
        if y[i] < y[i - 1] and y[i] < y[i + 1]:
            out.append(i)
    return out


def parse_args():
    ap = argparse.ArgumentParser(description="Curva BOP M_J=0 para Rb*-KRb")
    ap.add_argument("--reuse", action="store_true",
                    help="reutiliza plots/archive/bop_curve_MJ0.npz en vez de rebarrer R")
    ap.add_argument("--step", type=float, default=5.0, help="paso base en R [a0]")
    ap.add_argument("--ymin", type=float, default=-120.0,
                    help="límite inferior de la VENTANA DE ENERGÍA [GHz]")
    ap.add_argument("--ymax", type=float, default=25.0,
                    help="límite superior de la ventana de energía [GHz]")
    ap.add_argument("--ap-threshold", type=float, default=1.0e5,
                    help="|A_p| [a0^3] por encima del cual un punto se marca "
                         "como dominado por la resonancia de forma p")
    ap.add_argument("--exclude-resonance", action="store_true",
                    help="además de marcarlos, elimina de la curva los puntos "
                         "cuyo A_p se obtuvo interpolando a través del polo")
    ap.add_argument("--out", default="plots/archive/bop_curve_MJ0.png")
    return ap.parse_args()


def resonance_flags(R_array, ap_threshold):
    """
    (mask_ap, mask_gap) para cada R:
      mask_ap  : |A_p| supera el umbral -> zona dominada por la resonancia
      mask_gap : el R' remapeado cae en un hueco de malla de la tabla, es decir
                 A_p viene de interpolar a través de la divergencia. Estos
                 puntos NO son una lectura de la tabla.
    """
    m_ap, m_gap = [], []
    for R in R_array:
        try:
            ns = FERMI.n_star_of_l(3)
            _, A_p = FERMI.scattering.scattering_pair(R, ns, ns)
        except ValueError:
            A_p = 0.0
        m_ap.append(abs(A_p) > ap_threshold)
        m_gap.append(FERMI.bridges_gap(R))
    return np.array(m_ap), np.array(m_gap)


def main():
    args = parse_args()
    # El dominio lo fija el par MÁS restrictivo de toda la base (27s-27s),
    # calculado, no elegido a mano.
    R_lo, R_max = FERMI.domain_bounds()
    R_min = max(110.0, np.ceil(R_lo))
    print("=" * 84)
    print("CURVA BOP  M_J = 0   —   Rb*-KRb, manifold n=24 (l>=3) + 27s")
    print("=" * 84)
    print(f"  dim(bloque M_J=0) = {len(BLOCK)}")
    print(f"  dominio calculado sobre TODOS los pares (l1,l2): "
          f"R in [{R_lo:.2f}, {R_max:.2f}] a0")
    print(f"  el par más restrictivo es 27s-27s (el más ligado de la base)")
    print(f"  se traza desde R = {R_min:.1f} a0")
    print(f"  cero de energía: E(n=24) + KRb(N=0) = {E_MAN:.12e} E_h")

    npz = "plots/archive/bop_curve_MJ0.npz"
    if args.reuse and __import__("os").path.exists(npz):
        d = np.load(npz)
        res = {"R": d["R"], "E": d["E"], "index": d["index"], "context": d["context"],
               "overlap": np.ones_like(d["R"]), "refinements": [], "n_solves": 0}
        print(f"\n  (reutilizando {npz}, sin rebarrer R)")
        dt = 0.0
    else:
        t0 = time.perf_counter()
        res = trace_curve(
            solve, R_min, R_max, step=args.step,
            select_initial=select_lowest_manifold,
            overlap_threshold=0.7, max_refine=6, n_context=N_CONTEXT,
            progress=print,
        )
        dt = time.perf_counter() - t0
    R, E = res["R"], res["E"]
    y = (E - E_MAN) * GHZ

    print(f"\n  {len(R)} puntos trazados, {res['n_solves']} diagonalizaciones, {dt:.1f} s")
    print(f"  solapamiento |<v|v_prev>|²: min = {res['overlap'][1:].min():.4f}, "
          f"media = {res['overlap'][1:].mean():.4f}")
    print(f"  refinamientos por solapamiento bajo: {len(res['refinements'])}")
    if res["refinements"]:
        print("     R_desde   R_probado  solapam.  profundidad  paso_nuevo")
        for r in res["refinements"]:
            print(f"     {r['R_from']:9.3f} {r['R_try']:10.3f} {r['overlap']:9.4f} "
                  f"{r['depth']:11d} {r['new_step']:11.4f}")
    idx = res["index"]
    cambios = np.flatnonzero(np.diff(idx) != 0)
    print(f"\n  índice ordenado del estado seguido: {idx[0]} -> {idx[-1]}, "
          f"{len(cambios)} cambios")
    if len(cambios):
        print("     cambia en R = " + ", ".join(f"{R[i+1]:.1f}" for i in cambios[:20]))
        print("     (un cambio de índice = cruce genuino entre especies sigma_v de M_J=0)")

    # -------- características cuantitativas --------
    print("\n" + "-" * 84)
    print("CARACTERÍSTICAS DE LA CURVA (calculadas de nuestros números)")
    print("-" * 84)
    mins = local_minima(R, y)
    maxs = [i for i in range(1, len(y) - 1) if y[i] > y[i - 1] and y[i] > y[i + 1]]
    print(f"  {len(mins)} mínimos locales, {len(maxs)} máximos locales\n")
    print("    #  |  R_min [a0] |  E-E_man [GHz] | prof. vs máx. derecho [GHz] | ΔR al anterior")
    print("   ----|-------------|----------------|-----------------------------|---------------")
    prev_R = None
    seps = []
    for j, i in enumerate(mins):
        nxt = [m for m in maxs if m > i]
        depth = (y[nxt[0]] - y[i]) if nxt else float("nan")
        dR = (R[i] - prev_R) if prev_R is not None else float("nan")
        if prev_R is not None:
            seps.append(R[i] - prev_R)
        prev_R = R[i]
        print(f"   {j:4d}| {R[i]:11.2f} | {y[i]:14.4f} | {depth:27.4f} | "
              + (f"{dR:14.2f}" if np.isfinite(dR) else "           ---"))
    if len(seps) >= 2:
        seps = np.array(seps)
        print(f"\n   separación entre mínimos: media {seps.mean():.2f} a0, "
              f"sigma {seps.std():.2f} a0, min {seps.min():.2f}, max {seps.max():.2f}")
        print(f"   coeficiente de variación = {seps.std()/seps.mean():.3f}  "
              f"({'aprox. regular' if seps.std()/seps.mean() < 0.25 else 'NO regular'})")

    print(f"\n  Rango de la curva: [{y.min():.3f}, {y.max():.3f}] GHz")
    print(f"  Mínimo global: {y.min():.4f} GHz en R = {R[np.argmin(y)]:.2f} a0")

    # afirmación del texto del paper: desplazada > 20 GHz para R <~ 1200 a0
    frac = float(np.mean(np.abs(y) > 20.0))
    below = R[np.abs(y) <= 20.0]
    print(f"\n  Texto del paper: 'shifted more than 20 GHz for R <~ 1200 a0' (M_J=0)")
    print(f"    fracción de nuestro dominio con |E-E_man| > 20 GHz: {frac*100:.2f}%")
    if len(below):
        print(f"    puntos con |E-E_man| <= 20 GHz: {len(below)} de {len(R)}, "
              f"en R in [{below.min():.1f}, {below.max():.1f}] a0")
    else:
        print("    ningún punto del dominio queda por debajo de 20 GHz")

    # ---------------- ventana de energía y marcado de la resonancia -------
    ctx = res["context"]
    Y = (ctx - E_MAN) * GHZ                 # índice fijo = curva adiabática
    m_ap, m_gap = resonance_flags(R, args.ap_threshold)
    print("\n" + "-" * 84)
    print("VENTANA DE ENERGÍA Y RESONANCIA DE FORMA p")
    print("-" * 84)
    print(f"  ventana: [{args.ymin}, {args.ymax}] GHz")
    print(f"  umbral |A_p| = {args.ap_threshold:.1e} a0^3 -> {m_ap.sum()} de {len(R)} "
          f"puntos marcados como dominados por la resonancia")
    print(f"  huecos de malla en la tabla: {FERMI.large_gaps()}")
    print(f"  puntos cuyo A_p se interpola A TRAVÉS del polo: {m_gap.sum()}"
          + (f"  -> R = {', '.join(f'{r:.1f}' for r in R[m_gap])}" if m_gap.any() else ""))
    if args.exclude_resonance and m_gap.any():
        print("  --exclude-resonance: esos puntos se eliminan de la curva dibujada")

    # curva adiabática de carácter manifold: la identificamos en el borde
    w_edge, V_edge = solve(float(R[-1]))
    k_man = next(k for k in range(len(w_edge))
                 if float(np.sum(V_edge[:, k][IS_MANIFOLD] ** 2)) > 0.5)
    print(f"  curva adiabática más baja de carácter manifold en R_max: k = {k_man}")

    keep = ~m_gap if args.exclude_resonance else np.ones(len(R), bool)

    # ---------------- plot --------------------------------------------
    fig, axes = plt.subplots(2, 1, figsize=(9.5, 9.5), sharex=True)
    for ax, ylim, ttl in (
        (axes[0], (float(np.nanmin(Y[:, k_man])) * 1.05, 100.0), "Rango completo"),
        (axes[1], (args.ymin, args.ymax), "Ventana de energía (comparable a la Fig. 2)"),
    ):
        for lo, hi in FERMI.large_gaps():
            pass
        if m_gap.any():
            ax.axvspan(R[m_gap].min() - 2.5, R[m_gap].max() + 2.5,
                       color="gold", alpha=0.30, zorder=0,
                       label="$A_p$ interpolado a través del polo")
        for c in range(Y.shape[1]):
            ax.plot(R, Y[:, c], color="0.78", lw=0.6, zorder=1)
        ax.plot(R[keep], Y[keep, k_man], color="crimson", lw=3.0, zorder=4,
                label=f"adiabática k={k_man} (carácter manifold $n=24$)")
        ax.plot(R, y, color="darkblue", lw=1.5, ls="--", zorder=3,
                label="seguida por solapamiento (diabática)")
        if m_ap.any():
            ax.plot(R[m_ap], Y[m_ap, k_man], "o", ms=5, mfc="none",
                    mec="darkorange", mew=1.6, zorder=5,
                    label=fr"$|A_p| > 10^{{{np.log10(args.ap_threshold):.0f}}}\,a_0^3$")
        ax.axhline(0, color="k", lw=0.8, ls=":", zorder=2)
        ax.axhline(-20, color="seagreen", lw=1.2, ls="--", zorder=2,
                   label="$-20$ GHz (texto del paper)")
        ax.set_ylim(*ylim); ax.set_xlim(R[0], R[-1]); ax.grid(alpha=0.25)
        ax.set_ylabel(r"$E - E_{n=24}$  [GHz]"); ax.set_title(ttl, fontsize=10)
        ax.legend(loc="lower right", fontsize=8)
    for N in range(7):
        axes[1].axhline(-63.40 + 1.114 * N * (N + 1), color="tab:orange", lw=0.7, alpha=0.6)
        axes[1].axhline(1.114 * N * (N + 1), color="tab:blue", lw=0.7, alpha=0.6)
    axes[1].set_xlabel(r"$R$  [$a_0$]")
    fig.suptitle("Rb*-KRb, curvas BOP $M_J=0$   "
                 r"$H=H_a+B\mathbf{N}^2-\mathbf{d}\cdot\mathbf{F}_{ryd}+V_{Fermi}$"
                 "\nnaranja: umbrales 27s+KRb(N) · azul: manifold n=24+KRb(N)", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.945])
    fig.savefig(args.out, dpi=150)
    print(f"\n  PNG guardado en {args.out}")
    np.savez("plots/archive/bop_curve_MJ0.npz", R=R, E=E, y=y, index=idx, context=ctx)
    print("  datos en plots/archive/bop_curve_MJ0.npz")
    print("=" * 84)


if __name__ == "__main__":
    main()
