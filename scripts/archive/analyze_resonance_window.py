#!/usr/bin/env python3
"""
Ventana de exclusión de la resonancia de forma p, y análisis de la curva BOP
M_J=0 fuera de ella.

Hace tres cosas, todas sobre curvas YA calculadas (no rebarre R):

  1. Compara `plots/archive/bop_curve_MJ0.npz` (p_interpolation="inverse", por defecto)
     con `plots/archive/bop_curve_MJ0_linear.npz` (método antiguo) para verificar que el
     cambio de interpolación NO altera nada fuera de la resonancia.
  2. Define la ventana de exclusión con `ScatteringLengths.p_resonance_window()`
     — centro = polo del interpolante de 1/A_p, anchura = FWHM de |A_p|,
     margen = margin_factor × FWHM — y dice qué puntos caen dentro.
  3. Repite las características de la curva (mínimos locales, profundidades,
     "desplazada > 20 GHz para R <~ 1200 a0") EXCLUYENDO esa ventana, y genera
     la figura final con la región sombreada y rotulada.

Por qué se excluye: el pseudopotencial de Fermi de rango cero DIVERGE en la
resonancia de forma p (A_p -> inf). Es una limitación conocida del modelo, no un
fallo numérico; la vía de corrección (Omont 1977) no está implementada por
decisión explícita. Ver docs/analysis_interpolacion_polo_Ap.md.

Ejecutar:  poetry run python scripts/analyze_resonance_window.py
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_bop_curve as bop   # noqa: E402  (reutiliza base, H y utilidades)

NPZ_INVERSE = "plots/archive/bop_curve_MJ0.npz"
NPZ_LINEAR = "plots/archive/bop_curve_MJ0_linear.npz"
DOC = "docs/analysis_interpolacion_polo_Ap.md"
RULE = "=" * 84


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--margin-factor", type=float, default=2.0,
                    help="semianchura de la ventana en unidades de FWHM")
    ap.add_argument("--ymin", type=float, default=-120.0)
    ap.add_argument("--ymax", type=float, default=25.0)
    ap.add_argument("--out", default="plots/archive/bop_curve_MJ0_excluded.png")
    return ap.parse_args()


# ------------------------------------------------------------------ utilidades
def segment_features(R, y, label):
    """Mínimos/máximos locales estrictos dentro de UN segmento contiguo."""
    mins = [i for i in range(1, len(y) - 1) if y[i] < y[i - 1] and y[i] < y[i + 1]]
    maxs = [i for i in range(1, len(y) - 1) if y[i] > y[i - 1] and y[i] > y[i + 1]]
    rows = []
    for i in mins:
        nxt = [m for m in maxs if m > i]
        depth = (y[nxt[0]] - y[i]) if nxt else float("nan")
        rows.append({"segmento": label, "R": R[i], "E": y[i], "depth": depth})
    return rows


def k_manifold_at(R_edge):
    """Índice de la curva adiabática más baja con carácter de manifold."""
    w, V = bop.solve(float(R_edge))
    return next(k for k in range(len(w))
                if float(np.sum(V[:, k][bop.IS_MANIFOLD] ** 2)) > 0.5)


# ------------------------------------------------------------------ principal
def main():
    args = parse_args()
    print(RULE)
    print("VENTANA DE EXCLUSIÓN DE LA RESONANCIA DE FORMA p  —  curva BOP M_J=0")
    print(RULE)

    # ---------------------------------------------------------- 1. ventana
    win = bop.SYSTEM.p_resonance_window(margin_factor=args.margin_factor)
    HARTREE_MEV = 27211.386245988
    print("\n[1] CRITERIO DE LA VENTANA (todo calculado de rvsAP.dat)\n")
    print("    centro  = cero del interpolante lineal de 1/A_p en eps=k^2/2")
    print("              entre los dos nodos donde A_p cambia de signo")
    print(f"              eps_polo = {win['eps_pole']:.8f} E_h "
          f"= {win['eps_pole']*HARTREE_MEV:.3f} meV")
    print(f"              R_polo   = {win['R_pole']:.3f} a0")
    print("    anchura = FWHM de |A_p|: nodos con |A_p| >= max|A_p|/2")
    print(f"              R in [{win['R_fwhm_lo']:.3f}, {win['R_fwhm_hi']:.3f}] a0"
          f"  ->  FWHM = {win['fwhm']:.3f} a0")
    print(f"    margen  = +/- {args.margin_factor:g} x FWHM "
          f"= +/- {args.margin_factor*win['fwhm']:.3f} a0")
    print(f"\n    VENTANA EXCLUIDA: R in [{win['R_lo']:.2f}, {win['R_hi']:.2f}] a0"
          f"   (anchura {win['R_hi']-win['R_lo']:.2f} a0)")
    print("\n    Por qué NO se ensancha más: la ventana se define sobre A_p, que es")
    print("    donde está la divergencia. La cola energética que la resonancia deja")
    print("    en las curvas se trata aparte, en [2].")

    # ------------------------------------------------------------ datos
    d_inv, d_lin = np.load(NPZ_INVERSE), np.load(NPZ_LINEAR)
    R, R_lin = d_inv["R"], d_lin["R"]
    if not (len(R) == len(R_lin) and np.allclose(R, R_lin)):
        raise SystemExit("las dos mallas en R no coinciden: comparación no directa")
    keep = (R < win["R_lo"]) | (R > win["R_hi"])
    k_man = k_manifold_at(R[-1])
    Y_inv = (d_inv["context"] - bop.E_MAN) * bop.GHZ
    Y_lin = (d_lin["context"] - bop.E_MAN) * bop.GHZ
    y_inv, y_lin = d_inv["y"], d_lin["y"]

    # ------------------------------------- 2. ¿basta con ±margin×FWHM?
    print("\n" + "-" * 84)
    print("[2] DIAGNÓSTICO DEL MARGEN, Y LA COLA BUTTERFLY QUE SOBREVIVE")
    print("-" * 84)
    print("    Qué queda fuera de la ventana según el factor de margen. |Δ| es")
    print("    inverse vs linear sobre las 60 adiabáticas.\n")
    print("    factor | ventana [a0]           | n fuera | min E_k=%d fuera [GHz] | max|Δ| fuera"
          % k_man)
    print("    -------|------------------------|---------|----------------------|-------------")
    Dall = np.abs(Y_inv - Y_lin)
    for f in (1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0):
        lo, hi = win["R_pole"] - f * win["fwhm"], win["R_pole"] + f * win["fwhm"]
        kp = (R < lo) | (R > hi)
        mark = "  <-- adoptado" if abs(f - args.margin_factor) < 1e-9 else ""
        print(f"    {f:6g} | [{lo:8.2f}, {hi:8.2f}] | {kp.sum():7d} | "
              f"{Y_inv[kp, k_man].min():20.2f} | {Dall[kp].max():9.3e}{mark}")

    # dónde vuelve la k_man a caber en la ventana de energía de la figura
    right = R > win["R_hi"]
    off = right & (np.abs(Y_inv[:, k_man]) > abs(args.ymin))
    R_tail = float(R[off].max()) if off.any() else float(win["R_hi"])
    tail = (R > win["R_hi"]) & (R <= R_tail)
    print(f"\n    La ventana quita la DIVERGENCIA, no su cola: justo fuera del borde")
    print(f"    derecho la k={k_man} vale {Y_inv[R > win['R_hi'], k_man][0]:.1f} GHz. "
          f"Sólo vuelve a caber en")
    print(f"    la ventana de energía de la Fig. 2 (|E| <= {abs(args.ymin):.0f} GHz) "
          f"a partir de R = {R_tail:.1f} a0.")
    print(f"    -> BANDA DE COLA: R in ({win['R_hi']:.2f}, {R_tail:.2f}] a0, "
          f"{tail.sum()} puntos.")
    print("       Se dibuja aparte y se excluye TAMBIÉN de la estadística del")
    print("       criterio de 20 GHz, porque ahí la curva k=%d ya no es la del" % k_man)
    print("       manifold sino el estado butterfly que genera la resonancia.")
    keep_strict = keep & ~tail

    print("\n" + "-" * 84)
    print("[3] p_interpolation='inverse' FRENTE A 'linear'  (misma malla de R)")
    print("-" * 84)
    print(f"    {len(R)} puntos; {keep.sum()} fuera de la ventana, "
          f"{(~keep).sum()} dentro")
    print(f"    curva adiabática de carácter manifold en R_max: k = {k_man}\n")
    print("    curva                     | max|Δ| FUERA | mediana|Δ| FUERA | max|Δ| DENTRO")
    print("    --------------------------|--------------|------------------|--------------")
    for name, a, b in (("adiabática k=%d" % k_man, Y_inv[:, k_man], Y_lin[:, k_man]),
                       ("diabática (solapamiento)", y_inv, y_lin)):
        da = np.abs(a - b)
        print(f"    {name:25s} | {da[keep].max():9.3e} GHz | "
              f"{np.median(da[keep]):11.3e} GHz | {da[~keep].max():9.3e} GHz")

    # todas las curvas del contexto, no sólo la k=6
    D_all = np.abs(Y_inv - Y_lin)
    print(f"\n    sobre las {Y_inv.shape[1]} curvas adiabáticas a la vez:")
    print(f"      max|Δ| FUERA de la ventana  = {D_all[keep].max():.3e} GHz")
    print(f"      max|Δ| DENTRO de la ventana = {D_all[~keep].max():.3e} GHz")
    i, j = np.unravel_index(np.argmax(np.where(keep[:, None], D_all, -1)), D_all.shape)
    print(f"      el máximo de fuera está en R = {R[i]:.1f} a0, curva k={j}")
    worst = np.argsort(-D_all[keep].max(axis=1))[:5]
    Rk = R[keep]
    print("      los 5 R con mayor discrepancia fuera de la ventana:")
    for w_ in worst:
        print(f"        R = {Rk[w_]:8.2f} a0   max sobre curvas = "
              f"{D_all[keep][w_].max():.3e} GHz")

    # ---------------------------------- 4. características fuera de ventana
    y6 = Y_inv[:, k_man]
    left = R < win["R_lo"]
    right = R > win["R_hi"]
    print("\n" + "-" * 84)
    print(f"[4] CARACTERÍSTICAS DE LA ADIABÁTICA k={k_man}, EXCLUYENDO LA VENTANA")
    print("-" * 84)
    print("    Los extremos locales se buscan por separado en cada segmento, para")
    print("    que el hueco de la exclusión no cree mínimos ni máximos falsos.")
    print("    † = el mínimo cae en la banda de cola butterfly, no es comparable.\n")
    rows = (segment_features(R[left], y6[left], f"izq  R<{win['R_lo']:.1f}")
            + segment_features(R[right], y6[right], f"der  R>{win['R_hi']:.1f}"))
    print("     #  | segmento              |  R_min [a0] |  E-E_man [GHz] | prof. vs máx. der. [GHz] | ΔR al anterior")
    print("    ----|-----------------------|-------------|----------------|--------------------------|---------------")
    prev = None
    prev_seg = None
    seps = []
    for n, r in enumerate(rows):
        same = prev is not None and r["segmento"] == prev_seg
        dR = (r["R"] - prev) if prev is not None else float("nan")
        if same:
            seps.append(dR)
        flag = "†" if win["R_hi"] < r["R"] <= R_tail else " "
        prev, prev_seg = r["R"], r["segmento"]
        dep = f"{r['depth']:24.4f}" if np.isfinite(r["depth"]) else " " * 20 + "---"
        sep_txt = (f"{dR:13.2f}{'*' if not same else ' '}"
                   if np.isfinite(dR) else "           ---")
        print(f"    {n:3d}{flag}| {r['segmento']:21s} | {r['R']:11.2f} | "
              f"{r['E']:14.4f} | {dep} | " + sep_txt)
    print("    (* = salto que atraviesa la ventana excluida; no cuenta como "
          "separación física)")
    if len(seps) >= 2:
        s = np.array(seps)
        print(f"\n    separación entre mínimos CONTIGUOS (sin los saltos marcados *): "
              f"media {s.mean():.2f} a0, sigma {s.std():.2f} a0, "
              f"min {s.min():.2f}, max {s.max():.2f}")
        cv = s.std() / s.mean()
        print(f"    coeficiente de variación = {cv:.3f}  "
              f"({'aprox. regular' if cv < 0.25 else 'NO regular'})")
    elif seps:
        print(f"\n    sólo queda {len(seps)} separación entre mínimos contiguos "
              f"({seps[0]:.2f} a0): con la ventana fuera no hay estadística "
              "de regularidad que sostener")

    print(f"\n    rango de la curva fuera de la ventana: "
          f"[{y6[keep].min():.3f}, {y6[keep].max():.3f}] GHz")
    print(f"    mínimo global fuera de la ventana: {y6[keep].min():.4f} GHz "
          f"en R = {R[keep][np.argmin(y6[keep])]:.2f} a0  (es cola butterfly)")
    print(f"    mínimo global fuera de ventana Y de cola: {y6[keep_strict].min():.4f} "
          f"GHz en R = {R[keep_strict][np.argmin(y6[keep_strict])]:.2f} a0")
    print(f"    (dentro de la ventana la curva llega a {y6[~keep].min():.1f} GHz, "
          "valor SIN significado físico)")

    print("\n    Texto del paper: 'shifted more than 20 GHz for R <~ 1200 a0' (M_J=0)")
    print("    (nuestro dominio llega sólo a R = %.1f a0)" % R[-1])
    for name, yy in ((f"adiabática k={k_man}", y6), ("diabática", y_inv)):
        f_all = float(np.mean(np.abs(yy) > 20.0))
        f_out = float(np.mean(np.abs(yy[keep]) > 20.0))
        f_str = float(np.mean(np.abs(yy[keep_strict]) > 20.0))
        below = R[keep_strict][np.abs(yy[keep_strict]) <= 20.0]
        print(f"      {name:20s}: |E-E_man| > 20 GHz en")
        print(f"      {'':20s}    {f_all*100:6.2f}% sin excluir nada")
        print(f"      {'':20s}    {f_out*100:6.2f}% excluida la ventana")
        print(f"      {'':20s}    {f_str*100:6.2f}% excluidas ventana + cola butterfly")
        if len(below):
            print(f"      {'':20s}  los {len(below)} puntos por debajo de 20 GHz "
                  f"están en R in [{below.min():.1f}, {below.max():.1f}] a0")

    # ------------------------------------------------------------ 5. figura
    ymask = np.where(keep[:, None], Y_inv, np.nan)
    fig, axes = plt.subplots(2, 1, figsize=(9.5, 9.5), sharex=True)
    band_label = ("región excluida: divergencia conocida del pseudopotencial de\n"
                  "rango cero en la resonancia de forma p — ver\n" + DOC)
    for panel, (ax, ylim, ttl, masked) in enumerate((
        (axes[0], (float(np.nanmin(Y_inv[:, k_man])) * 1.05, 100.0),
         "Rango completo — se MUESTRA la divergencia que luego se excluye", False),
        (axes[1], (args.ymin, args.ymax),
         "Ventana de energía comparable a la Fig. 2 — región excluida en gris", True),
    )):
        ax.axvspan(win["R_lo"], win["R_hi"], color="0.45", alpha=0.32, zorder=0,
                   label=band_label if panel == 0 else None)
        ax.axvspan(win["R_hi"], R_tail, color="0.45", alpha=0.13, hatch="///",
                   ec="0.6", lw=0.0, zorder=0,
                   label=(f"cola butterfly, no comparable "
                          f"($R\\leq{R_tail:.0f}\\,a_0$)" if panel == 0 else None))
        ax.axvline(win["R_pole"], color="0.35", lw=1.0, ls=":", zorder=1,
                   label=(f"polo de $A_p$, $R={win['R_pole']:.1f}\\,a_0$"
                          if panel == 0 else None))
        src = ymask if masked else Y_inv
        for c in range(Y_inv.shape[1]):
            ax.plot(R, src[:, c], color="0.78", lw=0.6, zorder=1)
        ax.plot(R, src[:, k_man], color="crimson", lw=3.0, zorder=4,
                label=f"adiabática k={k_man} (carácter manifold $n=24$)")
        ax.plot(R, np.where(keep, y_inv, np.nan) if masked else y_inv,
                color="darkblue", lw=1.5, ls="--", zorder=3,
                label="seguida por solapamiento (diabática)")
        ax.axhline(0, color="k", lw=0.8, ls=":", zorder=2)
        ax.axhline(-20, color="seagreen", lw=1.2, ls="--", zorder=2,
                   label="$-20$ GHz (texto del paper)")
        ax.set_ylim(*ylim)
        ax.set_xlim(R[0], R[-1])
        ax.grid(alpha=0.25)
        ax.set_ylabel(r"$E - E_{n=24}$  [GHz]")
        ax.set_title(ttl, fontsize=10)
        ax.legend(loc="lower right", fontsize=7.5)
    for N in range(7):
        axes[1].axhline(-63.40 + 1.114 * N * (N + 1), color="tab:orange", lw=0.7, alpha=0.6)
        axes[1].axhline(1.114 * N * (N + 1), color="tab:blue", lw=0.7, alpha=0.6)
    axes[1].annotate(
        "REGIÓN EXCLUIDA\ndivergencia conocida del pseudopotencial\n"
        "de rango cero en la resonancia de forma p\n"
        f"(rayado: cola butterfly, tampoco comparable)\nver {DOC}",
        xy=(win["R_pole"], 0.80 * args.ymin), xycoords="data",
        xytext=(R[0] + 25, 0.78 * args.ymin), textcoords="data",
        fontsize=7.5, ha="left", va="center", color="0.10",
        bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="0.5", alpha=0.92),
        arrowprops=dict(arrowstyle="->", color="0.25", lw=1.1))
    axes[1].set_xlabel(r"$R$  [$a_0$]")
    fig.suptitle(
        "Rb*-KRb, curvas BOP $M_J=0$   "
        r"$H=H_a+B\mathbf{N}^2-\mathbf{d}\cdot\mathbf{F}_{ryd}+V_{Fermi}$"
        f"\ninterpolación de $1/A_p$ · ventana excluida "
        f"$R\\in[{win['R_lo']:.1f},\\,{win['R_hi']:.1f}]\\,a_0$ "
        f"= $R_{{polo}}\\pm{args.margin_factor:g}\\times$FWHM$({win['fwhm']:.1f}\\,a_0)$",
        fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig(args.out, dpi=150)
    print(f"\n    PNG guardado en {args.out}")
    print(RULE)


if __name__ == "__main__":
    main()
