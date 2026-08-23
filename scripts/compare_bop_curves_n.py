#!/usr/bin/env python3
"""
Figura comparativa de curvas BOP polares para varios n.

Fase A del plan docs/PLAN_figuras_publicacion.md. NO recalcula nada: lee los
.npz que ya produjo scripts/compute_bop_curve.py para cada n (mismo formato,
mismo cero de energía relativo a cada manifold) y los superpone en un único
eje, en unidades reducidas R/(2n^2) para poder comparar la forma del pozo
entre manifolds de tamaño muy distinto.

USO
---
    poetry run python scripts/compare_bop_curves_n.py --n 24 25 26 27
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from trimero.systems.polar_molecule import MOLECULES, get_molecule


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n", type=int, nargs="+", default=[24, 25, 26, 27])
    ap.add_argument("--molecule", choices=tuple(MOLECULES), default="krb")
    ap.add_argument("--mj", type=int, default=0)
    ap.add_argument("--npz-dir", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--ymin", type=float, default=-25.0)
    ap.add_argument("--ymax", type=float, default=1.0)
    return ap.parse_args()


def local_minima(y):
    return [i for i in range(1, len(y) - 1) if y[i] < y[i - 1] and y[i] < y[i + 1]]


def main():
    args = parse_args()
    molecule = get_molecule(args.molecule)
    root = f"plots/rb_{molecule.key}_polar"
    args.npz_dir = args.npz_dir or os.path.join(root, "data")
    args.out = args.out or os.path.join(root, "figures", f"fig_compare_n_MJ{args.mj}.png")
    colors = plt.cm.viridis(np.linspace(0.1, 0.85, len(args.n)))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.6))

    summary = {}
    for n, c in zip(args.n, colors):
        npz = os.path.join(args.npz_dir, f"fig1_ad_MJ{args.mj}_n{n}.npz")
        if not os.path.exists(npz):
            raise FileNotFoundError(
                f"{npz} no existe. Corre antes: "
                f"scripts/compute_bop_curve.py --n-manifold {n} --mj {args.mj}")
        d = np.load(npz)
        if "molecule" in d.files and str(d["molecule"]) != molecule.key:
            raise ValueError(f"{npz} pertenece a {d['molecule']}, no a {molecule.key}")
        R, E = d["R"], d["E"]
        i_deep = int(np.nanargmin(E))
        mins = local_minima(E)
        summary[n] = {"deep_E": float(E[i_deep]), "deep_R": float(R[i_deep]),
                      "n_min": len(mins), "E_end": float(E[-1])}

        ax1.plot(R, E, color=c, lw=1.8, label=f"n={n}")
        ax2.plot(R / (2 * n * n), E, color=c, lw=1.8, label=f"n={n}")

    for ax, xlabel in ((ax1, r"$R$  [$a_0$]"), (ax2, r"$R / 2n^2$")):
        ax.axhline(0.0, color="k", lw=0.9, ls=":", zorder=0)
        ax.set_xlabel(xlabel)
        ax.set_ylim(args.ymin, args.ymax)
        ax.grid(alpha=0.25)
        ax.legend(fontsize=9, frameon=False)
    ax1.set_ylabel(r"$V(R) = E - E_{n,\,l\geq3}$  [GHz]")
    ax1.set_title("escala absoluta")
    ax2.set_title("escala reducida ($R/2n^2$)")

    fig.suptitle(
        rf"Rb*-{molecule.label}, curvas BOP $M_J=%d$ para varios $n$   —   "
        r"$H_{ad} = H_A + H_{mol}$  (Ec. 1, Aguilera-Fernández et al. 2015)"
        % args.mj, fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    fig.savefig(args.out, dpi=150)
    print(f"PNG guardado en {args.out}\n")

    print(f"{'n':>4}  {'pozo mas profundo (GHz)':>24}  {'R pozo (a0)':>12}  "
          f"{'minimos locales':>16}  {'E(R=1800) (GHz)':>17}")
    for n in args.n:
        s = summary[n]
        print(f"{n:>4}  {s['deep_E']:>24.4f}  {s['deep_R']:>12.1f}  "
              f"{s['n_min']:>16d}  {s['E_end']:>17.4f}")
    return summary


if __name__ == "__main__":
    main()
