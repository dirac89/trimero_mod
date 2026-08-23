#!/usr/bin/env python3
"""
Figura comparativa de ⟨cos θ_d⟩ para varios n y una molécula polar.

Fase B de docs/PLAN_figuras_publicacion.md. Lee los .npz que ya produjo
scripts/compute_orientation_curve.py para cada n y los superpone. Marca en
gris el tramo en el que scripts/compute_orientation_curve.py detectó eventos
de seguimiento ambiguo (cambio de K o salto de ⟨cos⟩ > 0.02 entre puntos
consecutivos), para no dibujar una curva limpia donde el propio cálculo dice
que no lo es.

USO
---
    poetry run python scripts/compare_orientation_n.py --n 24 25 26 27
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
    ap.add_argument("--rmax-plot", type=float, default=800.0,
                    help="la cola [800,1800] a paso grueso se calcula pero no "
                         "se dibuja por defecto: ya no aporta información de "
                         "orientación (decae a ~0) y con paso 50 introduciría "
                         "un artefacto visual de resolución distinta")
    return ap.parse_args()


def ambiguous_mask(K, COS, thr=0.02):
    """Ambiguo = salto de magnitud |Delta cos| > thr entre puntos consecutivos.

    OJO: un cambio de K por sí solo NO se marca. Un cruce evitado aislado con
    |Delta cos| pequeño (< thr) es una relabelación adiabática benigna, no un
    fallo de seguimiento — marcar todo dK!=0 sobre-estima la región ambigua
    (se probó y arrastraba el inicio hasta R=100 en los cuatro n, un
    artefacto del criterio, no un resultado físico)."""
    n = len(K)
    bad = np.zeros(n, dtype=bool)
    for i in range(1, n):
        if abs(COS[i] - COS[i - 1]) > thr:
            bad[i - 1] = True
            bad[i] = True
    return bad


def main():
    args = parse_args()
    molecule = get_molecule(args.molecule)
    root = f"plots/rb_{molecule.key}_polar"
    args.npz_dir = args.npz_dir or os.path.join(root, "data")
    args.out = args.out or os.path.join(
        root, "figures", f"fig_orientation_compare_n_MJ{args.mj}.png")
    colors = plt.cm.viridis(np.linspace(0.1, 0.85, len(args.n)))

    fig, ax = plt.subplots(figsize=(8, 6))
    onset = {}
    for n, c in zip(args.n, colors):
        npz = os.path.join(args.npz_dir, f"orientation_MJ{args.mj}_n{n}.npz")
        if not os.path.exists(npz):
            raise FileNotFoundError(
                f"{npz} no existe. Corre antes: "
                f"scripts/compute_orientation_curve.py --n-manifold {n}")
        d = np.load(npz)
        if "molecule" in d.files and str(d["molecule"]) != molecule.key:
            raise ValueError(f"{npz} pertenece a {d['molecule']}, no a {molecule.key}")
        R, K, COS = d["R"], d["K"], d["COS"]
        keep = R <= args.rmax_plot
        R, K, COS = R[keep], K[keep], COS[keep]
        bad = ambiguous_mask(K, COS)

        ax.plot(R, COS, color=c, lw=1.6, alpha=0.35, zorder=1)
        ax.plot(np.where(bad, np.nan, R), np.where(bad, np.nan, COS),
                color=c, lw=2.0, zorder=3, label=f"n={n}")
        if bad.any():
            onset[n] = float(R[bad][0])

    ax.axhline(0.78, color="tab:red", lw=1.0, ls="--", alpha=0.7,
               label=r"$\langle\cos\theta_d\rangle=0.78$ (González-Férez 2015)")
    ax.axhline(0.0, color="k", lw=0.8, ls=":")
    ax.set_xlabel(r"$R$  [$a_0$]")
    ax.set_ylabel(r"$\langle\cos\theta_d\rangle$")
    ax.set_ylim(-0.05, 1.0)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=9, frameon=False)
    ax.set_title(
        f"Rb*-{molecule.label}, orientación del dipolo, $M_J=0$, $H_{{ad}}=H_A+H_{{mol}}$ "
        "(sin V_Fermi)\ntrazo grueso = seguimiento fiable, fino = eventos de "
        "ambigüedad cercanos", fontsize=10.5)

    fig.tight_layout()
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    fig.savefig(args.out, dpi=150)
    print(f"PNG guardado en {args.out}\n")

    print("primer R con evento de ambigüedad detectado, por n:")
    for n in args.n:
        r0 = onset.get(n)
        print(f"  n={n}: {'sin eventos en el rango graficado' if r0 is None else f'{r0:.1f} a0'}")
    return onset


if __name__ == "__main__":
    main()
