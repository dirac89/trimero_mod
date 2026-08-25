#!/usr/bin/env python3
"""
Figura comparativa de ⟨cos θ_d⟩ para varios n y una molécula polar.

Fase B de docs/PLAN_figuras_publicacion.md. Lee los .npz que ya produjo
scripts/compute_orientation_curve.py para cada n y los superpone. Marca en
gris el tramo en el que scripts/compute_orientation_curve.py detectó eventos
de seguimiento ambiguo (cambio de K o salto de ⟨cos⟩ > 0.02 entre puntos
consecutivos), para no dibujar una curva limpia donde el propio cálculo dice
que no lo es.

`--molecule` es OBLIGATORIO: el default silencioso `krb` es lo que hizo que las
Fases A y B salieran con KRb creyéndose RbCs. Además se comprueba contra el
campo `molecule` de cada .npz, y los .npz sin ese campo (esquema anterior al
registro de moléculas) se rechazan en vez de suponer que son de la molécula
pedida. Ver `docs/analysis_faseA_curvas_bop_varios_n.md` §0.

USO
---
    poetry run python scripts/compare_orientation_n.py --molecule rbcs --n 24 25 26 27
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
    ap.add_argument("--molecule", choices=tuple(MOLECULES), required=True,
                    help="molécula polar (OBLIGATORIO, sin valor por defecto)")
    ap.add_argument("--mj", type=int, default=0)
    ap.add_argument("--npz-dir", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--ymin", type=float, default=-1.02,
                    help="con KRb bastaba -0.05; con RbCs la curva se vuelca a "
                         "valores negativos (dipolo antiparalelo) y recortarla "
                         "escondería el resultado")
    ap.add_argument("--ymax", type=float, default=1.02)
    ap.add_argument("--rmax-plot", type=float, default=800.0,
                    help="R máximo dibujado. El default 800 viene de la ronda de "
                         "KRb, donde la cola se calculaba a paso 50 y ⟨cos⟩ ya "
                         "había decaído. Con RbCs NO vale: sigue orientada más "
                         "allá de 800 a0 (Fase A §6.1), así que hay que subirlo "
                         "hasta donde llegue la malla fina")
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

    # Dos paneles, misma decisión que la Fase A: el absoluto por sí solo no
    # sostiene la afirmación de escala; el reducido R/2n^2 sí.
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(13.2, 5.8), sharey=True)
    onset = {}
    grids = []
    flip = {}
    for n, c in zip(args.n, colors):
        npz = os.path.join(args.npz_dir, f"orientation_MJ{args.mj}_n{n}.npz")
        if not os.path.exists(npz):
            raise FileNotFoundError(
                f"{npz} no existe. Corre antes: "
                f"scripts/compute_orientation_curve.py "
                f"--molecule {molecule.key} --n-manifold {n}")
        d = np.load(npz)
        if "molecule" not in d.files:
            raise ValueError(
                f"{npz} no declara su molécula (esquema anterior al registro de "
                f"moléculas). No se puede suponer que sea {molecule.key}: "
                f"recalcúlalo con compute_orientation_curve.py "
                f"--molecule {molecule.key} --n-manifold {n}.")
        if str(d["molecule"]) != molecule.key:
            raise ValueError(f"{npz} pertenece a {d['molecule']}, no a {molecule.key}")
        R, K, COS = d["R"], d["K"], d["COS"]
        # Malla homogénea entre n: los .npz de RbCs del primer intento mezclaban
        # dos mallas distintas (unos [100,1800] con cola a paso 50, otro
        # [400,1800] a paso 5), y superponerlas dibuja curvas que no son
        # comparables punto a punto.
        grid = (float(R[0]), float(R[-1]), len(R))
        if grids and grid != grids[-1][1]:
            raise ValueError(
                f"malla incompatible: n={n} tiene R=[{grid[0]:.0f},{grid[1]:.0f}] "
                f"con {grid[2]} puntos, y n={grids[-1][0]} tiene "
                f"R=[{grids[-1][1][0]:.0f},{grids[-1][1][1]:.0f}] con "
                f"{grids[-1][1][2]} puntos. Recalcula con los mismos parámetros.")
        grids.append((n, grid))
        keep = R <= args.rmax_plot
        R, K, COS = R[keep], K[keep], COS[keep]
        bad = ambiguous_mask(K, COS)

        ax.plot(R, COS, color=c, lw=1.6, alpha=0.35, zorder=1)
        ax.plot(np.where(bad, np.nan, R), np.where(bad, np.nan, COS),
                color=c, lw=2.0, zorder=3, label=f"n={n}")
        ax2.plot(R / (2 * n * n), COS, color=c, lw=1.6, alpha=0.35, zorder=1)
        ax2.plot(np.where(bad, np.nan, R / (2 * n * n)), np.where(bad, np.nan, COS),
                 color=c, lw=2.0, zorder=3, label=f"n={n}")
        if bad.any():
            onset[n] = float(R[bad][0])
        zero = np.flatnonzero((COS[:-1] > 0) & (COS[1:] <= 0))
        if len(zero):
            flip[n] = float(R[zero[0] + 1])

    for a in (ax, ax2):
        a.axhline(0.78, color="tab:red", lw=1.0, ls="--", alpha=0.7,
                  label=r"$\langle\cos\theta_d\rangle=0.78$ (González-Férez 2015)")
        a.axhline(0.0, color="k", lw=0.8, ls=":")
        a.set_ylim(args.ymin, args.ymax)
        a.grid(alpha=0.25)
    ax.set_xlabel(r"$R$  [$a_0$]")
    ax2.set_xlabel(r"$R / 2n^2$")
    ax.set_ylabel(r"$\langle\cos\theta_d\rangle$")
    ax.set_title("escala absoluta", fontsize=10)
    ax2.set_title(r"escala reducida ($R/2n^2$)", fontsize=10)
    ax.legend(fontsize=9, frameon=False, loc="lower left")
    fig.suptitle(
        f"Rb*-{molecule.label}, orientación del dipolo, $M_J=0$, $H_{{ad}}=H_A+H_{{mol}}$ "
        "(sin V_Fermi)\n"
        f"{molecule.label}:  $B$ = {molecule.B_ghz:.6f} GHz,  "
        f"$d$ = {molecule.dipole_debye:.3f} D   —   "
        "trazo grueso = seguimiento fiable, fino = eventos de "
        "ambigüedad cercanos", fontsize=10.5)

    fig.tight_layout(rect=[0, 0, 1, 0.90])
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    fig.savefig(args.out, dpi=150)
    print(f"PNG guardado en {args.out}\n")

    print("primer R con evento de ambigüedad detectado, por n:")
    for n in args.n:
        r0 = onset.get(n)
        print(f"  n={n}: {'sin eventos en el rango graficado' if r0 is None else f'{r0:.1f} a0'}")
    print("\nvuelco de signo de <cos theta_d> (dipolo paralelo -> antiparalelo):")
    for n in args.n:
        r0 = flip.get(n)
        if r0 is None:
            print(f"  n={n}: no se observa en el rango graficado")
        else:
            print(f"  n={n}: R = {r0:.0f} a0   R/2n^2 = {r0/(2*n*n):.4f}")
    return onset


if __name__ == "__main__":
    main()
