#!/usr/bin/env python
"""
Curvas de potencial adiabático del trímero Rydberg LINEAL SIMÉTRICO en campo DC.

Aguilera-Fernández, Schmelcher & González-Férez, J. Phys. B 49, 124002 (2016).
Reproduce las Figs. 3 (F=0) y 4/5 (F ≠ 0): Rb(5s)–Rb*(n=35, l≥3)–Rb(5s) con los
dos átomos neutros en θ=0 y θ=π a la misma distancia R del core.

    poetry run python scripts/compute_trimer_curves.py                # Σ, todas las F
    poetry run python scripts/compute_trimer_curves.py --symmetry Pi
    poetry run python scripts/compute_trimer_curves.py --fields 0 500 --s-wave-only

Escribe `plots/rb_neutral_perturber/trimer_lineal_<SIM>_n<n>.npz` y su PNG. El .npz lleva `R`, `E`
(GHz, relativa al manifold Rb(n,l≥3) sin campo), `W` (peso de manifold de cada
autovector) y `fields_V_per_m`.

⚠️ Este es el sistema del PERTURBADOR NEUTRO, no el polar Rb*-KRb. Ver
`.claude/CLAUDE.md`: son dos sistemas físicos distintos y no se mezclan.
"""

import argparse
from pathlib import Path

import numpy as np

from trimero.systems.rb_neutral_perturber.linear_trimer import (
    SymmetricLinearTrimer,
    field_au,
)

SYMMETRY_M = {"Sigma": 0, "Pi": 1, "Delta": 2}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--n-manifold", type=int, default=35,
                   help="n del manifold degenerado (por defecto 35, el del paper)")
    p.add_argument("--symmetry", choices=sorted(SYMMETRY_M), default="Sigma",
                   help="simetría molecular: Σ=m_l 0, Π=|m_l| 1, Δ=|m_l| 2")
    p.add_argument("--fields", type=float, nargs="+", default=[0.0, 100.0, 300.0, 500.0],
                   help="intensidades de campo DC en V/m")
    p.add_argument("--rmin", type=float, default=1000.0)
    p.add_argument("--rmax", type=float, default=None,
                   help="por defecto el último nodo de la tabla (2448 a₀)")
    p.add_argument("--s-wave-only", action="store_true",
                   help="apaga la onda p (comparación con la Fig. 3(a))")
    p.add_argument("--dimer", action="store_true",
                   help="un solo perturbador en θ=0: la curva de referencia")
    p.add_argument("--weight", type=float, default=0.5,
                   help="peso de manifold mínimo para dibujar una curva")
    p.add_argument("--outdir", type=Path, default=Path("plots/rb_neutral_perturber"),
                   help="directorio de salida (por defecto plots/rb_neutral_perturber)")
    p.add_argument("--no-plot", action="store_true")
    return p.parse_args()


def main():
    args = parse_args()
    m_l = SYMMETRY_M[args.symmetry]

    trimer = SymmetricLinearTrimer(
        n_manifold=args.n_manifold,
        p_wave=not args.s_wave_only,
        n_perturbers=1 if args.dimer else 2,
    )
    R = trimer.table_R(R_min=args.rmin, R_max=args.rmax)
    if R.size == 0:
        raise SystemExit(f"ningún nodo de la tabla en [{args.rmin}, {args.rmax}]")

    print(f"n={args.n_manifold}  {args.symmetry} (m_l={m_l})  "
          f"{'dímero' if args.dimer else 'trímero simétrico'}  "
          f"{'sólo onda s' if args.s_wave_only else 'ondas s+p'}")
    print(f"{len(R)} nodos nativos de la tabla, R = {R[0]:.0f} .. {R[-1]:.0f} a₀ "
          f"(sin interpolar A_s ni A_p)")

    E, W = [], []
    for F in args.fields:
        e, w = trimer.spectra(R, field_au(F), m_l)
        E.append(e)
        W.append(w)
        vis = np.where(w > args.weight, e, np.nan)
        print(f"  F = {F:6.1f} V/m ({field_au(F):.4e} u.a.):  "
              f"mín(manifold) = {np.nanmin(vis):8.3f} GHz")
    E, W = np.array(E), np.array(W)

    args.outdir.mkdir(parents=True, exist_ok=True)
    tag = (f"trimer_lineal_{args.symmetry}_n{args.n_manifold}"
           f"{'_swave' if args.s_wave_only else ''}{'_dimer' if args.dimer else ''}")
    npz = args.outdir / f"{tag}.npz"
    np.savez_compressed(npz, R=R, E=E, W=W,
                        fields_V_per_m=np.array(args.fields),
                        n_manifold=args.n_manifold, m_l=m_l,
                        s_wave_only=args.s_wave_only, dimer=args.dimer)
    print(f"-> {npz}")

    if not args.no_plot:
        plot(args, trimer, R, E, W, tag)


def plot(args, trimer, R, E, W, tag):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(args.fields)
    fig, axes = plt.subplots(n, 1, figsize=(7.0, 2.6 * n), sharex=True, squeeze=False)
    for k, F in enumerate(args.fields):
        ax = axes[k][0]
        vis = np.where(W[k] > args.weight, E[k], np.nan)
        ax.plot(R / 1000.0, vis, lw=0.8, color="C0")
        # Los niveles vecinos (38s, 37p, 36d) con su propia etiqueta: son los
        # estados Rb(5s)Rb(n'l')Rb(5s) que el paper dibuja aparte.
        other = np.where(W[k] <= args.weight, E[k], np.nan)
        ax.plot(R / 1000.0, other, lw=0.8, color="C3", ls="--")
        ax.axhline(0.0, color="0.7", lw=0.6)
        ax.set_ylabel(r"$\epsilon(R,R)$ [GHz]")
        ax.text(0.02, 0.06, f"F = {F:g} V/m", transform=ax.transAxes)
        lo = np.nanmin(vis)
        ax.set_ylim(max(lo * 1.05, -40.0), 2.0)
    axes[-1][0].set_xlabel(r"$R$ [$10^3\,a_0$]")
    axes[0][0].set_title(
        f"Trímero lineal simétrico Rb(5s)Rb({args.n_manifold},l≥3)Rb(5s) — "
        f"{args.symmetry}"
        + ("  (sólo onda s)" if args.s_wave_only else "  (ondas s+p)")
    )
    fig.tight_layout()
    png = args.outdir / f"{tag}.png"
    fig.savefig(png, dpi=150)
    print(f"-> {png}")


if __name__ == "__main__":
    main()
