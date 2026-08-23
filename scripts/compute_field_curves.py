#!/usr/bin/env python3
"""Curvas BOP de Rb*+molécula polar para varios campos DC paralelos a Z."""

import argparse
from pathlib import Path
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from trimero.systems.polar_molecule import MOLECULES, get_molecule
from trimero.systems.polar_rydberg import GHZ_PER_HARTREE, PolarBOPSystem
from trimero.systems.rb_krb_polar.rb_defects import DELTA0_NS_PAPER


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--molecule", choices=tuple(MOLECULES), default="rbcs")
    parser.add_argument("--n-manifold", type=int, default=25)
    parser.add_argument("--n-max", type=int, default=6)
    parser.add_argument("--mj", type=int, default=0)
    parser.add_argument("--fields", type=float, nargs="+", default=[0, 100, 300, 500])
    parser.add_argument("--rmin", type=float, default=400.0)
    parser.add_argument("--rmax", type=float, default=1800.0)
    parser.add_argument("--step", type=float, default=5.0)
    parser.add_argument("--weight", type=float, default=0.5)
    parser.add_argument("--keep", type=int, default=120)
    parser.add_argument("--out-root", type=Path, default=None)
    parser.add_argument("--ymin", type=float, default=-70.0)
    parser.add_argument("--ymax", type=float, default=5.0)
    return parser.parse_args()


def main():
    args = parse_args()
    molecule = get_molecule(args.molecule)
    root = args.out_root or Path(f"plots/rb_{molecule.key}_polar")
    data_dir, figures_dir = root / "data", root / "figures"
    data_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    system = PolarBOPSystem(
        molecule=molecule, n_manifold=args.n_manifold, N_max=args.n_max,
        delta0_ns=DELTA0_NS_PAPER,
    )
    R = np.arange(args.rmin, args.rmax + 1e-9, args.step)
    mask = system.is_manifold(args.mj)
    field_matrices = {
        field: system.external_field_matrix(args.mj, field) for field in args.fields
    }
    results = {}

    for field in args.fields:
        started = time.perf_counter()
        reference_path = data_dir / f"fig1_ad_MJ{args.mj}_n{args.n_manifold}.npz"
        if field == 0.0 and reference_path.exists():
            with np.load(reference_path) as reference:
                compatible = (
                    np.array_equal(reference["R"], R)
                    and str(reference["molecule"]) == molecule.key
                    and int(reference["N_max"]) == args.n_max
                    and reference["spectrum"].shape[1] >= min(args.keep, len(system.block(args.mj)))
                )
                if compatible:
                    keep = min(args.keep, reference["spectrum"].shape[1])
                    results[field] = {
                        "R": reference["R"].copy(), "E": reference["E"].copy(),
                        "K": reference["K"].copy(), "W": reference["W"].copy(),
                        "spectrum": reference["spectrum"][:, :keep].copy(),
                    }
            if field in results:
                print(f"F=0: reutilizando {reference_path}")
                continue
        E = np.full(len(R), np.nan)
        K = np.full(len(R), -1, dtype=int)
        W = np.zeros(len(R))
        keep = min(args.keep, len(system.block(args.mj)))
        spectrum = np.empty((len(R), keep))
        Vext = field_matrices[field]
        for i, radius in enumerate(R):
            H = system.hamiltonian(float(radius), args.mj)
            if field != 0.0:
                H = H + Vext
            values, vectors = np.linalg.eigh(H)
            spectrum[i] = (values[:keep] - system.E_manifold) * GHZ_PER_HARTREE
            for k, value in enumerate(values):
                weight = float(np.sum(vectors[:, k][mask] ** 2))
                if weight > args.weight:
                    E[i] = (value - system.E_manifold) * GHZ_PER_HARTREE
                    K[i], W[i] = k, weight
                    break
            if i % 40 == 0:
                print(f"F={field:g} V/m: R={radius:.0f} ({i+1}/{len(R)})", flush=True)
        results[field] = dict(R=R, E=E, K=K, W=W, spectrum=spectrum)
        tag = f"field_F{field:g}_MJ{args.mj}_n{args.n_manifold}_Nmax{args.n_max}.npz"
        np.savez(
            data_dir / tag, **results[field], molecule=molecule.key,
            field_v_per_m=field, n_manifold=args.n_manifold, N_max=args.n_max,
            M_J=args.mj, character_weight=args.weight, schema_version=1,
        )
        print(f"F={field:g}: {time.perf_counter()-started:.1f}s, "
              f"E=[{np.nanmin(E):.4f},{np.nanmax(E):.4f}] GHz")

    ncols = 2
    nrows = int(np.ceil(len(args.fields) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(12.5, 4.8 * nrows),
                             sharex=True, sharey=True, squeeze=False)
    colors = plt.cm.plasma(np.linspace(0.08, 0.88, len(args.fields)))
    for ax, field, color in zip(axes.flat, args.fields, colors):
        data = results[field]
        visible = np.any(
            (data["spectrum"] >= args.ymin) & (data["spectrum"] <= args.ymax), axis=0
        )
        ax.plot(R, data["spectrum"][:, visible], color="0.72", lw=0.45, alpha=0.8)
        ax.plot(R, data["E"], color=color, lw=2.4,
                label="curva principal de carácter manifold")
        ax.axhline(0.0, color="k", lw=0.8, ls=":")
        ax.set(xlim=(args.rmin, args.rmax), ylim=(args.ymin, args.ymax),
               title=f"F = {field:g} V/m")
        ax.grid(alpha=0.25)
        ax.legend(frameon=False, fontsize=8, loc="lower right")
    for ax in axes.flat[len(args.fields):]:
        ax.set_visible(False)
    for ax in axes[-1]:
        if ax.get_visible():
            ax.set_xlabel(r"$R$ [$a_0$]")
    for row in axes:
        row[0].set_ylabel(r"$E-E_{n,l\geq3}$ [GHz]")
    fig.suptitle(f"Rb*-{molecule.label}, n={args.n_manifold}, $M_J={args.mj}$, "
                 f"$N_{{max}}={args.n_max}$ — campo DC paralelo a Z",
                 fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    png = figures_dir / f"field_comparison_MJ{args.mj}_n{args.n_manifold}_Nmax{args.n_max}.png"
    fig.savefig(png, dpi=150)
    print(f"Figura: {png}")
    return results


if __name__ == "__main__":
    main()
