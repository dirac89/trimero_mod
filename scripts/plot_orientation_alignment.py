#!/usr/bin/env python3
"""Figura tipo Fig. 3: orientación y alineamiento de Rb*+RbCs."""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from compute_orientation_curve import cos_theta_matrix, cos2_theta_matrix
from trimero.systems.polar_molecule import get_molecule
from trimero.systems.polar_rydberg import PolarBOPSystem
from trimero.systems.rb_krb_polar.rb_defects import DELTA0_NS_PAPER


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--n", type=int, nargs=2, default=[25, 29])
    p.add_argument("--molecule", default="rbcs", choices=("krb", "rbcs"))
    p.add_argument("--mj", type=int, default=0)
    p.add_argument("--n-max", type=int, default=6)
    p.add_argument("--rmin", type=float, default=400.0)
    p.add_argument("--rmax", type=float, default=1800.0)
    p.add_argument("--npz-dir", type=Path, default=None)
    p.add_argument("--out", type=Path, default=None)
    return p.parse_args()


def isolated_rotor_curves(molecule, R, sign, n_rot=30):
    """Estado fundamental del rotor en el campo uniforme E=sign/R²."""
    template = PolarBOPSystem(molecule=molecule, N_max=1,
                              delta0_ns=DELTA0_NS_PAPER)
    N = np.arange(n_rot + 1)
    C = np.zeros((n_rot + 1, n_rot + 1))
    for n in N:
        if n < n_rot:
            C[n, n + 1] = template.hmol.cos_theta_element(n, 0, n + 1, 0)
            C[n + 1, n] = C[n, n + 1]
    C2 = C @ C
    orientation = np.empty_like(R)
    alignment = np.empty_like(R)
    for i, radius in enumerate(R):
        H = np.diag(molecule.B_au * N * (N + 1)) - molecule.d_au * sign * C / radius**2
        _, vectors = np.linalg.eigh(H)
        v = vectors[:, 0]
        orientation[i] = v @ C @ v
        alignment[i] = v @ C2 @ v
    return orientation, alignment


def main():
    args = parse_args()
    molecule = get_molecule(args.molecule)
    root = Path(f"plots/rb_{molecule.key}_polar")
    data_dir = args.npz_dir or root / "data"
    output = args.out or root / "figures" / "fig_orientation_alignment_n25_n29_MJ0.png"

    datasets = []
    for n in args.n:
        path = data_dir / f"orientation_MJ{args.mj}_n{n}.npz"
        with np.load(path) as data:
            if "COS2" not in data.files:
                raise ValueError(f"{path} no contiene COS2; vuelve a calcularlo")
            keep = (data["R"] >= args.rmin) & (data["R"] <= args.rmax)
            datasets.append((n, data["R"][keep], data["COS"][keep], data["COS2"][keep]))

    Rref = np.linspace(args.rmin, args.rmax, 500)
    parallel = isolated_rotor_curves(molecule, Rref, +1.0)
    antiparallel = isolated_rotor_curves(molecule, Rref, -1.0)

    fig, axes = plt.subplots(1, 2, figsize=(12.4, 5.1), sharex=True)
    colors = ("#315b9a", "#e69f00")
    widths = (1.5, 2.8)
    for (n, R, cos, cos2), color, width in zip(datasets, colors, widths):
        axes[0].plot(R, cos, color=color, lw=width,
                     label=rf"Rb($n={n},\,l\geq3$)–{molecule.label}")
        axes[1].plot(R, cos2, color=color, lw=width,
                     label=rf"Rb($n={n},\,l\geq3$)–{molecule.label}")
    for ax, index in zip(axes, (0, 1)):
        ax.plot(Rref, parallel[index], color="#b12a90", lw=2.0, ls=(0, (1, 2)),
                label=r"rotor, $\mathbf{E}=+\hat Z/R^2$")
        ax.plot(Rref, antiparallel[index], color="0.25", lw=1.5, ls=(0, (1, 3)),
                label=r"rotor, $\mathbf{E}=-\hat Z/R^2$")
        ax.set_xlim(args.rmin, args.rmax)
        ax.set_xlabel(r"$R$ [$a_0$]")
        ax.grid(color="0.88", lw=0.7)
        ax.legend(frameon=False, fontsize=8.5, loc="best")
    axes[0].axhline(0, color="0.25", lw=0.7)
    axes[0].set_ylabel(r"orientación $\langle\cos\theta_d\rangle$")
    axes[0].set_ylim(-1.0, 1.0)
    axes[1].set_ylabel(r"alineamiento $\langle\cos^2\theta_d\rangle$")
    axes[1].set_ylim(0.0, 1.0)
    axes[0].text(0.97, 0.95, "(a)", transform=axes[0].transAxes, ha="right", va="top")
    axes[1].text(0.97, 0.95, "(b)", transform=axes[1].transAxes, ha="right", va="top")
    fig.suptitle(rf"Rb*+{molecule.label}: orientación y alineamiento, $M_J={args.mj}$, "
                 rf"$N_{{max}}={args.n_max}$")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)
    print(f"Figura: {output}")


if __name__ == "__main__":
    main()
