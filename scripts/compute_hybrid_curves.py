#!/usr/bin/env python3
"""
Curvas adiabáticas del sistema HÍBRIDO Rb*(n=35) - Rb(5s, θ=π) - RbCs(θ=0).

    poetry run python scripts/compute_hybrid_curves.py

H(R1, R2) = H_A(n=35) + H_mol^KRb(R2) + V_Fermi^π(R1)

  H_A        : manifold Rb(n=35, l>=3) hidrogenoide (defectos nulos para l>=3).
  H_mol^KRb  : rotor rígido del RbCs (B=490.17 MHz, d=1.225 D) en el campo del
               ion y del electrón Rydberg, sobre el eje +Z (θ=0).
  V_Fermi^π  : pseudopotencial de contacto s+p del perturbador NEUTRO de Rb(5s)
               en θ=π, con la fase de paridad (-1)^{l+l'} verificada por fuerza
               bruta (tests/systems/hybrid_neutral_polar/test_transformacion_theta_pi.py).

El buen número cuántico es M_J = m_l + M_N (test_conservacion_mj.py).

La parte de Fermi sólo depende de R1: se evalúa UNA vez por R1 y se suma en
cada R2 con el mismo orden de términos que `HybridNeutralPolar.hamiltonian`,
así que los números son idénticos a los del método.

Fases 1-2 del plan híbrido: primera producción física tras el arreglo de
wigner_3j (hermiticidad 3e-19 relativo; ver docs/analysis_wigner3j_orden_canonico.md).
"""
import argparse
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from trimero.systems.hybrid_neutral_polar import (
    B_RBCS_MHZ,
    D_RBCS_DEBYE,
    GHZ_PER_HARTREE,
    HybridNeutralPolar,
)

RULE = "=" * 92


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-manifold", type=int, default=35)
    ap.add_argument("--n-max", type=int, default=2,
                    help="N_max del rotor RbCs (los tests de límite usan 2)")
    ap.add_argument("--mj", type=int, default=0)
    ap.add_argument("--r1", type=float, nargs="+", default=[600.0, 900.0, 1100.0])
    ap.add_argument("--rmin", type=float, default=500.0)
    ap.add_argument("--rmax", type=float, default=1500.0)
    ap.add_argument("--step", type=float, default=25.0)
    ap.add_argument("--keep", type=int, default=40,
                    help="autovalores más bajos guardados por punto")
    return ap.parse_args()


def local_minima(y):
    return [i for i in range(1, len(y) - 1) if y[i] < y[i - 1] and y[i] < y[i + 1]]


def main():
    args = parse_args()
    mol = HybridNeutralPolar(n_manifold=args.n_manifold, N_max=args.n_max)
    out_dir = Path("plots") / "hybrid_neutral_polar"
    out_dir.mkdir(parents=True, exist_ok=True)

    print(RULE)
    print("Fases 1-2 — Rb*(n=35)-Rb(5s,θ=π)-RbCs(θ=0): curvas E0(R2) a R1 fijo")
    print(RULE)
    print(f"  d(RbCs) = {mol.d_debye} D = {mol.d_au:.12f} ea0   "
          f"B(RbCs) = {mol.B_mhz} MHz = {mol.B_au:.6e} Eh")
    print(f"  bloque M_J={args.mj}: dim = {len(mol.block(args.mj))}   "
          f"N_max = {args.n_max}")
    print(f"  cero de energías: E(n=35,l>=3) = {mol.E_manifold:.12e} Eh")

    R2 = np.arange(args.rmin, args.rmax + 1e-9, args.step)
    diag = mol.rydberg_diagonal(args.mj)
    blk = mol.block(args.mj)
    keep = min(args.keep, len(blk))

    # sanidad: hermiticidad del ensamblado en el primer punto del primer barrido
    Vf0 = mol.fermi_pi_matrix(blk, args.r1[0])
    H0 = np.diag(diag) + mol.hmol.build(blk, float(R2[0])) + Vf0
    asym = float(np.linalg.norm(H0 - H0.T))
    print(f"\n  sanidad ||H-H^T||_F en (R1={args.r1[0]:.0f}, R2={R2[0]:.0f}): "
          f"{asym:.3e} Eh")

    curves = {}
    for R1 in args.r1:
        print(f"\n  --- R1 = {R1:.0f} a0 (perturbador neutro en theta=pi) ---")
        Vf = Vf0 if R1 == args.r1[0] else mol.fermi_pi_matrix(blk, float(R1))
        E0 = np.full(len(R2), np.nan)
        SP = np.empty((len(R2), keep))
        t0 = time.perf_counter()
        for i, r2 in enumerate(R2):
            H = np.diag(diag) + mol.hmol.build(blk, float(r2)) + Vf
            w = np.linalg.eigvalsh(H)
            E0[i] = (w[0] - mol.E_manifold) * GHZ_PER_HARTREE
            SP[i] = (w[:keep] - mol.E_manifold) * GHZ_PER_HARTREE
            if i % 10 == 0:
                print(f"    R2 = {r2:8.1f} ({i+1}/{len(R2)})  "
                      f"E0 = {E0[i]:+9.4f} GHz", flush=True)
        dt = time.perf_counter() - t0
        mins = local_minima(E0)
        i_deep = int(np.argmin(E0))
        print(f"    {len(R2)} puntos en {dt:.1f} s ({dt/len(R2):.2f} s/punto)")
        print(f"    E0: rango [{E0.min():+.4f}, {E0.max():+.4f}] GHz; "
              f"mínimo global {E0[i_deep]:+.4f} GHz en R2 = {R2[i_deep]:.0f} a0")
        print(f"    mínimos locales de E0: {len(mins)} en "
              f"{[int(R2[i]) for i in mins]} a0")
        curves[R1] = {"R2": R2, "E0": E0, "spectrum": SP}
        np.savez(out_dir / f"hybrid_curves_R1{int(R1)}_n{args.n_manifold}.npz",
                 R2=R2, E0=E0, spectrum=SP, R1=R1,
                 M_J=args.mj, N_max=args.n_max)
        print(f"    datos en {out_dir}/hybrid_curves_R1{int(R1)}_n{args.n_manifold}.npz")

    # ------------------------------------------------------------- figura
    fig, ax = plt.subplots(figsize=(7.4, 5.4))
    colors = {600.0: "tab:blue", 900.0: "tab:red", 1100.0: "tab:green"}
    for R1, d in curves.items():
        vis = np.any((d["spectrum"] > d["E0"].min() - 15.0)
                     & (d["spectrum"] < 20.0), axis=0)
        ax.plot(d["R2"], d["spectrum"][:, vis], color="0.85", lw=0.4, zorder=1)
        ax.plot(d["R2"], d["E0"], lw=2.2, color=colors.get(R1),
                label=f"$R_1 = {R1:.0f}\\,a_0$ (neutro $\\theta=\\pi$)")
    ax.axhline(0.0, color="k", lw=0.8, ls=":")
    ax.set_xlabel(r"$R_2$  [$a_0$]  (distancia al RbCs, $\theta=0$)")
    ax.set_ylabel(r"$E_0 - E_{n=35,\,l\geq3}$  [GHz]")
    ax.set_title(f"Híbrido Rb*(n={args.n_manifold})-Rb-RbCs: autovalor más bajo "
                 f"del bloque $M_J={args.mj}$, $N_{{max}}={args.n_max}$")
    ax.grid(alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    png = out_dir / f"hybrid_curves_MJ{args.mj}_n{args.n_manifold}.png"
    fig.savefig(png, dpi=150)
    print(f"\n  PNG en {png}")
    print(RULE)


if __name__ == "__main__":
    main()
