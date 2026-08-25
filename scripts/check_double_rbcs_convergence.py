#!/usr/bin/env python3
"""Convergencia rotacional puntual de Rb*+RbCs+RbCs."""

import argparse
import time

import numpy as np

from trimero.systems.double_polar_rydberg import GHZ_PER_HARTREE, RbTwoRbCsSystem
from trimero.systems.rb_krb_polar.rb_defects import DELTA0_NS_PAPER


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--n-manifold", type=int, default=20)
    p.add_argument("--n-max", type=int, nargs="+", default=[3, 4])
    p.add_argument("--geometry", choices=("symmetric", "unilateral"),
                   default="symmetric")
    p.add_argument("--separation", type=float, default=300.0)
    p.add_argument("--r", type=float, nargs="+", default=[400, 800, 1400])
    p.add_argument("--mj", type=int, nargs="+", default=[0, 1])
    p.add_argument("--fields", type=float, nargs="+", default=[0, 500])
    p.add_argument("--weight", type=float, default=0.5)
    p.add_argument("--k", type=int, default=120)
    p.add_argument("--max-k", type=int, default=480)
    p.add_argument("--sigma-ghz", type=float, default=-30.0)
    return p.parse_args()


def characteristic(system, radius, mj, geometry, separation, field, args):
    seed_energy, seed_vector = system.manifold_seed(
        radius, mj, geometry, separation, field
    )
    sigma = (seed_energy - system.E_manifold) * GHZ_PER_HARTREE
    k = args.k
    while True:
        values, vectors = system.solve_near(
            radius, mj, geometry, separation, field, k=k, sigma_ghz=sigma,
        )
        weights = np.sum(vectors[system.is_manifold(mj), :] ** 2, axis=0)
        overlaps = (seed_vector @ vectors) ** 2
        selected = int(np.argmax(overlaps))
        if weights[selected] > args.weight:
            break
        if k >= args.max_k:
            raise RuntimeError("ningún estado supera el peso de manifold pedido")
        k = min(2 * k, args.max_k)
    vector = vectors[:, selected]
    return np.array([
        (values[selected] - system.E_manifold) * GHZ_PER_HARTREE,
        vector @ (system.orientation_matrix(1, mj) @ vector),
        vector @ (system.orientation_matrix(2, mj) @ vector),
        weights[selected],
    ])


def main():
    args = parse_args()
    results = {}
    for nmax in args.n_max:
        system = RbTwoRbCsSystem(
            n_manifold=args.n_manifold, N_max=nmax,
            delta0_ns=DELTA0_NS_PAPER,
        )
        for mj in args.mj:
            print(f"Nmax={nmax}, MJ={mj}, dim={len(system.block(mj))}")
            for field in args.fields:
                for radius in args.r:
                    started = time.perf_counter()
                    value = characteristic(
                        system, radius, mj, args.geometry,
                        args.separation, field, args,
                    )
                    results[(nmax, mj, field, radius)] = value
                    print(
                        f"  F={field:g}, R={radius:g}: E={value[0]:.6f} GHz, "
                        f"cos=({value[1]:+.6f},{value[2]:+.6f}), "
                        f"W={value[3]:.4f}, t={time.perf_counter()-started:.1f}s"
                    )
    print("\nDiferencias entre truncaciones consecutivas:")
    for low, high in zip(args.n_max, args.n_max[1:]):
        for mj in args.mj:
            for field in args.fields:
                for radius in args.r:
                    delta = np.abs(
                        results[(high, mj, field, radius)][:3]
                        - results[(low, mj, field, radius)][:3]
                    )
                    ok = delta[0] < 0.01 and np.all(delta[1:] < 0.01)
                    print(
                        f"  {low}->{high}, MJ={mj}, F={field:g}, R={radius:g}: "
                        f"dE={delta[0]:.6f} GHz, dcos=({delta[1]:.6f},"
                        f"{delta[2]:.6f}) {'PASA' if ok else 'NO CONVERGE'}"
                    )


if __name__ == "__main__":
    main()
