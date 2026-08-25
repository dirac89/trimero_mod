#!/usr/bin/env python3
"""Convergencia de la base pendular contraída para Rb*+RbCs+RbCs."""

import argparse
import time

from trimero.systems.double_polar_rydberg import ContractedRbTwoRbCsSystem
from trimero.systems.rb_krb_polar.rb_defects import DELTA0_NS_PAPER


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--n-manifold", type=int, default=20)
    p.add_argument("--primitive-n-max", type=int, default=8)
    p.add_argument("--rotor-keep", type=int, nargs="+", default=[1, 2, 3, 4])
    p.add_argument("--geometry", choices=("symmetric", "unilateral"),
                   default="symmetric")
    p.add_argument("--separation", type=float, default=300.0)
    p.add_argument("--r", type=float, default=800.0)
    p.add_argument("--mj", type=int, default=0)
    p.add_argument("--field", type=float, default=0.0)
    p.add_argument("--k", type=int, default=24)
    args = p.parse_args()

    previous = None
    for keep in args.rotor_keep:
        system = ContractedRbTwoRbCsSystem(
            n_manifold=args.n_manifold,
            primitive_N_max=args.primitive_n_max,
            rotor_keep=keep,
            delta0_ns=DELTA0_NS_PAPER,
        )
        point = system.point(
            args.r, args.mj, args.geometry, args.separation, args.field
        )
        started = time.perf_counter()
        result = system.characteristic(point, k=args.k)
        elapsed = time.perf_counter() - started
        print(
            f"keep={keep}, dim={result['dimension']}, E={result['E']:.6f} GHz, "
            f"cos=({result['COS1']:+.6f},{result['COS2']:+.6f}), "
            f"W={result['W']:.4f}, t={elapsed:.1f}s"
        )
        if previous is not None:
            print(
                f"  Δ: dE={abs(result['E']-previous['E']):.6f} GHz, "
                f"dcos=({abs(result['COS1']-previous['COS1']):.6f},"
                f"{abs(result['COS2']-previous['COS2']):.6f})"
            )
        previous = result


if __name__ == "__main__":
    main()
