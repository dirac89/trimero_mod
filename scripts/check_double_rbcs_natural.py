#!/usr/bin/env python3
"""Comprueba la contracción natural de Rb*+RbCs+RbCs en un punto."""

import argparse
import time

import numpy as np

from trimero.systems.double_polar_rydberg import (
    ContractedRbTwoRbCsSystem,
    GHZ_PER_HARTREE,
    RbTwoRbCsSystem,
)
from trimero.systems.rb_krb_polar.rb_defects import DELTA0_NS_PAPER


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-manifold", type=int, default=20)
    parser.add_argument("--source-n-max", type=int, default=6)
    parser.add_argument("--target-n-max", type=int, default=8)
    parser.add_argument("--natural-keep", type=int, nargs="+", default=[2, 3, 4])
    parser.add_argument("--geometry", choices=("symmetric", "unilateral"),
                        default="symmetric")
    parser.add_argument("--separation", type=float, default=300.0)
    parser.add_argument("--r", type=float, default=800.0)
    parser.add_argument("--mj", type=int, default=0)
    parser.add_argument("--field", type=float, default=0.0)
    parser.add_argument("--k", type=int, default=40)
    parser.add_argument("--source-k", type=int, default=120)
    parser.add_argument("--weight", type=float, default=0.5)
    parser.add_argument("--source-cache", type=str)
    parser.add_argument("--save-source", type=str)
    return parser.parse_args()


def source_state(system, args):
    block = system.block(args.mj)
    if args.source_cache:
        cached = np.load(args.source_cache)
        vector = cached["vector"]
        if vector.shape != (len(block),):
            raise ValueError("el vector cacheado no coincide con la base fuente")
        if "energy_ghz" in cached:
            energy = float(cached["energy_ghz"])
        elif "E" in cached:
            # Compatibilidad con los primeros benchmarks: E era la energía
            # absoluta en Hartree, no el desplazamiento en GHz.
            energy = (float(cached["E"]) - system.E_manifold) * GHZ_PER_HARTREE
        else:
            raise ValueError("la caché no contiene energy_ghz ni E")
        weight = (
            float(cached["weight"])
            if "weight" in cached
            else float(np.sum(vector[system.is_manifold(args.mj)] ** 2))
        )
        return block, vector, energy, weight

    seed_energy, seed_vector = system.manifold_seed(
        args.r, args.mj, args.geometry, args.separation, args.field
    )
    sigma = (seed_energy - system.E_manifold) * GHZ_PER_HARTREE
    values, vectors = system.solve_near(
        args.r, args.mj, args.geometry, args.separation, args.field,
        k=args.source_k, sigma_ghz=sigma,
    )
    weights = np.sum(vectors[system.is_manifold(args.mj), :] ** 2, axis=0)
    overlaps = (seed_vector @ vectors) ** 2
    eligible = np.flatnonzero(weights > args.weight)
    if not len(eligible):
        raise RuntimeError("ningún autoestado fuente supera el peso solicitado")
    selected = eligible[np.argmax(overlaps[eligible])]
    energy = (values[selected] - system.E_manifold) * GHZ_PER_HARTREE
    return block, vectors[:, selected], float(energy), float(weights[selected])


def main():
    args = parse_args()
    primitive = RbTwoRbCsSystem(
        n_manifold=args.n_manifold, N_max=args.source_n_max,
        delta0_ns=DELTA0_NS_PAPER,
    )
    started = time.perf_counter()
    block, vector, source_energy, source_weight = source_state(primitive, args)
    if args.save_source:
        np.savez_compressed(
            args.save_source, vector=vector, energy_ghz=source_energy,
            weight=source_weight, n_manifold=args.n_manifold,
            source_n_max=args.source_n_max, geometry=args.geometry,
            separation=args.separation, R=args.r, M_J=args.mj,
            field_v_per_m=args.field,
        )
    print(
        f"fuente Nmax={args.source_n_max}, dim={len(block)}, "
        f"E={source_energy:.9f} GHz, W={source_weight:.6f}, "
        f"t={time.perf_counter()-started:.1f}s"
    )

    previous = None
    for keep in args.natural_keep:
        contracted = ContractedRbTwoRbCsSystem(
            n_manifold=args.n_manifold,
            primitive_N_max=args.target_n_max,
            rotor_keep=keep,
            delta0_ns=DELTA0_NS_PAPER,
        )
        point = contracted.natural_point(
            block, vector, args.r, args.mj, args.geometry,
            args.separation, args.field, natural_keep=keep,
        )
        started = time.perf_counter()
        result = contracted.characteristic(point, k=args.k)
        elapsed = time.perf_counter() - started
        print(
            f"keep={keep}, dim={result['dimension']}, E={result['E']:.9f} GHz, "
            f"cos=({result['COS1']:+.6f},{result['COS2']:+.6f}), "
            f"W={result['W']:.6f}, t={elapsed:.1f}s, "
            f"|E-Esource|={abs(result['E']-source_energy):.6f} GHz"
        )
        if previous is not None:
            delta_e = abs(result["E"] - previous["E"])
            delta_cos = max(
                abs(result["COS1"] - previous["COS1"]),
                abs(result["COS2"] - previous["COS2"]),
            )
            print(
                f"  k anterior→{keep}: dE={delta_e:.6f} GHz, "
                f"dcos_max={delta_cos:.6f} "
                f"{'PASA' if delta_e < 0.01 and delta_cos < 0.01 else 'NO CONVERGE'}"
            )
        previous = result


if __name__ == "__main__":
    main()
