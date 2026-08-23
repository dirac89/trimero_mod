#!/usr/bin/env python3
"""Comprueba convergencia rotacional de Rb*+molécula polar en varios radios."""

import argparse

from trimero.systems.polar_molecule import MOLECULES, get_molecule
from trimero.systems.polar_rydberg import PolarBOPSystem
from trimero.systems.rb_krb_polar.rb_defects import DELTA0_NS_PAPER


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--molecule", choices=tuple(MOLECULES), default="rbcs")
    parser.add_argument("--n-manifold", type=int, default=25)
    parser.add_argument("--mj", type=int, default=0)
    parser.add_argument("--n-max", type=int, nargs="+", default=[4, 6, 8])
    parser.add_argument("--r", type=float, nargs="+", default=[500.0, 1000.0, 1500.0])
    parser.add_argument("--weight", type=float, default=0.5)
    parser.add_argument("--tolerance-ghz", type=float, default=0.1)
    return parser.parse_args()


def main():
    args = parse_args()
    molecule = get_molecule(args.molecule)
    rows = {}
    for cutoff in args.n_max:
        system = PolarBOPSystem(
            molecule=molecule, n_manifold=args.n_manifold, N_max=cutoff,
            delta0_ns=DELTA0_NS_PAPER,
        )
        rows[cutoff] = [system.character_curve(r, args.mj, args.weight)[0] for r in args.r]
        print(f"N_max={cutoff}: " + ", ".join(
            f"E({r:g})={energy:.8f} GHz" for r, energy in zip(args.r, rows[cutoff])))
    previous, latest = args.n_max[-2:]
    errors = [abs(a - b) for a, b in zip(rows[previous], rows[latest])]
    worst = max(errors)
    print(f"max |E(N_max={latest})-E(N_max={previous})| = {worst:.6g} GHz")
    if worst > args.tolerance_ghz:
        raise SystemExit(
            f"NO CONVERGE dentro de {args.tolerance_ghz:g} GHz; aumenta N_max"
        )
    return rows


if __name__ == "__main__":
    main()
