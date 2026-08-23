#!/usr/bin/env python
"""
Esquemas de geometría molecular desde la línea de órdenes.

Envoltorio fino de `trimero.visualization.geometry_diagram`: no contiene
lógica de dibujo ni física, sólo traduce argumentos.

Cada cuerpo se declara con `--body`, repetible, en formato
`clave=valor,clave=valor`:

    kind=neutral_atom|polar_molecule   (por defecto neutral_atom)
    label=Rb                           texto junto a la esfera
    R=900                              distancia en a₀ (proporción relativa)
    theta=180                          ángulo polar respecto a +Z, 0..180
    phi=0                              ángulo azimutal, opcional
    dipole=0                           dibuja d⃗ con ese ángulo polar, opcional
    color=#c2543a                      opcional

Ejemplos:

    # El híbrido de este proyecto: Rb neutro en θ=π, RbCs polar en θ=0
    poetry run python scripts/draw_geometry.py \\
        --body kind=neutral_atom,label=Rb,R=600,theta=180 \\
        --body kind=polar_molecule,label=RbCs,R=900,theta=0,dipole=0 \\
        --electron r=1200,theta=35 \\
        --title "Rb*(35,l>=3) + Rb + RbCs" \\
        -o plots/geometry/hybrid_neutral_polar_esquema.png --formats png pdf

    # Atajo con el convenio R1/R2 de HybridNeutralPolar
    poetry run python scripts/draw_geometry.py --hybrid 600 900 -o plots/geometry/h.png

    # Configuración no colineal cualquiera
    poetry run python scripts/draw_geometry.py \\
        --body label=A,R=900,theta=60 --body label=B,R=1200,theta=140,phi=180 \\
        -o /tmp/no_colineal.png
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")   # utilidad de script: nunca abre ventana

from trimero.visualization.geometry_diagram import (  # noqa: E402
    Body,
    Electron,
    draw_from_system_config,
    draw_geometry,
)

# Nombres cortos admitidos en --body / --electron -> campos de las dataclases.
FIELD_ALIASES = {"theta": "theta_deg", "phi": "phi_deg", "dipole": "dipole_deg"}
FLOAT_FIELDS = {"R", "r", "theta_deg", "phi_deg", "dipole_deg", "size"}


def _parse_spec(text: str, what: str) -> dict:
    """`clave=valor,clave=valor` -> dict con los tipos ya convertidos."""
    spec = {}
    for chunk in text.split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        if "=" not in chunk:
            raise argparse.ArgumentTypeError(
                f"{what}: '{chunk}' no tiene la forma clave=valor")
        key, value = chunk.split("=", 1)
        key = FIELD_ALIASES.get(key.strip(), key.strip())
        value = value.strip()
        if key in FLOAT_FIELDS:
            try:
                spec[key] = float(value)
            except ValueError:
                raise argparse.ArgumentTypeError(
                    f"{what}: '{key}' esperaba un número, no {value!r}")
        else:
            spec[key] = value
    return spec


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--body", action="append", default=[], metavar="SPEC",
                   help="un cuerpo, repetible. Ver ejemplos arriba.")
    p.add_argument("--hybrid", nargs=2, type=float, metavar=("R1", "R2"),
                   help="atajo: R1 = neutro en θ=180°, R2 = polar en θ=0°, "
                        "el convenio de HybridNeutralPolar")
    p.add_argument("--electron", metavar="SPEC",
                   help="el electrón Rydberg, p. ej. r=1200,theta=35")
    p.add_argument("--ion-label", default=r"Rb$^+$")
    p.add_argument("--title", default=None)
    p.add_argument("-o", "--output", type=Path, required=True,
                   help="ruta de salida; la extensión la fija --formats")
    p.add_argument("--formats", nargs="+", default=["png"])
    p.add_argument("--projection", choices=("xz", "oblique"), default="xz")
    p.add_argument("--radius-scale", choices=("proportional", "relative"),
                   default="proportional")
    p.add_argument("--dpi", type=int, default=200)
    p.add_argument("--figsize", nargs=2, type=float, default=[6.4, 6.4])
    p.add_argument("--no-distance-labels", action="store_true")
    return p.parse_args()


def main():
    args = parse_args()
    if not args.body and args.hybrid is None:
        raise SystemExit("nada que dibujar: usa --body (repetible) o --hybrid")

    electron = None
    if args.electron:
        electron = Electron(**_parse_spec(args.electron, "--electron"))

    common = dict(
        electron=electron,
        ion_label=args.ion_label,
        title=args.title,
        output_path=args.output,
        formats=tuple(args.formats),
        projection=args.projection,
        radius_scale=args.radius_scale,
        figsize=tuple(args.figsize),
        dpi=args.dpi,
        show_distance_labels=not args.no_distance_labels,
    )

    if args.hybrid is not None and not args.body:
        R1, R2 = args.hybrid
        _, _, paths = draw_from_system_config({"R1": R1, "R2": R2}, **common)
    else:
        bodies = [Body(**_parse_spec(spec, "--body")) for spec in args.body]
        if args.hybrid is not None:
            R1, R2 = args.hybrid
            bodies = [Body(kind="neutral_atom", label="Rb", R=R1, theta_deg=180.0),
                      Body(kind="polar_molecule", label="RbCs", R=R2,
                           theta_deg=0.0, dipole_deg=0.0)] + bodies
        _, _, paths = draw_geometry(bodies, **common)

    for path in paths:
        print(f"-> {path}")


if __name__ == "__main__":
    main()
