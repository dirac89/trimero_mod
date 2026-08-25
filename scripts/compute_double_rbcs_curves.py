#!/usr/bin/env python3
"""BOP y orientación de Rb*+RbCs+RbCs en dos geometrías colineales."""

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from trimero.systems.double_polar_rydberg import GHZ_PER_HARTREE, RbTwoRbCsSystem
from trimero.systems.rb_krb_polar.rb_defects import DELTA0_NS_PAPER


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--geometry", choices=("symmetric", "unilateral", "both"),
                   default="both")
    p.add_argument("--separation", type=float, default=300.0)
    p.add_argument("--n-manifold", type=int, default=20)
    p.add_argument("--n-max", type=int, default=3)
    p.add_argument("--mj", type=int, nargs="+", default=[0, 1])
    p.add_argument("--fields", type=float, nargs="+", default=[0, 100, 300, 500])
    p.add_argument("--rmin", type=float, default=200.0)
    p.add_argument("--rmax", type=float, default=1800.0)
    p.add_argument("--step", type=float, default=10.0)
    p.add_argument("--weight", type=float, default=0.5)
    p.add_argument("--overlap", type=float, default=0.7)
    p.add_argument("--min-substep", type=float, default=0.5,
                   help="paso mínimo de la continuación adaptativa [a0]")
    p.add_argument("--bisect-max-depth", type=int, default=6,
                   help="número máximo de bisecciones por hueco")
    p.add_argument("--sigma-ghz", type=float, default=-30.0)
    p.add_argument("--k", type=int, default=80)
    p.add_argument("--max-k", type=int, default=320)
    p.add_argument("--context", type=int, default=60)
    p.add_argument("--workers", type=int, default=1,
                   help="casos (geometría, MJ, campo) ejecutados en paralelo")
    p.add_argument("--out-root", type=Path,
                   default=Path("plots/rb_rbcs_rbcs_polar"))
    p.add_argument("--reuse", action="store_true")
    p.add_argument("--no-plot", action="store_true")
    p.add_argument("--ymin", type=float, default=None,
                   help="límite inferior común; automático si se omite")
    p.add_argument("--ymax", type=float, default=20.0)
    return p.parse_args()


def dataset_path(root, geometry, n, nmax, mj, field, separation):
    dtag = f"_D{separation:g}" if geometry == "unilateral" else ""
    return root / "data" / (
        f"bop_{geometry}{dtag}_n{n}_Nmax{nmax}_MJ{mj}_F{field:g}.npz"
    )


def compatible(path, args, geometry, mj, field, R):
    try:
        with np.load(path) as data:
            return (
                np.array_equal(data["R"], R)
                and str(data["geometry"]) == geometry
                and int(data["n_manifold"]) == args.n_manifold
                and int(data["N_max"]) == args.n_max
                and int(data["M_J"]) == mj
                and float(data["field_v_per_m"]) == field
                and float(data["separation_a0"]) == args.separation
                and int(data["schema_version"]) == 1
            )
    except (KeyError, OSError, ValueError):
        return False


def solve_with_k(system, radius, mj, geometry, separation, field, sigma, k):
    return system.solve_near(
        radius, mj, geometry, separation, field, k=k, sigma_ghz=sigma
    )


def _continuation_attempt(
    system, radius, mj, geometry, args, field, sigma, mask,
    reference_vector=None, seed_vector=None,
):
    """Resuelve un radio conservando intacta la duplicación existente de k."""
    k = args.k
    while True:
        values, vectors = solve_with_k(
            system, radius, mj, geometry, args.separation, field, sigma, k,
        )
        weights = np.sum(vectors[mask, :] ** 2, axis=0)
        eligible = np.flatnonzero(weights > args.weight)
        reference = seed_vector if reference_vector is None else reference_vector
        overlaps = (reference @ vectors) ** 2
        selected = (
            int(eligible[np.argmax(overlaps[eligible])])
            if len(eligible) else int(np.argmax(overlaps))
        )
        quality = float(overlaps[selected])
        if weights[selected] > args.weight:
            if quality >= args.overlap or k >= args.max_k:
                break
        if k >= args.max_k:
            raise RuntimeError(
                f"no se localizó la rama en R={radius:g}, MJ={mj}, F={field:g}"
            )
        k = min(2 * k, args.max_k)

    relative = (values - system.E_manifold) * GHZ_PER_HARTREE
    return {
        "values": values,
        "vectors": vectors,
        "weights": weights,
        "selected": selected,
        "quality": quality,
        "vector": vectors[:, selected],
        "relative": relative,
        "energy": float(relative[selected]),
        "k": k,
    }


def _bisect_continuation(
    system, start_radius, start_vector, start_sigma, target_radius,
    mj, geometry, args, field, mask, depth=0, initial_result=None,
):
    """Cierra un hueco con puntos temporales que no forman parte del grid."""
    target = initial_result or _continuation_attempt(
        system, target_radius, mj, geometry, args, field, start_sigma, mask,
        reference_vector=start_vector,
    )
    if target["quality"] >= args.overlap:
        return target

    half_step = abs(target_radius - start_radius) / 2.0
    if depth >= args.bisect_max_depth or half_step < args.min_substep:
        return target

    midpoint = (start_radius + target_radius) / 2.0
    middle = _bisect_continuation(
        system, start_radius, start_vector, start_sigma, midpoint,
        mj, geometry, args, field, mask, depth=depth + 1,
    )
    if middle["quality"] < args.overlap:
        return target

    return _bisect_continuation(
        system, midpoint, middle["vector"], middle["energy"], target_radius,
        mj, geometry, args, field, mask, depth=depth + 1,
    )


def sweep(system, R, mj, geometry, args, field):
    """Siembra en Rmax y sigue hacia dentro por máximo solapamiento."""
    mask = system.is_manifold(mj)
    C1 = system.orientation_matrix(1, mj)
    C2 = system.orientation_matrix(2, mj)
    order = np.arange(len(R) - 1, -1, -1)
    E = np.full(len(R), np.nan)
    W = np.full(len(R), np.nan)
    O = np.full(len(R), np.nan)
    COS1 = np.full(len(R), np.nan)
    COS2 = np.full(len(R), np.nan)
    KUSED = np.zeros(len(R), dtype=int)
    spectrum = np.full((len(R), args.context), np.nan)
    warnings = []
    previous = None
    previous_radius = None
    last_good_vector = None
    last_good_radius = None
    last_good_sigma = None
    sigma = args.sigma_ghz
    seed_vector = None
    started = time.perf_counter()

    for count, idx in enumerate(order, 1):
        radius = float(R[idx])
        if previous is None:
            seed_energy, seed_vector = system.manifold_seed(
                radius, mj, geometry, args.separation, field
            )
            sigma = (seed_energy - system.E_manifold) * GHZ_PER_HARTREE
        result = _continuation_attempt(
            system, radius, mj, geometry, args, field, sigma, mask,
            reference_vector=previous, seed_vector=seed_vector,
        )
        if (
            previous is not None
            and result["quality"] < args.overlap
            and result["k"] == args.max_k
            and last_good_vector is not None
        ):
            result = _bisect_continuation(
                system, last_good_radius, last_good_vector, last_good_sigma,
                radius, mj, geometry, args, field, mask,
                initial_result=(
                    result if previous_radius == last_good_radius else None
                ),
            )

        vector = result["vector"]
        relative = result["relative"]
        selected = result["selected"]
        E[idx] = result["energy"]
        W[idx] = result["weights"][selected]
        O[idx] = result["quality"]
        COS1[idx] = float(vector @ (C1 @ vector))
        COS2[idx] = float(vector @ (C2 @ vector))
        KUSED[idx] = result["k"]
        take = min(args.context, len(relative))
        spectrum[idx, :take] = relative[:take]
        if O[idx] < args.overlap or W[idx] < args.weight:
            warnings.append(
                (radius, float(O[idx]), float(W[idx]), int(result["k"]))
            )
        previous = vector
        previous_radius = radius
        sigma = float(E[idx])
        if O[idx] >= args.overlap:
            last_good_vector = vector
            last_good_radius = radius
            last_good_sigma = sigma
        if count == 1 or count % 20 == 0 or count == len(R):
            print(
                f"{geometry} MJ={mj} F={field:g}: R={radius:.0f} "
                f"({count}/{len(R)}), E={E[idx]:.3f} GHz, "
                f"k={result['k']}", flush=True
            )

    elapsed = time.perf_counter() - started
    print(f"  terminado en {elapsed:.1f}s; avisos de continuidad: {len(warnings)}")
    return dict(
        R=R, E=E, W=W, overlap=O, COS1=COS1, COS2=COS2,
        k_used=KUSED, spectrum=spectrum,
        warning_R=np.array([w[0] for w in warnings]),
        warning_overlap=np.array([w[1] for w in warnings]),
        warning_weight=np.array([w[2] for w in warnings]),
        elapsed_seconds=elapsed,
    )


def save(path, data, args, geometry, mj, field, dim):
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(
        path, **data, geometry=geometry, separation_a0=args.separation,
        molecule="rbcs", n_manifold=args.n_manifold, N_max=args.n_max,
        M_J=mj, field_v_per_m=field, character_weight=args.weight,
        overlap_threshold=args.overlap, min_substep_a0=args.min_substep,
        bisect_max_depth=args.bisect_max_depth, dimension=dim, schema_version=1,
    )


def plot_geometry(results, geometry, args):
    nfields, nmj = len(args.fields), len(args.mj)
    fig, axes = plt.subplots(2 * nmj, nfields, figsize=(4.0 * nfields, 3.2 * 2 * nmj),
                             sharex=True, squeeze=False)
    for col, field in enumerate(args.fields):
        for row, mj in enumerate(args.mj):
            data = results[(mj, field)]
            ax_e, ax_o = axes[2 * row, col], axes[2 * row + 1, col]
            ax_e.plot(data["R"], data["spectrum"], color="0.78", lw=0.45)
            ax_e.plot(data["R"], data["E"], color="#315b9a", lw=2.2)
            uncertain = data["overlap"] < args.overlap
            ax_e.scatter(
                data["R"][uncertain], data["E"][uncertain], s=13,
                marker="x", color="#b2182b", linewidths=0.8, zorder=4,
                label="solapamiento bajo",
            )
            ax_e.axhline(0, color="0.25", lw=0.7, ls=":")
            ax_e.set_ylim(args.ymin, args.ymax)
            ax_e.set_title(f"F={field:g} V/m")
            ax_e.set_ylabel(rf"$M_J={mj}$  $E-E_n$ [GHz]")
            ax_o.plot(data["R"], data["COS1"], color="#315b9a", lw=1.8,
                      label=r"$\langle\cos\theta_1\rangle$")
            ax_o.plot(data["R"], data["COS2"], color="#e69f00", lw=1.8, ls="--",
                      label=r"$\langle\cos\theta_2\rangle$")
            ax_o.axhline(0, color="0.25", lw=0.7)
            ax_o.set_ylim(-1.02, 1.02)
            ax_o.set_ylabel("orientación")
            if row == 0 and col == 0:
                ax_o.legend(frameon=False, fontsize=8)
                ax_e.legend(frameon=False, fontsize=7, loc="lower right")
            for ax in (ax_e, ax_o):
                ax.grid(color="0.9", lw=0.6)
        axes[-1, col].set_xlabel(r"$R$ [$a_0$]")
    geometry_label = (
        r"RbCs($-R$)–Rb*–RbCs($+R$)"
        if geometry == "symmetric"
        else rf"Rb*–RbCs($R$)–RbCs($R+{args.separation:g}$)"
    )
    fig.suptitle(
        rf"{geometry_label}, $n={args.n_manifold}$, "
        rf"$N_{{max}}={args.n_max}$ (exploratorio)",
        fontsize=14,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out = args.out_root / "figures" / (
        f"comparison_{geometry}_n{args.n_manifold}_Nmax{args.n_max}.png"
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=170)
    plt.close(fig)
    print(f"Figura: {out}")


def compute_case(args, geometry, mj, field, R):
    """Ejecuta y persiste un caso independiente; apto para multiprocessing."""
    system = RbTwoRbCsSystem(
        n_manifold=args.n_manifold, N_max=args.n_max,
        delta0_ns=DELTA0_NS_PAPER,
    )
    path = dataset_path(
        args.out_root, geometry, args.n_manifold, args.n_max,
        mj, field, args.separation,
    )
    if args.reuse and path.exists() and compatible(
        path, args, geometry, mj, field, R
    ):
        return geometry, mj, field, path, True
    data = sweep(system, R, mj, geometry, args, field)
    save(path, data, args, geometry, mj, field, len(system.block(mj)))
    return geometry, mj, field, path, False


def main():
    args = parse_args()
    if args.rmax <= args.rmin or args.step <= 0:
        raise SystemExit("se requiere rmax>rmin y step>0")
    if args.min_substep <= 0 or args.bisect_max_depth < 0:
        raise SystemExit("min-substep debe ser >0 y bisect-max-depth >=0")
    geometries = ("symmetric", "unilateral") if args.geometry == "both" else (args.geometry,)
    R = np.arange(args.rmin, args.rmax + 1e-9, args.step)
    if args.workers < 1:
        raise SystemExit("workers debe ser >=1")
    tasks = [
        (geometry, mj, field)
        for geometry in geometries for mj in args.mj for field in args.fields
    ]
    completed = []
    if args.workers == 1:
        completed = [compute_case(args, *task, R) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            futures = {
                executor.submit(compute_case, args, *task, R): task for task in tasks
            }
            for future in as_completed(futures):
                result = future.result()
                completed.append(result)
                geometry, mj, field, path, reused = result
                action = "reutilizado" if reused else "calculado"
                print(f"Caso {action}: {geometry}, MJ={mj}, F={field:g}: {path}")

    if not args.no_plot:
        if args.ymin is None:
            minima = []
            for geometry in geometries:
                for mj in args.mj:
                    for field in args.fields:
                        path = dataset_path(
                            args.out_root, geometry, args.n_manifold, args.n_max,
                            mj, field, args.separation,
                        )
                        with np.load(path) as saved:
                            minima.append(float(np.nanmin(saved["E"])))
            args.ymin = 10.0 * np.floor(min(minima) / 10.0)
        for geometry in geometries:
            results = {}
            for mj in args.mj:
                for field in args.fields:
                    path = dataset_path(
                        args.out_root, geometry, args.n_manifold, args.n_max,
                        mj, field, args.separation,
                    )
                    with np.load(path) as saved:
                        results[(mj, field)] = {
                            key: saved[key] for key in saved.files
                        }
            plot_geometry(results, geometry, args)


if __name__ == "__main__":
    main()
