#!/usr/bin/env python3
"""
Curvas BOP de Rb* unido a KRb o RbCs: script polar de producción.

Sustituye a los doce `run_*.py` / `analyze_*.py` de las rondas de exploración,
que están en `scripts/archive/`. Lo que aquí se calcula es lo que
`docs/STATUS.md` declara vigente, y nada más.

MODELO FÍSICO
-------------
`--system polar` (el único vigente) es el Hamiltoniano de Aguilera-Fernández,
Sadeghpour, Schmelcher & González-Férez, J. Phys.: Conf. Ser. 635, 012023
(2015), Ec. 1:

    H_ad(R) = H_A + H_mol
            = diag(E_ryd) + [B·N² - d·F_ion(R) - d·F_elec(R)]

La molécula es un DIPOLO PUNTUAL en el campo del Rydberg. **No hay pseudopotencial de
Fermi**: el de contacto modela un perturbador NEUTRO (la otra línea, la de
Aguilera-Fernández 2016), y meterlo aquí fue la premisa equivocada que
`docs/analysis_fig1_carga_dipolo_sin_fermi.md` §1 corrige.

Consecuencia práctica: no hay remapeo k(R), ni ventana de exclusión de la
resonancia de onda p, ni tope de dominio en 2n²a₀. H_A y H_mol se evalúan con
la función de onda hidrogenoide, definida en todo R, así que el rango sale
entero y continuo.

Se llama a `BOPSystem.hamiltonian(R, M_J, fermi=False)`, que es literalmente
`diag(rydberg_diagonal) + hmol.build(blk, R)` — la misma expresión que el
`hamiltonian_ad()` del script archivado con el que se produjeron los números
verificados. Mismos números bit a bit, sin reimplementar el Hamiltoniano.

⚠️ `fermi=False` NO significa «este sistema tiene V_Fermi y lo apago». Significa
que el pseudopotencial no forma parte de este modelo. Que el flag exista es
deuda técnica documentada en `docs/STATUS.md`: `BOPSystem` sigue construyendo un
`FermiPseudopotential` en `__post_init__` aunque el lado polar no lo use.

IDENTIFICACIÓN DE LA CURVA
--------------------------
Por CARÁCTER (peso de manifold > `--weight`), no por índice fijo. Sigue siendo
necesario aunque no haya V_Fermi: los estados de (n+1)d y (n+2)p producen cruces
evitados con las del manifold, y un índice fijo cambiaría de objeto por el
camino. Ver `docs/analysis_base_correcta_3_vecinos.md` §5.1.

QUÉ MOLÉCULA POLAR
------------------
`--molecule` es OBLIGATORIO y no tiene valor por defecto. Antes lo tenía
(`krb`), y esa es exactamente la razón por la que la Fase A del plan de figuras
se calculó con KRb creyendo que era RbCs: un default silencioso no aparece en la
línea de órdenes que uno copia al documento, así que nada delata la molécula
usada. Ahora hay que decirla siempre, y queda escrita en el `.npz`
(`molecule`, `B_hz`, `d_debye`) y en el título de la figura.

Cada molécula se instancia con SU clase de sistema (`SYSTEM_FOR_MOLECULE`):
`rbcs` → `RbRbCsPolarSystem`, `krb` → `PolarBOPSystem(molecule=KRB)`.

USO
---
    # Fig. 1 de Aguilera-Fernández 2015, Rb*-KRb (los números verificados)
    poetry run python scripts/compute_bop_curve.py --molecule krb --n-manifold 25 --mj 0 1

    # Rb*-RbCs, un manifold distinto, sin figura
    poetry run python scripts/compute_bop_curve.py --molecule rbcs --n-manifold 24 --mj 0 --no-plot
"""
import argparse
import os
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from trimero.systems.polar_molecule import MOLECULES, get_molecule
from trimero.systems.polar_rydberg import GHZ_PER_HARTREE as GHZ, PolarBOPSystem
from trimero.systems.rb_krb_polar.bop_system import BOPSystem
from trimero.systems.rb_krb_polar.rb_defects import DELTA0_NS_PAPER
from trimero.systems.rb_rbcs_polar import RbRbCsPolarSystem

RULE = "=" * 92
N_KEEP = 250          # autovalores guardados por punto, para el fondo de la figura

# Clase de sistema por molécula. `PolarBOPSystem` es genérico y valdría para las
# dos, pero la clase específica es la que fija sus constantes sin que quien
# llama pueda contradecirlas por descuido, que es el fallo que se corrige aquí.
SYSTEM_FOR_MOLECULE = {"rbcs": RbRbCsPolarSystem}


def build_system(molecule, **kwargs):
    """Instancia el sistema polar de esta molécula, con su clase específica."""
    cls = SYSTEM_FOR_MOLECULE.get(molecule.key)
    if cls is None:
        return PolarBOPSystem(molecule=molecule, **kwargs)
    system = cls(**kwargs)          # la clase específica fija molecule ella misma
    assert system.molecule is molecule
    return system


def parse_args():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--molecule", choices=tuple(MOLECULES), required=True,
                    help="molécula polar (OBLIGATORIO, sin valor por defecto): "
                         + " o ".join(sorted(MOLECULES)))
    ap.add_argument("--n-manifold", type=int, default=25,
                    help="n del manifold cuasi-degenerado (l >= 3)")
    ap.add_argument("--n-max", type=int, default=6,
                    help="corte rotacional N_max (RbCs n=25 validado con 6)")
    ap.add_argument("--mj", type=int, nargs="+", default=[0, 1],
                    help="bloques M_J a calcular")
    ap.add_argument("--rmin", type=float, default=400.0)
    ap.add_argument("--rmax", type=float, default=1800.0)
    ap.add_argument("--step", type=float, default=5.0)
    ap.add_argument("--weight", type=float, default=0.5,
                    help="peso mínimo de manifold para aceptar la curva")
    ap.add_argument("--npz-dir", default=None,
                    help="directorio de datos; por defecto depende de --molecule")
    ap.add_argument("--reuse", action="store_true",
                    help="reutiliza el .npz si existe, en vez de rebarrer R")
    ap.add_argument("--out", default=None,
                    help="PNG; por defecto plots/rb_<mol>_polar/figures/")
    ap.add_argument("--no-plot", action="store_true")
    ap.add_argument("--ymin", type=float, default=-25.0)
    ap.add_argument("--ymax", type=float, default=1.0)
    return ap.parse_args()


# ------------------------------------------------------------ Hamiltoniano
def character_curve(sysm, M_J, R, weight=0.5):
    """
    (E - E_manifold [GHz], k, peso, todos los autovalores) de la curva
    adiabática más baja con CARÁCTER de manifold en este R.
    """
    if isinstance(sysm, PolarBOPSystem):
        H = sysm.hamiltonian(R, M_J)
    else:  # compatibilidad con la regresión del BOPSystem histórico
        H = sysm.hamiltonian(R, M_J, fermi=False)
    w, V = np.linalg.eigh(H)
    mask = sysm.is_manifold(M_J)
    for k in range(len(w)):
        wm = float(np.sum(V[:, k][mask] ** 2))
        if wm > weight:
            return (w[k] - sysm.E_manifold) * GHZ, k, wm, w
    return float("nan"), -1, 0.0, w


def local_minima(y):
    return [i for i in range(1, len(y) - 1) if y[i] < y[i - 1] and y[i] < y[i + 1]]


# ------------------------------------------------------------------ barrido
def sweep(sysm, M_J, R, weight):
    blk = sysm.block(M_J)
    print(f"\n  M_J={M_J}: dim(bloque) = {len(blk)}   {len(R)} puntos")
    t0 = time.perf_counter()
    E = np.full(len(R), np.nan)
    K = np.full(len(R), -1, dtype=int)
    W = np.zeros(len(R))
    keep = min(N_KEEP, len(blk))
    SP = np.empty((len(R), keep))
    for i, r in enumerate(R):
        E[i], K[i], W[i], w_all = character_curve(sysm, M_J, float(r), weight)
        SP[i] = w_all[:keep]
        if i % 40 == 0:
            print(f"    R = {r:8.2f} ({i+1}/{len(R)})", flush=True)
    dt = time.perf_counter() - t0
    print(f"    {len(R)} diagonalizaciones en {dt:.1f} s ({dt/len(R):.2f} s/punto)")
    return {"R": R, "E": E, "K": K, "W": W,
            "spectrum": (SP - sysm.E_manifold) * GHZ}


def report(sysm, M_J, d, rmax):
    Rr, y = d["R"], d["E"]
    i_deep = int(np.nanargmin(y))
    mins = local_minima(y)
    print(f"\n    M_J = {M_J}")
    print(f"      pozo MÁS PROFUNDO   : {y[i_deep]:9.4f} GHz en R = {Rr[i_deep]:.1f} a0")
    print(f"      E(R = {Rr[-1]:.0f} a0)     : {y[-1]:9.4f} GHz   "
          f"(k = {d['K'][-1]}, peso = {d['W'][-1]:.4f})")
    print(f"      mínimos locales     : {len(mins)}")
    print(f"      rango de la curva   : [{np.nanmin(y):.4f}, {np.nanmax(y):.4f}] GHz")
    for lvl in (-10.0, -5.0, -2.0, -1.0):
        idx = np.flatnonzero(y > lvl)
        cross = float(Rr[idx[0]]) if len(idx) else float("nan")
        print(f"      cruza {lvl:6.1f} GHz subiendo en R = {cross:8.1f} a0")
    return {"deep_E": float(y[i_deep]), "deep_R": float(Rr[i_deep]),
            "E_end": float(y[-1]), "n_min": len(mins)}


def make_plot(sysm, args, data, thr, out):
    n = sysm.n_manifold
    mj = list(args.mj)
    fig, axes = plt.subplots(1, len(mj), figsize=(6.2 * len(mj), 5.6),
                             sharey=True, squeeze=False)
    for ax, M_J in zip(axes[0], mj):
        d = data[M_J]
        Rr, y, SP = d["R"], d["E"], d["spectrum"]
        vis = np.any((SP > args.ymin - 5.0) & (SP < args.ymax + 5.0), axis=0)
        ax.plot(Rr, SP[:, vis], color="0.62", lw=0.5, alpha=0.85, zorder=1)
        ax.plot([], [], color="0.62", lw=0.5,
                label=f"resto del bloque ({int(vis.sum())} curvas en la ventana)")
        ax.plot(Rr, y, color="0.05", lw=2.4, zorder=4,
                label=f"más baja de carácter manifold (>{args.weight:.0%} de peso)")
        for N, c in ((5, "tab:orange"), (6, "tab:green")):
            ax.axhline(thr[N], color=c, lw=1.3, ls="--", alpha=0.9, zorder=3,
                       label=f"$\\Delta E_{{{sysm.n_s}s}}+{N*(N+1)}B$ "
                             f"($N={N}$) = {thr[N]:.2f} GHz")
        ax.axhline(0.0, color="k", lw=0.9, ls=":", zorder=2)
        ax.axvline(2.0 * n * n, color="tab:purple", lw=1.0, ls="-.", alpha=0.5,
                   zorder=2, label=f"$2n^2 = {2*n*n}\\,a_0$ (referencia de escala)")
        ax.set_xlim(args.rmin, args.rmax)
        ax.set_ylim(args.ymin, args.ymax)
        ax.grid(alpha=0.25)
        ax.set_xlabel(r"$R$  [$a_0$]")
        ax.set_title(f"$M_J = {M_J}$", fontsize=12)
    axes[0][0].set_ylabel(
        rf"$V(R) = E - [E_{{n={n},\,l\geq3}} + E_{{{sysm.molecule.label}}}(N=0)]$  [GHz]")
    fig.suptitle(
        f"Rb*-{sysm.molecule.label}, curvas BOP del manifold $n={n}$   —   "
        r"$H_{ad} = H_A + H_{mol}$  (Ec. 1 de Aguilera-Fernández et al. 2015)"
        "\nSIN pseudopotencial de Fermi: no hay remapeo $k(R)$, ni ventana de "
        "resonancia, ni tope de dominio\n"
        f"base: manifold $({n},l\\geq3)$ + ${sysm.levels[2]}d$ + "
        f"${sysm.levels[1]}p$ + ${sysm.levels[0]}s$",
        fontsize=10.5)
    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=8.5,
               frameon=False, bbox_to_anchor=(0.5, 0.005))
    fig.tight_layout(rect=[0, 0.11, 1, 0.88])
    out_dir = os.path.dirname(out)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    fig.savefig(out, dpi=150)
    print(f"\n  PNG guardado en {out}")


def main():
    args = parse_args()
    molecule = get_molecule(args.molecule)
    root = f"plots/rb_{molecule.key}_polar"
    args.npz_dir = args.npz_dir or f"{root}/data"
    sysm = build_system(molecule, n_manifold=args.n_manifold,
                        N_max=args.n_max,
                        delta0_ns=DELTA0_NS_PAPER)
    n = sysm.n_manifold
    L = {0: "s", 1: "p", 2: "d"}

    print(RULE)
    print(f"Curvas BOP Rb*-{molecule.label}   —   H_ad = H_A + H_mol   (SIN pseudopotencial "
          "de Fermi)")
    print(RULE)
    print(f"\n  molécula: {molecule.label} (key={molecule.key})   "
          f"B = {molecule.B_ghz:.6f} GHz   d = {molecule.dipole_debye:.3f} D   "
          f"clase: {type(sysm).__name__}")
    print(f"\n  manifold n={n} (l={sysm.l_min}..{sysm.l_max}) + "
          + " + ".join(f"{sysm.levels[l]}{L[l]}" for l in (2, 1, 0)))
    print(f"  delta0_ns = {DELTA0_NS_PAPER}   cero de energía: "
          f"E(n={n}, l>=3) + {molecule.label}(N=0) = {sysm.E_manifold:.12e} E_h")
    print(f"  rango: R ∈ [{args.rmin:.0f}, {args.rmax:.0f}] a0, paso "
          f"{args.step:.0f} a0 — completo, sin recortes")

    thr = sysm.thresholds_ns_ghz((5, 6))
    print(f"\n  umbrales asintóticos {sysm.n_s}s + {molecule.label}(N):")
    for N in (5, 6):
        print(f"    ΔE({sysm.n_s}s) + {N*(N+1)}B (N={N}) = {thr[N]:9.4f} GHz")

    R = np.arange(args.rmin, args.rmax + 1e-9, args.step)
    data, summary = {}, {}
    for M_J in args.mj:
        npz = os.path.join(args.npz_dir, f"fig1_ad_MJ{M_J}_n{n}.npz")
        if args.reuse and os.path.exists(npz):
            d = np.load(npz)
            if "molecule" not in d.files or str(d["molecule"]) != molecule.key:
                raise ValueError(
                    f"{npz} no contiene metadatos compatibles con {molecule.key}; "
                    "recalcula sin --reuse"
                )
            data[M_J] = {k: d[k] for k in d.files}
            print(f"\n  M_J={M_J}: reutilizando {npz}")
        else:
            data[M_J] = sweep(sysm, M_J, R, args.weight)
            os.makedirs(args.npz_dir, exist_ok=True)
            np.savez(npz, **data[M_J], molecule=molecule.key,
                     B_hz=molecule.rotational_constant_hz,
                     d_debye=molecule.dipole_debye,
                     n_manifold=n, N_max=sysm.N_max, M_J=M_J,
                     character_weight=args.weight, schema_version=1)
            print(f"    datos en {npz}")

    print("\n" + "-" * 92)
    print("FORMA DE LA CURVA")
    print("-" * 92)
    for M_J in args.mj:
        summary[M_J] = report(sysm, M_J, data[M_J], args.rmax)

    if not args.no_plot:
        out = args.out or os.path.join(
            root, "figures",
            "fig1_ad_" + "_".join(f"MJ{m}" for m in args.mj) + f"_n{n}.png")
        make_plot(sysm, args, data, thr, out)

    print(RULE)
    return summary


if __name__ == "__main__":
    main()
