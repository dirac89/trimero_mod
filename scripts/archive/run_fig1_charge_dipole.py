#!/usr/bin/env python3
"""
Fig. 1 de Aguilera-Fernández et al. 2015: curvas BOP de Rb*-KRb con

    H_ad = H_A + H_mol          (su Ec. 1)

y NADA MÁS. Sin pseudopotencial de Fermi.

    poetry run python scripts/run_fig1_charge_dipole.py

POR QUÉ NO HAY V_Fermi AQUÍ
---------------------------
Ese paper trata KRb como un DIPOLO PUNTUAL que siente el campo eléctrico del
Rydberg, no como un centro de dispersión de contacto para el electrón. H_A es
el Rydberg libre (defectos cuánticos de Marinescu-Sadeghpour-Dalgarno) y H_mol
es puramente carga-dipolo, B·N² - d·F_ryd.

El pseudopotencial de `hamiltonians/fermi_krb.py` pertenece a otra línea
(perturbador NEUTRO, Aguilera-Fernández 2016). Sigue siendo código válido para
ese otro problema, pero NO entra en esta comparación: por eso este script NO lo
importa ni lo instancia, y arma el Hamiltoniano a mano en `hamiltonian_ad()` en
vez de llamar a `BOPSystem.hamiltonian(fermi=False)`. Se ve de un vistazo qué
términos hay.

CONSECUENCIA: H_A y H_mol no tienen restricción de dominio en R — se evalúan con
la función de onda hidrogenoide, que decae pero está definida en todo R. No hay
remapeo k(R), ni ventana de exclusión de la resonancia de onda p, ni cola
butterfly, ni tope en 2n²a₀. El rango [400, 1800] a₀ se calcula entero y directo.

Ver docs/analysis_fig1_carga_dipolo_sin_fermi.md.
"""
import argparse
import os
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from trimero.systems.rb_krb_polar.charge_dipole import B_KRB_GHZ
from trimero.systems.rb_krb_polar.bop_system import GHZ_PER_HARTREE as GHZ, BOPSystem
from trimero.systems.rb_krb_polar.rb_defects import DELTA0_NS_PAPER as D0

RULE = "=" * 92
N_KEEP = 250          # autovalores guardados por punto (para el fondo del plot)


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-manifold", type=int, default=25)
    ap.add_argument("--mj", type=int, nargs="+", default=[0, 1])
    ap.add_argument("--rmin", type=float, default=400.0)
    ap.add_argument("--rmax", type=float, default=1800.0)
    ap.add_argument("--step", type=float, default=5.0)
    ap.add_argument("--ymin", type=float, default=-25.0)
    ap.add_argument("--ymax", type=float, default=1.0)
    ap.add_argument("--reuse", action="store_true")
    ap.add_argument("--out", default="plots/fig1_ad_MJ0_MJ1_n25.png")
    return ap.parse_args()


# --------------------------------------------------------- Hamiltoniano
def hamiltonian_ad(sysm, M_J, R):
    """
    H_ad(R) = H_A + H_mol, escrito término a término.

    H_A   : diagonal, energías Rydberg del manifold y de los tres vecinos.
    H_mol : B·N² - d·F_ion(R) - d·F_elec(R), construido por
            ChargeDipoleHamiltonian.

    No hay un tercer término. Esto no es `fermi=False`: es que aquí el
    pseudopotencial no forma parte del modelo.
    """
    blk = sysm.block(M_J)
    return np.diag(sysm.rydberg_diagonal(M_J)) + sysm.hmol.build(blk, R)


def character_curve_ad(sysm, M_J, R, weight=0.5):
    """
    (E-E_manifold [GHz], k, peso, todos los autovalores) de la curva adiabática
    más baja con CARÁCTER de manifold.

    Sigue haciendo falta el criterio de carácter aunque no haya V_Fermi: los
    estados de (n+1)d y (n+2)p producen cruces evitados con las curvas del
    manifold, así que un índice fijo cambiaría de objeto por el camino.
    Ver docs/analysis_base_correcta_3_vecinos.md §5.1.
    """
    w, V = np.linalg.eigh(hamiltonian_ad(sysm, M_J, R))
    mask = sysm.is_manifold(M_J)
    for k in range(len(w)):
        if float(np.sum(V[:, k][mask] ** 2)) > weight:
            return ((w[k] - sysm.E_manifold) * GHZ, k,
                    float(np.sum(V[:, k][mask] ** 2)), w)
    return float("nan"), -1, 0.0, w


# --------------------------------------------------------- utilidades
def local_minima(y):
    return [i for i in range(1, len(y) - 1) if y[i] < y[i - 1] and y[i] < y[i + 1]]


def local_maxima(y):
    return [i for i in range(1, len(y) - 1) if y[i] > y[i - 1] and y[i] > y[i + 1]]


def main():
    args = parse_args()
    sysm = BOPSystem(n_manifold=args.n_manifold, delta0_ns=D0)
    n = sysm.n_manifold
    L = {0: "s", 1: "p", 2: "d"}
    print(RULE)
    print(f"Fig. 1 — Rb*-KRb, H_ad = H_A + H_mol (SIN pseudopotencial de Fermi)")
    print(RULE)
    print(f"\n  manifold n={n} (l>={sysm.l_min}..{sysm.l_max}) + "
          + " + ".join(f"{sysm.levels[l]}{L[l]}" for l in (2, 1, 0)))
    print(f"  delta0_ns = {D0}   cero de energía: E(n={n}, l>=3) + KRb(N=0) "
          f"= {sysm.E_manifold:.12e} E_h")
    print(f"  rango: R ∈ [{args.rmin}, {args.rmax}] a0, paso {args.step} a0  "
          "— completo, sin recortes")
    print("  NO hay: remapeo k(R), ventana de resonancia p, cola butterfly ni "
          "tope en 2n²a0.")
    print(f"  (2n² = {2*n*n} a0 se marca en la figura sólo como referencia de "
          "escala.)")

    thr = {N: sysm.delta_E_ghz(0) + B_KRB_GHZ * N * (N + 1) for N in (5, 6)}
    print(f"\n  umbrales asintóticos {sysm.n_s}s + KRb(N):")
    print(f"    ΔE({sysm.n_s}s) + 30B (N=5) = {thr[5]:9.4f} GHz")
    print(f"    ΔE({sysm.n_s}s) + 42B (N=6) = {thr[6]:9.4f} GHz")

    R = np.arange(args.rmin, args.rmax + 1e-9, args.step)
    data = {}
    for M_J in args.mj:
        npz = f"plots/fig1_ad_MJ{M_J}_n{n}.npz"
        if args.reuse and os.path.exists(npz):
            d = np.load(npz)
            data[M_J] = {k: d[k] for k in d.files}
            print(f"\n  M_J={M_J}: reutilizando {npz}")
            continue
        blk = sysm.block(M_J)
        print(f"\n  M_J={M_J}: dim(bloque) = {len(blk)}   "
              f"{len(R)} puntos")
        t0 = time.perf_counter()
        E = np.full(len(R), np.nan)
        K = np.full(len(R), -1, dtype=int)
        W = np.zeros(len(R))
        keep = min(N_KEEP, len(blk))
        SP = np.empty((len(R), keep))
        for i, r in enumerate(R):
            E[i], K[i], W[i], w_all = character_curve_ad(sysm, M_J, float(r))
            SP[i] = w_all[:keep]
            if i % 40 == 0:
                print(f"    R = {r:8.2f} ({i+1}/{len(R)})", flush=True)
        dt = time.perf_counter() - t0
        print(f"    {len(R)} diagonalizaciones en {dt:.1f} s "
              f"({dt/len(R):.2f} s/punto)")
        data[M_J] = {"R": R, "E": E, "K": K, "W": W,
                     "spectrum": (SP - sysm.E_manifold) * GHZ}
        np.savez(npz, **data[M_J])
        print(f"    datos en {npz}")

    # ------------------------------------------- [3] límite en R = R_max
    print("\n" + "-" * 92)
    print(f"[3] LÍMITE EN EL BORDE SUPERIOR R = {args.rmax:.0f} a0")
    print("-" * 92)
    print("    Sin dominio truncado, la curva del manifold debe tender al cero")
    print("    de energía (manifold + KRb(N=0)) cuando R -> infinito.\n")
    print("     M_J | E(R_max) [GHz] |  k  | peso manifold | E(R_max-100) [GHz] "
          "| pendiente [GHz/100 a0]")
    print("    -----|----------------|-----|---------------|--------------------"
          "|----------------------")
    for M_J in args.mj:
        d = data[M_J]
        Rr, Ee = d["R"], d["E"]
        i_end = len(Rr) - 1
        i_100 = int(np.argmin(np.abs(Rr - (args.rmax - 100.0))))
        slope = (Ee[i_end] - Ee[i_100]) / ((Rr[i_end] - Rr[i_100]) / 100.0)
        print(f"    {M_J:4d} | {Ee[i_end]:14.5f} | {d['K'][i_end]:3d} | "
              f"{d['W'][i_end]:13.4f} | {Ee[i_100]:18.5f} | {slope:+21.5f}")

    # ------------------------------------------- [5] forma de la curva
    print("\n" + "-" * 92)
    print("[5] FORMA DE LA CURVA — números para comparar con la figura real")
    print("-" * 92)
    stats = {}
    for M_J in args.mj:
        d = data[M_J]
        Rr, y = d["R"], d["E"]
        mins, maxs = local_minima(y), local_maxima(y)
        i_deep = int(np.nanargmin(y))
        # dónde la curva sube por encima de umbrales de "casi cero"
        def cross_up(level):
            idx = np.flatnonzero(y > level)
            return float(Rr[idx[0]]) if len(idx) else float("nan")
        stats[M_J] = {"deep_R": float(Rr[i_deep]), "deep_E": float(y[i_deep]),
                      "n_min": len(mins), "n_max": len(maxs),
                      "mins": [(float(Rr[i]), float(y[i])) for i in mins]}
        print(f"\n    M_J = {M_J}")
        print(f"      pozo MÁS PROFUNDO: {y[i_deep]:.4f} GHz en R = {Rr[i_deep]:.1f} a0")
        print(f"      mínimos locales: {len(mins)}   máximos locales: {len(maxs)}"
              f"   -> ~{len(mins)} oscilaciones visibles")
        print(f"      rango de la curva: [{np.nanmin(y):.4f}, {np.nanmax(y):.4f}] GHz")
        for lvl in (-10.0, -5.0, -2.0, -1.0):
            print(f"      cruza {lvl:6.1f} GHz (subiendo) en R = {cross_up(lvl):8.1f} a0")
        print("      los 8 mínimos más externos:")
        for Rw, Ew in stats[M_J]["mins"][-8:]:
            print(f"        R = {Rw:7.1f} a0   E = {Ew:9.4f} GHz   "
                  f"R/2n² = {Rw/(2*n*n):.3f}")

    # ------------------------------------------- [4] figura
    fig, axes = plt.subplots(1, len(args.mj), figsize=(6.2 * len(args.mj), 5.6),
                             sharey=True, squeeze=False)
    for ax, M_J in zip(axes[0], args.mj):
        d = data[M_J]
        Rr, y, SP = d["R"], d["E"], d["spectrum"]
        # fondo: TODAS las curvas guardadas, sin filtrar
        vis = np.any((SP > args.ymin - 5.0) & (SP < args.ymax + 5.0), axis=0)
        ax.plot(Rr, SP[:, vis], color="0.62", lw=0.5, alpha=0.85, zorder=1)
        ax.plot([], [], color="0.62", lw=0.5,
                label=f"resto del bloque ({int(vis.sum())} curvas en la ventana)")
        ax.plot(Rr, y, color="0.05", lw=2.4, zorder=4,
                label="más baja de carácter manifold (>50 % de peso)")
        for N, c in ((5, "tab:orange"), (6, "tab:green")):
            ax.axhline(thr[N], color=c, lw=1.3, ls="--", alpha=0.9, zorder=3,
                       label=f"$\\Delta E_{{{sysm.n_s}s}}+{N*(N+1)}B$ "
                             f"($N={N}$) = {thr[N]:.2f} GHz")
        ax.axhline(0.0, color="k", lw=0.9, ls=":", zorder=2)
        ax.axvline(2.0 * n * n, color="tab:purple", lw=1.0, ls="-.", alpha=0.5,
                   zorder=2, label=f"$2n^2 = {2*n*n}\\,a_0$ (referencia)")
        ax.set_xlim(args.rmin, args.rmax)
        ax.set_ylim(args.ymin, args.ymax)
        ax.grid(alpha=0.25)
        ax.set_xlabel(r"$R$  [$a_0$]")
        ax.set_title(f"$M_J = {M_J}$", fontsize=12)
    axes[0][0].set_ylabel(r"$V(R) = E - [E_{n=%d,\,l\geq3} + E_{KRb}(N=0)]$  [GHz]"
                          % n)
    fig.suptitle(
        f"Rb*-KRb, curvas BOP del manifold $n={n}$   —   "
        r"$H_{ad} = H_A + H_{mol}$  (Ec. 1 de Aguilera-Fernández et al. 2015)"
        "\nSIN pseudopotencial de Fermi: no hay remapeo $k(R)$, ni ventana de "
        "resonancia, ni tope de dominio — el rango 400-1800 $a_0$ es directo\n"
        f"base: manifold $({n},l\\geq3)$ + ${sysm.levels[2]}d$ + "
        f"${sysm.levels[1]}p$ + ${sysm.levels[0]}s$",
        fontsize=10.5)
    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, fontsize=8.5,
               frameon=False, bbox_to_anchor=(0.5, 0.005))
    fig.tight_layout(rect=[0, 0.11, 1, 0.88])
    fig.savefig(args.out, dpi=150)
    print(f"\n  PNG guardado en {args.out}")
    print(RULE)


if __name__ == "__main__":
    main()
