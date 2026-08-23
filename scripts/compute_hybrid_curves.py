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

IDENTIFICACIÓN DE LA CURVA (corregido tras docs/analysis_hibrido_caracter_E0.md)
-------------------------------
Por CARÁCTER: autoestado más bajo con peso > --weight en el manifold
(n=35, l>=3), como en scripts/compute_bop_curve.py del sistema polar. El
autovalor más bajo sin más es el UMBRAL 37p (−102.360 GHz), no ligadura.

Estado de REFERENCIA (v3): "el más bajo con carácter" sobre el Hamiltoniano
COMPLETO no define un objeto continuo en R1 — V_Fermi^π sumerge niveles de
carácter manifold hasta cientos de GHz por debajo (p.ej. −286.9 GHz en
R1=800), igual que ocurre en el módulo neutro validado
(SymmetricLinearTrimer da −329.2 GHz en R=800; las figuras del paper 2016
recortan en −40 GHz). Esos niveles sumergidos son un canal aparte. La curva
de ligadura del manifold se define aquí como el estado que continúa
adiabáticamente al autoestado MÁS BAJO CON CARÁCTER de H_A+H_mol (sin
Fermi) evaluado en R2 = R2[0]; a partir de ahí, seguimiento adiabático por
máximo solapamiento de autovector entre puntos consecutivos de R2.

El seguimiento NO se interrumpe cuando el estado pierde carácter: se
registra el peso en cada punto y se avisa explícitamente. Física esperada:
al subir R2 la curva cruza el umbral 38s (−20.267 GHz) y intercambia
carácter con ese canal; el tramo por debajo de −20.3 GHz es ligadura de
manifold limpia. El criterio puramente puntual se guarda siempre
(K_puntual/E_puntual) para comparar.

La parte de Fermi sólo depende de R1: se evalúa UNA vez por R1 y se suma en
cada R2 con el mismo orden de términos que `HybridNeutralPolar.hamiltonian`,
así que los números son idénticos a los del método.
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
NEIGH_LABELS = ("38s", "37p", "36d")   # l = 0, 1, 2


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-manifold", type=int, default=35)
    ap.add_argument("--n-max", type=int, default=6,
                    help="N_max del rotor RbCs (convergido: max|d(6-4)| "
                         "= 0.013 GHz en R1=900, R2 in {500,1000,1500})")
    ap.add_argument("--mj", type=int, default=0)
    ap.add_argument("--r1", type=float, nargs="+", default=[600.0, 900.0, 1100.0])
    ap.add_argument("--rmin", type=float, default=500.0)
    ap.add_argument("--rmax", type=float, default=1500.0)
    ap.add_argument("--step", type=float, default=25.0)
    ap.add_argument("--weight", type=float, default=0.5,
                    help="peso mínimo de manifold para aceptar la curva")
    ap.add_argument("--keep", type=int, default=40,
                    help="autovalores más bajos guardados por punto")
    return ap.parse_args()


def pesos_autovector(estados, v):
    """(peso manifold, pesos vecinos [38s,37p,36d]) de un autovector."""
    man = 0.0
    neigh = [0.0, 0.0, 0.0]
    for s, c in zip(estados, v):
        p = c * c
        if s[0] >= 3:
            man += p
        else:
            neigh[s[0]] += p
    return man, neigh


def local_minima(y):
    return [i for i in range(1, len(y) - 1) if y[i] < y[i - 1] and y[i] < y[i + 1]]


def main():
    args = parse_args()
    mol = HybridNeutralPolar(n_manifold=args.n_manifold, N_max=args.n_max)
    out_dir = Path("plots") / "hybrid_neutral_polar"
    out_dir.mkdir(parents=True, exist_ok=True)

    print(RULE)
    print("Rb*(n=35)-Rb(5s,theta=pi)-RbCs(theta=0): curvas de CARÁCTER manifold")
    print(RULE)
    print(f"  d(RbCs) = {mol.d_debye} D = {mol.d_au:.12f} ea0   "
          f"B(RbCs) = {mol.B_mhz} MHz = {mol.B_au:.6e} Eh")
    print(f"  bloque M_J={args.mj}: dim = {len(mol.block(args.mj))}   "
          f"N_max = {args.n_max}   umbral de carácter = {args.weight:.0%}")
    print(f"  cero de energías: E(n={args.n_manifold},l>=3) = "
          f"{mol.E_manifold:.12e} Eh")

    R2 = np.arange(args.rmin, args.rmax + 1e-9, args.step)
    diag = mol.rydberg_diagonal(args.mj)
    blk = mol.block(args.mj)
    estados = blk.states
    keep = min(args.keep, len(blk))

    # sanidad: hermiticidad del ensamblado en el primer punto del primer barrido
    Vf0 = mol.fermi_pi_matrix(blk, args.r1[0])
    H0 = np.diag(diag) + mol.hmol.build(blk, float(R2[0])) + Vf0
    print(f"\n  sanidad ||H-H^T||_F en (R1={args.r1[0]:.0f}, R2={R2[0]:.0f}): "
          f"{float(np.linalg.norm(H0 - H0.T)):.3e} Eh")
    del H0

    # ------------------------------------------------ estado de referencia:
    # más bajo con carácter de H_A+H_mol SIN Fermi en R2[0]; todas las curvas
    # R1 continúan ESTE objeto (continuidad en la fuerza de V_Fermi)
    w_r, V_r = np.linalg.eigh(np.diag(diag) + mol.hmol.build(blk, float(R2[0])))
    mask_man = np.array([s[0] >= 3 for s in estados], dtype=bool)
    pm_ref = np.sum(V_r[mask_man] ** 2, axis=0)
    ok_ref = np.flatnonzero(pm_ref > args.weight)
    jref = int(ok_ref[int(np.argmin(w_r[ok_ref]))])
    v_ref = V_r[:, jref].copy()
    E_ref = (w_r[jref] - mol.E_manifold) * GHZ_PER_HARTREE
    print(f"  referencia sin Fermi en R2={R2[0]:.0f}: k={jref}, "
          f"E={E_ref:+.4f} GHz, peso manifold={pm_ref[jref]:.4f}")
    del V_r

    curves = {}
    for R1 in args.r1:
        print(f"\n  --- R1 = {R1:.0f} a0 (perturbador neutro en theta=pi) ---")
        Vf = Vf0 if R1 == args.r1[0] else mol.fermi_pi_matrix(blk, float(R1))
        E = np.full(len(R2), np.nan)          # curva con seguimiento adiabático
        K = np.full(len(R2), -1, dtype=int)   # índice seguido
        Wman = np.zeros(len(R2))
        Wneigh = np.zeros((len(R2), 3))
        OVL = np.ones(len(R2))                # solapamiento con el punto anterior
        K_puntual = np.full(len(R2), -1, dtype=int)  # criterio puntual (comparación)
        E_puntual = np.full(len(R2), np.nan)
        SP = np.empty((len(R2), keep))
        t0 = time.perf_counter()
        mask_man = np.array([s[0] >= 3 for s in estados], dtype=bool)
        masks_neigh = {li: np.array([s[0] == li for s in estados], dtype=bool)
                       for li in (0, 1, 2)}
        v_prev = None
        n_sin_caracter = 0
        for i, r2 in enumerate(R2):
            w, V = np.linalg.eigh(np.diag(diag) + mol.hmol.build(blk, float(r2)) + Vf)
            SP[i] = (w[:keep] - mol.E_manifold) * GHZ_PER_HARTREE

            pm = np.sum(V[mask_man] ** 2, axis=0)
            pn = np.column_stack([np.sum(V[masks_neigh[li]] ** 2, axis=0)
                                  for li in (0, 1, 2)])

            # criterio puntual: autoestado más bajo con carácter de manifold
            ok = np.flatnonzero(pm > args.weight)
            kp = int(ok[int(np.argmin(w[ok]))]) if len(ok) else -1
            K_puntual[i] = kp
            if kp >= 0:
                E_puntual[i] = (w[kp] - mol.E_manifold) * GHZ_PER_HARTREE

            # seguimiento adiabático: en R2[0], continuación del estado de
            # referencia sin Fermi; después, máximo solapamiento con el punto
            # anterior
            if i == 0:
                ovs = np.abs(V.T @ v_ref)
                j = int(np.argmax(ovs))
                ovl = float(ovs[j])
            else:
                ovs = np.abs(V.T @ v_prev)
                j = int(np.argmax(ovs))
                ovl = float(ovs[j])

            # el seguimiento por solapamiento NO se interrumpe si el estado
            # pierde carácter: se registra y se avisa (intercambio con el
            # canal 38s, umbral a −20.267 GHz)
            if j < 0:
                E[i] = np.nan
                print(f"    !! R2={r2:.0f}: seguimiento roto", flush=True)
            else:
                E[i] = (w[j] - mol.E_manifold) * GHZ_PER_HARTREE
                K[i], Wman[i], OVL[i], Wneigh[i] = j, float(pm[j]), ovl, pn[j]
                if pm[j] <= args.weight:
                    n_sin_caracter += 1
                    print(f"    !! R2={r2:.0f}: el estado seguido cae al "
                          f"{pm[j]:.1%} de manifold "
                          f"(vecinos: 38s={pn[j][0]:.2f}, 37p={pn[j][1]:.2f}, "
                          f"36d={pn[j][2]:.2f})", flush=True)
            v_prev = V[:, j].copy() if j >= 0 else None
            if i % 10 == 0:
                print(f"    R2 = {r2:8.1f} ({i+1}/{len(R2)})  "
                      f"E = {E[i]:+9.4f} GHz  peso={Wman[i]:.3f}  "
                      f"ovl={OVL[i]:.4f}", flush=True)
        dt = time.perf_counter() - t0

        diverge = np.flatnonzero(
            (K != K_puntual) & (K_puntual >= 0) & ~np.isnan(E_puntual))
        mins = local_minima(E)
        i_deep = int(np.nanargmin(E))
        print(f"    {len(R2)} puntos en {dt:.1f} s ({dt/len(R2):.2f} s/punto)")
        print(f"    puntos donde el estado seguido cae <= {args.weight:.0%} "
              f"de manifold: {n_sin_caracter}")
        print(f"    seguimiento vs criterio puntual: difieren en {len(diverge)}"
              f"{': ' + str([int(R2[i]) for i in diverge[:8]]) if len(diverge) else ''}")
        print(f"    solapamiento mínimo entre puntos consecutivos: "
              f"{OVL.min():.6f}")
        print(f"    E: rango [{np.nanmin(E):+.4f}, {np.nanmax(E):+.4f}] GHz; "
              f"mínimo global {E[i_deep]:+.4f} GHz en R2 = {R2[i_deep]:.0f} a0")
        print(f"    mínimos locales: {len(mins)} en {[int(R2[i]) for i in mins]} a0")
        print(f"    pesos vecinos en el extremo R2={R2[-1]:.0f}: "
              + ", ".join(f"{lb}={wn:.4f}"
                          for lb, wn in zip(NEIGH_LABELS, Wneigh[-1])))
        curves[R1] = {"R2": R2, "E": E, "K": K, "Wman": Wman,
                      "Wneigh": Wneigh, "overlap": OVL, "spectrum": SP}
        base = f"hybrid_caracter_R1{int(R1)}_n{args.n_manifold}_Nmax{args.n_max}"
        np.savez(out_dir / f"{base}.npz",
                 R2=R2, E=E, K=K, Wman=Wman, Wneigh=Wneigh, overlap=OVL,
                 spectrum=SP, K_puntual=K_puntual, E_puntual=E_puntual,
                 R1=R1, M_J=args.mj, N_max=args.n_max, weight=args.weight)
        print(f"    datos en {out_dir}/{base}.npz")

    # ------------------------------------------------------------- figura
    fig, ax = plt.subplots(figsize=(7.4, 5.4))
    colors = {600.0: "tab:blue", 900.0: "tab:red", 1100.0: "tab:green"}
    allE = np.concatenate([d["E"] for d in curves.values()])
    lo, hi = np.nanmin(allE) - 8.0, np.nanmax(allE) + 8.0
    for R1, d in curves.items():
        vis = np.any((d["spectrum"] > lo) & (d["spectrum"] < hi), axis=0)
        ax.plot(d["R2"], d["spectrum"][:, vis], color="0.85", lw=0.4, zorder=1)
        ax.plot(d["R2"], d["E"], lw=2.2, color=colors.get(R1),
                label=f"$R_1 = {R1:.0f}\\,a_0$ (neutro $\\theta=\\pi$)")
    ax.axhline(0.0, color="k", lw=0.8, ls=":")
    ax.set_xlabel(r"$R_2$  [$a_0$]  (distancia al RbCs, $\theta=0$)")
    ax.set_ylabel(r"$E_{\mathrm{car}} - E_{n=35,\,l\geq3}$  [GHz]")
    ax.set_title(f"Híbrido Rb*(n={args.n_manifold})-Rb-RbCs: curva más baja de "
                 f"carácter manifold\n"
                 f"(>{args.weight:.0%} de peso, seguimiento adiabático, "
                 f"$M_J={args.mj}$, $N_{{max}}={args.n_max}$)")
    ax.grid(alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    png = out_dir / f"hybrid_curves_caracter_MJ{args.mj}_n{args.n_manifold}.png"
    fig.savefig(png, dpi=150)
    print(f"\n  PNG en {png}")
    print(RULE)


if __name__ == "__main__":
    main()
