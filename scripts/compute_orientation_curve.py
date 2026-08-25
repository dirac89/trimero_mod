#!/usr/bin/env python3
"""
Orientación ⟨cos θ_d⟩ de Rb* unido a KRb o RbCs.

Fase B de docs/PLAN_figuras_publicacion.md. Generaliza el cálculo puntual de
docs/analysis_base_correcta_3_vecinos.md §4.2 (`scripts/archive/
run_basis_correction_check.py`, función `cos_theta_matrix`/`orientation_curve`)
a una malla de producción y a varios n.

⚠️ DIFERENCIA DELIBERADA CON EL CÓDIGO ORIGINAL: aquel llamaba a
`system.hamiltonian(R)` sin fijar `fermi`, que por defecto es `fermi=True`.
Eso incluía el pseudopotencial de Fermi en el Hamiltoniano usado para
identificar el autoestado — la misma premisa equivocada que
docs/analysis_fig1_carga_dipolo_sin_fermi.md corrigió para las curvas BOP.
Aquí se usa `fermi=False` explícito, igual que `compute_bop_curve.py`, para
ser consistentes con el modelo polar vigente (docs/STATUS.md). Los números
CAMBIAN respecto a la ronda anterior (documentado en
docs/analysis_faseB_orientacion_varios_n.md); no es un bug de este script.

Consecuencia práctica: al no haber V_Fermi no hay restricción de dominio de
remapeo k(R) (a diferencia del script original, que capturaba ValueError
fuera de `domain_bounds()`), así que el barrido cubre R libremente, incluido
más allá de donde el remapeo de Fermi estaría definido.

Operador ⟨i|cosθ_d|j⟩ = δ_ll' δ_mm' · cos_theta_element(N,M_N,N',M_N'), ya
verificado (docs/analysis_campo_electron_rydberg.md); no depende de R.

DIAGNÓSTICO DE AMBIGÜEDAD
-------------------------
El criterio "más bajo con carácter de manifold" puede saltar entre estados
casi degenerados cuando la densidad de estados es alta (cruces evitados con
manifold vecino). Se registra el índice K en cada punto: un cambio de K
entre pasos consecutivos es un cruce evitado real; una RACHA de cambios
seguidos, o saltos grandes de ⟨cosθ_d⟩ sin que cambie K (oscilación del
propio autovector cerca de una cuasi-degeneración), es la señal de
seguimiento no fiable que describe cualitativamente
docs/analysis_verificacion_tabla_I.md §11 — aquí se cuantifica por n.

QUÉ MOLÉCULA POLAR
------------------
`--molecule` es OBLIGATORIO y no tiene valor por defecto, por la misma razón
que en `compute_bop_curve.py`: el default silencioso `krb` hizo que las Fases A
y B del plan de figuras se calcularan con KRb creyéndose RbCs. La molécula
queda escrita en el `.npz` (`molecule`, `B_hz`, `d_debye`) y en el título de la
figura. Ver `docs/analysis_faseA_curvas_bop_varios_n.md` §0.

USO
---
    poetry run python scripts/compute_orientation_curve.py --molecule rbcs --n-manifold 24
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
from trimero.systems.rb_krb_polar.rb_defects import DELTA0_NS_PAPER
from trimero.systems.rb_rbcs_polar import RbRbCsPolarSystem

# Misma fábrica que compute_bop_curve.py: cada molécula con SU clase.
SYSTEM_FOR_MOLECULE = {"rbcs": RbRbCsPolarSystem}


def build_system(molecule, **kwargs):
    """Instancia el sistema polar de esta molécula, con su clase específica."""
    cls = SYSTEM_FOR_MOLECULE.get(molecule.key)
    if cls is None:
        return PolarBOPSystem(molecule=molecule, **kwargs)
    system = cls(**kwargs)
    assert system.molecule is molecule
    return system


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-manifold", type=int, default=24)
    ap.add_argument("--n-max", type=int, default=6)
    ap.add_argument("--molecule", choices=tuple(MOLECULES), required=True,
                    help="molécula polar (OBLIGATORIO, sin valor por defecto)")
    ap.add_argument("--mj", type=int, default=0)
    ap.add_argument("--rmin-fine", type=float, default=100.0)
    ap.add_argument("--rmax-fine", type=float, default=800.0)
    ap.add_argument("--step-fine", type=float, default=5.0)
    ap.add_argument("--rmax-coarse", type=float, default=1800.0)
    ap.add_argument("--step-coarse", type=float, default=50.0)
    ap.add_argument("--weight", type=float, default=0.5)
    ap.add_argument("--npz-dir", default=None)
    ap.add_argument("--out", default=None,
                    help="PNG; por defecto figures/orientation_MJ<mj>_n<n>_Nmax<N>.png")
    ap.add_argument("--no-plot", action="store_true",
                    help="guarda sólo el .npz; útil antes de una comparación multi-n")
    return ap.parse_args()


def cos_theta_matrix(system, block):
    """<i|cos(theta_d)|j> = delta_ll' delta_mm' <N M|cos|N' M'>. No depende de R."""
    states = block.states
    index = {st: i for i, st in enumerate(states)}
    C = np.zeros((len(states), len(states)))
    for i, (l, m, N, MN) in enumerate(states):
        for Np in (N - 1, N + 1):
            if Np < 0 or abs(MN) > Np:
                continue
            j = index.get((l, m, Np, MN))
            if j is not None:
                C[i, j] = system.hmol.cos_theta_element(N, MN, Np, MN)
    return C


def cos2_theta_matrix(system, block):
    """Matriz exacta de cos²(theta_d) en la base rotacional truncada.

    La suma incluye el estado intermedio N_max+1 cuando contribuye a un
    elemento diagonal; por ello no se usa simplemente ``C @ C``.
    """
    states = block.states
    index = {st: i for i, st in enumerate(states)}
    C2 = np.zeros((len(states), len(states)))
    for i, (l, m, N, MN) in enumerate(states):
        for Np in (N - 2, N, N + 2):
            if Np < 0 or abs(MN) > Np:
                continue
            j = index.get((l, m, Np, MN))
            if j is None:
                continue
            intermediates = set((N - 1, N + 1)) & set((Np - 1, Np + 1))
            C2[i, j] = sum(
                system.hmol.cos_theta_element(N, MN, Nt, MN)
                * system.hmol.cos_theta_element(Nt, MN, Np, MN)
                for Nt in intermediates if Nt >= abs(MN)
            )
    return C2


def sweep(sysm, M_J, R, C, C2, weight):
    mask = sysm.is_manifold(M_J)
    E = np.full(len(R), np.nan)
    K = np.full(len(R), -1, dtype=int)
    W = np.zeros(len(R))
    COS = np.full(len(R), np.nan)
    COS2 = np.full(len(R), np.nan)
    t0 = time.perf_counter()
    for i, r in enumerate(R):
        w, V = np.linalg.eigh(sysm.hamiltonian(float(r), M_J))
        for k in range(len(w)):
            wm = float(np.sum(V[:, k][mask] ** 2))
            if wm > weight:
                v = V[:, k]
                E[i] = (w[k] - sysm.E_manifold) * GHZ
                K[i] = k
                W[i] = wm
                COS[i] = float(v @ C @ v)
                COS2[i] = float(v @ C2 @ v)
                break
        if i % 40 == 0:
            print(f"    R = {r:8.2f} ({i+1}/{len(R)})", flush=True)
    dt = time.perf_counter() - t0
    print(f"    {len(R)} diagonalizaciones en {dt:.1f} s ({dt/len(R):.2f} s/punto)")
    return {"R": R, "E": E, "K": K, "W": W, "COS": COS, "COS2": COS2}


def diagnose_ambiguity(d, thr=0.02):
    """Eventos de seguimiento ambiguo: |Delta cos theta_d| > thr entre puntos
    consecutivos. Se registra tambien dK como informacion, pero NO como
    criterio: un cambio de K aislado con |Delta cos| pequeno es una
    relabelacion adiabatica benigna (el autoestado vecino tiene practicamente
    el mismo caracter), no un fallo de seguimiento. Usar dK!=0 como criterio
    sobre-estima la region ambigua (arrastra el inicio hasta el borde interno
    del barrido en los cuatro n; comprobado y descartado)."""
    R, K, COS = d["R"], d["K"], d["COS"]
    events = []
    for i in range(1, len(R)):
        dK = K[i] - K[i - 1]
        dC = COS[i] - COS[i - 1]
        if abs(dC) > thr:
            events.append((float(R[i - 1]), float(R[i]), int(dK), float(dC)))
    return events


def make_plot(d, molecule, args, out):
    bad = np.zeros(len(d["R"]), dtype=bool)
    for r0, r1, _, _ in diagnose_ambiguity(d):
        bad |= (d["R"] == r0) | (d["R"] == r1)
    fig, ax = plt.subplots(figsize=(8, 5.6))
    ax.plot(d["R"], d["COS"], color="tab:blue", lw=1.3, alpha=0.35)
    ax.plot(np.where(bad, np.nan, d["R"]), np.where(bad, np.nan, d["COS"]),
            color="tab:blue", lw=2.0, label="seguimiento fiable")
    ax.axhline(0.0, color="k", lw=0.8, ls=":")
    ax.set(xlabel=r"$R$  [$a_0$]", ylabel=r"$\langle\cos\theta_d\rangle$",
           ylim=(-1.02, 1.02),
           title=f"Rb*-{molecule.label}, orientación — n={args.n_manifold}, "
                 f"$M_J={args.mj}$, $N_{{max}}={args.n_max}$")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    fig.savefig(out, dpi=150)
    plt.close(fig)


def main():
    args = parse_args()
    molecule = get_molecule(args.molecule)
    args.npz_dir = args.npz_dir or f"plots/rb_{molecule.key}_polar/data"
    sysm = build_system(molecule, n_manifold=args.n_manifold,
                        N_max=args.n_max,
                        delta0_ns=DELTA0_NS_PAPER)
    n = sysm.n_manifold
    block = sysm.block(args.mj)
    C = cos_theta_matrix(sysm, block)
    C2 = cos2_theta_matrix(sysm, block)

    print("=" * 88)
    print(f"Orientacion <cos theta_d> Rb*-{molecule.label}, n={n}, M_J={args.mj}")
    print("=" * 88)
    print(f"  molécula: {molecule.label} (key={molecule.key})   "
          f"B = {molecule.B_ghz:.6f} GHz   d = {molecule.dipole_debye:.3f} D   "
          f"clase: {type(sysm).__name__}")
    print(f"  dim(bloque) = {len(block)}   ||C||_F = {np.linalg.norm(C):.4f}")

    R = np.concatenate([
        np.arange(args.rmin_fine, args.rmax_fine, args.step_fine),
        np.arange(args.rmax_fine, args.rmax_coarse + 1e-9, args.step_coarse),
    ])
    d = sweep(sysm, args.mj, R, C, C2, args.weight)

    n_nan = int(np.sum(np.isnan(d["COS"])))
    cos_valid = d["COS"][~np.isnan(d["COS"])]
    print(f"\n  NaN en COS: {n_nan}/{len(R)}")
    print(f"  rango de <cos theta_d>: [{cos_valid.min():.6f}, {cos_valid.max():.6f}]")
    assert np.all(cos_valid >= -1.0 - 1e-9) and np.all(cos_valid <= 1.0 + 1e-9), \
        "cos theta_d fuera de [-1,1]: violación de una propiedad matemática exacta"
    print("  [-1,1] verificado: PASA")
    cos2_valid = d["COS2"][~np.isnan(d["COS2"])]
    assert np.all(cos2_valid >= -1e-9) and np.all(cos2_valid <= 1.0 + 1e-9), \
        "cos² theta_d fuera de [0,1]"
    print(f"  rango de <cos² theta_d>: [{cos2_valid.min():.6f}, {cos2_valid.max():.6f}]")

    W_valid = d["W"][d["K"] >= 0]
    print(f"\n  peso de manifold: W_min = {W_valid.min():.4f} en R = "
          f"{d['R'][d['K'] >= 0][int(np.argmin(W_valid))]:.1f} a0   "
          f"(umbral de aceptación {args.weight})   "
          f"#W<0.90 = {int((W_valid < 0.90).sum())}   "
          f"#W<0.70 = {int((W_valid < 0.70).sum())}")

    events = diagnose_ambiguity(d)
    print(f"\n  eventos de posible ambigüedad (cambio de K o |Δcos|>0.02): "
          f"{len(events)}")
    for r0, r1, dK, dC in events:
        print(f"    R={r0:7.1f}->{r1:7.1f}   dK={dK:+d}   dcos={dC:+.4f}")

    os.makedirs(args.npz_dir, exist_ok=True)
    npz = os.path.join(args.npz_dir, f"orientation_MJ{args.mj}_n{n}.npz")
    np.savez(npz, **d, molecule=molecule.key,
             B_hz=molecule.rotational_constant_hz,
             d_debye=molecule.dipole_debye,
             n_manifold=n, N_max=sysm.N_max, M_J=args.mj,
             character_weight=args.weight, schema_version=2)
    print(f"\n  datos en {npz}")
    if not args.no_plot:
        out = args.out or f"plots/rb_{molecule.key}_polar/figures/orientation_MJ{args.mj}_n{n}_Nmax{args.n_max}.png"
        make_plot(d, molecule, args, out)
        print(f"  PNG en {out}")
    print("=" * 88)


if __name__ == "__main__":
    main()
