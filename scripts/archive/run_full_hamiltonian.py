#!/usr/bin/env python3
"""
Escaneo en R del Hamiltoniano COMPLETO H = H_a + B·N² - d·F_ryd(R,r),
con los DOS términos de F_ryd (Ec. 4):

    F_ryd = e·R/R³            (ion Rb⁺)
          + e·(r-R)/|r-R|³    (electrón Rydberg, expansión Ec. A.6-A.10)

y comparación de la magnitud de ambos términos.

Ejecutar:  poetry run python src/run_full_hamiltonian_scan.py
"""

import time

import numpy as np

from trimero.systems.rb_krb_polar.charge_dipole import (
    B_KRB_AU,
    D_KRB_AU,
    HZ_PER_HARTREE,
    ChargeDipoleHamiltonian,
    RydbergElectronField,
    rydberg_diagonal,
)
from trimero.basis.quantum import CoupledBasis
from trimero.basis.radial import RadialBasis

EH_TO_GHZ = HZ_PER_HARTREE / 1.0e9
EH_TO_MHZ = HZ_PER_HARTREE / 1.0e6


def main():
    basis = CoupledBasis(N_max=6, manifold_l_min=3)
    block = basis.get_block(0)
    radial = RadialBasis()
    efield = RydbergElectronField(radial)

    h_ion = ChargeDipoleHamiltonian()                          # sólo ion
    h_all = ChargeDipoleHamiltonian(electron_field=efield)     # ion + electrón

    E_ryd = rydberg_diagonal(block)
    E_rot = h_all.rotational_diagonal(block)
    ref = np.sort(E_ryd + E_rot)                               # espectro desacoplado

    print("=" * 78)
    print("H = H_a + B·N² - d·F_ryd(R,r)   —   bloque M_J = 0   —   AMBOS términos")
    print("=" * 78)
    print(f"  dim = {len(block)};  B = {B_KRB_AU:.6e} E_h;  d = {D_KRB_AU:.6e} e·a0")
    print(f"  manifold Rydberg: n=24 (l=3..23) + 27s;  <r>_{{24,l}} ≈ 588-864 a0")
    print(f"  punto de retorno externo clásico 2n² = {2*24**2} a0")

    # ------------------------------------------------------------------
    # 1) Magnitud relativa de los dos términos
    # ------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("1) MAGNITUD DE CADA TÉRMINO DEL ACOPLAMIENTO  (V = H_mol - B·N² diagonal)")
    print("-" * 78)
    print("      R    |  ||V_ion||_F  |  ||V_ele||_F  | ratio ele/ion |  max|V_ion| |  max|V_ele|")
    print("   --------|---------------|---------------|---------------|-------------|------------")
    for R in (600.0, 900.0, 1200.0, 1500.0, 2000.0, 3000.0):
        H_ion = h_ion.build(block, R)
        H_all = h_all.build(block, R)
        diagrot = np.diag(E_rot)
        V_ion = H_ion - diagrot
        V_ele = H_all - H_ion                 # el término del electrón, aislado
        n_ion, n_ele = np.linalg.norm(V_ion), np.linalg.norm(V_ele)
        print(f"   {R:8.0f}| {n_ion:.7e} | {n_ele:.7e} | {n_ele/n_ion:13.3f} | "
              f"{np.abs(V_ion).max():.4e} | {np.abs(V_ele).max():.4e}")

    print("\n  Los dos términos se CANCELAN parcialmente (el monopolo del electrón")
    print("  anula el campo del ion), así que el cociente de normas de arriba NO mide")
    print("  quién domina. Lo que el ion no puede hacer en absoluto es mezclar l:")
    print("\n      R    | max|V_ion| (l'=l) | max|V_ele| (l'=l) | max|V_ele| (l'≠l)")
    print("   --------|-------------------|-------------------|------------------")
    l_of = np.array([st[0] for st in block.states])
    same_l = l_of[:, None] == l_of[None, :]
    for R in (900.0, 1500.0, 3000.0):
        H_ion = h_ion.build(block, R)
        H_all = h_all.build(block, R)
        V_ion = H_ion - np.diag(E_rot)
        V_ele = H_all - H_ion
        print(f"   {R:8.0f}| {np.abs(V_ion[same_l]).max():17.4e} | "
              f"{np.abs(V_ele[same_l]).max():17.4e} | {np.abs(V_ele[~same_l]).max():.4e}")
    print("\n  El acoplamiento del ion es estrictamente diagonal en l; el del electrón")
    print("  conecta l distintos DENTRO del manifold cuasi-degenerado n=24, que es")
    print("  de donde sale el efecto grande.")

    # ------------------------------------------------------------------
    # 2) Escaneo del espectro con el Hamiltoniano completo
    # ------------------------------------------------------------------
    print("\n" + "-" * 78)
    print("2) ESCANEO EN R DEL ESPECTRO COMPLETO (ion + electrón)")
    print("-" * 78)
    print("  Desplazamiento máximo del espectro respecto del caso desacoplado,")
    print("  con SÓLO el ion y con AMBOS términos. La última columna responde a")
    print("  'cuál domina y por cuánto'.\n")
    print("      R    | sólo ion [MHz] | ion+electrón [MHz] | cociente | shift·R⁶ (completo)")
    print("   --------|----------------|--------------------|----------|--------------------")
    for R in (600.0, 900.0, 1200.0, 1500.0, 2000.0, 3000.0, 6000.0, 12000.0):
        ev_i = np.sort(np.linalg.eigvalsh(np.diag(E_ryd) + h_ion.build(block, R)))
        ev_a = np.sort(np.linalg.eigvalsh(np.diag(E_ryd) + h_all.build(block, R)))
        di = np.abs(ev_i - ref).max()
        da = np.abs(ev_a - ref).max()
        print(f"   {R:8.0f}| {di*EH_TO_MHZ:14.4f} | {da*EH_TO_MHZ:18.4f} | "
              f"{da/di:8.2f} | {da*R**6:.3e}")

    # ------------------------------------------------------------------
    # 3) Espectro a R típico
    # ------------------------------------------------------------------
    R = 1500.0
    print("\n" + "-" * 78)
    print(f"3) ESPECTRO A R = {R:g} a0 (bloque M_J=0, {len(block)} autovalores)")
    print("-" * 78)
    t0 = time.perf_counter()
    H = np.diag(E_ryd) + h_all.build(block, R)
    ev_all = np.sort(np.linalg.eigvalsh(H))
    ev_ion = np.sort(np.linalg.eigvalsh(np.diag(E_ryd) + h_ion.build(block, R)))
    print(f"  construido y diagonalizado en {time.perf_counter()-t0:.2f} s")
    asym = np.linalg.norm(H - H.T) / np.linalg.norm(H)
    print(f"  asimetría relativa ||H-H.T||/||H|| = {asym:.2e}  "
          f"(con el término del electrón ya no es bit a bit: hay redondeo en las")
    print(f"   sumas sobre k y en las integrales radiales)")
    print("\n   k |   E completo [E_h]   |  E sólo ion [E_h]    | dif [MHz] | E-E0 [GHz]")
    print("  ---|----------------------|----------------------|-----------|-----------")
    for k in list(range(8)):
        print(f"  {k:3d}| {ev_all[k]:.15e} | {ev_ion[k]:.15e} | "
              f"{(ev_all[k]-ev_ion[k])*EH_TO_MHZ:9.3f} | {(ev_all[k]-ev_all[0])*EH_TO_GHZ:10.4f}")
    print("=" * 78)


if __name__ == "__main__":
    main()
