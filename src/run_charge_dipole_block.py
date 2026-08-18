#!/usr/bin/env python3
"""
Construcción y diagonalización de H = H_a + H_mol para un bloque M_J.

    H_mol = B·N² - d·F_ion(R),   F_ion = e·R/R³        (término del ion Rb⁺)

El término del campo del electrón Rydberg (e·(r-R)/|r-R|³, Ec. A.6-A.10)
NO está incluido: implementación parcial a propósito.

Ejecutar:  poetry run python src/run_charge_dipole_block.py
"""

import numpy as np

from charge_dipole import (
    B_KRB_AU,
    D_KRB_AU,
    HZ_PER_HARTREE,
    ChargeDipoleHamiltonian,
    rydberg_diagonal,
)
from quantum_basis import CoupledBasis

EH_TO_GHZ = HZ_PER_HARTREE / 1.0e9


def total_hamiltonian(h, block, R):
    """H = H_a (diagonal Rydberg) + H_mol (B·N² - d·F_ion)."""
    return np.diag(rydberg_diagonal(block)) + h.build(block, R)


def main():
    basis = CoupledBasis(N_max=6, manifold_l_min=3)
    block = basis.get_block(0)
    h = ChargeDipoleHamiltonian()

    print("=" * 72)
    print("H = H_a + B·N² - d·F_ion(R)   —   bloque M_J = 0")
    print("=" * 72)
    print(f"  dim(bloque M_J=0) = {len(block)}")
    print(f"  B = {B_KRB_AU:.6e} E_h ({B_KRB_AU * EH_TO_GHZ:.4f} GHz)")
    print(f"  d = {D_KRB_AU:.6e} e·a0")
    print("  término del electrón Rydberg: NO incluido (ronda siguiente)")

    # ---- R intermedio --------------------------------------------------
    R = 1500.0
    H = total_hamiltonian(h, block, R)
    ev = np.linalg.eigvalsh(H)

    print(f"\n--- R = {R:g} a0 -------------------------------------------------")
    print(f"  |F_ion| = 1/R²  = {1.0 / R**2:.6e} a.u.")
    print(f"  d/R²            = {D_KRB_AU / R**2:.6e} E_h")
    print(f"  ||H||_F         = {np.linalg.norm(H):.9e}")
    print(f"  simétrica       = {np.array_equal(H, H.T)}")
    print(f"  autovalores: {len(ev)}, todos reales = {np.all(np.isreal(ev))}")
    print(f"  E_min = {ev.min():.15e} E_h")
    print(f"  E_max = {ev.max():.15e} E_h")
    print("\n  10 autovalores más bajos [E_h]         [GHz rel. a E_min]")
    for k in range(10):
        print(f"    {k:2d}  {ev[k]:.15e}   {(ev[k] - ev[0]) * EH_TO_GHZ:12.6f}")

    # ---- repetición del test 1 con esta implementación parcial ---------
    ref = np.sort(rydberg_diagonal(block) + h.rotational_diagonal(block))

    print("\n--- Test 1 (límite R grande) con la implementación parcial ------")
    print("     R [a0]   |  max|E - (E_Ryd + B N(N+1))| [E_h] | shift·R⁴ [E_h·a0⁴]")
    print("  ------------|-------------------------------------|-------------------")
    for R_scan in (500.0, 1000.0, 1500.0, 3000.0, 10000.0):
        evs = np.sort(np.linalg.eigvalsh(total_hamiltonian(h, block, R_scan)))
        dmax = np.abs(evs - ref).max()
        print(f"   {R_scan:10.1f} | {dmax:35.6e} | {dmax * R_scan**4:.6e}")

    R_big = 10000.0
    evb = np.sort(np.linalg.eigvalsh(total_hamiltonian(h, block, R_big)))
    dmax = np.abs(evb - ref).max()
    print(f"\n  R = {R_big:g} a0:  Δmax = {dmax:.6e} E_h   (tolerancia 1e-10)")
    print(f"  {'✓ PASA' if dmax < 1e-10 else '✗ FALLA'}: se recupera E_Ryd + B·N(N+1)")
    print("\n  El producto shift·R⁴ es aprox. constante: confirma que el shift es")
    print("  de SEGUNDO orden en un acoplamiento que cae como 1/R², como debe.")
    print("=" * 72)


if __name__ == "__main__":
    main()
