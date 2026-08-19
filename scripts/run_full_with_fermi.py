#!/usr/bin/env python3
"""
Espectro del Hamiltoniano COMPLETO  H = H_a + H_mol + V_Fermi  para el bloque M_J=0.

    H_a      energías Rydberg (manifold n=24, l>=3, + 27s)
    H_mol    B·N² - d·F_ryd(R,r)   (ion Rb⁺ + electrón Rydberg)
    V_Fermi  pseudopotencial de Fermi s+p, con remapeo k(R) de las tablas n=35

Ejecutar:  poetry run python scripts/run_full_with_fermi.py
"""

import numpy as np

from trimero.basis.quantum import CoupledBasis
from trimero.basis.radial import RadialBasis
from trimero.hamiltonians.charge_dipole import (
    HZ_PER_HARTREE,
    ChargeDipoleHamiltonian,
    RydbergElectronField,
    rydberg_diagonal,
)
from trimero.hamiltonians.fermi_krb import (
    FermiPseudopotential,
    ScatteringLengths,
    n_star_of_l,
)

GHZ = HZ_PER_HARTREE / 1.0e9


def main():
    basis = CoupledBasis(N_max=6, manifold_l_min=3)
    block = basis.get_block(0)
    radial = RadialBasis()
    scat = ScatteringLengths()
    h_mol = ChargeDipoleHamiltonian(electron_field=RydbergElectronField(radial))
    fermi = FermiPseudopotential(radial, scat)

    E_ryd = rydberg_diagonal(block)
    E0_manifold = -0.5 / 24.0**2

    print("=" * 78)
    print("H = H_a + H_mol + V_Fermi   —   bloque M_J = 0   —   dim", len(block))
    print("=" * 78)
    print("  H_mol: ion Rb⁺ + electrón Rydberg (ambos verificados)")
    print("  V_Fermi: s+p, perturbador puntual sobre el eje, remapeo k(R) desde n=35")
    print(f"  cero de energía en las tablas: E(n=24) = {E0_manifold:.12e} E_h")

    for R in (600.0, 900.0, 1100.0, 1500.0):
        print("\n" + "-" * 78)
        if not scat.in_domain(R, n_star_of_l(3)):
            k = scat.k_of_R(R, n_star_of_l(3))
            print(f"R = {R:g} a0 —— FUERA DE DOMINIO, no se calcula")
            print(f"  k(R) = {k}  (punto de retorno externo de n=24: 2n² = 1152 a0)")
            print(f"  R' remapeado = {scat.remap_R(R, n_star_of_l(3)):.1f} a0, "
                  f"tabla hasta 2448 a0")
            print("  El electrón no tiene momento clásicamente permitido a esta")
            print("  distancia: no hay longitud de dispersión definida. No se extrapola.")
            continue

        A_s, A_p = scat.scattering(R, n_star_of_l(3))
        H = np.diag(E_ryd) + h_mol.build(block, R) + fermi.build(block, R)
        ev = np.sort(np.linalg.eigvalsh(H))
        print(f"R = {R:g} a0    A_s = {A_s:+.6f} a0    A_p = {A_p:+.4f} a0³")
        print(f"  ||H-H.T||/||H|| = {np.linalg.norm(H-H.T)/np.linalg.norm(H):.2e}")
        print("\n    k |        E [E_h]         | E - E(n=24) [GHz] |  E - E_0 [GHz]")
        print("   ---|------------------------|-------------------|----------------")
        for j in range(10):
            print(f"   {j:3d}| {ev[j]:.15e} | {(ev[j]-E0_manifold)*GHZ:17.6f} | "
                  f"{(ev[j]-ev[0])*GHZ:14.6f}")
    print("=" * 78)


if __name__ == "__main__":
    main()
