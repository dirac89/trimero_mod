"""
PUNTO 4 de la ronda híbrida: tests de límites ANTES de mirar física nueva.

    a) R1→∞ (V_Fermi apagado EXACTO): el espectro debe coincidir BIT A BIT
       con rb_krb_polar puro reconstruido con d/B del RbCs.
    b) H_mol apagado EXACTO: el subespacio N=0, M_N=0 debe coincidir BIT A
       BIT con rb_neutral_perturber puro en n=35 (un perturbador en θ=π,
       cuyo espectro es el del dímero θ=0 del trímero conjugado por paridad).
    c) Hermiticidad de H_A + H_mol + V_Fermi(π).
    d) Control positivo: con AMBOS acoplamientos activos el espectro se
       aparta de ambos límites.

Las referencias se construyen con objetos FRESCOS dentro de cada test (no se
reutilizan los internos del híbrido): así la comparación valida el ensamblaje
completo, no un aliasing de instancias.

Si alguno falla: PARAR y reportar, no ajustar el test.

Ejecutar:  poetry run pytest tests/systems/hybrid_neutral_polar/test_limites_hibrido.py -s
"""

import time

import numpy as np
import pytest

from trimero.basis.quantum import CoupledBasis
from trimero.basis.radial import RadialBasis
from trimero.systems.hybrid_neutral_polar import HybridNeutralPolar
from trimero.systems.rb_krb_polar.charge_dipole import (
    ChargeDipoleHamiltonian,
    RydbergElectronField,
    rydberg_diagonal,
)
from trimero.systems.rb_neutral_perturber.linear_trimer import (
    SymmetricLinearTrimer,
)

N_MANIFOLD = 35
N_MAX = 2          # base reducida: suficiente para validar ensamblajes
R1 = 900.0         # Rb neutro (θ=π)
R2 = 800.0         # RbCs (θ=0)


@pytest.fixture(scope="module")
def hyb():
    return HybridNeutralPolar(n_manifold=N_MANIFOLD, N_max=N_MAX)


# ================================================================ TEST L-a
def test_la_limite_polar_bit_a_bit(hyb):
    """
    H(fermi=False) == diag(E_ryd) + H_mol^RbCs(R2), con la referencia
    montada desde cero con objetos frescos.
    """
    t0 = time.time()
    blk = hyb.block(0)
    H_hyb = hyb.hamiltonian(R1, R2, M_J=0, fermi=False)

    # Referencia independiente: base, radial, campo y H_mol nuevos.
    basis = CoupledBasis(N_max=N_MAX, manifold_l_min=3,
                         manifold_l_max=N_MANIFOLD - 1, neighbor_l=(0, 1, 2))
    radial = RadialBasis(n_manifold=N_MANIFOLD, l_min=3,
                         l_max=N_MANIFOLD - 1, neighbor_l=(0, 1, 2))
    hmol_ref = ChargeDipoleHamiltonian(
        B=hyb.B_au, d=hyb.d_au, electron_field=RydbergElectronField(radial))
    blk_ref = basis.get_block(0)
    assert blk_ref.states == blk.states, "las bases difieren"
    H_ref = (np.diag(rydberg_diagonal(blk_ref, n_manifold=N_MANIFOLD))
             + hmol_ref.build(blk_ref, R2))

    print(f"\n  L-a LÍMITE POLAR (R1={R1:g}→∞, R2={R2:g}), "
          f"bloque M_J=0 dim={len(blk)}  [{time.time()-t0:.1f}s]")
    print(f"    max|H_hyb - H_ref| = {np.max(np.abs(H_hyb - H_ref)):.3e}")
    assert np.array_equal(H_hyb, H_ref), "las matrices NO son bit a bit"

    w_hyb = np.linalg.eigvalsh(H_hyb)
    w_ref = np.linalg.eigvalsh(H_ref)
    print(f"    max|eigvalsh diff| = {np.max(np.abs(w_hyb - w_ref)):.3e}")
    assert np.array_equal(w_hyb, w_ref)
    print(f"    espectro polar puro reproducido BIT A BIT "
          f"(E[0]={w_hyb[0]:.12f} E_h)")


# ================================================================ TEST L-b
def test_lb_limite_neutro_bit_a_bit(hyb):
    """
    Subespacio N=0, M_N=0 del híbrido (mol=False) == dímero θ=π puro.

    Referencia: `SymmetricLinearTrimer` VALIDADO (radiales tabuladas,
    n_perturbers=1, F=0) conjugado por paridad. El conjugador se aplica como
    FLIP elemento a elemento con ±1.0 exactos —nunca D@M@D, que sumaría—.
    Si la matriz permutada coincide bit a bit, sus eigvalsh son idénticos.
    """
    t0 = time.time()
    m = 0
    blk = hyb.block(m)
    sub = [i for i, st in enumerate(blk.states) if st[2] == 0 and st[3] == 0]
    H_full = hyb.hamiltonian(R1, R2, M_J=m, mol=False)
    H_sub = H_full[np.ix_(sub, sub)]

    tri = SymmetricLinearTrimer(n_manifold=N_MANIFOLD,
                                radial_source="tabulated", n_perturbers=1)
    ls = tri.l_values(m)
    # Posición de cada l DENTRO del subespacio (no del bloque completo).
    sub_pos = {}
    for k, i in enumerate(sub):
        st = blk.states[i]
        assert st[3] == 0 and st[2] == 0
        sub_pos[st[0]] = k
    assert sorted(sub_pos) == ls, "el subespacio N=0 no cubre el bloque"

    # Matriz del híbrido reordenada al orden del trímero (permutación pura).
    orden = [sub_pos[l] for l in ls]
    H_sub_perm = H_sub[np.ix_(orden, orden)]

    # Dímero θ=0 validado y conjugación de paridad explícita.
    H_tri = tri.hamiltonian(R1, 0.0, m)
    signs = np.array([float((-1) ** l) for l in ls])
    H_pi_ref = signs[:, None] * H_tri * signs[None, :]

    print(f"\n  L-b LÍMITE NEUTRO (mol=False, R1={R1:g}), bloque m_l={m}, "
          f"dim={len(ls)}  [{time.time()-t0:.1f}s]")
    print(f"    max|H_sub - D·H_trimero·D| = "
          f"{np.max(np.abs(H_sub_perm - H_pi_ref)):.3e}")
    assert np.array_equal(H_sub_perm, H_pi_ref), (
        "las matrices NO son bit a bit")

    w = np.linalg.eigvalsh(H_sub_perm)
    w_ref = np.linalg.eigvalsh(H_pi_ref)
    assert np.array_equal(w, w_ref)
    ghz = (w - tri.E_manifold) * 6.579683920502e6
    print(f"    espectro neutro puro reproducido BIT A BIT")
    print(f"    E[0:3] rel. manifold [GHz] = {ghz[:3]}")
    # Control de no-vacuidad: el dímero θ=0 tiene que acoplar de verdad.
    assert np.max(np.abs(H_tri - np.diag(np.diag(H_tri)))) > 0.0


# ================================================================ TEST L-c
def test_lc_hermiticidad_completa(hyb):
    """||H - H^T||_F <= 1e-14 ||H||_F, convención del repo (con F_elec)."""
    print("\n  L-c HERMITICIDAD de H_A + H_mol + V_Fermi(π):")
    for (r1, r2, mj) in [(900.0, 800.0, 0), (600.0, 1200.0, 3),
                         (1100.0, 500.0, -2)]:
        H = hyb.hamiltonian(r1, r2, M_J=mj)
        asym = np.linalg.norm(H - H.T)
        rel = asym / np.linalg.norm(H)
        print(f"    R1={r1:g} R2={r2:g} M_J={mj:+d}: dim={H.shape[0]}, "
              f"||H-H^T||_F = {asym:.3e} (rel {rel:.1e})")
        assert asym <= 1e-14 * np.linalg.norm(H)


# ================================================================ TEST L-d
def test_ld_control_positivo_separa_limites(hyb):
    """
    Con ambos acoplamientos activos, el espectro del bloque M_J=0 se aparta
    de ambos límites por encima de 1e-8 E_h (~66 kHz). Se barren los tres
    R1 representativos de la ronda con R2 fijo.
    """
    print(f"\n  L-d CONTROL POSITIVO (ambos acoplamientos), R2={R2:g}:")
    print("    R1   | max|E_full-E_polar| | max|E_full-E_neutro| (E_h)")
    for r1 in (600.0, 900.0, 1100.0):
        w_full = hyb.eigvals(r1, R2, M_J=0)
        w_pol = hyb.eigvals(r1, R2, M_J=0, fermi=False)

        blk = hyb.block(0)
        sub = [i for i, st in enumerate(blk.states)
               if st[2] == 0 and st[3] == 0]
        w_neu = np.linalg.eigvalsh(
            hyb.hamiltonian(r1, R2, M_J=0, mol=False)[np.ix_(sub, sub)])
        # Ambos espectros viven en la misma base: la distancia entre tramas
        # ordenadas mide cuánto mueve el acoplamiento dipolar al límite
        # neutro. El polar es comparación directa elemento a elemento.
        n = min(len(w_full), len(w_neu))
        d_neu = np.max(np.abs(np.sort(w_full)[:n] - np.sort(w_neu)[:n]))
        d_pol = np.max(np.abs(w_full - w_pol))
        print(f"    {r1:5.0f} | {d_pol:.6e}          | {d_neu:.6e}")
        assert d_pol > 1e-8, (r1, d_pol)
        assert d_neu > 1e-8, (r1, d_neu)
