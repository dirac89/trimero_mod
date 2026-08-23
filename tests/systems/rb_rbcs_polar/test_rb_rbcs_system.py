import numpy as np

from trimero.systems.hybrid_neutral_polar import HybridNeutralPolar
from trimero.systems.polar_molecule import KRB, RBCS
from trimero.systems.polar_rydberg import PolarBOPSystem
from trimero.systems.rb_krb_polar.bop_system import BOPSystem
from trimero.systems.rb_krb_polar.rb_defects import DELTA0_NS_PAPER
from trimero.systems.rb_rbcs_polar import RbRbCsPolarSystem


def test_krb_generic_core_matches_legacy_polar_matrix_bit_for_bit():
    params = dict(n_manifold=24, N_max=1, delta0_ns=DELTA0_NS_PAPER)
    generic = PolarBOPSystem(molecule=KRB, **params)
    legacy = BOPSystem(**params)
    assert generic.block(0).states == legacy.block(0).states
    for radius in (500.0, 1000.0):
        assert np.array_equal(
            generic.hamiltonian(radius, 0),
            legacy.hamiltonian(radius, 0, fermi=False),
        )


def test_rbcs_wrapper_uses_catalogue_and_rotational_thresholds():
    system = RbRbCsPolarSystem(n_manifold=25, N_max=1)
    assert system.molecule is RBCS
    assert system.hmol.B == RBCS.B_au
    assert system.hmol.d == RBCS.d_au
    thresholds = system.thresholds_manifold_ghz((0, 1, 2))
    assert thresholds == {0: 0.0, 1: 2 * RBCS.B_ghz, 2: 6 * RBCS.B_ghz}


def test_rbcs_polar_core_matches_hybrid_polar_limit_bit_for_bit():
    params = dict(n_manifold=25, N_max=1)
    polar = RbRbCsPolarSystem(**params)
    hybrid = HybridNeutralPolar(**params, radial_source="hydrogenic")
    assert polar.block(0).states == hybrid.block(0).states
    assert np.array_equal(
        polar.hamiltonian(800.0, 0),
        hybrid.hamiltonian(900.0, 800.0, M_J=0, fermi=False),
    )


def test_hybrid_custom_constants_reach_charge_dipole_hamiltonian():
    hybrid = HybridNeutralPolar(
        n_manifold=25, N_max=0, radial_source="hydrogenic",
        d_debye=0.75, B_mhz=321.0,
    )
    assert hybrid.hmol.d == hybrid.d_au
    assert hybrid.hmol.B == hybrid.B_au
