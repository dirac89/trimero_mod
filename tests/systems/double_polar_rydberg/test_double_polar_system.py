import numpy as np
import pytest

from trimero.systems.double_polar_rydberg import (
    ContractedRbTwoRbCsSystem,
    RbTwoRbCsSystem,
    TwoRotorBasis,
)
from trimero.systems.polar_molecule import RBCS
from trimero.systems.polar_rydberg import PolarBOPSystem


@pytest.fixture(scope="module")
def system():
    return RbTwoRbCsSystem(n_manifold=6, N_max=1)


def test_basis_enumeration_and_mj():
    basis = TwoRotorBasis(1, range(4))
    rot = [(N, M) for N in range(2) for M in range(-N, N + 1)]
    for mj in (0, 1):
        block = basis.get_block(mj)
        brute = sum(
            m + M1 + M2 == mj
            for l in range(4) for m in range(-l, l + 1)
            for _, M1 in rot for _, M2 in rot
        )
        assert len(block) == brute
        assert all(st.m + st.M1 + st.M2 == mj for st in block.states)


@pytest.mark.parametrize("mj", [0, 1])
def test_every_component_is_hermitian(system, mj):
    for matrix in system.components(500.0, mj, "symmetric", 300.0, 100.0).values():
        assert np.max(np.abs((matrix - matrix.T).toarray())) < 1e-12


def test_zero_external_field_is_exact_zero(system):
    matrix = system.external_field_matrix(0.0, 0)
    assert matrix.nnz == 0


def test_collinear_dipole_dipole_known_matrix_element(system):
    block = system.block(0)
    index = {state: i for i, state in enumerate(block.states)}
    i = index[(3, 0, 0, 0, 0, 0)]
    j = index[(3, 0, 1, 0, 1, 0)]
    distance = 300.0
    got = system.dipole_dipole_matrix(distance, 0)[i, j]
    expected = -2.0 * system.molecule.d_au**2 / (3.0 * distance**3)
    assert got == pytest.approx(expected, abs=1e-15)


def test_disabling_vdd_removes_only_that_component():
    with_vdd = RbTwoRbCsSystem(n_manifold=6, N_max=1, include_dipole_dipole=True)
    without = RbTwoRbCsSystem(n_manifold=6, N_max=1, include_dipole_dipole=False)
    parts = with_vdd.components(500, 0, "unilateral", 300, 0)
    expected = sum((value for key, value in parts.items() if key != "dipole_dipole"))
    assert np.array_equal(
        without.hamiltonian(500, 0, "unilateral", 300, 0).toarray(),
        expected.toarray(),
    )


def test_one_rotor_subspace_reproduces_polar_system(system):
    polar = PolarBOPSystem(molecule=RBCS, n_manifold=6, N_max=1)
    double_block = system.block(0)
    keep = [i for i, st in enumerate(double_block.states) if st.N2 == st.M2 == 0]
    polar_index = {state: i for i, state in enumerate(polar.block(0).states)}
    order = [polar_index[(double_block.states[i].l, double_block.states[i].m,
                          double_block.states[i].N1, double_block.states[i].M1)]
             for i in keep]
    double_h = (
        system.atomic_rotational_matrix(0)
        + system.charge_dipole_matrix(1, 700, 0)
    ).toarray()[np.ix_(keep, keep)]
    polar_h = polar.hamiltonian(700, 0)[np.ix_(order, order)]
    assert np.allclose(double_h, polar_h, atol=1e-13)


def test_exchange_molecules_maps_symmetric_components(system):
    block = system.block(0)
    index = {st: i for i, st in enumerate(block.states)}
    permutation = np.array([
        index[(st.l, st.m, st.N2, st.M2, st.N1, st.M1)] for st in block.states
    ])
    h1 = system.charge_dipole_matrix(1, -500, 0).toarray()
    h2 = system.charge_dipole_matrix(2, -500, 0).toarray()
    assert np.allclose(h1[np.ix_(permutation, permutation)], h2, atol=1e-13)


def test_negative_axis_parity_relation(system):
    positive = system.charge_dipole_matrix(1, 500, 0).toarray()
    negative = system.charge_dipole_matrix(1, -500, 0).toarray()
    states = system.block(0).states
    for i, si in enumerate(states):
        for j, sj in enumerate(states):
            if si.l != sj.l or si.m != sj.m:
                assert negative[i, j] == pytest.approx(
                    (-1.0) ** (si.l + sj.l + 1) * positive[i, j], abs=1e-13
                )


def test_sparse_solver_matches_dense(system):
    H = system.hamiltonian(700, 0, "unilateral", 300, 100)
    dense = np.linalg.eigvalsh(H.toarray())
    values, _ = system.solve_near(700, 0, "unilateral", 300, 100, k=H.shape[0])
    assert np.allclose(values, dense, atol=1e-12)
    lowest, _ = system.solve_lowest(700, 0, "unilateral", 300, 100, k=12)
    assert np.allclose(lowest, dense[:12], atol=1e-12)


def test_orientation_bounds(system):
    values, vectors = system.solve_near(700, 0, "symmetric", 300, 0, k=20)
    for rotor in (1, 2):
        operator = system.orientation_matrix(rotor, 0)
        expectations = np.einsum("ij,ij->j", vectors, operator @ vectors)
        assert np.all(np.abs(expectations) <= 1.0 + 1e-12)


@pytest.mark.parametrize("geometry", ["symmetric", "unilateral"])
def test_untruncated_local_basis_reproduces_primitive_spectrum(geometry):
    primitive = RbTwoRbCsSystem(n_manifold=6, N_max=1)
    contracted = ContractedRbTwoRbCsSystem(
        n_manifold=6, primitive_N_max=1, rotor_keep=2,
    )
    point = contracted.point(700, 0, geometry, 300, 100)
    hp = primitive.hamiltonian(700, 0, geometry, 300, 100)
    hc = contracted.hamiltonian(point)
    assert hp.shape == hc.shape
    assert np.max(np.abs((hc - hc.T).toarray())) < 1e-12
    assert np.allclose(
        np.linalg.eigvalsh(hp.toarray()),
        np.linalg.eigvalsh(hc.toarray()),
        atol=2e-12,
    )


def test_contracted_characteristic_orientation_bounds():
    system = ContractedRbTwoRbCsSystem(
        n_manifold=6, primitive_N_max=2, rotor_keep=1,
    )
    result = system.characteristic(system.point(700, 0, "symmetric", 300, 0), k=12)
    assert -1 <= result["COS1"] <= 1
    assert -1 <= result["COS2"] <= 1
    assert 0 <= result["W"] <= 1


def test_complete_natural_basis_reproduces_primitive_spectrum():
    primitive = RbTwoRbCsSystem(n_manifold=6, N_max=1)
    H = primitive.hamiltonian(700, 0, "symmetric", 300, 0)
    _, vectors = np.linalg.eigh(H.toarray())
    contracted = ContractedRbTwoRbCsSystem(
        n_manifold=6, primitive_N_max=1, rotor_keep=2,
    )
    point = contracted.natural_point(
        primitive.block(0), vectors[:, 0], 700, 0, "symmetric", 300, 0,
        natural_keep=2,
    )
    natural_h = contracted.hamiltonian(point)
    assert natural_h.shape == H.shape
    assert np.allclose(
        np.linalg.eigvalsh(natural_h.toarray()),
        np.linalg.eigvalsh(H.toarray()),
        atol=2e-12,
    )


def test_transferred_natural_basis_updates_local_hamiltonian():
    source = RbTwoRbCsSystem(n_manifold=6, N_max=1)
    block = source.block(0)
    vector = np.zeros(len(block))
    vector[0] = 1.0
    contracted = ContractedRbTwoRbCsSystem(
        n_manifold=6, primitive_N_max=2, rotor_keep=2,
    )
    anchor = contracted.natural_point(block, vector, 700, natural_keep=2)
    transferred = contracted.transfer_point(anchor, 800)

    assert transferred.block is anchor.block
    for old, new in ((anchor.rotor1, transferred.rotor1),
                     (anchor.rotor2, transferred.rotor2)):
        for M in old.coefficients:
            assert np.array_equal(old.coefficients[M], new.coefficients[M])
    assert any(
        not np.allclose(anchor.rotor1.hamiltonians[M],
                        transferred.rotor1.hamiltonians[M])
        for M in anchor.rotor1.hamiltonians
    )
    matrix = contracted.hamiltonian(transferred)
    assert np.max(np.abs((matrix - matrix.T).data)) < 1e-12
