"""Pruebas aisladas de la continuación adaptativa de las BOP."""

import importlib.util
from pathlib import Path
from types import SimpleNamespace

import numpy as np


SCRIPT = Path(__file__).parents[3] / "scripts" / "compute_double_rbcs_curves.py"
SPEC = importlib.util.spec_from_file_location("compute_double_rbcs_curves", SCRIPT)
CURVES = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CURVES)


class RotatingBranchSystem:
    """Rama sintética cuyo autovector rota de forma controlada con R."""

    E_manifold = 0.0

    def __init__(self, angular_rate):
        self.angular_rate = angular_rate
        self.calls = []

    def state(self, radius):
        angle = self.angular_rate * (2.0 - radius)
        return np.array([np.cos(angle), np.sin(angle), 0.0])

    def is_manifold(self, mj):
        return np.ones(3, dtype=bool)

    def orientation_matrix(self, rotor, mj):
        return np.zeros((3, 3))

    def manifold_seed(self, radius, mj, geometry, separation, field):
        return 0.0, self.state(radius)

    def solve_near(
        self, radius, mj, geometry, separation, field, k, sigma_ghz,
    ):
        self.calls.append(float(radius))
        vectors = np.column_stack((self.state(radius), np.array([0.0, 0.0, 1.0])))
        values = np.array([radius, radius + 1.0]) / CURVES.GHZ_PER_HARTREE
        return values, vectors


def args(**overrides):
    values = dict(
        separation=300.0, weight=0.5, overlap=0.9, k=4, max_k=4,
        context=2, sigma_ghz=-30.0, min_substep=0.01,
        bisect_max_depth=6, catastrophic_overlap=0.15,
    )
    values.update(overrides)
    return SimpleNamespace(**values)


def comparable(result):
    return {
        key: value for key, value in result.items()
        if key != "elapsed_seconds"
    }


def test_sweep_is_identical_when_direct_overlap_passes():
    grid = np.array([0.0, 2.0])
    disabled_system = RotatingBranchSystem(angular_rate=0.1)
    adaptive_system = RotatingBranchSystem(angular_rate=0.1)

    disabled = CURVES.sweep(
        disabled_system, grid, 0, "symmetric",
        args(bisect_max_depth=0), 0.0,
    )
    adaptive = CURVES.sweep(
        adaptive_system, grid, 0, "symmetric", args(), 0.0,
    )

    assert disabled_system.calls == adaptive_system.calls == [2.0, 0.0]
    for key, expected in comparable(disabled).items():
        assert np.array_equal(comparable(adaptive)[key], expected, equal_nan=True)


def test_bisection_reduces_low_overlap_grid_points():
    grid = np.array([0.0, 2.0])
    # cos(theta)^2=0.8 para el hueco completo y ~0.947 para cada mitad.
    angular_rate = np.arccos(np.sqrt(0.8)) / 2.0
    disabled = CURVES.sweep(
        RotatingBranchSystem(angular_rate), grid, 0, "symmetric",
        args(bisect_max_depth=0), 0.0,
    )
    adaptive_system = RotatingBranchSystem(angular_rate)
    adaptive = CURVES.sweep(
        adaptive_system, grid, 0, "symmetric", args(), 0.0,
    )

    assert np.count_nonzero(disabled["overlap"] < 0.9) == 1
    assert np.count_nonzero(adaptive["overlap"] < 0.9) == 0
    assert 1.0 in adaptive_system.calls
    assert np.array_equal(adaptive["R"], grid)
    assert adaptive["E"].shape == disabled["E"].shape == grid.shape


class LostBranchSystem(RotatingBranchSystem):
    """Hace que la continuación elija una rama rota que la semilla recupera."""

    def __init__(self):
        super().__init__(angular_rate=0.0)

    def manifold_seed(self, radius, mj, geometry, separation, field):
        state = (
            np.array([1.0, 0.0, 0.0])
            if radius == 2.0 else np.array([0.0, 1.0, 0.0])
        )
        energy = 0.0 if radius == 2.0 else -10.0 / CURVES.GHZ_PER_HARTREE
        return energy, state

    def solve_near(
        self, radius, mj, geometry, separation, field, k, sigma_ghz,
    ):
        self.calls.append(float(radius))
        if radius == 2.0:
            vectors = np.eye(3)
            energies = np.array([0.0, 1.0, 2.0])
        else:
            correct = np.array([0.0, 1.0, 0.0])
            broken = np.array([np.sqrt(0.1), 0.0, np.sqrt(0.9)])
            distractor = np.array([0.0, 1.0, 0.0])
            vectors = np.column_stack((correct, broken, distractor))
            energies = np.array([-10.0, -1.0, 1.0])
        return energies / CURVES.GHZ_PER_HARTREE, vectors


def test_catastrophic_overlap_reseeds_lost_branch():
    result = CURVES.sweep(
        LostBranchSystem(), np.array([1.0, 2.0]), 0, "symmetric",
        args(overlap=0.7, catastrophic_overlap=0.15, bisect_max_depth=0),
        0.0,
    )

    assert result["resembrado"].tolist() == [True, False]
    assert result["W"][0] > 0.5
    assert result["overlap"][0] == 1.0
    assert result["E"][0] == -10.0
