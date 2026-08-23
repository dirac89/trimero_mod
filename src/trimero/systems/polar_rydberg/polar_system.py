"""Sistema BOP polar genérico, sin pseudopotencial de Fermi."""

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from trimero.basis.quantum import CoupledBasis, QuantumBasisBlock
from trimero.basis.radial import RadialBasis
from trimero.systems.polar_molecule import KRB, PolarMolecule, get_molecule
from trimero.systems.rb_krb_polar.charge_dipole import (
    HZ_PER_HARTREE,
    ChargeDipoleHamiltonian,
    RydbergElectronField,
    rydberg_diagonal,
)
from trimero.systems.rb_krb_polar.rb_defects import n_star_nl, neighbor_levels

__all__ = ["GHZ_PER_HARTREE", "PolarBOPSystem"]

GHZ_PER_HARTREE = HZ_PER_HARTREE / 1.0e9


@dataclass
class PolarBOPSystem:
    """Hamiltoniano ``H_A + B N² - d·F_ryd`` para un rotor polar."""

    molecule: PolarMolecule = KRB
    n_manifold: int = 24
    n_s: Optional[int] = None
    delta0_ns: Optional[float] = None
    N_max: int = 6
    l_min: int = 3
    neighbors: tuple = (0, 1, 2)
    include_electron_field: bool = True

    basis: CoupledBasis = field(init=False)
    radial: RadialBasis = field(init=False)
    hmol: ChargeDipoleHamiltonian = field(init=False)

    def __post_init__(self) -> None:
        self.molecule = get_molecule(self.molecule)
        levels = neighbor_levels(self.n_manifold)
        if self.n_s is not None:
            levels[0] = self.n_s
        self.neighbor_l = tuple(sorted(l for l in levels if l in self.neighbors))
        self.levels = {l: levels[l] for l in self.neighbor_l}
        self.n_s = self.levels.get(0)
        self.l_max = self.n_manifold - 1
        self.basis = CoupledBasis(
            N_max=self.N_max,
            manifold_l_min=self.l_min,
            manifold_l_max=self.l_max,
            neighbor_l=self.neighbor_l,
        )
        self.radial = RadialBasis(
            n_manifold=self.n_manifold,
            l_min=self.l_min,
            l_max=self.l_max,
            neighbor_l=self.neighbor_l,
        )
        field_op = RydbergElectronField(self.radial) if self.include_electron_field else None
        self.hmol = ChargeDipoleHamiltonian(
            B=self.molecule.B_au, d=self.molecule.d_au, electron_field=field_op
        )
        self._blocks = {}

    @property
    def E_manifold(self) -> float:
        return -0.5 / self.n_manifold**2

    @property
    def E_ns(self) -> float:
        return -0.5 / self.n_star_s() ** 2

    def n_star_s(self) -> float:
        if 0 not in self.levels:
            raise ValueError("la base no contiene un vecino ns")
        return n_star_nl(self.levels[0], 0, self.delta0_ns)

    def delta_E_ghz(self, l: int) -> float:
        if l not in self.levels:
            return 0.0
        nstar = n_star_nl(self.levels[l], l, self.delta0_ns)
        return (-0.5 / nstar**2 - self.E_manifold) * GHZ_PER_HARTREE

    def thresholds_level_ghz(self, l: int, N_values):
        delta = self.delta_E_ghz(l)
        return {int(N): delta + self.molecule.B_ghz * N * (N + 1) for N in N_values}

    def thresholds_ns_ghz(self, N_values):
        return self.thresholds_level_ghz(0, N_values)

    def thresholds_manifold_ghz(self, N_values):
        return {int(N): self.molecule.B_ghz * N * (N + 1) for N in N_values}

    def block(self, M_J: int = 0) -> QuantumBasisBlock:
        return self.basis.get_block(M_J)

    def rydberg_diagonal(self, M_J: int = 0) -> np.ndarray:
        key = ("E", M_J)
        if key not in self._blocks:
            self._blocks[key] = rydberg_diagonal(
                self.block(M_J), self.n_manifold, self.n_s, self.delta0_ns
            )
        return self._blocks[key]

    def is_manifold(self, M_J: int = 0) -> np.ndarray:
        key = ("m", M_J)
        if key not in self._blocks:
            self._blocks[key] = np.array([s[0] >= self.l_min for s in self.block(M_J).states])
        return self._blocks[key]

    def hamiltonian(self, R: float, M_J: int = 0) -> np.ndarray:
        block = self.block(M_J)
        return np.diag(self.rydberg_diagonal(M_J)) + self.hmol.build(block, R)

    def eigvals(self, R: float, M_J: int = 0) -> np.ndarray:
        return np.linalg.eigvalsh(self.hamiltonian(R, M_J))

    def solve(self, R: float, M_J: int = 0):
        return np.linalg.eigh(self.hamiltonian(R, M_J))

    def character_curve(self, R: float, M_J: int = 0, weight: float = 0.5):
        values, vectors = self.solve(R, M_J)
        mask = self.is_manifold(M_J)
        for k, value in enumerate(values):
            manifold_weight = float(np.sum(vectors[:, k][mask] ** 2))
            if manifold_weight > weight:
                return ((value - self.E_manifold) * GHZ_PER_HARTREE, k, manifold_weight)
        return float("nan"), -1, 0.0
