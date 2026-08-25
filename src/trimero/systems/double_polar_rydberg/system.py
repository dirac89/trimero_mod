"""Hamiltoniano disperso de Rb*+RbCs+RbCs en geometría axial."""

from dataclasses import dataclass, field

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix, diags
from scipy.sparse.linalg import eigsh

from trimero.basis.radial import RadialBasis
from trimero.mathlib.angular import gaunt
from trimero.systems.polar_molecule import RBCS, PolarMolecule, get_molecule
from trimero.systems.polar_rydberg import V_PER_M_PER_AU
from trimero.systems.rb_atom import Atom
from trimero.systems.rb_krb_polar.charge_dipole import (
    HZ_PER_HARTREE,
    ChargeDipoleHamiltonian,
    RydbergElectronField,
)
from trimero.systems.rb_krb_polar.rb_defects import (
    energy_rb,
    neighbor_levels,
)

from .basis import TwoRotorBasis, TwoRotorBasisBlock, TwoRotorState

GHZ_PER_HARTREE = HZ_PER_HARTREE / 1.0e9


@dataclass
class RbTwoRbCsSystem:
    """Dos rotores RbCs etiquetados, colineales con el núcleo Rydberg."""

    n_manifold: int = 20
    N_max: int = 3
    molecule: PolarMolecule = RBCS
    l_min: int = 3
    delta0_ns: float | None = None
    include_electron_field: bool = True
    include_dipole_dipole: bool = True

    basis: TwoRotorBasis = field(init=False)
    radial: RadialBasis = field(init=False)
    electron_field: RydbergElectronField | None = field(init=False)
    hmol: ChargeDipoleHamiltonian = field(init=False)

    def __post_init__(self):
        self.molecule = get_molecule(self.molecule)
        if self.n_manifold <= self.l_min:
            raise ValueError("n_manifold debe ser mayor que l_min")
        self.l_max = self.n_manifold - 1
        self.levels = neighbor_levels(self.n_manifold)
        self.l_values = tuple(sorted(set(range(self.l_min, self.l_max + 1)) | {0, 1, 2}))
        self.basis = TwoRotorBasis(self.N_max, self.l_values)
        self.radial = RadialBasis(
            n_manifold=self.n_manifold,
            l_min=self.l_min,
            l_max=self.l_max,
            neighbor_l=(0, 1, 2),
        )
        self.electron_field = (
            RydbergElectronField(self.radial) if self.include_electron_field else None
        )
        self.hmol = ChargeDipoleHamiltonian(
            B=self.molecule.B_au,
            d=self.molecule.d_au,
            electron_field=self.electron_field,
        )
        self._cache = {}

    @property
    def E_manifold(self) -> float:
        return -0.5 / self.n_manifold**2

    def block(self, M_J=0) -> TwoRotorBasisBlock:
        return self.basis.get_block(M_J)

    def is_manifold(self, M_J=0) -> np.ndarray:
        return np.fromiter(
            (state.l >= self.l_min for state in self.block(M_J).states), dtype=bool
        )

    def positions(self, R: float, geometry="symmetric", separation=300.0):
        if R <= 0 or separation <= 0:
            raise ValueError("R y separation deben ser positivos")
        if geometry == "symmetric":
            return (-float(R), float(R))
        if geometry == "unilateral":
            return (float(R), float(R + separation))
        raise ValueError("geometry debe ser 'symmetric' o 'unilateral'")

    @staticmethod
    def _sparse(entries, dim):
        if not entries:
            return csr_matrix((dim, dim), dtype=float)
        rows, cols, data = zip(*entries)
        return coo_matrix((data, (rows, cols)), shape=(dim, dim)).tocsr()

    @staticmethod
    def _rotor_c(Ni, Mi, Nj, Mj, q):
        if Mi != Mj + q:
            return 0.0
        return float(np.sqrt(4.0 * np.pi / 3.0) * gaunt(Ni, Mi, 1, q, Nj, Mj))

    def atomic_rotational_matrix(self, M_J=0):
        key = ("diag", M_J)
        if key not in self._cache:
            values = []
            energies = {}
            for st in self.block(M_J).states:
                if st.l not in energies:
                    n = self.levels.get(st.l, self.n_manifold)
                    energies[st.l] = energy_rb(n, st.l, self.delta0_ns)
                values.append(
                    energies[st.l]
                    + self.molecule.B_au
                    * (st.N1 * (st.N1 + 1) + st.N2 * (st.N2 + 1))
                )
            self._cache[key] = diags(values, format="csr")
        return self._cache[key]

    def charge_dipole_matrix(self, rotor: int, z: float, M_J=0):
        """Carga–dipolo de un rotor situado en z, sin energía rotacional."""
        if rotor not in (1, 2) or z == 0:
            raise ValueError("rotor debe ser 1 o 2 y z no puede ser cero")
        states = self.block(M_J).states
        index = {state: i for i, state in enumerate(states)}
        radius, sign = abs(float(z)), 1.0 if z > 0 else -1.0
        entries = []
        for i, si in enumerate(states):
            Ni, Mi = (si.N1, si.M1) if rotor == 1 else (si.N2, si.M2)
            for q in (-1, 0, 1):
                mj = si.m + q
                Mj = Mi - q
                for Nj in (Ni - 1, Ni + 1):
                    if Nj < 0 or Nj > self.N_max or abs(Mj) > Nj:
                        continue
                    for lj in self.l_values:
                        if abs(mj) > lj:
                            continue
                        if rotor == 1:
                            sj = TwoRotorState(lj, mj, Nj, Mj, si.N2, si.M2)
                        else:
                            sj = TwoRotorState(lj, mj, si.N1, si.M1, Nj, Mj)
                        j = index.get(sj)
                        if j is None:
                            continue
                        value = 0.0
                        if q == 0 and si.l == lj and si.m == mj:
                            value += -self.molecule.d_au * sign * self.hmol.cos_theta_element(
                                Ni, Mi, Nj, Mj
                            ) / radius**2
                        if self.electron_field is not None:
                            sti = (si.l, si.m, Ni, Mi)
                            stj = (lj, mj, Nj, Mj)
                            parity = 1.0 if sign > 0 else (-1.0) ** (si.l + lj + 1)
                            value += parity * self.hmol.electron_field_element(sti, stj, radius)
                        if value:
                            entries.append((i, j, value))
        return self._sparse(entries, len(states))

    def dipole_dipole_matrix(self, distance: float, M_J=0):
        if distance <= 0:
            raise ValueError("distance debe ser positiva")
        states = self.block(M_J).states
        index = {state: i for i, state in enumerate(states)}
        entries = []
        coefficient = {-1: -1.0, 0: -2.0, 1: -1.0}
        scale = self.molecule.d_au**2 / distance**3
        for i, si in enumerate(states):
            for q in (-1, 0, 1):
                M1j, M2j = si.M1 - q, si.M2 + q
                for N1j in (si.N1 - 1, si.N1 + 1):
                    if N1j < 0 or N1j > self.N_max or abs(M1j) > N1j:
                        continue
                    c1 = self._rotor_c(si.N1, si.M1, N1j, M1j, q)
                    for N2j in (si.N2 - 1, si.N2 + 1):
                        if N2j < 0 or N2j > self.N_max or abs(M2j) > N2j:
                            continue
                        c2 = self._rotor_c(si.N2, si.M2, N2j, M2j, -q)
                        sj = TwoRotorState(si.l, si.m, N1j, M1j, N2j, M2j)
                        j = index.get(sj)
                        value = scale * coefficient[q] * c1 * c2
                        if j is not None and value:
                            entries.append((i, j, value))
        return self._sparse(entries, len(states))

    def orientation_matrix(self, rotor: int, M_J=0):
        key = ("cos", rotor, M_J)
        if key not in self._cache:
            states = self.block(M_J).states
            index = {state: i for i, state in enumerate(states)}
            entries = []
            for i, st in enumerate(states):
                N, M = (st.N1, st.M1) if rotor == 1 else (st.N2, st.M2)
                for Nj in (N - 1, N + 1):
                    if Nj < 0 or Nj > self.N_max or abs(M) > Nj:
                        continue
                    sj = (
                        TwoRotorState(st.l, st.m, Nj, M, st.N2, st.M2)
                        if rotor == 1
                        else TwoRotorState(st.l, st.m, st.N1, st.M1, Nj, M)
                    )
                    j = index.get(sj)
                    if j is not None:
                        entries.append((i, j, self.hmol.cos_theta_element(N, M, Nj, M)))
            self._cache[key] = self._sparse(entries, len(states))
        return self._cache[key]

    def external_field_matrix(self, field_v_per_m: float, M_J=0):
        key = ("field", float(field_v_per_m), M_J)
        if key in self._cache:
            return self._cache[key]
        states = self.block(M_J).states
        dim = len(states)
        if field_v_per_m == 0:
            return csr_matrix((dim, dim), dtype=float)
        field_au = field_v_per_m / V_PER_M_PER_AU
        result = -self.molecule.d_au * field_au * (
            self.orientation_matrix(1, M_J) + self.orientation_matrix(2, M_J)
        )
        atom = Atom(self.n_manifold, self.l_min)
        index = {state: i for i, state in enumerate(states)}
        entries = []
        radial_cache = {}
        for i, st in enumerate(states):
            for lj in (st.l - 1, st.l + 1):
                sj = TwoRotorState(lj, st.m, st.N1, st.M1, st.N2, st.M2)
                j = index.get(sj)
                if j is None:
                    continue
                pair = (min(st.l, lj), max(st.l, lj))
                if pair not in radial_cache:
                    radial_cache[pair] = float(np.trapezoid(
                        self.radial.u(pair[0]) * self.radial.u(pair[1]) * self.radial.r,
                        self.radial.r,
                    ))
                value = atom.Vfield(st.l, lj, st.m, st.m, radial_cache[pair], field_au)
                if value:
                    entries.append((i, j, value))
        result = result + self._sparse(entries, dim)
        self._cache[key] = result.tocsr()
        return self._cache[key]

    def components(
        self, R: float, M_J=0, geometry="symmetric", separation=300.0,
        field_v_per_m=0.0,
    ):
        z1, z2 = self.positions(R, geometry, separation)
        distance = abs(z2 - z1)
        zero = csr_matrix(self.atomic_rotational_matrix(M_J).shape)
        return {
            "atomic_rotational": self.atomic_rotational_matrix(M_J),
            "charge_dipole_1": self.charge_dipole_matrix(1, z1, M_J),
            "charge_dipole_2": self.charge_dipole_matrix(2, z2, M_J),
            "dipole_dipole": (
                self.dipole_dipole_matrix(distance, M_J)
                if self.include_dipole_dipole else zero
            ),
            "external_field": self.external_field_matrix(field_v_per_m, M_J),
        }

    def hamiltonian(self, R: float, M_J=0, geometry="symmetric", separation=300.0,
                    field_v_per_m=0.0):
        parts = self.components(R, M_J, geometry, separation, field_v_per_m)
        return sum(parts.values(), start=csr_matrix(parts["atomic_rotational"].shape)).tocsr()

    def solve_near(self, R: float, M_J=0, geometry="symmetric", separation=300.0,
                   field_v_per_m=0.0, k=80, sigma_ghz=-30.0):
        H = self.hamiltonian(R, M_J, geometry, separation, field_v_per_m)
        if k >= H.shape[0] - 1:
            values, vectors = np.linalg.eigh(H.toarray())
        else:
            sigma = self.E_manifold + sigma_ghz / GHZ_PER_HARTREE
            values, vectors = eigsh(H, k=k, sigma=sigma, which="LM")
            order = np.argsort(values)
            values, vectors = values[order], vectors[:, order]
        return values, vectors

    def solve_lowest(self, R: float, M_J=0, geometry="symmetric", separation=300.0,
                     field_v_per_m=0.0, k=80):
        """Los k autovalores algebraicamente más bajos, para sembrar la BOP."""
        H = self.hamiltonian(R, M_J, geometry, separation, field_v_per_m)
        if k >= H.shape[0] - 1:
            values, vectors = np.linalg.eigh(H.toarray())
        else:
            values, vectors = eigsh(H, k=k, which="SA")
            order = np.argsort(values)
            values, vectors = values[order], vectors[:, order]
        return values, vectors

    def manifold_seed(self, R: float, M_J=0, geometry="symmetric", separation=300.0,
                      field_v_per_m=0.0):
        """Estado fundamental del Hamiltoniano proyectado sobre l>=l_min."""
        H = self.hamiltonian(R, M_J, geometry, separation, field_v_per_m)
        indices = np.flatnonzero(self.is_manifold(M_J))
        projected = H[indices][:, indices]
        if projected.shape[0] <= 2:
            values, vectors = np.linalg.eigh(projected.toarray())
            value, vector = values[0], vectors[:, 0]
        else:
            values, vectors = eigsh(projected, k=1, which="SA")
            value, vector = values[0], vectors[:, 0]
        embedded = np.zeros(H.shape[0])
        embedded[indices] = vector
        return float(value), embedded
