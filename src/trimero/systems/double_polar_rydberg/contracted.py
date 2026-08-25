"""Base pendular local contraída para el sistema de dos rotores."""

from dataclasses import dataclass
from typing import NamedTuple

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix, diags
from scipy.sparse.linalg import eigsh

from trimero.systems.polar_molecule import RBCS, PolarMolecule
from trimero.systems.polar_rydberg import V_PER_M_PER_AU
from trimero.systems.rb_atom import Atom
from trimero.systems.rb_krb_polar.rb_defects import energy_rb

from .system import GHZ_PER_HARTREE, RbTwoRbCsSystem


class ContractedState(NamedTuple):
    l: int
    m: int
    M1: int
    a1: int
    M2: int
    a2: int


@dataclass(frozen=True)
class LocalRotorBasis:
    """Autoestados pendulares separados por M; columnas en la base N."""

    energies: dict[int, np.ndarray]
    coefficients: dict[int, np.ndarray]
    N_values: dict[int, np.ndarray]
    hamiltonians: dict[int, np.ndarray]
    occupations: dict[int, np.ndarray] | None = None

    def modes(self):
        for M, values in self.energies.items():
            for alpha in range(len(values)):
                yield M, alpha


@dataclass(frozen=True)
class ContractedBlock:
    M_J: int
    states: tuple[ContractedState, ...]

    def __len__(self):
        return len(self.states)


@dataclass(frozen=True)
class ContractedPoint:
    R: float
    geometry: str
    separation: float
    field_v_per_m: float
    M_J: int
    z1: float
    z2: float
    rotor1: LocalRotorBasis
    rotor2: LocalRotorBasis
    block: ContractedBlock


class ContractedRbTwoRbCsSystem:
    """Hamiltoniano proyectado en estados pendulares de cada rotor."""

    def __init__(
        self, n_manifold=20, primitive_N_max=8, rotor_keep=2,
        molecule: PolarMolecule = RBCS, l_min=3, delta0_ns=None,
        include_electron_field=True, include_dipole_dipole=True,
    ):
        if primitive_N_max < 0 or rotor_keep < 1:
            raise ValueError("primitive_N_max>=0 y rotor_keep>=1")
        self.primitive_N_max = int(primitive_N_max)
        self.rotor_keep = int(rotor_keep)
        self.primitive = RbTwoRbCsSystem(
            n_manifold=n_manifold, N_max=primitive_N_max,
            molecule=molecule, l_min=l_min, delta0_ns=delta0_ns,
            include_electron_field=include_electron_field,
            include_dipole_dipole=include_dipole_dipole,
        )
        self.n_manifold = n_manifold
        self.l_min = l_min
        self.l_values = self.primitive.l_values
        self.molecule = self.primitive.molecule
        self.include_dipole_dipole = include_dipole_dipole
        self.E_manifold = self.primitive.E_manifold
        self._atomic = {}

    def _local_rotor(self, local_field_au: float) -> LocalRotorBasis:
        energies, coefficients, n_values = {}, {}, {}
        for M in range(-self.primitive_N_max, self.primitive_N_max + 1):
            N = np.arange(abs(M), self.primitive_N_max + 1)
            if not len(N):
                continue
            H = np.diag(self.molecule.B_au * N * (N + 1))
            for i, Ni in enumerate(N[:-1]):
                value = self.primitive.hmol.cos_theta_element(Ni, M, Ni + 1, M)
                H[i, i + 1] = H[i + 1, i] = -self.molecule.d_au * local_field_au * value
            values, vectors = np.linalg.eigh(H)
            keep = min(self.rotor_keep, len(values))
            energies[M] = values[:keep]
            coefficients[M] = vectors[:, :keep]
            n_values[M] = N
        hamiltonians = {M: np.diag(values) for M, values in energies.items()}
        return LocalRotorBasis(energies, coefficients, n_values, hamiltonians)

    def _natural_rotor(
        self, source_block, source_vector, rotor_number: int,
        local_field_au: float, keep: int,
    ) -> LocalRotorBasis:
        """Estados naturales por M, completados con modos pendulares locales."""
        local = self._local_rotor(local_field_au)
        source_nmax = max(
            (st.N1 if rotor_number == 1 else st.N2) for st in source_block.states
        )
        energies, coefficients, n_values, hamiltonians, occupations = {}, {}, {}, {}, {}
        for M in range(-self.primitive_N_max, self.primitive_N_max + 1):
            target_N = np.arange(abs(M), self.primitive_N_max + 1)
            dim = len(target_N)
            if not dim:
                continue
            environment = {}
            for index, st in enumerate(source_block.states):
                N, state_M = (st.N1, st.M1) if rotor_number == 1 else (st.N2, st.M2)
                if state_M != M:
                    continue
                env = (
                    (st.l, st.m, st.N2, st.M2)
                    if rotor_number == 1 else (st.l, st.m, st.N1, st.M1)
                )
                environment.setdefault(env, np.zeros(dim))[N - abs(M)] = source_vector[index]
            rho = np.zeros((dim, dim))
            for amplitude in environment.values():
                rho += np.outer(amplitude, amplitude)
            occ, natural = np.linalg.eigh(rho)
            order = np.argsort(occ)[::-1]
            occ, natural = occ[order], natural[:, order]
            selected = []
            selected_occ = []
            for column, value in zip(natural.T, occ):
                if value <= 1e-13 or len(selected) >= min(keep, dim):
                    break
                selected.append(column.copy())
                selected_occ.append(float(value))
            # Completar con modos locales de baja energía, ortogonalizados.
            for candidate in local.coefficients[M].T:
                if len(selected) >= min(keep, dim):
                    break
                vector = candidate.copy()
                for previous in selected:
                    vector -= previous * (previous @ vector)
                norm = np.linalg.norm(vector)
                if norm > 1e-10:
                    selected.append(vector / norm)
                    selected_occ.append(0.0)
            # rotor_keep puede ser menor que keep; completar con vectores canónicos.
            for candidate in np.eye(dim):
                if len(selected) >= min(keep, dim):
                    break
                vector = candidate.copy()
                for previous in selected:
                    vector -= previous * (previous @ vector)
                norm = np.linalg.norm(vector)
                if norm > 1e-10:
                    selected.append(vector / norm)
                    selected_occ.append(0.0)
            U = np.column_stack(selected)
            primitive_h = np.diag(self.molecule.B_au * target_N * (target_N + 1))
            for i, Ni in enumerate(target_N[:-1]):
                value = self.primitive.hmol.cos_theta_element(Ni, M, Ni + 1, M)
                primitive_h[i, i + 1] = primitive_h[i + 1, i] = (
                    -self.molecule.d_au * local_field_au * value
                )
            hamiltonians[M] = U.T @ primitive_h @ U
            energies[M] = np.diag(hamiltonians[M]).copy()
            coefficients[M] = U
            n_values[M] = target_N
            occupations[M] = np.array(selected_occ)
        return LocalRotorBasis(
            energies, coefficients, n_values, hamiltonians, occupations
        )

    def point(self, R, M_J=0, geometry="symmetric", separation=300.0,
              field_v_per_m=0.0) -> ContractedPoint:
        z1, z2 = self.primitive.positions(R, geometry, separation)
        external = field_v_per_m / V_PER_M_PER_AU
        rotor1 = self._local_rotor(np.sign(z1) / z1**2 + external)
        rotor2 = self._local_rotor(np.sign(z2) / z2**2 + external)
        modes1, modes2 = tuple(rotor1.modes()), tuple(rotor2.modes())
        by_M2 = {}
        for M2, a2 in modes2:
            by_M2.setdefault(M2, []).append(a2)
        states = []
        for l in self.l_values:
            for m in range(-l, l + 1):
                for M1, a1 in modes1:
                    M2 = M_J - m - M1
                    for a2 in by_M2.get(M2, ()):
                        states.append(ContractedState(l, m, M1, a1, M2, a2))
        block = ContractedBlock(int(M_J), tuple(states))
        return ContractedPoint(
            float(R), geometry, float(separation), float(field_v_per_m), int(M_J),
            z1, z2, rotor1, rotor2, block,
        )

    def natural_point(
        self, source_block, source_vector, R, M_J=0, geometry="symmetric",
        separation=300.0, field_v_per_m=0.0, natural_keep=None,
    ) -> ContractedPoint:
        """Punto contraído usando densidades reducidas de un estado fuente."""
        keep = self.rotor_keep if natural_keep is None else int(natural_keep)
        z1, z2 = self.primitive.positions(R, geometry, separation)
        external = field_v_per_m / V_PER_M_PER_AU
        rotor1 = self._natural_rotor(
            source_block, source_vector, 1, np.sign(z1) / z1**2 + external, keep
        )
        rotor2 = self._natural_rotor(
            source_block, source_vector, 2, np.sign(z2) / z2**2 + external, keep
        )
        modes1, modes2 = tuple(rotor1.modes()), tuple(rotor2.modes())
        by_M2 = {}
        for M2, a2 in modes2:
            by_M2.setdefault(M2, []).append(a2)
        states = []
        for l in self.l_values:
            for m in range(-l, l + 1):
                for M1, a1 in modes1:
                    M2 = M_J - m - M1
                    for a2 in by_M2.get(M2, ()):
                        states.append(ContractedState(l, m, M1, a1, M2, a2))
        return ContractedPoint(
            float(R), geometry, float(separation), float(field_v_per_m), int(M_J),
            z1, z2, rotor1, rotor2, ContractedBlock(int(M_J), tuple(states)),
        )

    def _transfer_rotor(
        self, reference: LocalRotorBasis, local_field_au: float,
    ) -> LocalRotorBasis:
        """Reutiliza el subespacio de un ancla y actualiza su Hamiltoniano local."""
        energies, coefficients, n_values, hamiltonians = {}, {}, {}, {}
        for M, U in reference.coefficients.items():
            N = reference.N_values[M]
            primitive_h = np.diag(self.molecule.B_au * N * (N + 1))
            for i, Ni in enumerate(N[:-1]):
                value = self.primitive.hmol.cos_theta_element(
                    int(Ni), M, int(Ni + 1), M
                )
                primitive_h[i, i + 1] = primitive_h[i + 1, i] = (
                    -self.molecule.d_au * local_field_au * value
                )
            projected = U.T @ primitive_h @ U
            hamiltonians[M] = projected
            energies[M] = np.diag(projected).copy()
            coefficients[M] = U.copy()
            n_values[M] = N.copy()
        occupations = None
        if reference.occupations is not None:
            occupations = {
                M: values.copy() for M, values in reference.occupations.items()
            }
        return LocalRotorBasis(
            energies, coefficients, n_values, hamiltonians, occupations
        )

    def transfer_point(
        self, reference: ContractedPoint, R, M_J=None, geometry=None,
        separation=None, field_v_per_m=None,
    ) -> ContractedPoint:
        """Transfiere los subespacios rotacionales de un ancla a otro punto."""
        M_J = reference.M_J if M_J is None else int(M_J)
        geometry = reference.geometry if geometry is None else geometry
        separation = reference.separation if separation is None else float(separation)
        field_v_per_m = (
            reference.field_v_per_m
            if field_v_per_m is None else float(field_v_per_m)
        )
        if M_J != reference.M_J:
            raise ValueError("la transferencia entre bloques M_J no está definida")
        z1, z2 = self.primitive.positions(R, geometry, separation)
        external = field_v_per_m / V_PER_M_PER_AU
        rotor1 = self._transfer_rotor(
            reference.rotor1, np.sign(z1) / z1**2 + external
        )
        rotor2 = self._transfer_rotor(
            reference.rotor2, np.sign(z2) / z2**2 + external
        )
        return ContractedPoint(
            float(R), geometry, separation, field_v_per_m, M_J, z1, z2,
            rotor1, rotor2, reference.block,
        )

    @staticmethod
    def _sparse(entries, dim):
        if not entries:
            return csr_matrix((dim, dim), dtype=float)
        rows, cols, data = zip(*entries)
        return coo_matrix((data, (rows, cols)), shape=(dim, dim)).tocsr()

    def _rotor_element(self, rotor: LocalRotorBasis, Mi, ai, Mj, aj, q):
        if Mi != Mj + q or Mi not in rotor.energies or Mj not in rotor.energies:
            return 0.0
        Ni, Nj = rotor.N_values[Mi], rotor.N_values[Mj]
        ci = rotor.coefficients[Mi][:, ai]
        cj = rotor.coefficients[Mj][:, aj]
        total = 0.0
        for ii, n_i in enumerate(Ni):
            for jj, n_j in enumerate(Nj):
                total += ci[ii] * cj[jj] * self.primitive._rotor_c(
                    int(n_i), Mi, int(n_j), Mj, q
                )
        return float(total)

    def _atomic_energy(self, l):
        if l not in self._atomic:
            n = self.primitive.levels.get(l, self.n_manifold)
            self._atomic[l] = energy_rb(n, l, self.primitive.delta0_ns)
        return self._atomic[l]

    def hamiltonian(self, point: ContractedPoint):
        states = point.block.states
        index = {state: i for i, state in enumerate(states)}
        diagonal = [self._atomic_energy(st.l) for st in states]
        entries = []
        field_au = point.field_v_per_m / V_PER_M_PER_AU
        atom = Atom(self.n_manifold, self.l_min)
        radial_cache = {}

        # Hamiltonianos locales de ambos rotores; no tienen por qué ser
        # diagonales en una base natural.
        for i, st in enumerate(states):
            for beta in range(len(point.rotor1.energies[st.M1])):
                sj = ContractedState(st.l, st.m, st.M1, beta, st.M2, st.a2)
                j = index.get(sj)
                value = point.rotor1.hamiltonians[st.M1][st.a1, beta]
                if j is not None and value:
                    entries.append((i, j, value))
            for beta in range(len(point.rotor2.energies[st.M2])):
                sj = ContractedState(st.l, st.m, st.M1, st.a1, st.M2, beta)
                j = index.get(sj)
                value = point.rotor2.hamiltonians[st.M2][st.a2, beta]
                if j is not None and value:
                    entries.append((i, j, value))

        # Stark electrónico, identidad sobre ambos modos pendulares.
        if field_au:
            for i, st in enumerate(states):
                for lj in (st.l - 1, st.l + 1):
                    sj = ContractedState(lj, st.m, st.M1, st.a1, st.M2, st.a2)
                    j = index.get(sj)
                    if j is None:
                        continue
                    pair = (min(st.l, lj), max(st.l, lj))
                    if pair not in radial_cache:
                        radial_cache[pair] = float(np.trapezoid(
                            self.primitive.radial.u(pair[0])
                            * self.primitive.radial.u(pair[1])
                            * self.primitive.radial.r,
                            self.primitive.radial.r,
                        ))
                    value = atom.Vfield(st.l, lj, st.m, st.m, radial_cache[pair], field_au)
                    if value:
                        entries.append((i, j, value))

        # Campo del electrón Rydberg sobre cada rotor.
        if self.primitive.electron_field is not None:
            for rotor_number, (rotor, z) in enumerate(
                ((point.rotor1, point.z1), (point.rotor2, point.z2)), start=1
            ):
                radius = abs(z)
                negative = z < 0
                modes = tuple(rotor.modes())
                for i, st in enumerate(states):
                    Mi, ai = (st.M1, st.a1) if rotor_number == 1 else (st.M2, st.a2)
                    for q in (-1, 0, 1):
                        mj, Mj = st.m + q, Mi - q
                        for aj in (a for M, a in modes if M == Mj):
                            rot = self._rotor_element(rotor, Mi, ai, Mj, aj, q)
                            if not rot:
                                continue
                            for lj in self.l_values:
                                if abs(mj) > lj:
                                    continue
                                sj = (
                                    ContractedState(lj, mj, Mj, aj, st.M2, st.a2)
                                    if rotor_number == 1
                                    else ContractedState(lj, mj, st.M1, st.a1, Mj, aj)
                                )
                                j = index.get(sj)
                                if j is None:
                                    continue
                                electronic = self.primitive.electron_field.field_element(
                                    st.l, st.m, lj, mj, -q, radius
                                )
                                parity = (-1.0) ** (st.l + lj + 1) if negative else 1.0
                                value = -self.molecule.d_au * (-1.0 if q else 1.0) \
                                    * rot * electronic * parity
                                if value:
                                    entries.append((i, j, value))

        # Interacción dipolo--dipolo en el eje global Z.
        if self.include_dipole_dipole:
            distance = abs(point.z2 - point.z1)
            factor = self.molecule.d_au**2 / distance**3
            qfactor = {-1: -1.0, 0: -2.0, 1: -1.0}
            modes1, modes2 = tuple(point.rotor1.modes()), tuple(point.rotor2.modes())
            for i, st in enumerate(states):
                for q in (-1, 0, 1):
                    M1j, M2j = st.M1 - q, st.M2 + q
                    for a1j in (a for M, a in modes1 if M == M1j):
                        c1 = self._rotor_element(
                            point.rotor1, st.M1, st.a1, M1j, a1j, q
                        )
                        for a2j in (a for M, a in modes2 if M == M2j):
                            c2 = self._rotor_element(
                                point.rotor2, st.M2, st.a2, M2j, a2j, -q
                            )
                            sj = ContractedState(st.l, st.m, M1j, a1j, M2j, a2j)
                            j = index.get(sj)
                            value = factor * qfactor[q] * c1 * c2
                            if j is not None and value:
                                entries.append((i, j, value))
        return diags(diagonal, format="csr") + self._sparse(entries, len(states))

    def orientation_matrix(self, point: ContractedPoint, rotor_number: int):
        rotor = point.rotor1 if rotor_number == 1 else point.rotor2
        states, entries = point.block.states, []
        index = {state: i for i, state in enumerate(states)}
        for i, st in enumerate(states):
            M, alpha = (st.M1, st.a1) if rotor_number == 1 else (st.M2, st.a2)
            for beta in range(len(rotor.energies[M])):
                sj = (
                    ContractedState(st.l, st.m, M, beta, st.M2, st.a2)
                    if rotor_number == 1
                    else ContractedState(st.l, st.m, st.M1, st.a1, M, beta)
                )
                j = index.get(sj)
                value = self._rotor_element(rotor, M, alpha, M, beta, 0)
                if j is not None and value:
                    entries.append((i, j, value))
        return self._sparse(entries, len(states))

    def is_manifold(self, point: ContractedPoint):
        return np.fromiter((st.l >= self.l_min for st in point.block.states), dtype=bool)

    def solve_seed(self, point: ContractedPoint, k=40):
        H = self.hamiltonian(point)
        mask = self.is_manifold(point)
        indices = np.flatnonzero(mask)
        projected = H[indices][:, indices]
        value_p, vector_p = eigsh(projected, k=1, which="SA")
        seed = np.zeros(H.shape[0])
        seed[indices] = vector_p[:, 0]
        values, vectors = eigsh(H, k=min(k, H.shape[0] - 2), sigma=value_p[0], which="LM")
        order = np.argsort(values)
        values, vectors = values[order], vectors[:, order]
        selected = int(np.argmax((seed @ vectors) ** 2))
        return values, vectors, selected

    def characteristic(self, point: ContractedPoint, k=40):
        values, vectors, selected = self.solve_seed(point, k)
        vector = vectors[:, selected]
        return {
            "E": (values[selected] - self.E_manifold) * GHZ_PER_HARTREE,
            "W": float(np.sum(vector[self.is_manifold(point)] ** 2)),
            "COS1": float(vector @ (self.orientation_matrix(point, 1) @ vector)),
            "COS2": float(vector @ (self.orientation_matrix(point, 2) @ vector)),
            "dimension": len(point.block),
        }
