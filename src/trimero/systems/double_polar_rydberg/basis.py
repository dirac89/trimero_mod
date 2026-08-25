"""Base producto electrónica × rotor 1 × rotor 2 bloqueada por M_J."""

from dataclasses import dataclass
from typing import Iterator, NamedTuple


class TwoRotorState(NamedTuple):
    l: int
    m: int
    N1: int
    M1: int
    N2: int
    M2: int


@dataclass(frozen=True)
class TwoRotorBasisBlock:
    M_J: int
    states: tuple[TwoRotorState, ...]

    def __len__(self) -> int:
        return len(self.states)

    @property
    def dim(self) -> int:
        return len(self.states)


class TwoRotorBasis:
    """Enumera sólo el bloque solicitado para evitar la base producto completa."""

    def __init__(self, N_max: int, l_values):
        if N_max < 0:
            raise ValueError("N_max debe ser no negativo")
        self.N_max = int(N_max)
        self.l_values = tuple(sorted(set(int(l) for l in l_values)))
        if not self.l_values or self.l_values[0] < 0:
            raise ValueError("l_values debe contener enteros no negativos")
        self._blocks = {}

    def rotor_states(self) -> Iterator[tuple[int, int]]:
        for N in range(self.N_max + 1):
            for M in range(-N, N + 1):
                yield N, M

    def get_block(self, M_J: int) -> TwoRotorBasisBlock:
        M_J = int(M_J)
        if M_J not in self._blocks:
            rot = tuple(self.rotor_states())
            states = []
            for l in self.l_values:
                for m in range(-l, l + 1):
                    for N1, M1 in rot:
                        M2 = M_J - m - M1
                        for N2 in range(abs(M2), self.N_max + 1):
                            states.append(TwoRotorState(l, m, N1, M1, N2, M2))
            self._blocks[M_J] = TwoRotorBasisBlock(M_J, tuple(states))
        return self._blocks[M_J]
