"""
Generación de base cuántica acoplada {(l, m_l, N, M_N)} para Rb*-KRb.
Bloqueo por M_J = m_l + M_N para construcción eficiente de matriz Hamiltoniana.

Base científica: González-Férez, Sadeghpour & Schmelcher,
                New J. Phys. 17, 013021 (2015)
"""

from typing import List, Tuple, Dict, Optional
from collections import defaultdict
import numpy as np


class QuantumBasisBlock:
    """
    Representación de un bloque M_J de la base Rb*-KRb.

    Cada bloque contiene todos los estados {(l, m_l, N, M_N)}
    que satisfacen m_l + M_N = M_J (constante para el bloque).
    """

    def __init__(self, M_J: int, states: List[Tuple[int, int, int, int]]):
        """
        Args:
            M_J (int): valor del proyector M_J = m_l + M_N
            states (List[Tuple]): lista de (l, m_l, N, M_N) en este bloque
        """
        self.M_J = M_J
        self.states = states
        self.dim = len(states)

        # Índice local: para mapear estado (l,m_l,N,M_N) → índice en matriz 0..dim-1
        self._state_to_index = {state: i for i, state in enumerate(states)}

    def state_index(self, l: int, m_l: int, N: int, M_N: int) -> int:
        """
        Retorna índice lineal del estado dentro de este bloque.

        Args:
            l, m_l, N, M_N: números cuánticos

        Returns:
            int: índice 0..dim-1

        Raises:
            KeyError: si el estado no está en este bloque
        """
        state = (l, m_l, N, M_N)
        return self._state_to_index[state]

    def __len__(self) -> int:
        """Dimensión del bloque."""
        return self.dim

    def __repr__(self) -> str:
        return f"QuantumBasisBlock(M_J={self.M_J}, dim={self.dim})"


class CoupledBasis:
    """
    Base cuántica acoplada completa para Rb*-KRb con bloqueo por M_J.

    Contiene 59 bloques (M_J = -29 a 29), cada uno independiente.
    """

    def __init__(self, N_max: int = 6, manifold_l_min: int = 3,
                 manifold_l_max: int = 23, neighbor_l=(0, 1, 2)):
        """
        Args:
            N_max (int): número máximo de cuantos rotacionales de KRb (típico: 6)
            manifold_l_min (int): mínimo orbital angular en manifold (típico: 3)
            manifold_l_max (int): máximo orbital angular del manifold, = n-1.
                Por defecto 23 (manifold n=24); para n=25 hay que pasar 24.
                Es el único sitio donde el n del manifold entra en la base: el
                resto de la física lo fija `rydberg_diagonal(n_manifold=...)`.
            neighbor_l: valores de l de los niveles vecinos individuales.
                **(0, 1, 2)** = (n+3)s, (n+2)p, (n+1)d, que es la base del
                paper (Aguilera-Fernández et al. 2015, arXiv:1507.07972).
                `(0,)` reproduce la base INCOMPLETA usada antes de leer ese
                texto, sólo manifold + (n+3)s; se conserva para poder medir
                cuánto cambia añadir los otros dos, no para producción.

        Ojo: cada l identifica unívocamente un nivel, porque los tres vecinos
        tienen l = 0, 1, 2 y el manifold empieza en l = 3. Por eso el estado
        sigue siendo (l, m_l, N, M_N) sin necesidad de llevar n en la tupla.
        """
        self.N_max = N_max
        self.manifold_l_min = manifold_l_min
        self.manifold_l_max = manifold_l_max
        self.neighbor_l = tuple(sorted(neighbor_l))
        if any(l >= manifold_l_min for l in self.neighbor_l):
            raise ValueError(
                f"los vecinos {self.neighbor_l} deben tener l < manifold_l_min="
                f"{manifold_l_min}: si no, l ya no identifica el nivel")

        # Generar todos los estados sin restricción M_J
        self.all_states = self._enumerate_all_states()

        # Particionar por M_J
        self.blocks = self._partition_by_mj()

    def _enumerate_all_states(self) -> List[Tuple[int, int, int, int]]:
        """
        Enumera TODOS los estados {(l, m_l, N, M_N)} de la base.

        Restricciones:
        - l ∈ {manifold_l_min..manifold_l_max} (manifold Rb; excluye l=0,1,2
          por defecto cuántico). Con los valores por defecto, {3..23} = n=24.
        - l ∈ neighbor_l, uno por nivel vecino individual: l=0 es (n+3)s,
          l=1 es (n+2)p, l=2 es (n+1)d
        - m_l ∈ [-l, l]
        - N ∈ [0, N_max]
        - M_N ∈ [-N, N]

        Returns:
            List[(l, m_l, N, M_N)]: todos los estados, sin filtro M_J
        """
        states = []
        l_all = (list(range(self.manifold_l_min, self.manifold_l_max + 1))
                 + list(self.neighbor_l))
        for l in l_all:
            for m_l in range(-l, l + 1):
                for N in range(0, self.N_max + 1):
                    for M_N in range(-N, N + 1):
                        states.append((l, m_l, N, M_N))
        return states

    def _partition_by_mj(self) -> Dict[int, QuantumBasisBlock]:
        """
        Particiona todos los estados por M_J = m_l + M_N.

        Returns:
            Dict[M_J, QuantumBasisBlock]: bloques indexados por M_J
        """
        blocks_dict = defaultdict(list)

        for state in self.all_states:
            l, m_l, N, M_N = state
            M_J = m_l + M_N
            blocks_dict[M_J].append(state)

        # Convertir a QuantumBasisBlock
        blocks = {}
        for M_J in sorted(blocks_dict.keys()):
            blocks[M_J] = QuantumBasisBlock(M_J, blocks_dict[M_J])

        return blocks

    def get_block(self, M_J: int) -> QuantumBasisBlock:
        """
        Retorna el bloque correspondiente a M_J.

        Args:
            M_J (int): valor del proyector

        Returns:
            QuantumBasisBlock

        Raises:
            KeyError: si M_J no existe
        """
        return self.blocks[M_J]

    def block_M_J_values(self) -> List[int]:
        """Retorna lista de valores M_J disponibles, ordenada."""
        return sorted(self.blocks.keys())

    def total_dimension(self) -> int:
        """Dimensión total (suma de todos los bloques)."""
        return sum(len(b) for b in self.blocks.values())

    def num_blocks(self) -> int:
        """Número de bloques no vacíos."""
        return len(self.blocks)

    def print_statistics(self):
        """Imprime estadísticas sobre la base y bloques."""
        print("=" * 70)
        print("ESTADÍSTICAS DE BASE CUÁNTICA Rb*-KRb")
        print("=" * 70)

        print(f"\nParámetros:")
        print(f"  N_max (rotor KRb): {self.N_max}")
        print(f"  Manifold l_min: {self.manifold_l_min}")
        print(f"  Total de estados (sin bloqueo): {len(self.all_states)}")
        print(f"  Número de bloques M_J: {self.num_blocks()}")
        print(f"  Rango M_J: [{min(self.block_M_J_values())}, {max(self.block_M_J_values())}]")

        print(f"\nDimensiones por bloque:")
        print(f"  M_J  | Dim")
        print(f"  -----|-------")
        for M_J in self.block_M_J_values():
            dim = len(self.blocks[M_J])
            print(f"  {M_J:4d} | {dim:6d}")

        print(f"\n  Suma de todos: {self.total_dimension()}")
        print(f"  Esperado (568 × 49): {568 * 49}")
        print(f"  Coincide: {'✓ SÍ' if self.total_dimension() == 568 * 49 else '✗ NO'}")

        print("\n" + "=" * 70)

    def block_M_J_0_detail(self):
        """Imprime desglose detallado del bloque M_J=0 por l."""
        if 0 not in self.blocks:
            print("Bloque M_J=0 no existe")
            return

        block = self.blocks[0]
        print("=" * 70)
        print("BLOQUE M_J=0 (DESGLOSE POR l)")
        print("=" * 70)

        by_l = defaultdict(list)
        for state in block.states:
            l, m_l, N, M_N = state
            by_l[l].append(state)

        print(f"\n  l  | #estados")
        print(f"  ---|----------")
        total_check = 0
        for l in sorted(by_l.keys()):
            count = len(by_l[l])
            total_check += count
            print(f"  {l:3d} | {count:8d}")

        print(f"  ---|----------")
        print(f" TOT | {total_check:8d}")
        print(f"\n  Esperado: 1016")
        print(f"  Coincide: {'✓ SÍ' if total_check == 1016 else '✗ NO'}")


def test_basis():
    """Test y verificación de la clase CoupledBasis."""
    print("Inicializando base cuántica...")
    basis = CoupledBasis(N_max=6, manifold_l_min=3)

    # Estadísticas
    basis.print_statistics()

    # Detalle M_J=0
    basis.block_M_J_0_detail()

    # Test de acceso a estado específico
    print("\n" + "=" * 70)
    print("TEST DE ACCESO A ESTADO")
    print("=" * 70)

    block_0 = basis.get_block(0)
    print(f"\nBloque M_J=0: {block_0}")

    # Ejemplo: estado (l=5, m_l=-3, N=3, M_N=0) está en M_J=0
    # porque m_l + M_N = -3 + 0 = -3 ≠ 0, así que NO está
    # Probemos (l=5, m_l=3, N=3, M_N=-3): m_l+M_N = 3-3 = 0 ✓
    l, m_l, N, M_N = 5, 3, 3, -3
    M_J_test = m_l + M_N
    if M_J_test == 0:
        try:
            idx = block_0.state_index(l, m_l, N, M_N)
            print(f"\nEstado ({l}, {m_l}, {N}, {M_N}): M_J={M_J_test}")
            print(f"  Encontrado en índice {idx} del bloque M_J=0 ✓")
        except KeyError:
            print(f"\nEstado ({l}, {m_l}, {N}, {M_N}): M_J={M_J_test}")
            print(f"  ERROR: no encontrado en bloque (bug en code)")
    else:
        print(f"\nEstado ({l}, {m_l}, {N}, {M_N}): M_J={M_J_test} ≠ 0")
        print("  (Test saltado: estado no está en bloque M_J=0)")

    print("\n" + "=" * 70)
    print("Test completado.")
    print("=" * 70)


if __name__ == "__main__":
    test_basis()
