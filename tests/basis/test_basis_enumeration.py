#!/usr/bin/env python3
"""
Enumeración exhaustiva de estados para base Rb*-KRb con bloqueo por M_J.
Verifica contra números proporcionados por usuario:
- Bloque M_J=0: 1016 estados
- Bloques totales: 59
- Suma: 27.832 = 568 × 49
"""

import numpy as np
from typing import List, Tuple
from collections import defaultdict


def enumerate_basis_states(N_max: int = 6) -> List[Tuple[int, int, int, int]]:
    """
    Enumera TODOS los estados {(l, m_l, N, M_N)} de la base Rb*-KRb sin restricción M_J.

    Restricciones:
    - l ∈ {3, 4, ..., 23} (manifold cuasi-degenerado n=24)
    - l = 0 (estado 27s)
    - m_l ∈ [-l, l]
    - N ∈ [0, N_max]
    - M_N ∈ [-N, N]

    Returns:
        List[(l, m_l, N, M_N)]: todos los estados, sin filtro M_J
    """
    states = []

    # Manifold principal: l ∈ {3..23}
    for l in range(3, 24):
        for m_l in range(-l, l + 1):
            for N in range(0, N_max + 1):
                for M_N in range(-N, N + 1):
                    states.append((l, m_l, N, M_N))

    # Estado 27s: l=0
    l = 0
    for m_l in range(-l, l + 1):  # m_l = 0 only
        for N in range(0, N_max + 1):
            for M_N in range(-N, N + 1):
                states.append((l, m_l, N, M_N))

    return states


def block_by_mj(states: List[Tuple[int, int, int, int]]) -> dict:
    """
    Particiona estados por M_J = m_l + M_N.

    Args:
        states: lista de (l, m_l, N, M_N)

    Returns:
        dict: {M_J: [list of states with that M_J]}
    """
    blocks = defaultdict(list)
    for state in states:
        l, m_l, N, M_N = state
        M_J = m_l + M_N
        blocks[M_J].append(state)
    return dict(blocks)


def print_statistics(blocks: dict):
    """Imprime estadísticas sobre bloques."""
    print("=" * 70)
    print("ESTADÍSTICAS DE BLOQUES M_J")
    print("=" * 70)

    M_J_values = sorted(blocks.keys())
    print(f"\nRango de M_J: {min(M_J_values)} a {max(M_J_values)}")
    print(f"Número de bloques no vacíos: {len(blocks)}")
    print()

    # Tabla de tamaños
    print(f"{'M_J':>4} | {'Dim':>6}")
    print("-" * 12)
    for M_J in M_J_values:
        dim = len(blocks[M_J])
        print(f"{M_J:4d} | {dim:6d}")

    print()
    print(f"Tamaño bloque M_J=0: {len(blocks[0])} (esperado: 1016)")
    print(f"Tamaño bloque M_J=1: {len(blocks.get(1, []))} (verificación)")
    print()

    total = sum(len(b) for b in blocks.values())
    expected_total = 568 * 49  # manifold × rotor KRb
    print(f"Suma de todos los bloques: {total}")
    print(f"Esperado (568 × 49): {expected_total}")
    print(f"Coincide: {'✓ SÍ' if total == expected_total else '✗ NO'}")
    print()


def verify_block_mj_0(blocks: dict):
    """Verifica en detalle el bloque M_J=0."""
    print("=" * 70)
    print("VERIFICACIÓN DETALLADA: BLOQUE M_J=0")
    print("=" * 70)

    block = blocks[0]
    print(f"\nTotal de estados en M_J=0: {len(block)}")

    # Desglose por l
    by_l = defaultdict(list)
    for state in block:
        l, m_l, N, M_N = state
        by_l[l].append(state)

    print(f"\nDesglose por l:")
    print(f"{'l':>3} | {'#estados':>8}")
    print("-" * 13)
    total_check = 0
    for l in sorted(by_l.keys()):
        count = len(by_l[l])
        total_check += count
        print(f"{l:3d} | {count:8d}")

    print(f"{'TOTAL':>3} | {total_check:8d}")
    print(f"\nConsistencia: {total_check} == {len(block)} : {'✓' if total_check == len(block) else '✗'}")


def _report_enumeration():
    """Enumeración + informe por pantalla. Helper, no es un test.

    Se conserva porque su salida documenta la estructura de bloques; las
    aserciones viven en los tests de abajo.
    """
    all_states = enumerate_basis_states(N_max=6)
    blocks = block_by_mj(all_states)
    print_statistics(blocks)
    verify_block_mj_0(blocks)
    return all_states, blocks


# ----------------------------------------------------------------------
# Los tres números que el docstring de este módulo declaraba como esperados
# se comprobaban imprimiendo ✓/✗, sin fallar nunca. Ahora son aserciones.
# ----------------------------------------------------------------------
EXPECTED_TOTAL = 568 * 49        # 27832 estados sin restricción de M_J
EXPECTED_N_BLOCKS = 59           # bloques distintos de M_J
EXPECTED_MJ0_SIZE = 1016         # estados en el bloque M_J = 0


def test_total_number_of_states():
    assert len(enumerate_basis_states(N_max=6)) == EXPECTED_TOTAL


def test_number_of_mj_blocks():
    blocks = block_by_mj(enumerate_basis_states(N_max=6))
    assert len(blocks) == EXPECTED_N_BLOCKS


def test_mj_zero_block_size():
    blocks = block_by_mj(enumerate_basis_states(N_max=6))
    assert len(blocks[0]) == EXPECTED_MJ0_SIZE


def test_blocks_partition_the_full_basis():
    """Los bloques deben ser una partición: ni pierden ni duplican estados."""
    states = enumerate_basis_states(N_max=6)
    blocks = block_by_mj(states)
    assert sum(len(b) for b in blocks.values()) == len(states)
    assert len({s for b in blocks.values() for s in b}) == len(states)


def test_matches_production_coupled_basis():
    """La enumeración de referencia y CoupledBasis deben coincidir.

    Es lo que este fichero verificaba de forma implícita: que la clase usada
    en producción enumera exactamente la misma base que el script.
    """
    from quantum_basis import CoupledBasis

    basis = CoupledBasis(N_max=6, manifold_l_min=3)
    reference = block_by_mj(enumerate_basis_states(N_max=6))

    assert basis.total_dimension() == EXPECTED_TOTAL
    assert basis.num_blocks() == len(reference)
    for m_j in basis.block_M_J_values():
        assert set(basis.get_block(m_j).states) == set(reference[m_j]), (
            f"el bloque M_J={m_j} difiere de la enumeración de referencia"
        )
