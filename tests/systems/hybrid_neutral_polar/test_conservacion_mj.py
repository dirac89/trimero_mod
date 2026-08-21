"""
PUNTO 3 de la ronda híbrida: ¿qué número cuántico se conserva de verdad?

En el sistema polar puro la simetría axial alrededor de R⃗ conserva
M_J = m_l + M_N. En el sistema neutro puro (sin campo) se conserva m_l.
Aquí conviven ambos acoplamientos y el perturbador de Fermi está SOBRE el
eje Z (θ=π), así que:

* H_A: diagonal, trivialmente diagonal en M_J.
* H_mol: conserva M_J por construcción del bloque (verificado aquí con los
  elementos "unchecked" entre bloques distintos).
* V_Fermi(θ=π): conserva m_l (reglas de selección Δm=0, |m|≤1) y es
  identidad sobre el rotor → conserva M_J. Verificación NO trivial porque
  el factor electrónico es distinto del de θ=0 (fase (−1)^{l₁+l₂}) y las
  reglas de selección podrían haberse roto al aplicar la fase.

La conclusión esperada —M_J sigue siendo el buen número— no se asume: cada
pieza se verifica numéricamente abajo, incluyendo que `matrix_element` /
`fermi_pi_element` LANZAN si alguien pide un elemento entre bloques
distintos (guardia anti-error de indexación).

Ejecutar:  poetry run pytest tests/systems/hybrid_neutral_polar/test_conservacion_mj.py -s
"""

import numpy as np
import pytest

from trimero.systems.hybrid_neutral_polar import HybridNeutralPolar

N_MANIFOLD = 35
R1 = 900.0     # Rb neutro (θ=π)
R2 = 800.0     # RbCs (θ=0)


@pytest.fixture(scope="module")
def hyb():
    return HybridNeutralPolar(n_manifold=N_MANIFOLD, N_max=2)


# ================================================================ TEST M1
def test_m1_hmol_conserva_mj(hyb):
    """
    Elementos de H_mol entre estados de bloques M_J DISTINTOS: deben ser
    cero exacto. Se usa `_matrix_element_unchecked`, que NO filtra por M_J,
    igual que se verificó el sistema polar en su ronda.
    """
    h = hyb.hmol
    MJ_values = [0, 1, -1, 2]
    worst = 0.0
    n_pairs = 0
    print("\n  H_mol entre bloques M_J distintos (elemento unchecked):")
    for a in range(len(MJ_values)):
        for b in range(a + 1, len(MJ_values)):
            sa = hyb.block(MJ_values[a]).states
            sb = hyb.block(MJ_values[b]).states
            mx = max(abs(h._matrix_element_unchecked(x, y, R2))
                     for x in sa for y in sb)
            worst = max(worst, mx)
            n_pairs += len(sa) * len(sb)
            print(f"    M_J={MJ_values[a]:+d} vs {MJ_values[b]:+d}: "
                  f"{len(sa)}x{len(sb)} = {len(sa)*len(sb)} pares, "
                  f"max|elem| = {mx:.1e}")
    print(f"  total {n_pairs} pares, peor max|elem| = {worst:.1e}")
    assert worst == 0.0


# ================================================================ TEST M2
def test_m2_vfermi_conserva_mj(hyb):
    """
    Elementos de V_Fermi(θ=π) entre bloques M_J distintos.

    El factor electrónico (`fermi_pi_element_unchecked`) ignora el rotor a
    propósito: puede ser NO nulo para pares con el mismo m_l pero distinto
    (N, M_N). El elemento FÍSICO añade la ortogonalidad del rotor; lo que
    debe ser cero EXACTO entre M_J distintos es:
      * el elemento completo cuando m_l difiere (selección Δm=0), y
      * el elemento completo siempre, tras imponer δ_{N N'}δ_{M_N M_N'},
    porque dos estados de bloques M_J distintos con el mismo m_l necesariamente
    difieren en M_N (y/o N), y el rotor es ortogonal.
    """
    MJ_values = [0, 1, -1, 2]
    worst_full = 0.0
    n_nonzero_electronic = 0
    n_pairs = 0
    print("\n  V_Fermi(θ=π) entre bloques M_J distintos:")
    for a in range(len(MJ_values)):
        for b in range(a + 1, len(MJ_values)):
            sa = hyb.block(MJ_values[a]).states
            sb = hyb.block(MJ_values[b]).states
            for x in sa:
                for y in sb:
                    # Elemento completo construido a mano: no se puede usar
                    # `fermi_pi_element` porque su guardia LANZA entre bloques
                    # distintos (así debe ser). La física a verificar es que
                    # factor electrónico × ortogonalidad rotacional == 0.
                    rotor_delta = 1.0 if (x[2] == y[2] and x[3] == y[3]) else 0.0
                    full = rotor_delta * hyb.fermi_pi_element_unchecked(x, y, R1)
                    worst_full = max(worst_full, abs(full))
                    if abs(hyb.fermi_pi_element_unchecked(x, y, R1)) > 0.0:
                        n_nonzero_electronic += 1
                    n_pairs += 1
    print(f"    {n_pairs} pares entre bloques distintos")
    print(f"    factores electrónicos no nulos encontrados: "
          f"{n_nonzero_electronic} (pares con m_l común)")
    print(f"    peor |elemento completo| entre bloques distintos = "
          f"{worst_full:.1e}")
    assert n_nonzero_electronic > 0, (
        "el barrido es vacuo: ningún par cruzado comparte m_l")
    assert worst_full == 0.0


# ================================================================ TEST M3
def test_m3_guardias_lanzan_valueerror(hyb):
    """Pedir un elemento entre M_J distintos debe LANZAR, no devolver 0."""
    s0 = hyb.block(0).states[0]
    s1 = hyb.block(1).states[0]
    with pytest.raises(ValueError):
        hyb.hmol.matrix_element(s0, s1, R2)
    with pytest.raises(ValueError):
        hyb.fermi_pi_element(s0, s1, R1)
    print("\n  Guardias OK: matrix_element y fermi_pi_element lanzan "
          "ValueError entre M_J distintos")


# ================================================================ TEST M4
def test_m4_soporte_no_vacio_y_estructura_del_bloque(hyb):
    """
    Control positivo: dentro de un mismo bloque SÍ hay elementos no nulos de
    ambos términos, y el bloque M_J=0 contiene exactamente los estados con
    m_l + M_N = 0.
    """
    blk = hyb.block(0)
    assert all(st[1] + st[3] == 0 for st in blk.states)

    Hmol = hyb.hmol.build(blk, R2)
    Vpi = hyb.fermi_pi_matrix(blk, R1)
    off_hmol = np.count_nonzero(np.triu(Hmol, 1))
    off_vpi = np.count_nonzero(np.triu(Vpi, 1))
    scale_v = np.max(np.abs(Vpi))
    ghz = scale_v * 6.579683920502e6
    print(f"\n  Bloque M_J=0, dim={len(blk)}:")
    print(f"    H_mol fuera de diagonal: {off_hmol} elementos no nulos")
    print(f"    V_Fermi(π) fuera de diagonal: {off_vpi}, "
          f"max|V| = {scale_v:.3e} E_h ({ghz:.1f} GHz)")
    assert off_hmol > 0
    assert off_vpi > 0
    assert np.max(np.abs(Hmol)) > 0.0

    # Y el soporte de V dentro del bloque: sólo estados con |m_l| <= 1.
    soporte = [st for st in blk.states if abs(st[1]) <= 1]
    idx_soporte = [blk.states.index(st) for st in soporte]
    fuera = Vpi.copy()
    fuera[np.ix_(idx_soporte, idx_soporte)] = 0.0
    print(f"    soporte |m_l|<=1: {len(soporte)} estados; "
          f"V fuera del soporte: {np.count_nonzero(fuera)} elementos")
    assert np.count_nonzero(fuera) == 0
