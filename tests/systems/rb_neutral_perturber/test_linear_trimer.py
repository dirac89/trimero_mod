"""
Tests analíticos del trímero lineal simétrico en campo DC.

Aguilera-Fernández, Schmelcher & González-Férez, J. Phys. B 49, 124002 (2016).

Orden de la batería (cada uno cierra una premisa del anterior):

  T1  el buen número cuántico es m_l, no M_J
  T2  V_Fermi y F·z conservan m_l por separado -> H bloque-diagonal en m_l
  T3  la geometría de DOS perturbadores es el factor (−1)^{l₁+l₂}, verificado
      por fuerza bruta contra ψ y ∇ψ numéricos, no contra la fórmula cerrada
  T4  con F=0 el espectro depende sólo de R: bloques ±m_l degenerados, y para
      |m_l| ≥ 2 no depende ni de R
  T5  con F=0 los autovectores del bloque Σ son gerade o ungerade puros
  T6  el campo rompe esa paridad
  T7  límite A_s = A_p = 0: los autovalores son las energías atómicas
  T8  límite hidrogenoide puro: escalera Stark lineal (3/2)·n·q·F exacta
  T9  con n*=35 el remapeo k(R) es la identidad y la tabla se lee nativa
"""

import numpy as np
import pytest
from scipy.special import sph_harm_y

from trimero.systems.rb_atom import Atom
from trimero.systems.rb_neutral_perturber.fermi_krb import (
    ScatteringLengths,
    hydrogenic_R,
)
from trimero.systems.rb_neutral_perturber.linear_trimer import (
    HARTREE_TO_GHZ,
    SymmetricLinearTrimer,
    field_au,
)

N = 35
R_TEST = 1500.0


@pytest.fixture(scope="module")
def trimer():
    return SymmetricLinearTrimer(n_manifold=N)


# ------------------------------------------------------------------ T1
def test_t1_el_buen_numero_cuantico_es_m_l(trimer):
    """
    Sin rotor (N_max=0) M_J colapsa a m_l: cada bloque tiene un único m_l.

    Es la premisa de todo lo demás. En Rb*-KRb el bloque M_J mezclaba varios
    m_l porque M_N podía compensarlos; aquí M_N ≡ 0.
    """
    for M_J in trimer.basis.block_M_J_values():
        block = trimer.basis.get_block(M_J)
        assert {(N_, MN) for (_l, _m, N_, MN) in block.states} == {(0, 0)}
        assert {m for (_l, m, _N, _MN) in block.states} == {M_J}

    # Dimensiones: el bloque m_l=0 tiene un estado por l (0..n-1).
    assert trimer.l_values(0) == list(range(N))
    assert trimer.l_values(1) == list(range(1, N))
    assert trimer.l_values(5) == list(range(5, N))
    assert trimer.basis.total_dimension() == N * N


# ------------------------------------------------------------------ T2
def test_t2_los_tres_terminos_conservan_m_l(trimer):
    """V_Fermi y F·z son ambos Δm_l = 0, luego H es bloque-diagonal en m_l."""
    atom = Atom(N, 0)
    for l1 in (0, 3, 7, 20):
        for l2 in (1, 4, 8, 21):
            for m1 in (-1, 0, 1, 2):
                for m2 in (-1, 0, 1, 2):
                    if m1 == m2:
                        continue
                    if m1 <= l1 and m2 <= l2:
                        assert trimer.pseudo.electron_element(l1, m1, l2, m2, R_TEST) == 0.0
                        assert atom.Angular_dc_field(l1, l2, m1, m2) == 0.0

    # Y el pseudopotencial se anula idénticamente para |m_l| >= 2:
    # sobre el eje sólo sobreviven Y_{l0} y el gradiente de Y_{l,±1}.
    for m in (2, 3, 5):
        for l1 in (5, 9):
            for l2 in (5, 11):
                assert trimer.pseudo.electron_element(l1, m, l2, m, R_TEST) == 0.0


# ------------------------------------------------------------------ T3
def _psi(n_eff, l, m, xyz):
    """ψ_{nlm}(r⃗) hidrogenoide en cartesianas."""
    x, y, z = xyz
    r = np.sqrt(x * x + y * y + z * z)
    theta = np.arccos(np.clip(z / r, -1.0, 1.0))
    phi = np.arctan2(y, x)
    return hydrogenic_R(n_eff, l, r) * sph_harm_y(l, m, theta, phi)


def _grad_psi(n_eff, l, m, xyz, h=1e-3):
    """∇ψ por diferencias centradas en cartesianas (el eje es singular en θ)."""
    g = []
    for k in range(3):
        p = list(xyz)
        p[k] += h
        plus = _psi(n_eff, l, m, p)
        p[k] -= 2 * h
        minus = _psi(n_eff, l, m, p)
        g.append((plus - minus) / (2 * h))
    return np.array(g)


@pytest.mark.parametrize("m", [0, 1])
@pytest.mark.parametrize("l1,l2", [(4, 6), (4, 7), (5, 5), (3, 8)])
def test_t3_dos_perturbadores_dan_el_factor_de_paridad(trimer, m, l1, l2):
    """
    ⟨l₁m|V(θ=π)|l₂m⟩ = (−1)^{l₁+l₂} ⟨l₁m|V(θ=0)|l₂m⟩.

    Se comprueba evaluando el pseudopotencial DIRECTAMENTE sobre ψ y ∇ψ
    numéricos en R⃗ = (0,0,±R): no se reutiliza ninguna forma cerrada, así que
    valida a la vez la geometría y las fórmulas de `fermi_krb`.
    """
    A_s, A_p = trimer.scattering.scattering(R_TEST, float(N))

    def V(sign):
        p = (0.0, 0.0, sign * R_TEST)
        v = 2 * np.pi * A_s * np.conj(_psi(N, l1, m, p)) * _psi(N, l2, m, p)
        g1 = np.conj(_grad_psi(N, l1, m, p))
        g2 = _grad_psi(N, l2, m, p)
        v += 6 * np.pi * A_p * np.dot(g1, g2)
        return v

    v_plus, v_minus = V(+1.0), V(-1.0)
    assert abs(v_plus.imag) < 1e-12 * max(abs(v_plus), 1e-30)
    expected = (-1.0) ** (l1 + l2) * v_plus.real
    assert v_minus.real == pytest.approx(expected, rel=1e-6, abs=1e-18)

    # ...y la forma cerrada de `fermi_krb` reproduce el término θ=0.
    closed = trimer.pseudo.electron_element(l1, m, l2, m, R_TEST)
    assert closed == pytest.approx(v_plus.real, rel=2e-4, abs=1e-18)

    # De ahí sale el factor de paridad que usa el módulo.
    total = trimer.parity_factor(l1, l2) * closed
    assert total == pytest.approx(v_plus.real + v_minus.real, rel=2e-4, abs=1e-18)


def test_t3b_el_pseudopotencial_total_se_anula_para_l_impar(trimer):
    """l₁+l₂ impar -> V(θ=0) y V(θ=π) se cancelan EXACTAMENTE."""
    V = trimer.pseudopotential_matrix(R_TEST, 0)
    ls = trimer.l_values(0)
    for i, l1 in enumerate(ls):
        for j, l2 in enumerate(ls):
            if (l1 + l2) % 2 == 1:
                assert V[i, j] == 0.0
    assert np.any(V != 0.0)


# ------------------------------------------------------------------ T4
def test_t4_sin_campo_el_espectro_solo_depende_de_R(trimer):
    """
    Los bloques +m_l y −m_l son degenerados, y |m_l| ≥ 2 no siente el
    perturbador: su espectro es constante en R, igual a las energías atómicas.
    """
    for R in (1200.0, 1800.0, 2400.0):
        np.testing.assert_allclose(
            trimer.energies(R, 0.0, +1), trimer.energies(R, 0.0, -1), rtol=1e-12
        )

    E_atom = np.sort((trimer.diagonal(3) - trimer.E_manifold) * HARTREE_TO_GHZ)
    for R in (1200.0, 2400.0):
        np.testing.assert_allclose(trimer.energies(R, 0.0, 3), E_atom, atol=1e-10)

    # ...mientras que Σ y Π sí dependen de R (es la ligadura ULRM).
    assert not np.allclose(trimer.energies(1200.0, 0.0, 0),
                           trimer.energies(1800.0, 0.0, 0), atol=1e-6)


# ------------------------------------------------------------------ T5
def _parity_purity(vec, ls):
    """Peso en l par menos peso en l impar, en valor absoluto: 1 = puro."""
    w = np.abs(vec) ** 2
    even = sum(w[i] for i, l in enumerate(ls) if l % 2 == 0)
    return abs(2.0 * even - w.sum()) / w.sum()


def test_t5_sin_campo_los_autovectores_son_gerade_o_ungerade(trimer):
    """
    Con F=0 el bloque Σ se parte en l pares (gerade) y l impares (ungerade):
    es la afirmación que el paper cita de su Ref. [15].
    """
    ls = trimer.l_values(0)
    checked = 0
    for R in (1200.0, 1800.0, 2400.0):
        vals, vecs = np.linalg.eigh(trimer.hamiltonian(R, 0.0, 0))
        # V tiene rango bajo: la mayor parte del manifold queda EXACTAMENTE
        # degenerada en E_manifold, y dentro de ese subespacio `eigh` devuelve
        # una base arbitraria que mezcla paridades sin que eso signifique nada.
        # La afirmación física es sobre los estados que sí se separan.
        gap = np.diff(np.sort(vals))
        isolated = [k for k in range(len(vals))
                    if min(gap[max(k - 1, 0)], gap[min(k, len(gap) - 1)]) > 1e-14]
        assert len(isolated) >= 6
        for k in isolated:
            assert _parity_purity(vecs[:, k], ls) == pytest.approx(1.0, abs=1e-8)
            checked += 1
    assert checked >= 18


# ------------------------------------------------------------------ T6
def test_t6_el_campo_rompe_la_paridad(trimer):
    """F·z tiene Δl=±1: acopla gerade con ungerade (§III.A del paper)."""
    ls = trimer.l_values(0)
    _, vecs = np.linalg.eigh(trimer.hamiltonian(1800.0, field_au(500.0), 0))
    purity = [_parity_purity(vecs[:, k], ls) for k in range(vecs.shape[1])]
    assert min(purity) < 0.5


# ------------------------------------------------------------------ T7
def test_t7_limite_de_acoplamiento_nulo():
    """A_s = A_p = 0 y F = 0: los autovalores son las energías atómicas."""
    sc = ScatteringLengths(enabled=False, n_star_table=float(N))
    t = SymmetricLinearTrimer(n_manifold=N, scattering=sc)
    for m in (0, 1):
        np.testing.assert_allclose(
            t.energies(1800.0, 0.0, m),
            np.sort((t.diagonal(m) - t.E_manifold) * HARTREE_TO_GHZ),
            atol=1e-10,
        )


# ------------------------------------------------------------------ T8
def test_t8_escalera_stark_lineal_del_manifold_hidrogenoide():
    """
    Manifold hidrogenoide puro (l = 0..n−1 degenerados) en campo F:
    autovalores exactos (3/2)·n·q·F con q = −(n−1), …, +(n−1) de dos en dos.

    Valida `Atom.Vfield` y el convenio de `exp_val_r.txt` a la vez: los
    elementos dipolares se ponen a mano con ⟨n,l|r|n,l+1⟩ = −(3/2)n√(n²−(l+1)²)
    y el resultado tiene que ser la escalera analítica.
    """
    atom, F = Atom(N, 0), 1e-9
    d = np.array([-1.5 * N * np.sqrt(N**2 - (l + 1) ** 2) for l in range(N - 1)])
    H = np.zeros((N, N))
    for l1 in range(N):
        for l2 in range(N):
            if abs(l1 - l2) == 1:
                H[l1, l2] = atom.Vfield(l1, l2, 0, 0, d[min(l1, l2)], F)
    got = np.sort(np.linalg.eigvalsh(H))
    want = np.sort(1.5 * N * np.arange(-(N - 1), N, 2) * F)
    np.testing.assert_allclose(got, want, rtol=1e-9, atol=1e-18)

    # El signo global de d es un gauge: cambiarlo no mueve el espectro.
    H2 = np.zeros((N, N))
    for l1 in range(N):
        for l2 in range(N):
            if abs(l1 - l2) == 1:
                H2[l1, l2] = atom.Vfield(l1, l2, 0, 0, -d[min(l1, l2)], F)
    np.testing.assert_allclose(np.sort(np.linalg.eigvalsh(H2)), want,
                               rtol=1e-9, atol=1e-18)


def test_t8b_exp_val_r_reproduce_los_dipolos_del_manifold(trimer):
    """
    De l ≥ 3 en adelante, `exp_val_r.txt` ES −(3/2)n√(n²−(l+1)²): el mismo
    convenio de signo que `hydrogenic_R`, que es el que usa V_Fermi.
    """
    for l in range(3, N - 1):
        want = -1.5 * N * np.sqrt(N**2 - (l + 1) ** 2)
        # El fichero sólo llega a 3.3e-4 relativo (peor caso l=7): residuo de
        # la integración numérica original. Por eso el módulo usa la forma
        # cerrada para l >= 3 y reserva el fichero para los tres vecinos.
        assert trimer.dipole_table[l] == pytest.approx(want, rel=5e-4)
        assert trimer.dipole[l] == want
    # Los tres vecinos sí salen del fichero, tal cual.
    np.testing.assert_array_equal(trimer.dipole[:3], trimer.dipole_table[:3])


# ------------------------------------------------------------------ T9
def test_t9_con_n_35_el_remapeo_es_la_identidad():
    """
    La tabla se generó para n=35: pedirle A_s/A_p con n*=35 la lee tal cual,
    sin remapeo ni interpolación, hasta su último nodo (2448 a₀).
    """
    sc = ScatteringLengths(n_star_table=float(N))
    for i in (0, 100, 213, 400, len(sc.R_table) - 1):
        R = float(sc.R_table[i])
        assert sc.remap_R(R, float(N)) == R
        A_s, A_p = sc.scattering(R, float(N))
        assert A_s == pytest.approx(sc.A_s_table[i], rel=1e-14)
        assert A_p == pytest.approx(sc.A_p_table[i], rel=1e-12)

    assert sc.R_table.max() == 2448.0
    assert 2.0 * N**2 == 2450  # punto de retorno clásico: la tabla llega justo
