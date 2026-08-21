"""
Anclas cuantitativas del trímero lineal simétrico contra el texto del paper.

Aguilera-Fernández, Schmelcher & González-Férez, J. Phys. B 49, 124002 (2016),
arXiv:1601.05049. El paper NO trae ninguna tabla numérica: todo lo comparable
son valores sueltos citados en el texto de §III.A, más las escalas de los ejes
de las Figs. 3-5. Cada test dice literalmente qué frase ancla.

Estos números fijan la física vigente de esta línea. Si uno se mueve, es un
cambio de física: para y repórtalo, igual que con los golden files.
"""

import numpy as np
import pytest

from trimero.systems.rb_neutral_perturber.fermi_krb import ScatteringLengths
from trimero.systems.rb_neutral_perturber.linear_trimer import (
    HARTREE_TO_GHZ,
    SymmetricLinearTrimer,
    field_au,
)

N = 35
pytestmark = pytest.mark.slow


@pytest.fixture(scope="module")
def trimer():
    return SymmetricLinearTrimer(n_manifold=N)


def _manifold_sorted(t, R, F_au, m_l, k=4):
    """Las k energías más bajas con peso de manifold > 0.5, en GHz."""
    E, W = t.spectrum(float(R), F_au, m_l)
    return np.sort(E[W > 0.5])[:k]


# ------------------------------------------------------------------ base
def test_energias_de_los_tres_vecinos(trimer):
    """
    «we include in the basis set the degenerate manifold Rb(n = 35, l ≥ 3) and
    the energetically closest neighboring Rydberg levels 38s, 37p and 36d»

    El 38s cae a −20.27 GHz del manifold, que es exactamente donde la Fig. 3(b)
    dibuja la curva Rb(5s)Rb(38s)Rb(5s) (eje del panel: −20 … 0 GHz).
    """
    rel = (trimer.diagonal(0) - trimer.E_manifold) * HARTREE_TO_GHZ
    ls = trimer.l_values(0)
    got = {l: rel[ls.index(l)] for l in (0, 1, 2, 3)}
    assert got[0] == pytest.approx(-20.267, abs=0.001)   # 38s
    assert got[1] == pytest.approx(-102.360, abs=0.001)  # 37p
    assert got[2] == pytest.approx(-54.019, abs=0.001)   # 36d
    assert got[3] == 0.0                                 # 35f: defecto despreciado


def test_resonancia_de_onda_p_de_la_tabla():
    """
    «The resonance of the p-wave scattering length at R ≈ 780 a0»

    La tabla `rvsAP.dat` sitúa el polo en 759.3 a₀, un 2.7 % por debajo del
    valor citado. Es una propiedad de LA TABLA, no del Hamiltoniano: se fija
    aquí para que la discrepancia quede registrada y no se redescubra.
    """
    w = ScatteringLengths(n_star_table=float(N)).p_resonance_window(float(N))
    assert w["R_pole"] == pytest.approx(759.3, abs=0.5)
    assert w["R_fwhm_lo"] == pytest.approx(747.0, abs=0.5)
    assert w["R_fwhm_hi"] == pytest.approx(771.0, abs=0.5)


# --------------------------------------------------------------- Σ, F = 0
def test_sigma_solo_onda_s_cabe_en_el_eje_de_la_fig_3a():
    """
    Fig. 3(a): sólo onda s, eje ε ∈ [−10, 0] GHz para R ∈ [500, 3000] a₀.
    El mínimo del cálculo es −9.87 GHz: apura el eje sin salirse.
    """
    t = SymmetricLinearTrimer(n_manifold=N, p_wave=False)
    R = t.table_R(R_min=500.0)
    E, W = t.spectra(R, 0.0, 0)
    depth = np.nanmin(np.where(W > 0.5, E, np.nan))
    assert depth == pytest.approx(-9.871, abs=0.01)
    assert -10.0 < depth < -9.0


def test_sigma_cruces_evitados_cerca_de_1500(trimer):
    """
    «The p-wave and s-wave dominated Σ-APCs suffer several avoided crossings
    close to the internuclear distance R ≈ 1500 a0»

    Salen tres contactos entre 1500 y 1560 a₀. Dos son cruces REALES (paridad
    opuesta, gap ~0) y el de en medio es el cruce evitado de verdad.
    """
    R = trimer.table_R(R_min=1350.0, R_max=1650.0)
    low = np.array([_manifold_sorted(trimer, r, 0.0, 0) for r in R])
    mins = []
    for k in range(3):
        g = low[:, k + 1] - low[:, k]
        j = int(np.argmin(g))
        mins.append((float(R[j]), float(g[j])))
    assert [r for r, _ in mins] == pytest.approx([1557.0, 1524.0, 1503.0], abs=6.0)
    assert mins[1][1] == pytest.approx(0.251, abs=0.01)   # el evitado
    assert mins[0][1] < 0.01 and mins[2][1] < 0.01        # los cruces reales


def test_sigma_pares_gerade_ungerade_degeneran_a_R_grande(trimer):
    """
    «At large separations ... the s-wave (p-wave) dominated Σ molecular states
    become degenerate and converge to the s-wave (p-wave) Σ-APC of the
    diatomic ULRM.»
    """
    e = _manifold_sorted(trimer, 2400.0, 0.0, 0)
    assert e[1] - e[0] == pytest.approx(0.074, abs=0.01)
    assert e[3] - e[2] == pytest.approx(0.001, abs=0.005)

    # ...y ese par coincide con el dímero de un solo perturbador en θ=0.
    dimer = SymmetricLinearTrimer(n_manifold=N, n_perturbers=1)
    ed = _manifold_sorted(dimer, 2400.0, 0.0, 0)
    assert ed[0] == pytest.approx(0.5 * (e[0] + e[1]), abs=0.15)


# --------------------------------------------------------------- Π, F = 0
def test_pi_cruce_libre_de_campo_en_1060(trimer):
    """
    «the crossing of the field-free APCs at R ≈ 1060 a0»

    Es un cruce REAL, no evitado: los dos estados tienen paridad l opuesta y
    sin campo no se acoplan. En malla fina el gap baja a 0.016 GHz en 1060.5 a₀.
    """
    R = np.arange(1050.0, 1075.0, 0.25)
    gaps = np.array([np.diff(_manifold_sorted(trimer, r, 0.0, 1, k=2))[0] for r in R])
    j = int(np.argmin(gaps))
    assert R[j] == pytest.approx(1060.5, abs=1.0)
    assert gaps[j] < 0.05

    ls = np.array(trimer.l_values(1))
    for r in (1056.0, 1068.0):
        E, W = trimer.spectrum(r, 0.0, 1)
        _, v = np.linalg.eigh(trimer.hamiltonian(r, 0.0, 1))
        idx = [i for i in np.argsort(E) if W[i] > 0.5][:2]
        even = sorted(float((np.abs(v[ls % 2 == 0, i]) ** 2).sum()) for i in idx)
        assert even == pytest.approx([0.0, 1.0], abs=1e-9)


def test_pi_minimo_en_1115_y_su_stark(trimer):
    """
    «the energy of the minimum appearing at R ≈ 1115 a0 for the lowest lying
    Π-APC is shifted 0.3 GHz from its field-free value for F = 500 V/m»

    Sale el mínimo en 1116 a₀ y un desplazamiento de −0.244 GHz. El paper no
    da el signo; aquí baja (más ligado), que es lo que hace la repulsión de
    niveles sobre un estado por debajo del manifold.
    """
    R = trimer.table_R(R_min=1050.0, R_max=1200.0)

    def lowest(F):
        return np.array([_manifold_sorted(trimer, r, field_au(F), 1, k=1)[0] for r in R])

    l0, l5 = lowest(0.0), lowest(500.0)
    j = int(np.argmin(l0))
    assert R[j] == pytest.approx(1116.0, abs=4.0)
    assert l0[j] == pytest.approx(-34.267, abs=0.01)
    assert l5.min() - l0.min() == pytest.approx(-0.244, abs=0.02)
    assert abs(l5.min() - l0.min()) == pytest.approx(0.3, abs=0.1)


# --------------------------------------------------------------- con campo
def test_el_campo_abre_el_cruce_pi(trimer):
    """
    «this crossing ... becomes an avoided crossing in the presence of the dc
    field». El gap crece monótonamente con F.
    """
    R = trimer.table_R(R_min=1020.0, R_max=1110.0)
    gaps = []
    for F in (0.0, 100.0, 300.0, 500.0):
        g = [np.diff(_manifold_sorted(trimer, r, field_au(F), 1, k=2))[0] for r in R]
        gaps.append(min(g))
    assert gaps == pytest.approx([0.867, 0.971, 1.581, 2.412], abs=0.02)
    assert np.all(np.diff(gaps) > 0)


def test_stark_cuadratico_del_38s(trimer):
    """
    «due to the quadratic Stark shift of the 38s Rydberg state, the APC of
    Rb(5s)Rb(38s)Rb(5s) is very weakly affected by the electric field»

    Los desplazamientos escalan como F² (razón 1 : 9 : 25) y a 500 V/m valen
    0.27 GHz, frente a los ±11 GHz lineales del manifold.
    """
    shifts = {}
    for F in (0.0, 100.0, 300.0, 500.0):
        E, W = trimer.spectrum(2100.0, field_au(F), 0)
        near = E[(W <= 0.5) & (np.abs(E + 20.28) < 3.0)]
        assert near.size == 1
        shifts[F] = float(near[0])
    base = shifts[0.0]
    assert base == pytest.approx(-20.281, abs=0.005)
    d = np.array([shifts[F] - base for F in (100.0, 300.0, 500.0)])
    assert d == pytest.approx([0.0110, 0.0978, 0.2723], abs=0.002)
    np.testing.assert_allclose(d / d[0], [1.0, 9.0, 25.0], rtol=0.02)


def test_stark_lineal_del_manifold_a_R_grande(trimer):
    """
    «At large internuclear distances ... all the APCs from the Rb(n = 35, l≥3)
    Rydberg manifold are shifted linearly in energy with the dc field strength,
    which corresponds with the Stark shift of the (n = 35, l ≥ 3) Rydberg
    manifold of an isolated Rb atom.»

    Para Π el abanico calculado coincide con ±(3/2)·n·(n−1−|m|)·F dentro del
    2.2 %: en 2448 a₀ el pseudopotencial todavía no es cero del todo y el 36d
    sigue empujando desde abajo. Σ se aparta más porque su bloque incluye
    además 38s y 37p.
    """
    for F in (100.0, 300.0, 500.0):
        E, W = trimer.spectrum(2448.0, field_au(F), 1)
        man = E[W > 0.5]
        edge = 1.5 * N * (N - 1 - 1) * field_au(F) * HARTREE_TO_GHZ
        assert man.max() == pytest.approx(+edge, rel=0.03)
        assert man.min() == pytest.approx(-edge, rel=0.03)
        # y el abanico es SIMÉTRICO y LINEAL en F, que es lo que afirma el texto
        assert man.max() + man.min() == pytest.approx(0.0, abs=0.04 * edge)


def test_el_convenio_de_signo_de_los_vecinos_es_gauge_sin_campo():
    """
    Sin campo, invertir el signo de una radial tabulada es una transformación
    diagonal ±1: el espectro NO se mueve. Con campo sí, porque los dipolos de
    `exp_val_r.txt` no cambian con él. Cuantifica la aproximación 2 de
    `linear_trimer`: 1.2 GHz en el peor punto de Σ a 500 V/m, 0.006 GHz de
    mediana; Π es insensible (0.002 GHz).
    """
    base = SymmetricLinearTrimer(n_manifold=N)
    flip = SymmetricLinearTrimer(n_manifold=N, neighbor_sign={0: -1.0, 1: -1.0, 2: -1.0})
    R = base.table_R(R_min=1000.0)

    def low(t, F, m):
        return np.array([_manifold_sorted(t, r, field_au(F), m) for r in R])

    np.testing.assert_allclose(low(base, 0.0, 0), low(flip, 0.0, 0), atol=1e-9)
    np.testing.assert_allclose(low(base, 0.0, 1), low(flip, 0.0, 1), atol=1e-9)

    d_sigma = np.abs(low(base, 500.0, 0) - low(flip, 500.0, 0))
    d_pi = np.abs(low(base, 500.0, 1) - low(flip, 500.0, 1))
    assert d_sigma.max() == pytest.approx(1.156, abs=0.05)
    assert np.median(d_sigma) < 0.02
    assert d_pi.max() < 0.01
