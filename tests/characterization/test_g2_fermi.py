"""G2: el potencial de Fermi no cambia de valor en ninguna de sus 4 ramas.

Este es el juez del paso 7 del refactor, donde las 8 ramas duplicadas de
Vs/Vp colapsan a una sola expresión vía WavefunctionSource. A partir de ese
paso, `Vsp()` deja de existir y la comparación pasa a `matrix_element`; los
valores esperados no cambian, sólo la llamada que los produce.
"""
from pathlib import Path

import numpy as np

from fermi_potentials import FermiPotentials

GOLDEN = Path(__file__).resolve().parent / "goldens" / "g2_fermi.npz"
ORDER = ["Ap1", "As1", "Dwave1", "Dwave2", "li", "lj", "mi", "mj",
         "n", "n2", "r1", "s", "theta1", "wave1", "wave2"]


def _cases():
    g = np.load(GOLDEN, allow_pickle=True)
    names = [str(x) for x in g["param_names"]]
    assert names == ORDER, f"orden de parámetros inesperado: {names}"
    return g, [dict(zip(names, row)) for row in g["params"]]


def _build(c):
    return FermiPotentials(
        c["s"], c["n"], c["n2"], int(c["li"]), int(c["lj"]),
        int(c["mi"]), int(c["mj"]), c["r1"], c["theta1"],
        c["As1"], c["wave1"], c["wave2"], c["Ap1"], c["Dwave1"], c["Dwave2"],
    )


def test_g2_vs_vp_vsp_bit_for_bit():
    g, cases = _cases()
    got_s = np.array([_build(c).Vs() for c in cases])
    got_p = np.array([_build(c).Vp() for c in cases])
    got_sp = np.array([_build(c).Vsp() for c in cases])

    assert np.array_equal(got_s, g["Vs"]), "Vs cambió"
    assert np.array_equal(got_p, g["Vp"]), "Vp cambió"
    assert np.array_equal(got_sp, g["Vsp"]), "Vsp cambió"


def test_g2_covers_all_four_branches():
    """Cada una de las 4 ramas li/lj debe aportar casos NO nulos.

    Sin esto, una rama podría estar cubierta sólo por ceros triviales y el
    golden no detectaría que se rompió.
    """
    g, cases = _cases()
    li = np.array([c["li"] for c in cases])
    lj = np.array([c["lj"] for c in cases])
    vs = g["Vs"]
    for hi_i, hi_j in [(True, True), (False, True), (True, False), (False, False)]:
        mask = ((li > 2) == hi_i) & ((lj > 2) == hi_j)
        assert mask.sum() > 0, f"rama li>2={hi_i}, lj>2={hi_j} sin casos"
        assert np.count_nonzero(vs[mask]) > 0, (
            f"rama li>2={hi_i}, lj>2={hi_j} sólo tiene ceros: el golden sería vacuo"
        )
