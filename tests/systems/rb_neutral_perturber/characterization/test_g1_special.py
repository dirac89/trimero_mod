"""G1: las funciones especiales puras no cambian de valor.

Tolerancia bit a bit: es aritmética determinista sin diagonalización de por
medio, así que cualquier diferencia es un cambio real, no ruido numérico.
"""
from pathlib import Path

import numpy as np

from trimero.mathlib.special import Spherical, DRnl, DOlm, DPhilm, hydrogenicR

GOLDEN = Path(__file__).resolve().parent / "goldens" / "g1_special.npz"


def _golden():
    return np.load(GOLDEN, allow_pickle=True)


def test_g1_radial_functions_bit_for_bit():
    g = _golden()
    got_h, got_d = [], []
    for (n, l) in g["nl"]:
        for r in g["r"]:
            got_h.append(hydrogenicR(int(n), int(l), 1.0, float(r)))
            got_d.append(DRnl(int(n), int(l), float(r)))

    assert np.array_equal(np.array(got_h), g["hydrogenicR"]), "hydrogenicR cambió"
    assert np.array_equal(np.array(got_d), g["DRnl"]), "DRnl cambió"


def test_g1_angular_functions_bit_for_bit():
    g = _golden()
    got_s, got_o, got_p = [], [], []
    for (l, m) in g["lm"]:
        for th in g["theta"]:
            got_s.append(Spherical(int(l), int(m), np.cos(float(th))))
            got_o.append(DOlm(int(l), int(m), float(th)))
            got_p.append(DPhilm(int(l), int(m), float(th)))

    assert np.array_equal(np.array(got_s), g["Spherical"]), "Spherical cambió"
    assert np.array_equal(np.array(got_o), g["DOlm"]), "DOlm cambió"
    assert np.array_equal(np.array(got_p), g["DPhilm"]), "DPhilm cambió"


def test_g1_golden_is_not_vacuous():
    """Un golden de ceros pasaría los tests anteriores sin vigilar nada."""
    g = _golden()
    for key in ("hydrogenicR", "DRnl", "Spherical", "DOlm", "DPhilm"):
        arr = g[key]
        assert np.all(np.isfinite(arr)), f"{key} contiene no finitos"
        assert np.count_nonzero(arr) > 0, f"{key} es idénticamente nulo"
