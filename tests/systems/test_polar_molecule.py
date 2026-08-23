import numpy as np
import pytest

from trimero.systems.polar_molecule import KRB, RBCS, PolarMolecule, get_molecule
from trimero.systems.rb_krb_polar.charge_dipole import (
    B_KRB_AU,
    D_KRB_AU,
    HZ_PER_HARTREE,
)


def test_catalogue_conversions_match_validated_constants():
    assert KRB.B_au == B_KRB_AU
    assert KRB.d_au == D_KRB_AU
    assert RBCS.B_au == 490.17e6 / HZ_PER_HARTREE
    assert np.isclose(RBCS.d_au, 0.481952126075, rtol=0.0, atol=1e-15)
    assert get_molecule("RbCs") is RBCS


@pytest.mark.parametrize(
    "kwargs", [
        {"rotational_constant_hz": 0.0, "dipole_debye": 1.0},
        {"rotational_constant_hz": 1.0, "dipole_debye": -1.0},
    ]
)
def test_invalid_molecular_constants_are_rejected(kwargs):
    with pytest.raises(ValueError):
        PolarMolecule("x", "X", **kwargs)
