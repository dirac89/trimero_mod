"""Configuración de producción del sistema polar Rb*+RbCs."""

from trimero.systems.polar_molecule import RBCS
from trimero.systems.polar_rydberg import PolarBOPSystem

__all__ = ["RBCS", "RbRbCsPolarSystem"]


class RbRbCsPolarSystem(PolarBOPSystem):
    def __init__(self, **kwargs):
        if "molecule" in kwargs:
            raise TypeError("RbRbCsPolarSystem fija molecule=RBCS")
        super().__init__(molecule=RBCS, **kwargs)
