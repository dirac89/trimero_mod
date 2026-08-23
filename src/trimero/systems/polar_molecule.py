"""Parámetros de rotores polares usados por los sistemas Rydberg."""

from dataclasses import dataclass

from trimero.systems.rb_krb_polar.charge_dipole import (
    DEBYE_TO_EA0,
    HZ_PER_HARTREE,
)

__all__ = ["PolarMolecule", "KRB", "RBCS", "MOLECULES", "get_molecule"]


@dataclass(frozen=True)
class PolarMolecule:
    """Constantes espectroscópicas de un rotor rígido polar."""

    key: str
    label: str
    rotational_constant_hz: float
    dipole_debye: float

    def __post_init__(self) -> None:
        if not self.key or not self.label:
            raise ValueError("key y label no pueden estar vacíos")
        if self.rotational_constant_hz <= 0.0:
            raise ValueError("la constante rotacional debe ser positiva")
        if self.dipole_debye < 0.0:
            raise ValueError("el momento dipolar no puede ser negativo")

    @property
    def B_au(self) -> float:
        return self.rotational_constant_hz / HZ_PER_HARTREE

    @property
    def B_ghz(self) -> float:
        return self.rotational_constant_hz / 1.0e9

    @property
    def d_au(self) -> float:
        return self.dipole_debye * DEBYE_TO_EA0


KRB = PolarMolecule("krb", "KRb", 1.114e9, 0.566)
RBCS = PolarMolecule("rbcs", "RbCs", 490.17e6, 1.225)
MOLECULES = {m.key: m for m in (KRB, RBCS)}


def get_molecule(value: str | PolarMolecule) -> PolarMolecule:
    if isinstance(value, PolarMolecule):
        return value
    try:
        return MOLECULES[value.lower()]
    except (AttributeError, KeyError) as exc:
        choices = ", ".join(sorted(MOLECULES))
        raise ValueError(f"molécula desconocida {value!r}; opciones: {choices}") from exc
