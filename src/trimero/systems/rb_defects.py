"""
Override LOCAL del defecto cuántico s de Rb para la validación de Rb*-KRb.

POR QUÉ EXISTE ESTE MÓDULO
--------------------------
`systems/atom.py` usa δ₀(ns) = 3.1311804 (Li, Mourachko, Noel & Gallagher,
PRA 67, 052502, 2003). Ese valor es correcto y es la fuente de verdad para
`hamiltonians/trimer.py` y `hamiltonians/fermi.py`, protegidos por los golden
files de caracterización. **No se toca.**

Sin embargo, la Tabla I de González-Férez, Sadeghpour & Schmelcher, NJP 17,
013021 (2015) NO se reproduce con ese valor: las series p, d y f concuerdan a
±0.004 GHz, pero los dos niveles s (27s y 28s) fallan por -0.299 y +0.265 GHz.
Ambos piden el MISMO desplazamiento constante de δ₀, +6.19e-4 (ver
`docs/analysis_verificacion_tabla_I.md`, §5), lo que descarta δ₂, la
dependencia en n, la masa reducida y la conversión de unidades.

El paper obtiene sus energías de su referencia [20], Marinescu, Sadeghpour &
Dalgarno, PRA 49, 982 (1994), que es un cálculo de POTENCIAL MODELO, no un
ajuste Rydberg-Ritz. Que la discrepancia esté confinada a l=0 —el canal más
penetrante en el core— encaja con esa diferencia de método.

⚠️ ESTO NO ES UNA CORRECCIÓN GENERAL DEL CÓDIGO. Es un valor específico para
poder comparar contra este paper. El comportamiento por defecto de
`Atom.E_Rb()` no cambia, y `delta0_ns=None` en todas las funciones de abajo
delega exactamente en él.
"""

from trimero.systems.atom import Atom

__all__ = [
    "DELTA0_NS_ATOM",
    "DELTA0_NS_PAPER",
    "DELTA2_NS",
    "mu_ns",
    "n_star_ns",
    "energy_rb",
]

# El de systems/atom.py (Li+ 2003). Comportamiento por defecto del código.
DELTA0_NS_ATOM = 3.1311804
# El que reproduce la Tabla I de González-Férez 2015 (vía Marinescu+ 1994).
DELTA0_NS_PAPER = 3.13180
# δ₂ tal como está en atom.py. Ojo: Li+ 2003 da 0.1784; la diferencia mueve
# 0.003 GHz, dos órdenes por debajo del desajuste que resuelve δ₀. Se mantiene
# el valor de atom.py para no introducir un segundo cambio a la vez.
DELTA2_NS = 0.1745312


def mu_ns(n: int, delta0_ns: float = DELTA0_NS_ATOM,
          delta2_ns: float = DELTA2_NS) -> float:
    """Defecto cuántico de la serie s: μ(n) = δ₀ + δ₂·(n-δ₀)⁻² (Rydberg-Ritz)."""
    return delta0_ns + delta2_ns * (n - delta0_ns) ** -2


def n_star_ns(n: int, delta0_ns: float = DELTA0_NS_ATOM,
              delta2_ns: float = DELTA2_NS) -> float:
    """n* efectivo de un estado ns."""
    return n - mu_ns(n, delta0_ns, delta2_ns)


def energy_rb(n: int, l: int, delta0_ns: float = None) -> float:
    """
    Energía del nivel Rydberg de Rb en E_h.

    Args:
        delta0_ns: si es None (por defecto) delega EXACTAMENTE en
            `Atom.E_Rb()`, sin cambio alguno de comportamiento. Si se pasa un
            valor, sólo afecta a l=0; p, d y l≥3 siguen viniendo de `Atom`.
    """
    if l != 0 or delta0_ns is None:
        return Atom(n, l).E_Rb()
    return -0.5 / n_star_ns(n, delta0_ns) ** 2
