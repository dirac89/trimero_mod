"""G3: las matrices intermedias del legado no cambian.

Reconstruye spV y las llamadas a Atom.Vfield interceptando el legado igual
que hace el generador, y compara bit a bit.
"""
from pathlib import Path

import numpy as np

from generate_goldens_g3 import capture  # noqa: F401  (mismo directorio)

GOLDEN = Path(__file__).resolve().parent / "goldens" / "g3_matrices.npz"


def test_g3_spv_and_vfield_bit_for_bit():
    g = np.load(GOLDEN)
    for dc, tag in [(0.0, "dc0"), (0.1, "dc0p1")]:
        spv, vf = capture(dc)
        assert np.array_equal(spv, g[f"spV_{tag}"]), f"spV cambió (dc={dc})"
        assert np.array_equal(vf, g[f"vfield_{tag}"]), f"Vfield cambió (dc={dc})"


def test_g3_golden_is_not_vacuous():
    g = np.load(GOLDEN)
    d0, d1 = g["spV_dc0"], g["spV_dc0p1"]
    assert np.count_nonzero(d0) > 0 and np.all(np.isfinite(d0))
    assert not np.array_equal(d0, d1), (
        "spV no depende del campo DC: el caso dc=0.1 no aporta cobertura"
    )
    assert not np.array_equal(d0[0], d0[1]), "las dos filas capturadas son iguales"


def test_g3_spv_symmetric_without_field():
    """Sin campo DC el Hamiltoniano SÍ es real simétrico."""
    g = np.load(GOLDEN)
    for row in g["spV_dc0"]:
        assert np.array_equal(row, row.T), "spV no es simétrica con dc=0"


def test_g3_spv_asymmetry_is_the_known_legacy_bug():
    """Con campo DC, spV es asimétrica. NO es ruido: es un bug del legado.

    En `Trimer_energies_field` la matriz `field` se construye con

        for i in range(n1 - 1):     # i llega sólo hasta n1-2
            for j in range(n1):     # j llega hasta n1-1

    así que el par (i=n1-1, j=n1-2) nunca se rellena mientras que
    (i=n1-2, j=n1-1) sí. El acoplamiento Stark entre las dos últimas capas l
    queda sólo en el triángulo superior, y como `np.linalg.eigh` lee por
    defecto el triángulo INFERIOR (UPLO='L'), ese acoplamiento se descarta.

    Este test congela el bug tal cual, para que arreglarlo sea una decisión
    consciente que obligue a regenerar el golden, no un cambio accidental.
    """
    g = np.load(GOLDEN)
    n1 = int(g["n1"])
    l_of = lambda k: int(np.floor(np.sqrt(k)))

    for row in g["spV_dc0p1"]:
        asym = np.argwhere(np.abs(row - row.T) > 0)
        assert len(asym) > 0, "la asimetría desapareció: ¿se arregló el bug?"
        blocks = {(l_of(r), l_of(c)) for r, c in asym}
        assert blocks == {(n1 - 2, n1 - 1), (n1 - 1, n1 - 2)}, (
            f"la asimetría ya no está confinada al borde del manifold: {blocks}"
        )
        # el elemento perdido vive arriba; abajo, que es lo que eigh lee, hay 0
        for r, c in asym:
            lo, hi = (r, c) if r > c else (c, r)
            assert row[lo, hi] == 0.0
            assert row[hi, lo] != 0.0
