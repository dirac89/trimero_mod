"""G4: la salida final del legado (autovalores por radio) no cambia.

Reejecuta `Trimer_energies_field` de punta a punta y compara los dos ficheros
que produce contra el golden. Es la verificación end-to-end del refactor: si
G1-G3 pasan pero G4 falla, el fallo está en el ensamblado, no en la física de
los componentes.

Tolerancia rtol=1e-12 en lugar de bit a bit porque `np.linalg.eigh` delega en
LAPACK y el último bit puede variar entre versiones de BLAS; exigir igualdad
exacta daría fallos espurios en otra máquina.

Marcado `slow`: cada caso tarda ~3.5 min. Para iterar rápido:
    poetry run pytest -m "not slow"
"""
import os
import sys
import tempfile
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[2]
GOLDEN = Path(__file__).resolve().parent / "goldens" / "g4_eigenvalues.npz"
RTOL = 1e-12


def _run_legacy(dc_field_au, n1=5):
    """Ejecuta el legado en un directorio limpio y devuelve sus dos salidas."""
    from trimer import Trimer_energies_field

    cwd = Path.cwd()
    with tempfile.TemporaryDirectory() as tmp:
        (Path(tmp) / "data").symlink_to(REPO / "data")
        os.chdir(tmp)
        try:
            Trimer_energies_field(n1, dc_field_au)
            return {
                unit: np.loadtxt(f"Trimer_R_sp_wave_N35_R_300_{unit}.dat")
                for unit in ("GHz", "au")
            }
        finally:
            os.chdir(cwd)


@pytest.mark.slow
@pytest.mark.parametrize("dc,tag", [(0.0, "dc0"), (0.1, "dc0p1")])
def test_g4_end_to_end(dc, tag):
    g = np.load(GOLDEN)
    got = _run_legacy(dc)
    for unit in ("GHz", "au"):
        want = g[f"{tag}_{unit}"]
        assert got[unit].shape == want.shape, f"cambió la forma de la salida ({unit})"
        np.testing.assert_allclose(
            got[unit], want, rtol=RTOL, atol=0.0,
            err_msg=f"la salida final cambió (dc={dc}, {unit})",
        )


def test_g4_golden_is_not_vacuous():
    """El golden debe contener autovalores reales, variados y finitos."""
    g = np.load(GOLDEN)
    for key in g.files:
        arr = g[key]
        assert np.all(np.isfinite(arr)), f"{key} contiene no finitos"
        assert arr.shape[0] > 100, f"{key} tiene muy pocas filas"
        # los autovalores (columnas 1+) deben variar con R, no ser constantes
        spread = np.ptp(arr[:, 1:], axis=0)
        assert np.count_nonzero(spread) > 0, f"{key}: autovalores constantes en R"


def test_g4_field_changes_the_spectrum():
    """dc=0.1 debe dar un espectro distinto de dc=0; si no, no aporta cobertura."""
    g = np.load(GOLDEN)
    assert not np.array_equal(g["dc0_au"], g["dc0p1_au"])


def test_g4_documents_the_unit_bug():
    """El fichero rotulado GHz contiene MHz (paso 4 del refactor).

    1 E_h = 6.5797e15 Hz = 6.5797e6 GHz, pero el factor aplicado es 6.5797e9,
    que es Hartree->MHz. Este test fija el comportamiento ACTUAL; el paso 4 lo
    corregirá y deberá actualizar tanto este test como el golden en GHz.
    """
    g = np.load(GOLDEN)
    au, ghz = g["dc0_au"], g["dc0_GHz"]
    mask = au[:, 1:] != 0.0
    ratio = np.median(ghz[:, 1:][mask] / au[:, 1:][mask])
    assert ratio == pytest.approx(6.579683920729e9, rel=1e-9), (
        f"el factor de conversión cambió: {ratio:.9e}"
    )
    hartree_to_ghz = 6.579683920502e15 / 1e9
    assert ratio / hartree_to_ghz == pytest.approx(1000.0, rel=1e-6), (
        "el factor ya no está 1000x desviado: ¿se aplicó el paso 4?"
    )
