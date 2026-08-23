"""
Regresión de `scripts/compute_bop_curve.py` contra los números ya verificados.

La referencia es `plots/rb_krb_polar/data/fig1_ad_MJ0_n25.npz`, producida por el script archivado
`scripts/archive/run_fig1_charge_dipole.py` y documentada en
`docs/analysis_fig1_carga_dipolo_sin_fermi.md` §5. Esos números están
verificados de forma independiente, así que aquí son la fuente de verdad: si
algo diverge, lo que está mal es el código, no la referencia.

Dos niveles:

  R1  el .npz de referencia sigue diciendo lo que dice el documento. Rápido, no
      diagonaliza nada. Protege contra que alguien regenere el .npz y cambie la
      referencia sin darse cuenta.

  R2  `compute_bop_curve.py` reproduce ese .npz. Marcado `slow`: son tres
      diagonalizaciones de 1113x1113 (~1.1 s cada una) más el montaje de
      `BOPSystem`. Se comparan con rtol=1e-12, no «aproximadamente»: el script
      nuevo llama a `BOPSystem.hamiltonian(R, M_J, fermi=False)`, que es la
      MISMA expresión que el `hamiltonian_ad()` del script archivado, así que
      debe salir bit a bit lo mismo.
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[3]
REF = REPO / "plots" / "rb_krb_polar" / "data" / "fig1_ad_MJ0_n25.npz"

# Números de docs/analysis_fig1_carga_dipolo_sin_fermi.md §5, n=25, M_J=0.
DOC_DEPTH_GHZ = -23.100
DOC_E_AT_1800_GHZ = -0.338
DOC_N_LOCAL_MINIMA = 8


def _load_script():
    """Importa scripts/compute_bop_curve.py sin instalarlo como paquete."""
    path = REPO / "scripts" / "compute_bop_curve.py"
    spec = importlib.util.spec_from_file_location("compute_bop_curve", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["compute_bop_curve"] = mod
    spec.loader.exec_module(mod)
    return mod


def test_r1_reference_npz_matches_the_documented_numbers():
    d = np.load(REF)
    R, E = d["R"], d["E"]

    assert (R[0], R[-1], len(R)) == (400.0, 1800.0, 281), "cambió la malla en R"

    mod = _load_script()
    n_min = len(mod.local_minima(E))

    assert round(float(np.nanmin(E)), 3) == DOC_DEPTH_GHZ
    assert round(float(E[-1]), 3) == DOC_E_AT_1800_GHZ
    assert n_min == DOC_N_LOCAL_MINIMA


@pytest.mark.slow
def test_r2_compute_bop_curve_reproduces_the_reference():
    from trimero.systems.rb_krb_polar.bop_system import BOPSystem
    from trimero.systems.rb_krb_polar.rb_defects import DELTA0_NS_PAPER

    mod = _load_script()
    d = np.load(REF)
    R, E = d["R"], d["E"]

    sysm = BOPSystem(n_manifold=25, delta0_ns=DELTA0_NS_PAPER)
    for R_probe in (400.0, 1000.0, 1800.0):
        i = int(np.flatnonzero(R == R_probe)[0])
        got, k, w, _ = mod.character_curve(sysm, 0, R_probe, weight=0.5)
        np.testing.assert_allclose(
            got, E[i], rtol=1e-12, atol=0.0,
            err_msg=f"la curva cambió en R = {R_probe} a0")
        assert k == int(d["K"][i]), f"cambió el índice de la curva en R = {R_probe}"
        np.testing.assert_allclose(w, d["W"][i], rtol=1e-12, atol=0.0)
