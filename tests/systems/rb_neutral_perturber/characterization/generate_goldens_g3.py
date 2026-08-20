"""
Genera el golden G3: las matrices intermedias del legado (paso 0 del refactor).

No modifica ni una línea del código legado. Intercepta en tiempo de ejecución
las dos fronteras por las que pasan las matrices:

  - `Atom.Vfield`      -> caracteriza la construcción de la matriz `field`
  - `np.linalg.eigh`   -> captura `spV` completa, justo antes de diagonalizar

Tras capturar las primeras filas aborta con una excepción propia, en lugar de
esperar a las 479 filas del bucle.

Uso:
    python tests/characterization/generate_goldens_g3.py
"""
import sys
import tempfile
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "src"))
GOLDENS = Path(__file__).resolve().parent / "goldens"

N1 = 5           # dim = 25; visita las 4 ramas i<3/i>2 x j<3/j>2
N_ROWS = 2       # filas 297 y 298 del bucle

from trimero.systems import rb_atom as atom  # noqa: E402
from trimero.systems.rb_neutral_perturber import trimer  # noqa: E402


class _EnoughRows(Exception):
    """Corta el bucle del legado cuando ya se capturó lo necesario."""


def capture(dc_field_au):
    spv_matrices = []
    vfield_calls = []

    orig_eigh = np.linalg.eigh
    orig_vfield = atom.Atom.Vfield

    def spy_eigh(m, *a, **kw):
        spv_matrices.append(np.array(m, copy=True))
        if len(spv_matrices) >= N_ROWS:
            raise _EnoughRows
        return orig_eigh(m, *a, **kw)

    def spy_vfield(self, li, lj, mi, mj, radial, strength):
        out = orig_vfield(self, li, lj, mi, mj, radial, strength)
        vfield_calls.append((li, lj, mi, mj, radial, strength, out))
        return out

    np.linalg.eigh = spy_eigh
    atom.Atom.Vfield = spy_vfield
    cwd = Path.cwd()
    try:
        with tempfile.TemporaryDirectory() as tmp:
            # el legado usa la ruta relativa "data/Wavefunction/"
            (Path(tmp) / "data").symlink_to(REPO / "data")
            import os
            os.chdir(tmp)
            try:
                trimer.Trimer_energies_field(N1, dc_field_au)
            except _EnoughRows:
                pass
    finally:
        np.linalg.eigh = orig_eigh
        atom.Atom.Vfield = orig_vfield
        import os
        os.chdir(cwd)

    return np.array(spv_matrices), np.array(vfield_calls, dtype=np.float64)


def _refuse_accidental_regeneration():
    """Los goldens son la referencia del refactor: regenerarlos desde el
    código ya refactorizado destruiría la red de seguridad sin avisar.

    Sobrescribir exige intención explícita:
        TRIMERO_REGENERATE_GOLDENS=1 python <este script>
    """
    import os
    existing = sorted(GOLDENS.glob("*.npz"))
    if existing and os.environ.get("TRIMERO_REGENERATE_GOLDENS") != "1":
        names = ", ".join(f.name for f in existing)
        raise SystemExit(
            f"Ya existen goldens ({names}).\n"
            "Regenerarlos invalida la referencia contra la que se verifica el "
            "refactor.\nSi de verdad quieres sobrescribirlos:\n"
            "    TRIMERO_REGENERATE_GOLDENS=1 python " + __file__
        )


if __name__ == "__main__":
    GOLDENS.mkdir(parents=True, exist_ok=True)
    _refuse_accidental_regeneration()
    out = {}
    for dc in (0.0, 0.1):
        tag = f"dc{dc:g}".replace(".", "p")
        spv, vf = capture(dc)
        out[f"spV_{tag}"] = spv
        out[f"vfield_{tag}"] = vf
        print(f"  dc={dc}: spV {spv.shape}, Vfield {vf.shape[0]} llamadas")
    np.savez(GOLDENS / "g3_matrices.npz", n1=N1, n_rows=N_ROWS, **out)
    print(f"G3 escrito -> goldens/g3_matrices.npz")
