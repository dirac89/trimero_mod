"""
Tests de la utilidad de dibujo `trimero.visualization.geometry_diagram`.

**Aquí no se comprueba física ninguna.** Sólo que el módulo produce ficheros,
que generaliza a un número y a unos ángulos arbitrarios de cuerpos, y que
falla con un mensaje útil cuando la entrada no tiene sentido.

Un test de imagen «se parece a la figura del paper» no existe a propósito: la
comparación visual la hace una persona mirando el PNG, no un assert.
"""

import subprocess
import sys
from pathlib import Path

import matplotlib
import pytest

matplotlib.use("Agg")

from trimero.visualization.geometry_diagram import (  # noqa: E402
    Body,
    Electron,
    draw_from_system_config,
    draw_geometry,
)

REPO = Path(__file__).resolve().parents[2]


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    import matplotlib.pyplot as plt
    plt.close("all")


def _nonempty(path: Path) -> bool:
    return path.exists() and path.stat().st_size > 0


# ------------------------------------------------------- casos requeridos
def test_un_solo_cuerpo(tmp_path):
    out = tmp_path / "uno.png"
    fig, ax, paths = draw_geometry(
        bodies=[Body(kind="neutral_atom", label="Rb", R=800, theta_deg=0)],
        output_path=out, title="Un cuerpo")
    assert paths == [out]
    assert _nonempty(out)


def test_dos_cuerpos_configuracion_actual_del_hibrido(tmp_path):
    """Rb neutro en θ=180°, RbCs polar en θ=0°: lo que dibuja el ejemplo."""
    out = tmp_path / "hibrido.png"
    _, _, paths = draw_geometry(
        bodies=[
            Body(kind="neutral_atom", label="Rb", R=600, theta_deg=180),
            Body(kind="polar_molecule", label="RbCs", R=900, theta_deg=0,
                 dipole_deg=0),
        ],
        electron=Electron(r=1200, theta_deg=35),
        output_path=out, title="Híbrido")
    assert _nonempty(paths[0])


@pytest.mark.parametrize("n", [3, 5, 8])
def test_generaliza_a_n_cuerpos_y_angulos_arbitrarios(tmp_path, n):
    """Número arbitrario de cuerpos y ángulos que no son 0 ni 180."""
    bodies = [
        Body(kind="polar_molecule" if i % 2 else "neutral_atom",
             label=f"X{i}", R=500 + 180 * i, theta_deg=180.0 * i / (n - 1))
        for i in range(n)
    ]
    out = tmp_path / f"n{n}.png"
    _, _, paths = draw_geometry(bodies, output_path=out)
    assert _nonempty(paths[0])


# ------------------------------------------------------------- formatos
def test_png_y_pdf(tmp_path):
    out = tmp_path / "esquema.png"
    _, _, paths = draw_geometry(
        bodies=[Body(label="A", R=700, theta_deg=0),
                Body(label="B", R=900, theta_deg=180)],
        output_path=out, formats=("png", "pdf"))
    assert [p.suffix for p in paths] == [".png", ".pdf"]
    assert all(_nonempty(p) for p in paths)


def test_crea_los_directorios_intermedios(tmp_path):
    out = tmp_path / "a" / "b" / "c.png"
    _, _, paths = draw_geometry([Body(label="A", R=1, theta_deg=0)],
                                output_path=out)
    assert _nonempty(paths[0])


def test_sin_output_no_escribe_nada(tmp_path):
    _, _, paths = draw_geometry([Body(label="A", R=1, theta_deg=0)])
    assert paths == []
    assert list(tmp_path.iterdir()) == []


# ------------------------------------------------------------- opciones
@pytest.mark.parametrize("scale", ["proportional", "relative"])
def test_ambas_escalas_de_radio(tmp_path, scale):
    out = tmp_path / f"{scale}.png"
    _, _, paths = draw_geometry(
        bodies=[Body(label="cerca", R=200, theta_deg=0),
                Body(label="lejos", R=2400, theta_deg=180)],
        output_path=out, radius_scale=scale)
    assert _nonempty(paths[0])


def test_proyeccion_oblicua_admite_phi(tmp_path):
    out = tmp_path / "obl.png"
    _, _, paths = draw_geometry(
        bodies=[Body(label="A", R=900, theta_deg=90, phi_deg=90),
                Body(label="B", R=900, theta_deg=90, phi_deg=0),
                Body(kind="polar_molecule", label="C", R=1200, theta_deg=0)],
        projection="oblique", output_path=out)
    assert _nonempty(paths[0])


def test_dibuja_sobre_un_eje_dado_sin_guardar():
    """Con `ax` externo compone paneles y no toca el disco."""
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.5))
    for a, R in zip(axes, (600, 1400)):
        _, got, paths = draw_geometry([Body(label="Rb", R=R, theta_deg=180)],
                                      ax=a, output_path="no/deberia/escribir.png")
        assert got is a
        assert paths == []
    assert not Path("no").exists()


# --------------------------------------------------------- config dict
def test_config_por_tipo(tmp_path):
    out = tmp_path / "cfg.png"
    _, _, paths = draw_from_system_config(
        {"neutro": {"R": 600, "theta_deg": 180, "label": "Rb"},
         "polar": {"R": 900, "theta_deg": 0, "label": "RbCs"}},
        output_path=out)
    assert _nonempty(paths[0])


def test_config_atajo_R1_R2_del_hibrido(tmp_path):
    """R1 -> neutro en θ=180°, R2 -> polar en θ=0°, sin traducir nada."""
    out = tmp_path / "h.png"
    fig, ax, paths = draw_from_system_config({"R1": 600, "R2": 900},
                                             output_path=out)
    assert _nonempty(paths[0])
    textos = {t.get_text() for t in ax.texts}
    assert {"Rb", "RbCs"} <= textos


def test_config_angulos_pisados(tmp_path):
    _, _, paths = draw_from_system_config(
        {"R1": 600, "R2": 900}, angles={"R1": 140, "R2": 20},
        output_path=tmp_path / "ang.png")
    assert _nonempty(paths[0])


# --------------------------------------------------------------- errores
def test_error_si_no_hay_nada_que_dibujar():
    with pytest.raises(ValueError, match="vacío"):
        draw_geometry(bodies=[])


def test_error_theta_fuera_de_rango():
    with pytest.raises(ValueError, match=r"\[0, 180\]"):
        Body(label="A", R=1, theta_deg=200)


def test_error_R_negativo():
    with pytest.raises(ValueError, match="R debe ser"):
        Body(label="A", R=-5)


def test_error_proyeccion_desconocida(tmp_path):
    with pytest.raises(ValueError, match="projection"):
        draw_geometry([Body(label="A", R=1, theta_deg=0)], projection="3d",
                      output_path=tmp_path / "x.png")


def test_error_cuerpo_invisible_en_xz_sugiere_oblique(tmp_path):
    with pytest.raises(ValueError, match="oblique"):
        draw_geometry([Body(label="A", R=900, theta_deg=90, phi_deg=90)],
                      output_path=tmp_path / "x.png")


def test_error_numero_suelto_con_clave_no_conocida(tmp_path):
    with pytest.raises(ValueError, match="abreviatura"):
        draw_from_system_config({"pepito": 900}, output_path=tmp_path / "x.png")


# ------------------------------------------------------------------ CLI
def test_cli_produce_png_y_pdf(tmp_path):
    out = tmp_path / "cli.png"
    env_src = str(REPO / "src")
    cmd = [sys.executable, str(REPO / "scripts" / "draw_geometry.py"),
           "--body", "kind=neutral_atom,label=Rb,R=600,theta=180",
           "--body", "kind=polar_molecule,label=RbCs,R=900,theta=0,dipole=0",
           "--electron", "r=1200,theta=35",
           "--title", "CLI", "-o", str(out), "--formats", "png", "pdf"]
    import os
    env = dict(os.environ)
    env["PYTHONPATH"] = env_src + os.pathsep + env.get("PYTHONPATH", "")
    res = subprocess.run(cmd, capture_output=True, text=True, env=env)
    assert res.returncode == 0, res.stderr
    assert _nonempty(out) and _nonempty(out.with_suffix(".pdf"))


def test_cli_atajo_hybrid(tmp_path):
    import os
    out = tmp_path / "h.png"
    env = dict(os.environ)
    env["PYTHONPATH"] = str(REPO / "src") + os.pathsep + env.get("PYTHONPATH", "")
    res = subprocess.run(
        [sys.executable, str(REPO / "scripts" / "draw_geometry.py"),
         "--hybrid", "600", "900", "-o", str(out)],
        capture_output=True, text=True, env=env)
    assert res.returncode == 0, res.stderr
    assert _nonempty(out)


def test_cli_sin_cuerpos_falla_con_mensaje(tmp_path):
    import os
    env = dict(os.environ)
    env["PYTHONPATH"] = str(REPO / "src") + os.pathsep + env.get("PYTHONPATH", "")
    res = subprocess.run(
        [sys.executable, str(REPO / "scripts" / "draw_geometry.py"),
         "-o", str(tmp_path / "x.png")],
        capture_output=True, text=True, env=env)
    assert res.returncode != 0
    assert "nada que dibujar" in res.stderr
