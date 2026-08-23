"""
Esquemas de geometría molecular: ion Rydberg en el origen + N perturbadores.

Dibuja figuras del estilo de las de Aguilera-Fernández et al. (moléculas de
Rydberg tri-, tetra- y penta-atómicas): el core iónico en el origen, el eje Z
vertical, uno o varios cuerpos (átomos neutros o moléculas polares) a distinta
distancia y ángulo, y opcionalmente el electrón Rydberg con su vector r⃗.

    from trimero.visualization.geometry_diagram import (
        Body, Electron, draw_geometry)

    draw_geometry(
        bodies=[
            Body(kind="polar_molecule", label="RbCs", R=900, theta_deg=0),
            Body(kind="neutral_atom",   label="Rb",   R=600, theta_deg=180),
        ],
        electron=Electron(r=1200, theta_deg=35),
        ion_label=r"Rb$^+$",
        output_path="plots/geometry/esquema.png",
        title="Híbrido neutro-polar",
    )

SIN FÍSICA
----------
El módulo **no importa nada del resto del repositorio**: sólo numpy y
matplotlib. Es deliberado — así sirve para cualquier proyecto y no puede
romperse al tocar un Hamiltoniano.

NO ESTÁ A ESCALA
----------------
`R` (y `Electron.r`) están en a₀ pero **sólo fijan proporciones relativas
entre cuerpos dentro del dibujo**, nunca un tamaño físico: los radios de las
esferas son constantes de estilo, no radios atómicos. Es la misma convención
que las figuras publicadas, que tampoco están a escala. Dos modos:

* `radius_scale="proportional"` (por defecto): longitud de flecha ∝ R,
  normalizada al mayor. Fiel a las razones entre distancias. Con razones
  extremas (R_min/R_max ≲ 0.15) el cuerpo más cercano queda pegado al ion.
* `radius_scale="relative"`: mapea [R_min, R_max] a [0.55, 1.0] del radio de
  dibujo. Preserva el ORDEN pero no las razones; útil cuando conviven un
  cuerpo a 200 a₀ y otro a 2400 a₀.

PROYECCIÓN
----------
La figura es plana. `theta_deg` es el ángulo polar respecto a **+Z** (0-180) y
`phi_deg` el azimutal en el plano XY (por defecto 0, sistema colineal/coplanar
en XZ). Dos proyecciones:

* `projection="xz"` (por defecto): pantalla = (X, Z), con
  X = R·sinθ·cosφ, Z = R·cosθ. φ entra sólo como escorzo. Es lo que dibujan
  las figuras del paper, que son coplanares. Ojo: un cuerpo con φ=90° se
  proyecta sobre el eje Z y puede solaparse con otro.
* `projection="oblique"`: añade el eje Y como diagonal recedente a 30°, de
  modo que φ ≠ 0 se ve de verdad. Para configuraciones no coplanares.

Ningún modo intenta ser una perspectiva realista: son esquemas.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Mapping, Optional, Sequence, Tuple, Union

import numpy as np

__all__ = [
    "Body",
    "Electron",
    "KIND_STYLE",
    "KIND_ALIASES",
    "draw_geometry",
    "draw_from_system_config",
]

# --------------------------------------------------------------- estilo
# Un color y un tamaño relativo por tipo de cuerpo. `Body.color` los pisa.
KIND_STYLE: Dict[str, Dict[str, object]] = {
    "neutral_atom":  {"color": "#3d6fb4", "size": 1.00},   # azul
    "polar_molecule": {"color": "#c2543a", "size": 1.15},  # teja
    "ion":           {"color": "#c9a227", "size": 1.20},   # oro
    "electron":      {"color": "#2b2b2b", "size": 0.42},   # gris muy oscuro
}

# Nombres admitidos en `draw_from_system_config`.
KIND_ALIASES: Dict[str, str] = {
    "neutro": "neutral_atom",
    "neutral": "neutral_atom",
    "neutral_atom": "neutral_atom",
    "atomo": "neutral_atom",
    "atom": "neutral_atom",
    "polar": "polar_molecule",
    "polar_molecule": "polar_molecule",
    "molecula": "polar_molecule",
    "molecule": "polar_molecule",
    "dipolo": "polar_molecule",
}

ARROW_COLOR = "#2e8b3d"      # verde de las flechas R⃗
ELECTRON_ARROW_COLOR = "#5a3d8c"
AXIS_COLOR = "#3a3a3a"


# --------------------------------------------------------------- datos
@dataclass
class Body:
    """
    Un perturbador del esquema.

    Args:
        kind: "neutral_atom" o "polar_molecule". Fija color y tamaño por
            defecto; cualquier otra cadena se acepta y cae al estilo de
            `neutral_atom` (para tipos nuevos, pásale `color`).
        label: texto libre junto a la esfera ("Rb", "RbCs", "KRb"...).
            Admite mathtext de matplotlib: r"Rb$^+$".
        R: distancia al origen en a₀. **Sólo proporciones relativas.**
        theta_deg: ángulo polar respecto a +Z, en [0, 180]. Cualquier valor,
            no sólo 0 y 180: la geometría no colineal ya está soportada.
        phi_deg: ángulo azimutal en el plano XY. 0 = plano XZ.
        color: color matplotlib; None → el del `kind`.
        size: factor multiplicativo sobre el radio de la esfera.
        distance_label: texto de la etiqueta de distancia junto a la flecha.
            None → se autogenera (R₁, R₂, …) en `draw_geometry`.
        dipole_deg: si no es None, dibuja una flechita d⃗ sobre el cuerpo con
            ese ángulo polar (convención de `ChargeDipoleHamiltonian`:
            θ_d = 0 es paralela a +Z). Extra opcional, para moléculas polares.
    """

    kind: str = "neutral_atom"
    label: str = ""
    R: float = 1.0
    theta_deg: float = 0.0
    phi_deg: float = 0.0
    color: Optional[str] = None
    size: float = 1.0
    distance_label: Optional[str] = None
    dipole_deg: Optional[float] = None

    def __post_init__(self):
        if self.R < 0:
            raise ValueError(f"R debe ser >= 0, no {self.R}")
        if not (0.0 <= self.theta_deg <= 180.0):
            raise ValueError(
                f"theta_deg debe estar en [0, 180] (ángulo POLAR respecto a "
                f"+Z), no {self.theta_deg}. Para el otro lado usa phi_deg=180."
            )


@dataclass
class Electron:
    """
    El electrón Rydberg, con su propio vector posición r⃗ desde el origen.

    Mismo esquema angular que `Body`. Se dibuja como punto pequeño; la flecha
    va en otro color para distinguirla de las R⃗ de los perturbadores.
    """

    r: float = 1.0
    theta_deg: float = 45.0
    phi_deg: float = 0.0
    label: str = r"e$^-$"
    distance_label: str = r"$\vec{r}$"
    color: Optional[str] = None
    size: float = 1.0

    def __post_init__(self):
        if self.r < 0:
            raise ValueError(f"r debe ser >= 0, no {self.r}")
        if not (0.0 <= self.theta_deg <= 180.0):
            raise ValueError(
                f"theta_deg debe estar en [0, 180], no {self.theta_deg}")


# ---------------------------------------------------------- geometría
def _unit_xyz(theta_deg: float, phi_deg: float) -> Tuple[float, float, float]:
    """Vector unitario cartesiano a partir de (θ, φ) en grados."""
    t, p = math.radians(theta_deg), math.radians(phi_deg)
    return (math.sin(t) * math.cos(p), math.sin(t) * math.sin(p), math.cos(t))


def _project(xyz: Tuple[float, float, float], projection: str) -> Tuple[float, float]:
    """(x, y, z) 3D -> (horizontal, vertical) de la figura."""
    x, y, z = xyz
    if projection == "xz":
        return x, z
    if projection == "oblique":
        # Eje Y recedente a 30°, escorzado a 0.5: esquema, no perspectiva.
        k = 0.5
        a = math.radians(30.0)
        return x + k * y * math.cos(a), z + k * y * math.sin(a)
    raise ValueError(f"projection debe ser 'xz' u 'oblique', no {projection!r}")


def _plot_radii(values: Sequence[float], radius_scale: str,
                span: float) -> np.ndarray:
    """
    Distancias de dibujo (en unidades de datos) para una lista de R.

    Ver la cabecera del módulo: `proportional` conserva las razones,
    `relative` sólo el orden.
    """
    v = np.asarray(values, dtype=float)
    if v.size == 0:
        return v
    vmax = float(v.max())
    if radius_scale == "proportional":
        if vmax <= 0.0:
            return np.full_like(v, 0.6 * span)
        return v / vmax * span
    if radius_scale == "relative":
        vmin = float(v.min())
        if v.size == 1 or math.isclose(vmax, vmin):
            return np.full_like(v, 0.8 * span)
        lo, hi = 0.55 * span, span
        return lo + (v - vmin) / (vmax - vmin) * (hi - lo)
    raise ValueError(
        f"radius_scale debe ser 'proportional' o 'relative', no {radius_scale!r}")


# ------------------------------------------------------------ dibujo
def _sphere(ax, x: float, y: float, radius: float, color: str,
            zorder: float = 3.0, n: int = 160) -> None:
    """
    Esfera con sombreado radial: un gradiente sencillo, no un render.

    Se construye una imagen RGBA cuadrada con el foco de luz arriba-izquierda
    y se recorta con un círculo. Es la forma más corta de conseguir el aspecto
    de las figuras de referencia sin dependencias extra.
    """
    from matplotlib.colors import to_rgb
    from matplotlib.patches import Circle

    u = np.linspace(-1.0, 1.0, n)
    X, Y = np.meshgrid(u, u)
    inside = X * X + Y * Y <= 1.0

    lx, ly = -0.42, 0.42                       # posición del foco
    d = np.sqrt((X - lx) ** 2 + (Y - ly) ** 2)

    base = np.asarray(to_rgb(color), dtype=float)
    shade = np.clip(1.18 - 0.62 * d, 0.16, 1.0)          # oscurece al alejarse
    rgb = base[None, None, :] * shade[..., None]
    highlight = np.clip(1.0 - d / 0.62, 0.0, 1.0) ** 2   # brillo especular
    rgb = rgb + (1.0 - rgb) * highlight[..., None] * 0.80

    img = np.dstack([np.clip(rgb, 0.0, 1.0), inside.astype(float)])
    im = ax.imshow(img, extent=(x - radius, x + radius, y - radius, y + radius),
                   origin="lower", zorder=zorder, interpolation="bilinear")
    im.set_clip_path(Circle((x, y), radius, transform=ax.transData))
    # Borde tenue: separa la esfera del fondo cuando el color es claro.
    ax.add_patch(Circle((x, y), radius, fill=False, lw=0.7,
                        edgecolor=(0, 0, 0, 0.35), zorder=zorder + 0.1))


def _arrow(ax, x0: float, y0: float, x1: float, y1: float, color: str,
           lw: float = 1.6, zorder: float = 2.0, ls: str = "-") -> None:
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0), zorder=zorder,
                arrowprops=dict(arrowstyle="-|>", color=color, lw=lw,
                                linestyle=ls, shrinkA=0, shrinkB=0,
                                mutation_scale=15))


def _axes_cross(ax, length: float, projection: str, fs: float) -> None:
    """Ejes Z (vertical) y X (horizontal) con flecha; Y sólo en oblicua."""
    # Tramos negativos como línea fina: los cuerpos en θ=180 caen ahí.
    ax.plot([-length, length], [0, 0], color=AXIS_COLOR, lw=0.7,
            alpha=0.45, zorder=1)
    ax.plot([0, 0], [-length, length], color=AXIS_COLOR, lw=0.7,
            alpha=0.45, zorder=1)
    _arrow(ax, 0, 0, 0, length, AXIS_COLOR, lw=1.3, zorder=1.5)
    _arrow(ax, 0, 0, length, 0, AXIS_COLOR, lw=1.3, zorder=1.5)
    ax.text(0.035 * length, length, "Z", fontsize=fs, ha="left", va="top",
            color=AXIS_COLOR)
    ax.text(length, 0.045 * length, "X", fontsize=fs, ha="right", va="bottom",
            color=AXIS_COLOR)
    if projection == "oblique":
        yx, yy = _project((0.0, length, 0.0), projection)
        _arrow(ax, 0, 0, yx, yy, AXIS_COLOR, lw=1.3, zorder=1.5)
        ax.text(yx, yy, " Y", fontsize=fs, ha="left", va="bottom",
                color=AXIS_COLOR)


def _place_label(ax, x: float, y: float, ux: float, uy: float,
                 offset: float, text: str, fontsize: float,
                 color: str = "black") -> None:
    """Etiqueta desplazada del centro en la dirección radial (ux, uy)."""
    if not text:
        return
    norm = math.hypot(ux, uy) or 1.0
    ux, uy = ux / norm, uy / norm
    ha = "left" if ux > 0.15 else ("right" if ux < -0.15 else "center")
    va = "bottom" if uy > 0.15 else ("top" if uy < -0.15 else "center")
    ax.text(x + ux * offset, y + uy * offset, text, fontsize=fontsize,
            ha=ha, va=va, color=color, zorder=6)


def _rotate(u: Tuple[float, float], deg: float) -> Tuple[float, float]:
    a = math.radians(deg)
    c, si = math.cos(a), math.sin(a)
    return (u[0] * c - u[1] * si, u[0] * si + u[1] * c)


def _label_direction(cx: float, cy: float, u: Tuple[float, float],
                     offset: float, obstacles: Sequence[Tuple[float, float]],
                     tilt: float = 42.0) -> Tuple[float, float]:
    """
    Dirección en la que colocar la etiqueta de un cuerpo.

    No se pone radialmente: en una configuración colineal (θ=0 o 180) eso deja
    el texto encima de la propia flecha y del d⃗. Se inclina `tilt` grados a un
    lado o al otro y se elige el lado cuyo punto de anclaje queda MÁS LEJOS de
    los obstáculos —los demás cuerpos, el ion, el electrón, las etiquetas de
    los ejes—. Se mide en distancia real y no en ángulo: dos cuerpos pueden
    apuntar casi igual y estar lejísimos, y al revés.
    """
    best, best_score = None, -math.inf
    for deg in (+tilt, -tilt):
        d = _rotate(u, deg)
        px, py = cx + d[0] * offset, cy + d[1] * offset
        score = min((math.hypot(px - ox, py - oy) for ox, oy in obstacles),
                    default=math.inf)
        if score > best_score:
            best, best_score = d, score
    return best


def _free_diagonal(directions: Sequence[Tuple[float, float]]
                   ) -> Tuple[float, float]:
    """La diagonal (±1,±1) más alejada de todas las direcciones ocupadas."""
    diagonals = [(-0.7071, -0.7071), (0.7071, -0.7071),
                 (-0.7071, 0.7071), (0.7071, 0.7071)]
    unit = []
    for dx, dy in directions:
        n = math.hypot(dx, dy)
        if n > 1e-9:
            unit.append((dx / n, dy / n))
    if not unit:
        return diagonals[0]
    # Mejor diagonal = la que maximiza la distancia angular mínima.
    return max(diagonals,
               key=lambda d: min(1.0 - (d[0] * u[0] + d[1] * u[1]) for u in unit))


# ----------------------------------------------------------- API
def draw_geometry(
    bodies: Iterable[Body],
    electron: Optional[Electron] = None,
    *,
    ion_label: str = r"Rb$^+$",
    ion_color: Optional[str] = None,
    title: Optional[str] = None,
    output_path: Optional[Union[str, Path]] = None,
    formats: Sequence[str] = ("png",),
    projection: str = "xz",
    radius_scale: str = "proportional",
    figsize: Tuple[float, float] = (6.4, 6.4),
    dpi: int = 200,
    body_size: float = 0.085,
    fontsize: float = 12.0,
    show_distance_labels: bool = True,
    ax=None,
):
    """
    Dibuja el esquema y, si se le da `output_path`, lo guarda.

    Args:
        bodies: los perturbadores. Número arbitrario, ángulos arbitrarios.
        electron: opcional; si se da, se dibuja con su flecha r⃗.
        ion_label / ion_color: la esfera del origen.
        title: título de la figura (None = sin título).
        output_path: ruta de salida. La extensión se IGNORA: los ficheros se
            escriben como `<raíz>.<fmt>` para cada fmt de `formats`. Los
            directorios intermedios se crean.
        formats: ("png",) por defecto; p. ej. ("png", "pdf").
        projection / radius_scale: ver la cabecera del módulo.
        body_size: radio de la esfera de un cuerpo, en unidades de datos
            (el radio de dibujo máximo es 1.0).
        show_distance_labels: etiquetas R₁, R₂, … junto a cada flecha.
        ax: si se pasa un eje de matplotlib, se dibuja ahí y NO se guarda ni
            se cierra la figura. Útil para componer paneles.

    Returns:
        (fig, ax, paths) — `paths` es la lista de ficheros escritos (vacía si
        no se pidió `output_path`).

    Raises:
        ValueError: si `bodies` está vacío y no hay electrón, o si
            `projection` / `radius_scale` / `formats` no son válidos.
    """
    # No se toca `matplotlib.use`: una librería no debe secuestrar el backend
    # del que la llama. `scripts/draw_geometry.py` fija Agg por su cuenta.
    import matplotlib.pyplot as plt

    bodies = list(bodies)
    if not bodies and electron is None:
        raise ValueError("nada que dibujar: `bodies` está vacío y no hay electrón")
    bad = [f for f in formats if not str(f).strip()]
    if bad or not list(formats):
        raise ValueError(f"`formats` inválido: {formats!r}")

    span = 1.0
    # El electrón entra en la MISMA normalización que los cuerpos: si no, su
    # flecha no sería comparable con las R⃗ y el esquema mentiría.
    all_R = [b.R for b in bodies] + ([electron.r] if electron else [])
    radii = _plot_radii(all_R, radius_scale, span)
    body_radii = radii[:len(bodies)]
    electron_radius = float(radii[-1]) if electron else None

    own_fig = ax is None
    if own_fig:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.figure

    ax.set_aspect("equal")
    ax.axis("off")

    axis_len = 1.30 * span
    _axes_cross(ax, axis_len, projection, fontsize)

    ion_r = body_size * float(KIND_STYLE["ion"]["size"])
    _sphere(ax, 0.0, 0.0, ion_r,
            ion_color or str(KIND_STYLE["ion"]["color"]), zorder=5)
    # La etiqueta del ion va al cuadrante diagonal más libre: con cuerpos en
    # θ=0 y θ=180 la esquina inferior-izquierda fija chocaría con una flecha.
    projected = [_project(_unit_xyz(b.theta_deg, b.phi_deg), projection)
                 for b in bodies]
    occupied = list(projected)
    if electron is not None:
        occupied.append(_project(
            _unit_xyz(electron.theta_deg, electron.phi_deg), projection))
    ix, iy = _free_diagonal(occupied)
    _place_label(ax, 0.0, 0.0, ix, iy, ion_r * 1.5, ion_label, fontsize)

    # --- pasada 1: dónde va cada cosa ------------------------------------
    placements = []
    for i, (b, rad) in enumerate(zip(bodies, body_radii), start=1):
        style = KIND_STYLE.get(b.kind, KIND_STYLE["neutral_atom"])
        sphere_r = body_size * float(style["size"]) * b.size

        px, py = projected[i - 1]
        norm = math.hypot(px, py)
        if norm < 1e-9:
            # (θ=90°, φ=±90°) apunta a ±Y: en la proyección XZ cae encima del
            # ion y el esquema mentiría. Se dice, no se coloca en cualquier sitio.
            raise ValueError(
                f"el cuerpo {i} ({b.label or b.kind!r}, θ={b.theta_deg}°, "
                f"φ={b.phi_deg}°) apunta fuera del plano de la proyección "
                f"{projection!r} y se proyectaría sobre el origen. "
                "Usa projection='oblique'."
            )
        ux, uy = px / norm, py / norm
        cx, cy = px * rad, py * rad          # el escorzo acorta la flecha

        dipole_dir = None
        if b.dipole_deg is not None:
            dxp, dyp = _project(_unit_xyz(b.dipole_deg, b.phi_deg), projection)
            dn = math.hypot(dxp, dyp)
            if dn > 1e-9:
                dipole_dir = (dxp / dn, dyp / dn)

        placements.append(dict(
            body=b, index=i, cx=cx, cy=cy, ux=ux, uy=uy,
            sphere_r=sphere_r, color=b.color or str(style["color"]),
            dipole_dir=dipole_dir))

    e_place = None
    if electron is not None:
        e_style = KIND_STYLE["electron"]
        e_r = body_size * float(e_style["size"]) * electron.size
        px, py = _project(_unit_xyz(electron.theta_deg, electron.phi_deg),
                          projection)
        norm = math.hypot(px, py)
        if norm < 1e-9:
            raise ValueError(
                f"el electrón (θ={electron.theta_deg}°, φ={electron.phi_deg}°) "
                f"se proyectaría sobre el origen en projection={projection!r}. "
                "Usa projection='oblique'."
            )
        e_place = dict(cx=px * electron_radius, cy=py * electron_radius,
                       ux=px / norm, uy=py / norm, sphere_r=e_r,
                       color=electron.color or str(e_style["color"]))

    # Puntos que ya llevan tinta y con los que una etiqueta no debe chocar.
    obstacles = [(0.0, 0.0)]                                   # el ion
    obstacles += [(p["cx"], p["cy"]) for p in placements]
    if e_place is not None:
        obstacles.append((e_place["cx"], e_place["cy"]))
    obstacles += [(0.035 * axis_len, axis_len),                # etiqueta «Z»
                  (axis_len, 0.045 * axis_len)]                # etiqueta «X»
    for p in placements:
        if p["dipole_dir"] is not None:
            dx, dy = p["dipole_dir"]
            obstacles.append((p["cx"] + dx * p["sphere_r"] * 1.5,
                              p["cy"] + dy * p["sphere_r"] * 1.5))
        if show_distance_labels:
            # El anclaje de la etiqueta R_i no depende de nada que se decida
            # después, así que se puede meter ya como obstáculo.
            t = max(math.hypot(p["cx"], p["cy"]) - p["sphere_r"], 0.0) * 0.55
            obstacles.append((p["ux"] * t - p["uy"] * 0.075,
                              p["uy"] * t + p["ux"] * 0.075))

    # --- pasada 2: dibujar -----------------------------------------------
    for p in placements:
        b, cx, cy = p["body"], p["cx"], p["cy"]
        ux, uy, sphere_r = p["ux"], p["uy"], p["sphere_r"]

        tip = max(math.hypot(cx, cy) - sphere_r, 0.0)
        _arrow(ax, 0.0, 0.0, ux * tip, uy * tip, ARROW_COLOR)
        _sphere(ax, cx, cy, sphere_r, p["color"], zorder=4)

        others = [o for o in obstacles
                  if math.hypot(o[0] - cx, o[1] - cy) > 1e-9]
        off = sphere_r * 1.75
        lx, ly = _label_direction(cx, cy, (ux, uy), off, others)
        _place_label(ax, cx, cy, lx, ly, off, b.label, fontsize)

        if show_distance_labels:
            text = b.distance_label
            if text is None:
                text = rf"$R_{{{p['index']}}}$" if len(bodies) > 1 else r"$R$"
            mx, my = ux * tip * 0.55, uy * tip * 0.55
            _place_label(ax, mx, my, -uy, ux, 0.075, text,
                         fontsize, color=ARROW_COLOR)

        if p["dipole_dir"] is not None:
            dx, dy = p["dipole_dir"]
            L = sphere_r
            _arrow(ax, cx - dx * L, cy - dy * L, cx + dx * L, cy + dy * L,
                   "#111111", lw=1.2, zorder=5.5)
            # El d⃗ se rotula de lado, no en la punta: la punta suele caer
            # sobre el eje Z (θ_d = 0 es el caso típico) o sobre la flecha R⃗.
            # Y se manda al lado contrario al de la etiqueta del cuerpo.
            perp = (-dy, dx)
            if perp[0] * lx + perp[1] * ly > 0.0:
                perp = (dy, -dx)
            _place_label(ax, cx + dx * L, cy + dy * L, perp[0], perp[1],
                         0.055, r"$\vec{d}$", fontsize * 0.85)

    if e_place is not None:
        cx, cy = e_place["cx"], e_place["cy"]
        ux, uy, e_r = e_place["ux"], e_place["uy"], e_place["sphere_r"]
        tip = max(math.hypot(cx, cy) - e_r, 0.0)
        _arrow(ax, 0.0, 0.0, ux * tip, uy * tip, ELECTRON_ARROW_COLOR,
               lw=1.4, ls="--")
        _sphere(ax, cx, cy, e_r, e_place["color"], zorder=4)
        others = [o for o in obstacles
                  if math.hypot(o[0] - cx, o[1] - cy) > 1e-9]
        off = e_r * 2.8
        lx, ly = _label_direction(cx, cy, (ux, uy), off, others)
        _place_label(ax, cx, cy, lx, ly, off, electron.label, fontsize)
        if show_distance_labels:
            _place_label(ax, ux * tip * 0.55, uy * tip * 0.55, -uy, ux, 0.075,
                         electron.distance_label, fontsize,
                         color=ELECTRON_ARROW_COLOR)

    lim = axis_len * 1.16
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    if title:
        ax.set_title(title, fontsize=fontsize * 1.1, pad=10)

    paths = []
    if output_path is not None and own_fig:
        stem = Path(output_path).with_suffix("")
        stem.parent.mkdir(parents=True, exist_ok=True)
        for fmt in formats:
            out = stem.with_suffix(f".{str(fmt).lstrip('.')}")
            fig.savefig(out, dpi=dpi, bbox_inches="tight",
                        facecolor="white")
            paths.append(out)
    return fig, ax, paths


def draw_from_system_config(
    config: Mapping[str, Union[Mapping[str, object], float, int]],
    **kwargs,
):
    """
    `draw_geometry` a partir de un diccionario, sin construir `Body` a mano.

    Dos formas admitidas, mezclables:

    1. **Por tipo**, valor = diccionario con los campos de `Body`::

           draw_from_system_config({
               "neutro": {"R": 600, "theta_deg": 180, "label": "Rb"},
               "polar":  {"R": 900, "theta_deg": 0,   "label": "RbCs"},
           })

       La clave fija el `kind` (ver `KIND_ALIASES`) y, si no se da `label`,
       también la etiqueta. Una clave desconocida se acepta: se dibuja como
       átomo neutro y su nombre se usa de etiqueta.

    2. **Abreviatura del híbrido**, valor numérico::

           draw_from_system_config({"R1": 600, "R2": 900})

       `R1` → perturbador NEUTRO en θ=180°, `R2` → molécula POLAR en θ=0°,
       que es exactamente el convenio de `HybridNeutralPolar.hamiltonian(R1,
       R2, ...)`. Así se llama con los mismos parámetros del sistema sin
       traducir nada. Los ángulos se pueden pisar con
       `angles={"R1": 180, "R2": 0}`.

    El resto de argumentos (`electron`, `title`, `output_path`, `formats`,
    `projection`, …) se pasan tal cual a `draw_geometry`.

    Raises:
        ValueError: si `config` está vacío o una entrada no es ni diccionario
            ni número.
    """
    angles = kwargs.pop("angles", None) or {}
    if not config:
        raise ValueError("`config` está vacío: no hay ningún cuerpo que dibujar")

    # Convenio del híbrido: R1 = neutro en θ=π, R2 = polar en θ=0.
    shorthand = {
        "r1": ("neutral_atom", 180.0, "Rb"),
        "r2": ("polar_molecule", 0.0, "RbCs"),
    }

    bodies = []
    for key, value in config.items():
        low = str(key).lower()
        if isinstance(value, Mapping):
            spec = dict(value)
            kind = str(spec.pop("kind", KIND_ALIASES.get(low, low)))
            spec.setdefault("label", str(key))
            spec.setdefault("theta_deg", angles.get(key, 0.0))
            bodies.append(Body(kind=KIND_ALIASES.get(kind, kind), **spec))
        elif isinstance(value, (int, float)) and not isinstance(value, bool):
            if low not in shorthand:
                raise ValueError(
                    f"la clave {key!r} lleva un número ({value}) pero no es una "
                    f"abreviatura conocida {sorted(shorthand)}. Usa un "
                    f"diccionario: {{{key!r}: {{'R': {value}, 'theta_deg': 0}}}}"
                )
            kind, theta, label = shorthand[low]
            bodies.append(Body(kind=kind, label=label, R=float(value),
                               theta_deg=float(angles.get(key, theta))))
        else:
            raise ValueError(
                f"la entrada {key!r} debe ser un diccionario de campos de "
                f"`Body` o un número; es {type(value).__name__}"
            )
    return draw_geometry(bodies, **kwargs)
