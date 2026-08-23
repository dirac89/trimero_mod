# Esquemas de geometría molecular (`geometry_diagram`)

**Fecha**: 2026-08-22
**Autor**: Javier Aguilera
**Relevancia**: Herramienta de figuras para manuscritos y presentaciones. **No
es un análisis físico** y no toca ningún resultado: se documenta aquí sólo por
la convención de `docs/`.
**Tipo**: tool

## Resumen

`src/trimero/visualization/geometry_diagram.py` dibuja esquemas del estilo de
los de Aguilera-Fernández et al. (moléculas de Rydberg tri-, tetra- y
penta-atómicas): ion Rydberg en el origen, eje Z vertical, N perturbadores
—átomos neutros o moléculas polares— a distintas distancias y ángulos, y
opcionalmente el electrón con su vector r⃗.

**Es una utilidad aislada.** No importa `basis`, `systems` ni `mathlib`: sólo
numpy y matplotlib. Se puede copiar fuera del repositorio y sigue funcionando.
No hay ninguna dependencia en la otra dirección tampoco — ningún módulo de
física la importa.

## Uso

```python
from trimero.visualization.geometry_diagram import Body, Electron, draw_geometry

draw_geometry(
    bodies=[
        Body(kind="neutral_atom",   label="Rb",   R=900, theta_deg=180),
        Body(kind="polar_molecule", label="RbCs", R=800, theta_deg=0, dipole_deg=0),
    ],
    electron=Electron(r=1225, theta_deg=40),
    ion_label=r"Rb$^+$",
    title="Híbrido neutro-polar",
    output_path="plots/geometry/esquema.png",
    formats=("png", "pdf"),
)
```

Desde consola, con `scripts/draw_geometry.py` (envoltorio fino, sin lógica):

```bash
poetry run python scripts/draw_geometry.py \
    --body kind=neutral_atom,label=Rb,R=900,theta=180 \
    --body kind=polar_molecule,label=RbCs,R=800,theta=0,dipole=0 \
    --electron r=1225,theta=40 \
    --title "Híbrido neutro-polar" \
    -o plots/geometry/hybrid_neutral_polar_esquema.png --formats png pdf
```

Y con el atajo del híbrido, que usa los **mismos parámetros** que
`HybridNeutralPolar.hamiltonian(R1, R2, …)` sin traducir nada — R1 es el
neutro en θ=π y R2 el polar en θ=0:

```python
from trimero.visualization.geometry_diagram import draw_from_system_config
draw_from_system_config({"R1": 900, "R2": 800}, output_path="plots/geometry/h.png")
```
```bash
poetry run python scripts/draw_geometry.py --hybrid 900 800 -o plots/geometry/h.png
```

`draw_from_system_config` también admite la forma por tipo, con cualquier campo
de `Body` dentro:

```python
draw_from_system_config({
    "neutro": {"R": 900, "theta_deg": 180, "label": "Rb"},
    "polar":  {"R": 800, "theta_deg": 0,   "label": "RbCs", "dipole_deg": 0},
})
```

## Tres decisiones que conviene conocer

1. **La figura NO está a escala** y no pretende estarlo, igual que las
   publicadas. `R` sólo fija proporciones **relativas entre cuerpos dentro del
   dibujo**; los radios de las esferas son constantes de estilo, no radios
   atómicos. Dos modos: `radius_scale="proportional"` (por defecto, longitud
   de flecha ∝ R) y `radius_scale="relative"` (mapea [R_min, R_max] a
   [0.55, 1.0] del radio de dibujo, para cuando conviven 200 a₀ y 2400 a₀ y el
   primero quedaría pegado al ion).

2. **Ángulos y número de cuerpos son arbitrarios desde el principio.**
   `theta_deg ∈ [0, 180]` es el ángulo polar respecto a +Z; no hay nada
   cableado a 0/π. `phi_deg` saca el sistema del plano XZ. Es deliberado: la
   generalización a geometrías no colineales (asimétrica, planar) no requerirá
   tocar el dibujo.

3. **Dos proyecciones.** `projection="xz"` (por defecto) proyecta sobre
   (X, Z), que es lo que dibujan las figuras coplanares del paper; φ entra
   sólo como escorzo. Un cuerpo en (θ=90°, φ=±90°) apuntaría a ±Y y caería
   sobre el ion, así que **se lanza `ValueError` sugiriendo
   `projection="oblique"`** en vez de colocarlo en cualquier sitio.
   `projection="oblique"` añade el eje Y como diagonal recedente a 30° y
   dibuja φ de verdad. Ninguna de las dos es una perspectiva realista.

## Detalles de implementación

* **Esferas**: gradiente radial con foco arriba-izquierda, construido como
  imagen RGBA y recortado con un círculo. Es la vía más corta al aspecto de
  las figuras de referencia sin dependencias extra.
* **Anticolisión de etiquetas**: cada etiqueta de cuerpo se inclina ±42°
  respecto a su flecha y se elige el lado cuyo anclaje queda más lejos de los
  obstáculos (los demás cuerpos, el ion, el electrón, las etiquetas «Z» y «X»,
  las etiquetas R_i, la etiqueta d⃗). Se mide en **distancia real**, no en
  ángulo: en la configuración colineal dos cuerpos apuntan casi igual y están
  lejísimos, y el criterio angular fallaba ahí.
* **`dipole_deg`** (opcional) dibuja una flechita d⃗ sobre el cuerpo, con el
  convenio de `ChargeDipoleHamiltonian`: θ_d = 0 es paralela a +Z. Es lo único
  del módulo que guiña a la física, y es sólo una etiqueta.
* **La librería NO llama a `matplotlib.use`**: secuestrar el backend de quien
  la importa sería incorrecto. `scripts/draw_geometry.py` fija `Agg` por su
  cuenta.
* **`ax=`**: si se le pasa un eje, dibuja ahí, no guarda y no cierra la figura.
  Para componer paneles.

## Ejemplo generado

`plots/geometry/hybrid_neutral_polar_esquema.{png,pdf}` — la configuración
vigente del híbrido: Rb neutro en θ=180° a R₁ = 900 a₀, RbCs polar en θ=0° a
R₂ = 800 a₀ (los valores de `tests/systems/hybrid_neutral_polar/test_limites_hibrido.py`),
ion Rb⁺ en el origen y el electrón a r ≈ n² = 1225 a₀.

## Tests

`tests/visualization/test_geometry_diagram.py`, **24 tests, ~6 s**. No
comprueban física ni comparan imágenes: sólo que se produce un fichero no
vacío, que generaliza a 1, 2, 3, 5 y 8 cuerpos con ángulos arbitrarios, que
PNG y PDF salen, y que los errores de entrada dan un mensaje útil (θ fuera de
rango, R negativo, proyección desconocida, cuerpo invisible en XZ). La
comparación visual con la figura de referencia la hace una persona mirando el
PNG, no un `assert`.

## Referencias

- Aguilera-Fernández, J., Schmelcher, P. & González-Férez, R. (2017),
  moléculas de Rydberg penta-atómicas — figura de geometría de referencia.
- [`analysis_trimero_lineal_campo_dc.md`](analysis_trimero_lineal_campo_dc.md)
  — el sistema lineal simétrico que estos esquemas ilustran.
