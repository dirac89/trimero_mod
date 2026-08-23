# Fase A — curvas BOP de Rb*-KRb comparadas entre varios n

**Fecha**: 2026-08-23 (título corregido el mismo día: "Rb*-RbCs" → "Rb*-KRb";
ver nota de corrección en `docs/PLAN_figuras_publicacion.md`. No cambia
ningún número de este documento, sólo la etiqueta del sistema — `rb_krb_polar`
corre con constantes de KRb, `B_KRB_GHZ=1.114`, `D_KRB_DEBYE=0.566`.)
**Autor**: Javier Aguilera (via Claude Code)
**Relevancia**: primera fase de `docs/PLAN_figuras_publicacion.md` (figuras de
publicación). Genera la comparación entre manifolds que sirve de base a las
Fases B-E (orientación, campo eléctrico, híbrido, conclusión Franck-Condon).
**Tipo**: analysis

## Resumen

Se generalizó `scripts/compute_bop_curve.py` (sin tocar su código: ya acepta
`--n-manifold` arbitrario) a n=24, 25, 26, 27 en `M_J=0`, mismo rango
R ∈ [400, 1800] a₀, paso 5 a₀, que la Fig. 1 ya validada. n=25 se reutilizó
del `.npz` existente sin recalcular. Se construyó un script nuevo,
`scripts/compare_bop_curves_n.py`, que superpone las cuatro curvas en escala
absoluta y en escala reducida R/2n². Los cuatro barridos son físicamente
sanos (sin NaN, peso de manifold ≥0.92 en todo el rango) y confirman la
tendencia de escala esperada: la estructura de la curva (el "codo" hacia la
asíntota) se desplaza a mayor R y el pozo más profundo dentro de la ventana
[400,1800] se hace menos profundo, al crecer n.

## Qué se corrió

Línea base (antes de cualquier cambio):

```
poetry run pytest -m "not slow"
133 passed, 16 deselected in 154.12s (0:02:34)
```

Coste computacional por n (dim del bloque M_J=0, tiempo post-warmup):

| n | dim(M_J=0) | s/punto | estimado barrido completo (281 pts) |
|---|---|---|---|
| 24 | 1064 | 1.00 | 4.7 min |
| 25 | 1113 | 1.16 | 5.4 min |
| 26 | 1162 | 1.49 | 7.0 min |
| 27 | 1211 | 1.64 | 7.7 min |

`domain_bounds()` (dominio del remapeo semiclásico de Fermi): 1090/1185/1285/
1388 a₀ para n=24/25/26/27. **No es una restricción real aquí**: el sistema
polar corre con `fermi=False` (`docs/STATUS.md`, sección "Rb*-KRb: el modelo
vigente"), y `H_A+H_mol` está definido en todo R con la función de onda
hidrogenoide — por eso el barrido cubre [400,1800] a₀ entero para los cuatro
n sin recortes, aunque ese rango exceda `domain_bounds()` en todos los casos.
Se comprobó explícitamente para no confundirlo con un límite físico.

Barrido de producción (`scripts/compute_bop_curve.py --n-manifold {n} --mj 0
--rmin 400 --rmax 1800 --step 5 --no-plot`, n=24,26,27; n=25 reutilizado):

```
=== n=24 ===
  M_J=0: dim(bloque) = 1064   281 puntos
    281 diagonalizaciones en 298.1 s (1.06 s/punto)
    pozo MÁS PROFUNDO   :  -25.5283 GHz en R = 400.0 a0
    E(R = 1800 a0)      :   -0.2072 GHz   (k = 55, peso = 1.0000)
    mínimos locales     : 7

=== n=25 (reutilizado, sin recalcular) ===
    pozo MÁS PROFUNDO   :  -23.1002 GHz en R = 400.0 a0   (= referencia de STATUS.md)
    E(R = 1800 a0)      :   -0.3376 GHz
    mínimos locales     : 8

=== n=26 ===
  M_J=0: dim(bloque) = 1162   281 puntos
    281 diagonalizaciones en 353.6 s (1.26 s/punto)
    pozo MÁS PROFUNDO   :  -20.9491 GHz en R = 400.0 a0
    E(R = 1800 a0)      :   -0.5742 GHz   (k = 55, peso = 1.0000)
    mínimos locales     : 9

=== n=27 ===
  M_J=0: dim(bloque) = 1211   281 puntos
    281 diagonalizaciones en 388.2 s (1.38 s/punto)
    pozo MÁS PROFUNDO   :  -19.0791 GHz en R = 400.0 a0
    E(R = 1800 a0)      :   -1.0293 GHz   (k = 54, peso = 1.0000)
    mínimos locales     : 10
```

Figura comparativa:

```
poetry run python scripts/compare_bop_curves_n.py --n 24 25 26 27
```

produce `plots/rb_krb_polar/figures/fig_compare_n_MJ0.png` (dos paneles: R absoluto,
y R/2n² reducido) y esta tabla:

```
   n   pozo mas profundo (GHz)   R pozo (a0)   minimos locales    E(R=1800) (GHz)
  24                  -25.5283         400.0                 7            -0.2072
  25                  -23.1002         400.0                 8            -0.3376
  26                  -20.9491         400.0                 9            -0.5742
  27                  -19.0791         400.0                10            -1.0293
```

## Verificación de sanidad numérica

Antes de interpretar nada: los cuatro barridos son limpios.

```
n=24: nan=0  W_min=0.9983  W_max=1.0000  K cambia entre 2 valores
n=25: nan=0  W_min=0.9985  W_max=1.0000  K cambia entre 3 valores
n=26: nan=0  W_min=0.9216  W_max=1.0000  K cambia entre 3 valores
n=27: nan=0  W_min=0.9984  W_max=1.0000  K cambia entre 2 valores
```

Sin NaN en ningún punto, peso de manifold siempre ≥0.92 (muy por encima del
umbral 0.5 de aceptación), y el índice k identificado por carácter cambia 2-3
veces en el barrido — consistente con los cruces evitados con (n+1)d/(n+2)p
que `docs/analysis_base_correcta_3_vecinos.md` §5.1 ya documentó: la curva
sigue el carácter, no un índice fijo, y eso es exactamente lo que se ve.

## Interpretación — la tendencia de escala, y una precisión importante

**Lo que el plan anticipaba ("pozos se desplazan a mayor R con n creciente,
profundidad decrece con n") se confirma, pero con un matiz que hay que dejar
explícito, no ocultar:**

1. **El "pozo más profundo"** tal como lo reporta `report()` (mínimo global
   de la curva en la ventana) está anclado en R=400 a₀ para los cuatro n. Eso
   **no** es el desplazamiento a mayor R que describe la tendencia física: es
   que R=400 a₀ es el borde inferior de una ventana ABSOLUTA fija, heredada
   de la Fig. 1 original, que no escala con n. Dentro de esa ventana, la
   curva para n grande simplemente sigue bajando hacia R=400 sin haber
   alcanzado aún su mínimo real (que estaría en R menor, fuera de rango).
   Reportar solo esa cifra como "el pozo se mueve" sería forzar una lectura
   que los números no sostienen — de ahí este apartado.

2. **La tendencia real y verificable está en el panel reducido de la figura**
   (`R/2n²`): las cuatro curvas COLAPSAN sobre la misma forma en la región
   R/2n² ≈ 0.9-1.1, donde ocurre la transición hacia la asíntota (el "codo").
   En escala absoluta ese codo está en R≈1100-1150 a₀ para n=24 y en
   R≈1450-1500 a₀ para n=27: se desplaza a mayor R con n, tal como predice el
   escalado ~2n² de un manifold de Rydberg. Esto es el criterio de sanidad
   correcto, no el mínimo global de una ventana fija.

3. **La profundidad SÍ decrece monótonamente con n** de forma robusta, tanto
   si se mide en R=400 a₀ fijo (25.53 → 23.10 → 20.95 → 19.08 GHz) como en
   términos de la anchura de la estructura oscilatoria interna (7→8→9→10
   mínimos locales: el manifold tiene más estados vecinos a mayor n, así que
   hay más cruces evitados en la misma ventana absoluta — coherente, no un
   artefacto).

4. **`E(R=1800 a₀)` decrece con n** en valor absoluto creciente en magnitud
   negativa (−0.21 → −0.34 → −0.57 → −1.03 GHz): a n mayor, R=1800 a₀ está
   relativamente MÁS cerca del manifold (porque 2n² crece), así que la curva
   todavía no se ha aplanado del todo en el borde derecho de la ventana fija.
   Coherente con el punto 2.

## Conclusiones y aplicación al proyecto

- El conjunto n=24,25,26,27 es válido y ya está calculado y verificado; sirve
  de base directa para la Fase B (⟨cosθ_d⟩) y para la comparación
  Franck-Condon de la Fase E.
- La figura de esta fase (`plots/rb_krb_polar/figures/fig_compare_n_MJ0.png`) debe
  presentarse con los DOS paneles, no solo el absoluto: el panel reducido es
  el que sostiene la afirmación de escala, el absoluto por sí solo induce a
  leer mal el "pozo más profundo".
- Ventana fija [400,1800] a₀ para los cuatro n es una decisión de estilo
  consistente con la Fig. 1 original, no una limitación del cálculo: el
  Hamiltoniano polar está definido en todo R (§ arriba, `domain_bounds()`).
  Si en una fase posterior se quisiera capturar el mínimo global real para
  cada n, haría falta extender `--rmin` por debajo de 400 a₀ para n>25 —no se
  hizo en esta ronda porque no formaba parte del alcance de la Fase A.

## Referencias

- `docs/STATUS.md` — modelo vigente Rb*-KRb, números de referencia n=25.
- `docs/analysis_base_correcta_3_vecinos.md` — identificación de la curva por
  carácter, no por índice.
- `docs/PLAN_figuras_publicacion.md` — plan maestro, Fase A.
