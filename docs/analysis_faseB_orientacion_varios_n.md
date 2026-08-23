# Fase B — orientación ⟨cosθ_d⟩ de Rb*-KRb comparada entre varios n

**Fecha**: 2026-08-23 (título corregido el mismo día: "Rb*-RbCs" → "Rb*-KRb";
ver nota de corrección en `docs/PLAN_figuras_publicacion.md`. No cambia
ningún número de este documento, sólo la etiqueta del sistema — `rb_krb_polar`
corre con constantes de KRb, `B_KRB_GHZ=1.114`, `D_KRB_DEBYE=0.566`.)
**Autor**: Javier Aguilera (via Claude Code)
**Relevancia**: segunda fase de `docs/PLAN_figuras_publicacion.md`. Generaliza
a n=24,25,26,27 el cálculo de orientación de una ronda anterior de esta
sesión, y corrige una dependencia que había quedado desactualizada.

## Resumen

Se recuperó el código de orientación de
`scripts/archive/run_basis_correction_check.py` (§4.2 de
`docs/analysis_base_correcta_3_vecinos.md`) y se reprodujo bit a bit para
n=24 como verificación de recuperación. Al hacerlo se encontró que ese código
depende de `BOPSystem.hamiltonian()` con su valor por defecto `fermi=True`,
es decir, **incluye el pseudopotencial de Fermi** — la misma premisa
equivocada que `docs/analysis_fig1_carga_dipolo_sin_fermi.md` corrigió para
las curvas BOP, pero que nunca se corrigió para este cálculo de orientación.
Se generalizó a n=24,25,26,27 con `fermi=False` explícito, consistente con
el modelo vigente (`docs/STATUS.md`, y con la Fase A). Los números cambian
de forma no trivial respecto a la ronda anterior — se documenta la magnitud
del cambio, no se oculta.

Se construyó además un diagnóstico cuantitativo de la región de seguimiento
ambiguo (en vez de la descripción cualitativa de la ronda anterior, sólo
para n=24), y se aplicó a los cuatro n: el inicio de la ambigüedad **escala
con 2n²**, igual que el "codo" de la curva BOP de la Fase A.

## 1. Recuperación del código — verificación bit a bit

`cos_theta_matrix()` (operador, no depende de R) y `orientation_curve()` se
copiaron sin modificar, mismo grid `Rs` que la ronda anterior, para n=24:

```
--- fermi=True (código original, sin tocar) ---
  R= 110.0  cos=+0.930173  W=0.999854
  R= 380.0  cos=+0.807119  W=0.999969
  R= 400.0  cos=+0.769989  W=0.999988
```

Coincide exactamente con la tabla "base correcta" de
`docs/analysis_base_correcta_3_vecinos.md` §4.2 (0.930173 / 0.807119 /
0.769989). **Recuperación confirmada bit a bit**, no se rehizo desde cero.

## 2. La dependencia que cambió: `fermi=True` → `fermi=False`

`orientation_curve()` llama `system.hamiltonian(float(R))` sin pasar
`fermi`, que por defecto en `BOPSystem.hamiltonian()` es `fermi=True`. Eso
NO es un error de la ronda anterior con la información que tenía entonces
(la corrección de premisa de `analysis_fig1_carga_dipolo_sin_fermi.md` es
posterior y sólo se aplicó a `compute_bop_curve.py`), pero sí es una
dependencia que **cambió** y que hay que corregir aquí para ser consistente
con la Fase A y con `docs/STATUS.md`.

Comparación directa, mismo n=24, mismo grid, única diferencia `fermi`:

| R [a₀] | ⟨cos⟩ fermi=True (ronda anterior) | ⟨cos⟩ fermi=False (vigente) | Δ |
|---|---|---|---|
| 110 | +0.930173 | +0.929287 | −0.0009 |
| 290 | +0.842667 | +0.833720 | −0.0089 |
| 380 | +0.807119 | +0.775509 | **−0.0316** |
| 400 | +0.769989 | +0.776910 | +0.0069 |
| 480 | +0.712722 | +0.700325 | −0.0124 |
| 520 | +0.684694 | +0.725774 | **+0.0411** |
| 560 | +0.670482 | +0.630452 | −0.0400 |

Diferencias de hasta 0.04 (varios % del rango de la curva), y algo más
importante que la magnitud: **con `fermi=True` la rama R∈[350,560] a₀ es
monótona decreciente**; con `fermi=False` (el modelo correcto) **no lo es**
— hay un mínimo local en R≈370 (0.7706) y un máximo local en R≈390 (0.7823)
antes de volver a bajar. Es decir, el pseudopotencial de Fermi, al estar mal
incluido, estaba SUAVIZANDO artificialmente una estructura que en el modelo
correcto (H_A+H_mol solo) ya es no monótona ahí. Se usa `fermi=False` de
aquí en adelante, consistente con la Fase A.

## 2 bis. Test de sanidad — acotación en [-1,1]

Propiedad matemática exacta del operador coseno; se verifica explícitamente
en `scripts/compute_orientation_curve.py` con un `assert`, no se asume:

```
n=24: rango <cos theta_d> = [-0.781820, 0.931180]   [-1,1] verificado: PASA
n=25: rango <cos theta_d> = [-0.770080, 0.933306]   [-1,1] verificado: PASA
n=26: rango <cos theta_d> = [-0.758023, 0.933760]   [-1,1] verificado: PASA
n=27: rango <cos theta_d> = [-0.745404, 0.933859]   [-1,1] verificado: PASA
```

Los valores negativos (hasta −0.78) aparecen sólo en el tramo grueso
R∈[800,1800] a₀, precisamente la región identificada como no fiable en §3 —
no se interpretan como resultado físico, ver §3.

## 3. Generalización a n=24,25,26,27

`scripts/compute_orientation_curve.py --n-manifold {n}`, M_J=0, `fermi=False`:
grid fino R∈[100,800) a₀ paso 5 (140 puntos) + cola gruesa R∈[800,1800] a₀
paso 50 (21 puntos), 161 puntos por n. Coste: 130–181 s por n (0.81–1.12
s/punto), consistente con la Fase A.

```
n=24: dim=1064  ||C||_F=17.7331   161 pts en 130.4 s (0.81 s/pt)
n=25: dim=1113  ||C||_F=18.1235   161 pts en 152.3 s (0.95 s/pt)
n=26: dim=1162  ||C||_F=18.5057   161 pts en 162.6 s (1.01 s/pt)
n=27: dim=1211  ||C||_F=18.8802   161 pts en 180.9 s (1.12 s/pt)
```

Cero NaN en los cuatro barridos.

## 4. Diagnóstico de ambigüedad — corregido, y por n

**Primer intento, descartado explícitamente**: marcar como "ambiguo" todo
punto donde el índice K cambia. Con ese criterio los cuatro n dan "ambiguo
desde R=100 a₀", que es un artefacto del criterio, no un resultado físico:
la mayoría de esos cambios de K son cruces evitados aislados con
|Δ⟨cos⟩| < 0.005 — una relabelación adiabática benigna entre dos estados de
carácter casi idéntico, exactamente lo que ya se documentó que pasa en
`docs/analysis_base_correcta_3_vecinos.md` §4.1 (el índice del estado de
manifold "sube de k=4 a k=52" simplemente porque hay más estados por debajo,
sin que cambie el observable). Se descartó ese criterio y no se reporta.

**Criterio usado**: |Δ⟨cosθ_d⟩| > 0.02 entre puntos consecutivos del grid
fino (paso 5 a₀) — un salto de esa magnitud en 5 a₀ sí indica que el
autoestado elegido cambió de carácter, no una relabelación cosmética.

| n | primer R con \|Δcos\|>0.02 | R/2n² en ese punto |
|---|---|---|
| 24 | 615 a₀ | 0.534 |
| 25 | 685 a₀ | 0.548 |
| 26 | 755 a₀ | 0.558 |
| 27 | ≥800 a₀ (no se observa en el tramo fino) | ≥0.549 |

**El inicio de la ambigüedad escala con 2n²** (razón 0.53–0.56 en los cuatro
casos), igual que el "codo" de la curva BOP identificado en la Fase A
(`docs/analysis_faseA_curvas_bop_varios_n.md` §Interpretación, punto 2). Es
el mismo fenómeno físico visto desde otro observable: al acercarse a esa
fracción de 2n² el manifold empieza a tener suficiente densidad de estados
casi degenerados como para que el criterio "más bajo con carácter de
manifold" dejar de identificar un único estado bien definido.

**El tramo grueso R∈[800,1800] a₀ (paso 50) no es fiable en ningún n**: los
saltos de ⟨cos⟩ ahí alcanzan 0.3–0.6 entre puntos consecutivos (ver
`plots/rb_krb_polar/data/orientation_MJ0_n*.npz`, columna `COS` para R≥800), muy
por encima de cualquier variación físicamente razonable en 50 a₀. Es
consistente con la región de máxima densidad de estados, exactamente donde
la Fase A también vio crecer el número de mínimos locales de la curva BOP.
Después de R≈1400–1450 a₀ el patrón vuelve a ser suave y monótono (se
acerca al "codo" exterior de la Fase A, donde el peso de manifold vuelve a
ser inequívoco) — pero con un paso de 50 a₀ no se puede afirmar que ese
tramo esté bien resuelto, sólo que el salto punto a punto vuelve a ser
pequeño. No se reporta ningún valor cuantitativo de esa cola; sólo se señala
la tendencia cualitativa.

## 5. Comparación con el paper — cruce de 0.78

Con `fermi=False`, el cruce hacia abajo de ⟨cosθ_d⟩=0.78 (González-Férez
2015) ocurre en:

| n | 24 | 25 | 26 | 27 |
|---|---|---|---|---|
| primer R con cos<0.78 | 365 a₀ | 365 a₀ | 395 a₀ | 390 a₀ |

Los cuatro caen en R≈365–395 a₀, consistente con el R≈390 a₀ citado en el
paper (para n=24, allí con `fermi=True`, cruce entre 380–400 a₀). El cambio
de modelo (§2) desplaza el cruce de n=24 unos 15–25 a₀ antes, pero sigue
dentro del mismo entorno del valor publicado. Como ya advertía
`docs/analysis_verificacion_tabla_I.md` §11, esto **no confirma ni refuta**
el 0.78 del paper de forma estricta (no se sabe a qué estado/R exactos se
refiere el valor publicado), pero la coincidencia de rango se mantiene tras
la corrección del modelo.

## 6. Figura

`plots/rb_krb_polar/figures/fig_orientation_compare_n_MJ0.png` — ⟨cosθ_d⟩ vs R para
los cuatro n, R∈[100,800] a₀ (el tramo grueso R>800 se calcula y se guarda en
el `.npz` pero **no se dibuja**: por §4 no es fiable, y dibujarlo con un
paso 10× más grueso que el resto introduciría un artefacto visual de
resolución, no información). Trazo grueso = tramo con seguimiento fiable
(criterio §4); trazo fino = a partir del primer evento de ambigüedad de cada
n. Las cuatro curvas colapsan bien para R≲400 a₀ y empiezan a divergir justo
donde cada una entra en su propia región de ambigüedad — visualmente
consistente con la tabla del §4.

## Conclusiones y aplicación al proyecto

- El conjunto n=24,25,26,27 es válido para ⟨cosθ_d⟩, igual que para las
  curvas BOP de la Fase A, con el mismo tipo de escalado en 2n² gobernando
  dónde deja de ser fiable el seguimiento adiabático simple.
- **Corrección importante para las fases siguientes**: cualquier cálculo que
  reutilice `BOPSystem.hamiltonian()` debe pasar `fermi=False` de forma
  explícita para el sistema polar vigente. El valor por defecto del
  parámetro (`fermi=True`) sigue siendo deuda técnica documentada en
  `docs/STATUS.md`, y ya ha producido dos rondas con resultados calculados
  bajo la premisa equivocada (BOP en su momento, orientación aquí).
- Para la Fase E (Franck-Condon comparativo), el rango R≲600–750 a₀ (según
  n) es el que tiene un ⟨cosθ_d⟩ bien definido; más allá, cualquier
  conclusión sobre orientación tendría que apoyarse en observables distintos
  o en una malla mucho más fina, no en esta curva.

## Referencias

- `docs/analysis_base_correcta_3_vecinos.md` §4.2 — código y resultado
  original (base correcta, `fermi=True` implícito).
- `docs/analysis_verificacion_tabla_I.md` §11 — primera comparación con
  González-Férez 2015, límites de esa comparación.
- `docs/analysis_fig1_carga_dipolo_sin_fermi.md` — corrección de premisa que
  este documento extiende al cálculo de orientación.
- `docs/analysis_faseA_curvas_bop_varios_n.md` — el mismo escalado en 2n²
  visto en la curva BOP.
- `docs/PLAN_figuras_publicacion.md` — plan maestro, Fase B.
