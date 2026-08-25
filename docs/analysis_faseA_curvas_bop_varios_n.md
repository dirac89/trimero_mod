# Fase A — curvas BOP de Rb*-RbCs comparadas entre varios n

**Fecha**: 2026-08-24 (recálculo completo con RbCs; ver §0)
**Autor**: Javier Aguilera (via Claude Code)
**Relevancia**: primera fase de `docs/PLAN_figuras_publicacion.md` (figuras de
publicación). Genera la comparación entre manifolds que sirve de base a las
Fases B-E (orientación, campo eléctrico, híbrido, conclusión Franck-Condon).
**Tipo**: analysis

---

## 0. Corrección: esta fase se calculó primero con la molécula equivocada

**Este documento no siempre dijo RbCs.** Su versión del 2026-08-23 reportaba
la Fase A como *Rb\*-KRb*, y sus números eran los de KRb. Queda escrito aquí,
en vez de sustituirlo en silencio, porque el mecanismo del error es reutilizable
y conviene que la próxima ronda lo reconozca.

Qué pasó, en orden:

1. El plan (`docs/PLAN_figuras_publicacion.md`) pedía el sistema polar con
   **constantes de RbCs** (d = 1.225 D, B = 490.17 MHz).
2. El único sistema polar que existía entonces era `rb_krb_polar.BOPSystem`,
   que fija **KRb** (d = 0.566 D, B = 1.114 GHz) dentro de
   `charge_dipole.py`. `scripts/compute_bop_curve.py` lo instanciaba
   directamente: `sysm = BOPSystem(n_manifold=args.n_manifold, ...)`, sin
   ningún parámetro de molécula.
3. Al detectar la discrepancia se resolvió **por la etiqueta y no por la
   física**: se anotó en el plan que "el sistema es Rb\*-KRb, no Rb\*-RbCs" y
   se renombraron los títulos. La nota decía textualmente *"No cambia ningún
   número: es una corrección de etiqueta, no de física"*. Era cierto sobre lo
   que había, y equivocado sobre lo que se quería: el objetivo era RbCs.
4. Los commits `ce328c5` y `0fa41a0` añadieron RbCs de verdad
   (`polar_molecule.py`, `polar_rydberg/polar_system.py`,
   `rb_rbcs_polar/RbRbCsPolarSystem`), y `compute_bop_curve.py` pasó a usar
   `PolarBOPSystem` con una opción `--molecule`… **cuyo valor por defecto era
   `krb`**. El hardcode se convirtió en un default silencioso, que es peor: la
   línea de órdenes que uno copia al documento (`--n-manifold 26 --mj 0`) no
   menciona ninguna molécula, así que nada delata cuál se usó.

La huella material del error todavía está en el repositorio: los `.npz` de
`plots/rb_krb_polar/data/` **no llevan campo `molecule`** (esquema anterior al
registro de moléculas), mientras que los de `plots/rb_rbcs_polar/data/` sí
(`molecule='rbcs'`, `B_hz=4.9017e8`, `d_debye=1.225`).

Qué se corrigió en esta ronda (2026-08-24):

- Fase A **recalculada entera** con `RbRbCsPolarSystem` para n=24,25,26,27.
- `--molecule` pasa a ser **obligatorio y sin valor por defecto** en
  `compute_bop_curve.py` y en `compare_bop_curves_n.py`.
- `compare_bop_curves_n.py` **rechaza** los `.npz` sin campo `molecule` en vez
  de suponer que son de la molécula pedida.
- La figura lleva ahora **B y d escritas en el título**: una figura que declara
  sus constantes se contradice sola si la molécula no es la que dice el pie.

Los resultados de KRb **no se han borrado ni se han sobrescrito**: siguen en
`plots/rb_krb_polar/` y son un cálculo válido de Rb\*-KRb, validado contra
Aguilera-Fernández 2015. Lo que dejan de ser es "la Fase A".

---

## 1. Resumen

Curvas BOP de Rb\*-RbCs en `M_J=0` para n = 24, 25, 26, 27, en R ∈ [400, 1800]
a₀ con paso 5 a₀ (281 puntos por n), base `(n,l≥3) + (n+1)d + (n+2)p + (n+3)s`,
`N_max=6`, sin pseudopotencial de Fermi. Los cuatro barridos son numéricamente
sanos (sin NaN, peso de manifold siempre por encima del umbral de aceptación) y
confirman la tendencia de escala esperada: la transición hacia la asíntota se
desplaza a mayor R con n y colapsa en escala reducida R/2n².

Frente a KRb, RbCs profundiza la curva un **factor ≈2.42-2.47** a R = 400 a₀,
y aparece un efecto nuevo que con KRb no existía: el peso de manifold cae hasta
**0.52** en puntos aislados (con KRb nunca bajaba de 0.92), lo que acerca el
criterio de identificación por carácter a su umbral. Se documenta en §5.

## 2. Punto de control antes del barrido

Antes de gastar los ~25 min del barrido completo se comprobó, en un único
punto (R = 400 a₀, n = 25, M_J = 0, mismo bloque de 1113 estados), que las dos
clases dan números **distintos** — si hubieran coincidido, algo seguiría mal:

```
             B (GHz)   d (D)     curva de carácter        traza(H)
  KRb        1.114000  0.566     -23.100158642 GHz        -8.875089136593534e-01 E_h   (k=53, peso 0.999680)
  RbCs       0.490170  1.225     -56.141407508 GHz        -8.899884251916879e-01 E_h   (k=50, peso 0.989595)

  Delta E(curva) = -33.041248866 GHz   (RbCs - KRb)
  cociente       = 2.430347
```

Doble lectura: (a) las clases no son intercambiables, el sanity check pasa;
(b) el valor de KRb, **−23.100158642 GHz**, es exactamente el número de
referencia de la Fig. 1 anclado en
`tests/systems/rb_krb_polar/test_regression_fig1.py` (−23.100 GHz), lo que
cierra la cadena de evidencia: la Fase A original era KRb, sin ambigüedad.

## 3. Qué se corrió

```
poetry run python scripts/compute_bop_curve.py --molecule rbcs \
    --n-manifold {24,25,26,27} --mj 0 --rmin 400 --rmax 1800 --step 5 --no-plot
```

Los cuatro n se lanzaron **en paralelo**, un proceso por n, con los hilos BLAS
limitados a 2 por proceso (`OMP_NUM_THREADS=2`) en una máquina de 8 núcleos:

| n | dim(M_J=0) | base | tiempo | s/punto |
|---|---|---|---|---|
| 24 | 1064 | (24,l≥3)+25d+26p+27s | 347.1 s | 1.24 |
| 25 | 1113 | (25,l≥3)+26d+27p+28s | 384.3 s | 1.37 |
| 26 | 1162 | (26,l≥3)+27d+28p+29s | 418.9 s | 1.49 |
| 27 | 1211 | (27,l≥3)+28d+29p+30s | 463.4 s | 1.65 |

Los s/punto **no son comparables** con los de la ronda de KRb (1.06-1.38 s):
allí los barridos corrieron de uno en uno con la máquina libre, aquí cuatro a
la vez. El coste por diagonalización no depende de la molécula — la dimensión
del bloque es idéntica, sólo cambian dos constantes del Hamiltoniano.

Umbrales asintóticos ns + RbCs(N), que salen del propio cálculo:

| n | ΔE(ns) + 30B (N=5) | ΔE(ns) + 42B (N=6) |
|---|---|---|
| 24 (27s) | −48.6954 GHz | −42.8133 GHz |
| 25 (28s) | −41.3586 GHz | −35.4766 GHz |
| 26 (29s) | −35.1119 GHz | −29.2299 GHz |
| 27 (30s) | −29.7601 GHz | −23.8780 GHz |

`domain_bounds()` (dominio del remapeo semiclásico de Fermi): 1090.16 /
1185.50 / 1284.82 / 1388.14 a₀ para n=24/25/26/27, con el borde inferior en
≈105-108 a₀. **No acota nada aquí**, y además no depende de la molécula: es
una propiedad del pseudopotencial de contacto Rb-e⁻ (de cómo se leen las
tablas de longitudes de dispersión), y el modelo polar corre sin V_Fermi.
`H_A+H_mol` usa la función de onda hidrogenoide, definida en todo R, así que
el barrido cubre [400,1800] a₀ entero para los cuatro n aunque ese rango
exceda `domain_bounds()` en todos los casos. Se comprueba explícitamente para
no confundirlo con un límite físico.

## 4. Resultados

```
    n    dim   E(400) GHz   pozo min GHz   R pozo   min. loc.   E(1800) GHz   2n^2
   24   1064     -61.6546       -61.6546      400           7       -1.2782   1152
   25   1113     -56.1414       -56.1414      400           8       -1.8377   1250
   26   1162     -51.2863       -51.2863      400           9       -2.7179   1352
   27   1211     -47.1173       -47.1173      400          10       -4.1929   1458
```

Comparación directa con la Fase A anterior (KRb), mismo n, mismo R, misma base:

| n | E(400) KRb | E(400) RbCs | razón | mín. locales KRb → RbCs |
|---|---|---|---|---|
| 24 | −25.5283 | **−61.6546** | 2.4151 | 7 → 7 |
| 25 | −23.1002 | **−56.1414** | 2.4303 | 8 → 8 |
| 26 | −20.9491 | **−51.2863** | 2.4481 | 9 → 9 |
| 27 | −19.0791 | **−47.1173** | 2.4696 | 10 → 10 |

El n=25 de RbCs (−56.1414 GHz, E(1800) = −1.8377 GHz, 8 mínimos) coincide con
`docs/analysis_rb_rbcs_curvas_n25.md`, que se calculó por separado — es una
comprobación cruzada independiente dentro del propio repositorio.

Los cuatro `.npz` recalculados salieron **bit a bit idénticos** a los que ya
había en `plots/rb_rbcs_polar/data/` de los commits `ce328c5`/`0fa41a0`
(`max|ΔE| = 0.000e+00 GHz`, `K` idéntico): el cálculo es reproducible y los
datos que ya estaban en el repositorio eran correctos. Lo que estaba mal era el
documento de la Fase A y la figura que lo acompañaba, no esos datos.

Figura:

```
poetry run python scripts/compare_bop_curves_n.py --molecule rbcs \
    --n 24 25 26 27 --mj 0 --ymin -65 --ymax 2
→ plots/rb_rbcs_polar/figures/fig_compare_n_MJ0.png
```

Dos notas sobre la figura:

- La ventana en Y **hay que ampliarla a −65 GHz**. El valor por defecto
  (−25 GHz), heredado de la Fig. 1 de KRb, recortaría las cuatro curvas de
  RbCs por debajo: el pozo de n=24 llega a −61.65 GHz.
- Vive en `plots/rb_rbcs_polar/figures/`, no en `plots/rb_krb_polar/figures/`.
  La ruta la fija la molécula (`plots/rb_<mol>_polar/`), y mezclar sistemas en
  un mismo directorio es justo lo que `.claude/CLAUDE.md` prohíbe. El
  `fig_compare_n_MJ0.png` que sigue en `plots/rb_krb_polar/figures/` es la
  figura de KRb: correcta como resultado de KRb, pero ya no es la Fase A.

## 5. Verificación de sanidad numérica

Metadatos de los cuatro `.npz` (todos): `molecule='rbcs'`, `B_hz=4.90170e+08`,
`d_debye=1.225`, `N_max=6`, `M_J=0`, `character_weight=0.5`.

```
    n  nan    W_min    W_max         K distintos    R[0]   R[-1]     N
   24    0   0.5631   1.0000     [51,52,53,54,55]     400    1800   281
   25    0   0.6116   1.0000  [50,51,52,53,54,55]     400    1800   281
   26    0   0.8098   0.9999  [50,51,52,53,54,55]     400    1800   281
   27    0   0.5221   0.9998  [50,51,52,53,54,55]     400    1800   281
```

Sin NaN en ningún punto, malla correcta (400→1800 a₀, 281 puntos), y peso de
manifold por encima del umbral 0.5 en todo el barrido. **Pero el margen es
mucho menor que con KRb, y esto es nuevo:**

| n | W_min RbCs | R(W_min) | #W<0.90 | #W<0.70 | W_min KRb |
|---|---|---|---|---|---|
| 24 | 0.5631 | 760 a₀ | 3 | 1 | 0.9983 |
| 25 | 0.6116 | 610 a₀ | 9 | 3 | 0.9985 |
| 26 | 0.8098 | 720 a₀ | 6 | 0 | 0.9216 |
| 27 | **0.5221** | 925 a₀ | 22 | 5 | 0.9984 |

Con KRb el peso no bajaba de 0.92 y el índice k cambiaba 2-3 veces; con RbCs
baja a 0.52 (n=27, R=925 a₀) y k recorre 5-6 valores. Es coherente con la
física, no un fallo numérico: el dipolo de RbCs es 2.16 veces mayor, así que
el acoplamiento carga-dipolo mezcla el manifold con los vecinos (n+1)d y
(n+2)p mucho más, y los cruces evitados son más anchos y más frecuentes.

**Limitación que hay que arrastrar a las fases siguientes**: en esos puntos
aislados el estado es una mezcla casi al 50 % y "la curva del manifold" es una
etiqueta discutible, no un hecho. Con n=27 hay 5 puntos por debajo de 0.70 y
el mínimo roza 0.52. No se ha forzado ni suavizado nada: los puntos están en
la curva tal como salen. Si en la Fase C (campo eléctrico) el peso cruzara
0.5, `character_curve` devolvería `NaN` en vez de un valor equivocado — el
código falla de forma visible, que es lo que se quiere, pero conviene vigilar
`W_min` en cada barrido nuevo en vez de darlo por bueno.

## 6. Interpretación física

**6.1. Por qué RbCs profundiza ≈2.43× y no 10×.** Los dos límites del rotor
en el campo del ion acotan la respuesta:

- campo fuerte (d·F ≫ B), E → −d·F: la razón tendería a
  d(RbCs)/d(KRb) = **2.1643**;
- campo débil (d·F ≪ B), E → −(dF)²/6B: la razón tendería a
  (d_R/d_K)²·(B_K/B_R) = **10.6458**.

Lo observado (2.4151 → 2.4696) está justo por encima del límite de campo
fuerte y a un orden de magnitud del perturbativo. El parámetro de régimen lo
confirma (`F_ion = 1/R²`, cota inferior porque ignora el campo del electrón):

| R (a₀) | d·F/B (KRb) | d·F/B (RbCs) |
|---|---|---|
| 400 | 8.22 | 40.43 |
| 800 | 2.06 | 10.11 |
| 1200 | 0.91 | 4.49 |
| 1800 | 0.41 | 2.00 |

A R = 400 a₀ RbCs está profundamente en régimen pendular (d·F/B ≈ 40): el
rotor está prácticamente orientado a lo largo del campo y la energía escala
como −d·F, de ahí que la razón se pegue a 2.164. El pequeño exceso sobre
2.164 va en la dirección correcta: la energía de punto cero pendular escala
como √(2B·d·F) y es relativamente mayor en KRb (B 2.27 veces mayor), lo que
hace la curva de KRb menos profunda de lo que predice el escalado lineal puro.
El número exacto sale de la diagonalización, no de esta estimación.

Consecuencia para las fases siguientes: **RbCs entra en régimen pendular a R
mucho mayores que KRb**. La transición d·F/B ≈ 1 cae en R ≈ 1150 a₀ para KRb y
en R ≈ 2550 a₀ para RbCs, fuera de la ventana. Toda la ventana [400,1800] a₀
es régimen orientado para RbCs, y sólo su mitad interior lo es para KRb. Esto
es directamente relevante para la Fase B (⟨cos θ_d⟩) y para la Fase C (competencia
entre el campo del ion y el campo DC externo).

**6.2. La tendencia de escala se confirma, con el mismo matiz que en KRb.**

1. El "pozo más profundo" que reporta `report()` está anclado en R = 400 a₀
   para los cuatro n. Eso **no** es el desplazamiento a mayor R de la
   tendencia física: R = 400 a₀ es el borde inferior de una ventana ABSOLUTA
   fija, heredada de la Fig. 1 original, que no escala con n. La curva para n
   grande sigue bajando hacia R = 400 sin haber alcanzado su mínimo real, que
   estaría a R menor, fuera de rango. Reportar sólo esa cifra como "el pozo se
   mueve" sería forzar una lectura que los números no sostienen.
2. La tendencia real está en el **panel reducido** (R/2n²): las cuatro curvas
   colapsan sobre la misma forma, con la subida brusca hacia la asíntota
   centrada en R/2n² ≈ 1.0 para los cuatro n. Midiendo el "codo" como el R
   donde la curva cruza el 10 % de su profundidad en 400 a₀:

   | n | R_codo (a₀) | R_codo/2n² | (KRb: R_codo) | (KRb: R/2n²) |
   |---|---|---|---|---|
   | 24 | 1445 | 1.2543 | 1390 | 1.2066 |
   | 25 | 1555 | 1.2440 | 1490 | 1.1920 |
   | 26 | 1665 | 1.2315 | 1595 | 1.1797 |
   | 27 | 1780 | 1.2209 | 1705 | 1.1694 |

   En escala absoluta el codo se desplaza ~110 a₀ por unidad de n; en escala
   reducida es prácticamente constante (1.22-1.25). Es el escalado ~2n² de un
   manifold de Rydberg, y es el criterio de sanidad correcto, no el mínimo
   global de una ventana fija. RbCs desplaza el codo sistemáticamente ~55-75 a₀
   hacia fuera respecto de KRb, coherente con §6.1: con dipolo mayor, el
   acoplamiento sigue siendo apreciable a R más grandes.
3. La **profundidad decrece monótonamente con n** (61.65 → 56.14 → 51.29 →
   47.12 GHz a R = 400 a₀ fijo), igual que con KRb.
4. El número de **mínimos locales crece con n** (7 → 8 → 9 → 10), idéntico a
   KRb: es estructura del manifold (más estados, más cruces evitados en la
   misma ventana absoluta), no de la molécula. Que este recuento sea el mismo
   para las dos moléculas es en sí una comprobación de consistencia.
5. **E(R=1800 a₀)** crece en magnitud con n (−1.28 → −1.84 → −2.72 → −4.19
   GHz): a n mayor, R = 1800 a₀ está relativamente más cerca del manifold
   (2n² crece), así que la curva aún no se ha aplanado en el borde derecho.
   Con RbCs estos valores son ~4-6 veces mayores en magnitud que con KRb
   (−0.21 → −1.03 GHz): la cola del potencial es mucho más larga.
6. **La curva baja por debajo del umbral asintótico ns + RbCs(N=5)** en los
   cuatro n (p. ej. n=24: −61.65 GHz frente a −48.70 GHz). Con KRb esto no
   pasaba en la misma ventana. Significa que en la región interior la curva
   del manifold está energéticamente por debajo de canales rotacionalmente
   excitados del ns vecino, y ahí es donde aparecen los pesos bajos de §5.

## 7. Conclusiones y aplicación al proyecto

- La Fase A queda **rehecha con RbCs** para n = 24, 25, 26, 27, verificada y
  reproducible. Es la base válida para las Fases B-E.
- **Las Fases B, C, D y E deben recalcularse con RbCs.** La Fase B
  (⟨cos θ_d⟩) está hoy documentada en `docs/analysis_faseB_orientacion_varios_n.md`
  con KRb, por el mismo default silencioso; hay datos de orientación de RbCs
  en `plots/rb_rbcs_polar/data/orientation_MJ0_n2*.npz`, pero el documento de
  la Fase B todavía no los usa. `scripts/compute_orientation_curve.py` y
  `scripts/compare_orientation_n.py` **siguen teniendo `--molecule` con
  default `krb`**: hay que quitárselo igual que se ha hecho aquí, antes de
  cerrar la Fase B.
- Presentar siempre la figura con los **dos paneles**: el reducido es el que
  sostiene la afirmación de escala; el absoluto por sí solo induce a leer mal
  el "pozo más profundo".
- Vigilar `W_min` en cada barrido nuevo (§5). Con RbCs el margen sobre el
  umbral de carácter es de 0.02 en el peor punto, no de 0.42 como con KRb.
- La ventana fija [400,1800] a₀ es una decisión de estilo consistente con la
  Fig. 1, no una limitación del cálculo. Para capturar el mínimo global real
  de cada n haría falta bajar `--rmin` por debajo de 400 a₀; no se hizo porque
  no está en el alcance de la Fase A.

## 8. Ficheros

- Datos: `plots/rb_rbcs_polar/data/fig1_ad_MJ0_n{24,25,26,27}.npz`
- Figura: `plots/rb_rbcs_polar/figures/fig_compare_n_MJ0.png`
- Código: `scripts/compute_bop_curve.py`, `scripts/compare_bop_curves_n.py`
  (ambos con `--molecule` obligatorio desde esta ronda)
- Sistema: `src/trimero/systems/rb_rbcs_polar/__init__.py` (`RbRbCsPolarSystem`)
  sobre `src/trimero/systems/polar_rydberg/polar_system.py`

## 9. Referencias

- `docs/PLAN_figuras_publicacion.md` — plan maestro, Fase A (y la nota de
  corrección del 2026-08-23 que este documento rectifica en §0).
- `docs/analysis_rb_rbcs_curvas_n25.md` — cálculo independiente de RbCs n=25,
  usado aquí como comprobación cruzada.
- `docs/analysis_base_correcta_3_vecinos.md` — identificación de la curva por
  carácter, no por índice.
- `docs/STATUS.md` — modelo polar vigente y números de referencia de KRb n=25.
- Aguilera-Fernández, Sadeghpour, Schmelcher & González-Férez, J. Phys.: Conf.
  Ser. **635**, 012023 (2015) — Ec. 1, el Hamiltoniano que se diagonaliza.
