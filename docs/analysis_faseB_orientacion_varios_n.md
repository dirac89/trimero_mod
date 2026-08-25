# Fase B — orientación ⟨cosθ_d⟩ de Rb*-RbCs comparada entre varios n

**Fecha**: 2026-08-24 (recálculo completo con RbCs; ver §0)
**Autor**: Javier Aguilera (via Claude Code)
**Relevancia**: segunda fase de `docs/PLAN_figuras_publicacion.md`. Da el
observable de orientación sobre el que se apoyan la Fase C (campo eléctrico
DC, que compite con el campo del Rydberg) y la Fase E (Franck-Condon
comparativo).
**Tipo**: analysis

---

## 0. Corrección: esta fase se calculó primero con la molécula equivocada

**Este documento no siempre dijo RbCs.** Su versión del 2026-08-23 reportaba
la Fase B como *Rb\*-KRb* y sus números eran los de KRb. Es el **mismo error
estructural** que ya documentó `docs/analysis_faseA_curvas_bop_varios_n.md`
§0, repetido en el otro par de scripts:
`scripts/compute_orientation_curve.py` y `scripts/compare_orientation_n.py`
tenían `--molecule` con `default="krb"`. Un default silencioso no aparece en
la línea de órdenes que uno copia al documento
(`--n-manifold 25`), así que nada delata la molécula usada. Queda escrito
aquí, en vez de sustituirlo en silencio, porque el mecanismo es reutilizable.

Al corregirlo aparecieron **dos problemas más**, que no estaban previstos y
que conviene dejar anotados porque los dos habrían producido un resultado
falso con muy buena pinta:

**(a) Los `.npz` de RbCs que ya existían no eran homogéneos entre sí.** Los de
n=24, 26 y 27 (commit `0fa41a0`) son `schema_version=1`, malla
R ∈ [100,1800] a₀ con la cola R>800 a paso 50, y **sin `COS2`**. El de n=25
(commit `ce328c5`) es `schema_version=2`, malla R ∈ [400,1800] a₀ a paso 5
uniforme y **con `COS2`**. Superponer las cuatro en una figura habría dibujado
curvas que no son comparables punto a punto: la de n=25 no cubre R<400 (donde
está la región más orientada) y las otras tres tienen la cola diez veces menos
resuelta. La comprobación de "¿ya están los datos?" no puede pararse en "el
fichero existe y dice `rbcs`".

**(b) El propio test de sanidad física estaba mal formulado, y falló.** Se
planteó como: *dado que RbCs entra en régimen pendular mucho antes que KRb
(Fase A §6.1), ⟨cosθ_d⟩ de RbCs debe ser MAYOR que el de KRb al mismo R*. Con
ese enunciado el test **falla** en R = 1200 y 1800 a₀. La razón no es que la
Fase A esté mal: es que ahí ⟨cosθ_d⟩ es **negativo** — el dipolo se ha dado la
vuelta y apunta antiparalelo al eje. "Más orientada" en esa región significa
más cerca de −1, no de +1. Con el criterio correcto, `|⟨cosθ_d⟩|`, el test
pasa en todos los R (§2). Investigar ese fallo es lo que destapó el resultado
del §6, que es el hallazgo principal de esta ronda.

Qué se corrigió en esta ronda (2026-08-24):

- `--molecule` **obligatorio y sin valor por defecto** en
  `compute_orientation_curve.py` y `compare_orientation_n.py`.
- `compare_orientation_n.py` **rechaza** los `.npz` sin campo `molecule`, y
  además **rechaza mallas incompatibles entre n** (salvaguarda nueva, por (a)).
- Fase B **recalculada entera** con `RbRbCsPolarSystem` para n=24,25,26,27 en
  una malla homogénea que es superset estricto de todo lo anterior.
- `compute_orientation_curve.py` informa de `W_min` en cada barrido, que es lo
  que la Fase A §7 pidió vigilar.
- La figura lleva **B y d escritas en el título** y ya no recorta los valores
  negativos (el `ylim` era `(-0.05, 1.0)`, que habría escondido el §6 entero).

Los resultados de KRb no se han borrado: siguen en `plots/rb_krb_polar/` y son
un cálculo válido de Rb\*-KRb. Lo que dejan de ser es "la Fase B". Nota
práctica: como esos `.npz` son de esquema antiguo y no declaran molécula, la
salvaguarda nueva los rechaza — regenerar la figura de KRb exigiría
recalcularlos con `--molecule krb`. Es deliberado: preferimos un fallo ruidoso
a una figura silenciosamente mal etiquetada.

---

## 1. Resumen

⟨cosθ_d⟩ y ⟨cos²θ_d⟩ de Rb\*-RbCs en `M_J=0` para n = 24, 25, 26, 27, en
R ∈ [100, 1800] a₀ a paso 5 a₀ uniforme (341 puntos por n), `N_max=6`,
criterio de carácter `peso > 0.5`, sin pseudopotencial de Fermi — los mismos
parámetros de la Fase A.

Tres resultados:

1. **RbCs está sistemáticamente más orientada que KRb** a igual R, como
   predecía la Fase A §6.1 (§2). En la región interior la diferencia es
   grande: a R = 400 a₀, ⟨cos⟩ = 0.885 frente a 0.777.
2. **El dipolo se da la vuelta** a R/2n² ≈ 0.87, de forma **continua y bien
   resuelta**, porque el campo del electrón Rydberg supera ahí al del ion
   (§6). Esto **corrige** al documento de la ronda de KRb, que descartaba los
   valores negativos como artefacto de resolución.
3. **El seguimiento del estado es mucho más ambiguo que con KRb** (§5): 68-109
   eventos de salto frente a 0-17, y empiezan mucho antes. La región limpia se
   reduce a R ≲ 400-500 a₀.

## 2. Punto de control antes de aceptar nada

Mismo n (25), mismo M_J, mismo bloque de 1113 estados, mismos R; única
diferencia, la molécula. Se contrasta además contra el rotor **aislado** en el
campo del ion E = ±1/R², que es el límite de orientación pura sin estructura
de Rydberg:

```
   R (a0) |        KRb   dF/B |       RbCs   dF/B | |RbCs|-|KRb| |  rot +1/R2  rot -1/R2  (RbCs)
      400 |   0.771388   8.22 |   0.886902  40.43 |    +0.115514 |   0.888694  -0.888694
      600 |   0.661299   3.65 |   0.818666  17.97 |    +0.157368 |   0.832804  -0.832804
      800 |   0.457831   2.06 |   0.630914  10.11 |    +0.173083 |   0.776520  -0.776520
     1200 |  -0.741734   0.91 |  -0.866655   4.49 |    +0.124921 |   0.661048  -0.661048
     1800 |  -0.380063   0.41 |  -0.714414   2.00 |    +0.334351 |   0.480763  -0.480763

  (k, peso) KRb : R=400:(k=53,W=1.000)  R=800:(k=54,W=0.999)  R=1800:(k=55,W=1.000)
  (k, peso) RbCs: R=400:(k=50,W=0.990)  R=800:(k=53,W=0.984)  R=1800:(k=55,W=1.000)

  [1] ¿números DISTINTOS en todos los R?                SI
  [2] ¿RbCs MAS orientada (|cos| mayor) en todos los R?  SI — sanidad física OK
```

Dos lecturas más allá de "son distintos":

- **RbCs está saturada a R = 400 a₀**: 0.886902 frente a 0.888694 del rotor
  aislado, el 99.8 % del límite. KRb, en cambio, da 0.771388 frente a 0.751789
  del rotor — la **excede** en un 2.6 %, porque el campo del electrón Rydberg
  suma al del ion y KRb, menos saturada, todavía responde a esa diferencia.
  RbCs ya no puede: está tan orientada que no le queda margen.
- El parámetro de régimen `d·F/B` (con F = 1/R², cota inferior) es 5× mayor
  para RbCs a todo R, exactamente el factor que anticipó la Fase A §6.1.

## 3. Qué se corrió

```
poetry run python scripts/compute_orientation_curve.py --molecule rbcs \
    --n-manifold {24,25,26,27} --mj 0 --n-max 6 --weight 0.5 \
    --rmin-fine 100 --rmax-fine 800 --step-fine 5 \
    --rmax-coarse 1800 --step-coarse 5
```

Malla **uniforme** R ∈ [100,1800] a₀ paso 5 a₀ (341 puntos), idéntica en los
cuatro n. Es superset estricto de las dos mallas heterogéneas del §0(a): del
n=25 antiguo (paso 5 en [400,1800]) y de los n=24/26/27 antiguos (paso 5 en
[100,800], paso 50 después). El paso 5 en toda la cola es lo que permite el
§6: con paso 50 no se puede distinguir un vuelco continuo de un salto.

| n | dim(M_J=0) | ‖C‖_F | tiempo | s/punto |
|---|---|---|---|---|
| 24 | 1064 | 17.7331 | 350.9 s | 1.03 |
| 25 | 1113 | 18.1235 | 395.7 s | 1.16 |
| 26 | 1162 | 18.5057 | 448.9 s | 1.32 |
| 27 | 1211 | 18.8802 | 495.1 s | 1.45 |

(Los cuatro en paralelo, `OMP_NUM_THREADS=2` por proceso, máquina de 8
núcleos; no comparables con los s/punto de la ronda de KRb, que corrió en
serie.)

**Reproducibilidad**: en los puntos donde la malla nueva y la vieja solapan,
los valores salen **bit a bit idénticos** (`max|Δcos| = 0.000e+00` en los
cuatro n; 161 puntos en común para n=24/26/27 y 281 para n=25). Los datos que
ya estaban en el repositorio eran correctos; lo que no servía era su
heterogeneidad y el documento que los describía.

## 4. Orientación en la región interior

⟨cosθ_d⟩ a tres radios de la región donde el seguimiento es limpio:

| n | RbCs(150) | KRb(150) | RbCs(300) | KRb(300) | RbCs(400) | KRb(400) |
|---|---|---|---|---|---|---|
| 24 | 0.941834 | 0.909347 | 0.914007 | 0.826971 | 0.885163 | 0.776910 |
| 25 | 0.941226 | 0.909315 | 0.912775 | 0.825306 | 0.886902 | 0.771388 |
| 26 | 0.942329 | 0.909239 | 0.910900 | 0.824101 | 0.884678 | 0.766770 |
| 27 | 0.911010 | 0.909249 | 0.890773 | 0.823425 | 0.883863 | 0.764121 |

Dos cosas: (i) la orientación de RbCs es **casi independiente de n** en la
región interior (0.911-0.942 a R=150 a₀), porque está saturada y ya no
distingue detalles del campo; (ii) la diferencia con KRb **crece con R**
(+0.03 a R=150, +0.09 a R=300, +0.11 a R=400): a R pequeño las dos están
orientadas y la ventaja de RbCs se nota poco; al debilitarse el campo, KRb se
desorienta antes.

Rangos globales (los cuatro n, malla completa):

```
    n   cos_min   cos_max  cos2_min  cos2_max
   24 -0.894826  0.946329  0.085222  0.900642
   25 -0.888975  0.946571  0.096947  0.901071
   26 -0.882820  0.946799  0.096609  0.901461
   27 -0.876304  0.946959  0.103078  0.901749
```

La acotación exacta `⟨cos⟩ ∈ [-1,1]` y `⟨cos²⟩ ∈ [0,1]` se verifica con
`assert` en los cuatro barridos: **PASA**. Cero NaN en los 4 × 341 puntos.

## 5. Ambigüedad de seguimiento — peor que con KRb, por n

Es lo que la Fase A §5 advirtió que había que vigilar (`W_min` bajó de 0.92
con KRb a 0.52 con RbCs), y se confirma en este observable. Criterio: evento =
`|Δ⟨cos⟩| > 0.02` entre puntos consecutivos del paso de 5 a₀ (el mismo de la
ronda de KRb; un cambio de índice K por sí solo NO se marca, ver §5.1).

Comparación **en la misma malla** (R < 800 a₀, paso 5, donde las dos moléculas
tienen datos comparables):

| n | KRb: #ev | KRb: 1er R | RbCs: #ev | RbCs: 1er R | mayor \|Δcos\| RbCs |
|---|---|---|---|---|---|
| 24 | 17 | 620 a₀ | 18 | **225 a₀** | 0.4105 |
| 25 | 10 | 690 a₀ | 20 | **230 a₀** | 0.3975 |
| 26 | 3 | 760 a₀ | 10 | **585 a₀** | 0.1710 |
| 27 | 0 | — | 24 | **150 a₀** | 0.3562 |

En la malla completa [100,1800] a₀ el recuento de RbCs es **68 / 69 / 78 /
109** eventos para n = 24/25/26/27. Peso de manifold en el mismo barrido:

| n | W_min | R(W_min) | #W<0.90 | #W<0.70 |
|---|---|---|---|---|
| 24 | 0.5631 | 760 a₀ | 5 | 1 |
| 25 | 0.5219 | 230 a₀ | 12 | 4 |
| 26 | 0.8098 | 720 a₀ | 6 | 0 |
| 27 | 0.5221 | 925 a₀ | 28 | 6 |

**Dónde deja de ser fiable, por tramos** (mayor |Δcos| observado en cada uno):

| n | R < 400 a₀ | 400 ≤ R < 600 | 600 ≤ R < 800 |
|---|---|---|---|
| 24 | 4 ev, max 0.0448 | 0 ev | 14 ev, max **0.4105** |
| 25 | 4 ev, max 0.0902 | 4 ev, max 0.0571 | 12 ev, max **0.3975** |
| 26 | 0 ev | 2 ev, max 0.1195 | 8 ev, max **0.1710** |
| 27 | 8 ev, max 0.0567 | 6 ev, max 0.1358 | 10 ev, max **0.3562** |

La lectura honesta: **R ≲ 400 a₀ es limpio** en los cuatro n (los eventos que
hay son de magnitud ≤0.09, relabelaciones menores sobre una curva suave), y la
degradación seria empieza en R ≈ 600 a₀, con saltos de hasta 0.41 — que no son
física, son el criterio de carácter cambiando de estado. La curva de la figura
marca ese tramo con trazo fino y no se debe leer punto a punto.

### 5.1. El escalado en 2n² de la ronda de KRb NO se reproduce

El documento de KRb reportaba que el inicio de la ambigüedad escalaba con 2n²
(razón 0.53-0.56 en los cuatro n) y lo interpretaba como el mismo fenómeno que
el "codo" de la Fase A. **Con RbCs eso no se sostiene.** El primer evento cae
en R/2n² = 0.195 / 0.184 / 0.433 / 0.103 para n = 24/25/26/27: dispersión de
un factor 4, sin patrón.

Con un criterio más robusto — primer evento **grande**, `|Δcos| > 0.10`, que
no puede ser una relabelación benigna — sale un patrón, pero **el contrario**
del que tenía KRb:

| n | RbCs: 1er \|Δcos\|>0.10 | R/2n² | KRb |
|---|---|---|---|
| 24 | 760 a₀ | 0.6597 | ninguno en R<800 |
| 25 | 610 a₀ | 0.4880 | ninguno en R<800 |
| 26 | 585 a₀ | 0.4327 | ninguno en R<800 |
| 27 | 470 a₀ | 0.3224 | ninguno en R<800 |

Decrece con n en vez de crecer, tanto en R absoluto como reducido. Y KRb no
tiene **ni un solo** evento de esa magnitud en R < 800 a₀, mientras RbCs los
tiene en los cuatro n. No se propone aquí una explicación cerrada: lo que se
puede afirmar con los datos es que a mayor n el manifold tiene más estados
vecinos y el dipolo grande de RbCs los mezcla más, así que la región limpia se
encoge al subir n. El escalado limpio en 2n² era una propiedad del caso débil
(KRb), no una ley general — conviene no seguir citándolo como tal.

## 6. El vuelco del dipolo: es física, no un fallo de seguimiento

⟨cosθ_d⟩ cambia de signo una vez en cada n. Con paso 5 a₀ se ve que la
transición es **continua y bien resuelta**:

```
  n=24: R  990->995   cos +0.0328 -> -0.0262  (|d|=0.0590)  K 54->54  W 0.966->0.971
  n=25: R 1080->1085  cos +0.0192 -> -0.0300  (|d|=0.0492)  K 54->54  W 0.686->0.941
  n=26: R 1170->1175  cos +0.0424 -> -0.0087  (|d|=0.0511)  K 53->53  W 0.964->0.966
  n=27: R 1270->1275  cos +0.0081 -> -0.0397  (|d|=0.0477)  K 53->53  W 0.962->0.962
```

Un solo cruce por n, el índice K **no cambia**, el peso de manifold es alto a
ambos lados, y el salto (≈0.05) es del orden de la variación suave de la curva
en 5 a₀. Nada de eso es compatible con un cambio de estado.

Y escala con 2n² de forma muy limpia:

| n | R del vuelco | R/2n² |
|---|---|---|
| 24 | 995 a₀ | 0.8637 |
| 25 | 1085 a₀ | 0.8680 |
| 26 | 1175 a₀ | 0.8691 |
| 27 | 1275 a₀ | 0.8745 |

**Verificación por un camino independiente.** El término de acoplamiento del
Hamiltoniano es `V = −d (F·cosθ_d)`, con `F = F_ion + F_elec`. Construyendo el
mismo `ChargeDipoleHamiltonian` con `B=0` se aísla `V`, y
`F_eff = −⟨V⟩ / (d ⟨cosθ_d⟩)` es el campo efectivo **con signo** que ve el
rotor. `F_ion = +1/R²` es positivo siempre, así que un `F_eff` negativo sólo
puede venir del electrón. Para n=25:

```
   R (a0)      <cos>   F_eff (u.a.)   F_ion=1/R^2  F_eff/F_ion    K       W
      400   0.886902   2.059880e-05  6.250000e-06       3.2958   50  0.9896
      800   0.630914   2.119998e-05  1.562500e-06      13.5680   53  0.9841
     1000   0.820716   1.711823e-05  1.000000e-06      17.1182   53  0.9610
     1050   0.373823   3.528165e-05  9.070295e-07      38.8980   53  0.9594
     1080   0.019182   4.793681e-04  8.573388e-07     559.1350   54  0.6857
     1085  -0.030024  -4.181689e-04  8.494553e-07    -492.2789   54  0.9407
     1100  -0.195874  -6.305231e-05  8.264463e-07     -76.2933   54  0.9365
     1200  -0.866655  -1.446044e-05  6.944444e-07     -20.8230   54  0.9647
     1800  -0.714414  -1.094026e-06  3.086420e-07      -3.5446   55  1.0000
```

`F_eff` cambia de signo **entre R = 1080 y 1085 a₀**, exactamente el mismo
intervalo en que ⟨cos⟩ cruza cero. (Los valores enormes en 1080/1085 son
artefacto de la definición: ⟨cos⟩ está en el denominador y pasa por cero; lo
que importa es el signo, no la magnitud, ahí.) Interpretación: el rotor sigue
al campo neto, y el campo neto se invierte cuando la molécula sale de la zona
donde domina el ion. Es el mismo fenómeno que `scripts/plot_orientation_alignment.py`
ya anticipaba al dibujar las dos ramas de rotor, `E = ±Ẑ/R²`, como referencia.

**Esto corrige a la ronda de KRb.** Aquel documento (§2bis y §4) decía: *"Los
valores negativos (hasta −0.78) aparecen sólo en el tramo grueso R∈[800,1800]
a₀, precisamente la región identificada como no fiable — no se interpretan
como resultado físico"*. Con paso 50 a₀ esa cautela era correcta: no había
resolución para afirmar nada. Con paso 5 a₀ se ve que el vuelco es real, y
además ocurre también en KRb (§2: −0.742 en R=1200 a₀ con peso 0.999). No era
un artefacto; era una región mal muestreada.

## 7. Cruce de 0.78 (González-Férez 2015)

| n | KRb | RbCs |
|---|---|---|
| 24 | 365 a₀ | **625 a₀** |
| 25 | 365 a₀ | **610 a₀** |
| 26 | 395 a₀ | **585 a₀** |
| 27 | 390 a₀ | **470 a₀** |

Con RbCs el cruce se desplaza 100-260 a₀ hacia fuera, que es la consecuencia
directa de §2 y §4: al estar más orientada, RbCs tarda más en caer por debajo
de 0.78. Coherente, pero conviene ser explícito sobre qué **no** significa: el
0.78 de González-Férez 2015 es un número de **KRb**, y la coincidencia de la
ronda de KRb con él (365-395 a₀ frente a ≈390 a₀ publicado) era la comparación
pertinente. Los 470-625 a₀ de RbCs **no comparan con nada publicado**; la línea
de 0.78 se mantiene en la figura como referencia visual heredada, no como
validación. Además, para n=24/25 el cruce cae en R ≈ 610-625 a₀, dentro de la
región que §5 marca como ya degradada, así que ese valor concreto tiene
incertidumbre de seguimiento; el de n=27 (470 a₀) es el más limpio de los
cuatro.

## 8. Figura

```
poetry run python scripts/compare_orientation_n.py --molecule rbcs \
    --n 24 25 26 27 --mj 0 --rmax-plot 1800
→ plots/rb_rbcs_polar/figures/fig_orientation_compare_n_MJ0.png
```

Dos paneles, misma decisión que la Fase A: el absoluto por sí solo no sostiene
la afirmación de escala, el reducido R/2n² sí — y aquí el colapso de las cuatro
curvas sobre el vuelco (R/2n² ≈ 0.87) es lo más visible de la figura.

Tres ajustes respecto de la versión de KRb, los tres necesarios:

- `--rmax-plot 1800` en vez del default 800. Ese default venía de que en KRb la
  cola se calculaba a paso 50 y ⟨cos⟩ ya había decaído; con RbCs la curva sigue
  orientada mucho más allá de 800 a₀ y ahí está el §6.
- `ylim` de (−0.05, 1.0) a (−1.02, 1.02). El anterior **habría recortado todos
  los valores negativos**, es decir, habría escondido el resultado principal.
- Título con B y d, para que la figura se contradiga sola si la molécula no es
  la que dice el pie.

Trazo grueso = seguimiento fiable; trazo fino = puntos con eventos de
ambigüedad cercanos (§5). El tramo 500-1000 a₀ sale mayoritariamente fino: es
correcto que así sea.

## 9. Conclusiones y aplicación al proyecto

- El conjunto n=24,25,26,27 queda calculado y verificado para RbCs, en malla
  homogénea, con `COS` y `COS2`. Es la base para las Fases C y E.
- **Región de confianza: R ≲ 400-500 a₀.** Ahí ⟨cosθ_d⟩ es limpio, saturado
  (0.88-0.94) y casi independiente de n. Cualquier conclusión cuantitativa de
  orientación en las fases siguientes debería apoyarse en ese tramo.
- **Para la Fase C (campo DC externo)**: RbCs está saturada por el campo del
  Rydberg a R ≲ 400 a₀ (99.8 % del límite del rotor), así que un campo externo
  de 100-500 V/m **no puede aumentar la orientación ahí** — no queda margen. Su
  efecto interesante estará en R ≳ 1000 a₀, cerca del vuelco del §6, donde el
  campo neto del Rydberg pasa por cero y un campo externo pequeño es
  comparativamente grande. Conviene diseñar el barrido de la Fase C alrededor
  de esa región, no de la interior.
- **Para la Fase E (Franck-Condon)**: el vuelco del §6 cae en R/2n² ≈ 0.87, es
  decir dentro de la ventana [400,1800] a₀ de la Fase A para los cuatro n. Un
  cambio de orientación del dipolo en mitad del pozo es relevante para los
  solapamientos vibracionales y no debería ignorarse.
- **No seguir citando el escalado en 2n² del inicio de la ambigüedad** (§5.1):
  era una propiedad del caso KRb, y con RbCs no se reproduce.

## 10. Ficheros

- Datos: `plots/rb_rbcs_polar/data/orientation_MJ0_n{24,25,26,27}.npz`
  (`schema_version=2`, con `COS2`, malla [100,1800] paso 5)
- Figura: `plots/rb_rbcs_polar/figures/fig_orientation_compare_n_MJ0.png`
  (+ las individuales `orientation_MJ0_n{n}_Nmax6.png`)
- Código: `scripts/compute_orientation_curve.py`,
  `scripts/compare_orientation_n.py` (ambos con `--molecule` obligatorio)
- Sistema: `RbRbCsPolarSystem` sobre `polar_rydberg/polar_system.py`

## 11. Referencias

- `docs/analysis_faseA_curvas_bop_varios_n.md` — §0 (el mismo error
  estructural, documentado primero ahí) y §6.1 (régimen pendular, que este
  documento verifica en el observable de orientación).
- `docs/analysis_base_correcta_3_vecinos.md` §4.2 — origen del cálculo de
  ⟨cosθ_d⟩.
- `docs/analysis_verificacion_tabla_I.md` §11 — límites de la comparación con
  González-Férez 2015.
- `docs/analysis_rb_rbcs_curvas_n25.md` — cálculo independiente de RbCs n=25.
- `docs/PLAN_figuras_publicacion.md` — plan maestro, Fase B.
- González-Férez, Schmelcher et al. (2015) — origen del valor 0.78, para KRb.
