# Interpolación consciente del polo para A_p — y una corrección a la ronda anterior

> ⚠️ **AVISO DE PREMISA (2026-08-20).** Este documento trata el
> **pseudopotencial de Fermi**, que **NO forma parte del Hamiltoniano de
> Aguilera-Fernández et al. 2015 / González-Férez et al. 2015** para Rb*-KRb.
> Su Ec. 1 es `H_ad = H_A + H_mol`, con KRb como **dipolo puntual**, no como
> centro de dispersión de contacto. Su contenido técnico **sigue siendo
> válido** —y sigue aplicando a la línea de perturbador NEUTRO
> (Aguilera-Fernández 2016)— pero **su uso para comparar con esos dos papers
> partía de una premisa equivocada**. Ver
> `docs/analysis_fig1_carga_dipolo_sin_fermi.md` §1.

**Fecha**: 2026-08-19
**Autor**: Javier Aguilera
**Relevancia**: Sustituye la interpolación lineal de `A_p` por una que respeta el
polo de la resonancia de forma p. En el proceso **invalida la conclusión
principal de `analysis_resonancia_onda_p.md`**.

## Resumen — y corrección

Se implementó la interpolación de `1/A_p` (opción `p_interpolation="inverse"`,
ahora por defecto) con justificación de rango efectivo y verificación empírica.

⚠️ **CORRECCIÓN A LA RONDA ANTERIOR.** Allí concluí que el pozo de −6472 GHz
«no es un resultado, es un artefacto de interpolación». **Eso era incorrecto.**
El punto R=565 a₀ sí caía en el hueco de malla, pero su vecino **R=566 a₀
remapea a R'=765.2, fuera del hueco**, y con un `A_p` leído directamente de la
tabla da **−7790 GHz** — más profundo aún. El pozo butterfly sobrevive con
datos tabulados reales; **no** depende de la interpolación.

Lo que realmente ocurre es lo que se sospechó al principio y descarté demasiado
pronto: **el pseudopotencial de Fermi de rango cero diverge en la resonancia**,
y ésa es la limitación física conocida, no un problema de nuestros datos.

## 1. Justificación: por qué 1/A_p es suave y A_p no

Para onda p, con `A_p(k) = -tan δ_p(k)/k³` (de modo que `A_p → a_p` cuando
`k→0`), la expansión de rango efectivo

```
k³ cot δ_p(k) = -1/a_p + ½ r_p k² + O(k⁴)
```

es **analítica en k²**. Como `1/A_p = -k³ cot δ_p`, se sigue directamente

```
1/A_p(k) = 1/a_p - ½ r_p k² - O(k⁴)
```

es decir, **1/A_p ES la función de rango efectivo**. Con forma Breit-Wigner
cerca de una resonancia aislada, `tan δ_p = (Γ/2)/(E_r - E)`:

```
A_p   = -(Γ/2) / [ k³ (E_r - E) ]      ->  POLO simple en E = E_r
1/A_p = -k³ (E_r - E) / (Γ/2)          ->  CERO simple, lineal en E
```

El polo de `A_p` es exactamente un cero simple de `1/A_p`. Interpolar
linealmente entre dos nodos que rodean un cero simple es correcto; hacerlo a
través de un polo no lo es.

### Verificación empírica sobre la propia tabla

Ajuste lineal en `ε = k²/2` sobre 20 nodos (R' = 723–792 a₀) que abarcan la
resonancia:

| función ajustada | R² |
|---|---|
| `A_p` vs ε | **0.2935** |
| `1/A_p` vs ε | **0.9883** |

El cero del ajuste global de `1/A_p` da el polo en **24.711 meV**, y la
interpolación entre los dos nodos que rodean el hueco lo sitúa en **24.732 meV**
— coinciden a 0.02 meV, y ambos mejoran la estimación anterior por punto medio
(24.82 meV). En nuestra coordenada, el polo está en **R = 562.77 a₀**.

## 2. Implementación

`ScatteringLengths(p_interpolation=...)`:

- `"inverse"` (**por defecto**): interpola `1/A_p` linealmente en la energía
  cinética local del electrón `ε = k²/2` —la variable en la que es analítica— y
  luego invierte. Devuelve `inf` si el interpolante cruza exactamente cero.
- `"linear"`: comportamiento anterior (interpola `A_p` en R'), conservado para
  poder comparar.

**A_s no se invierte nunca.** No tiene polo (`max|A_s| = 15.99 a₀`) pero **sí
cruza cero**: mínimo de Ramsauer en R' = 429 a₀, con `min|A_s| = 0.005693`. Ahí
`1/A_s` sería singular. Se interpola directamente en ambos modos.

Se añadió una comprobación en el constructor: si `rvsAP.dat` contuviera algún
`A_p = 0` exacto, lanza `ValueError` en vez de producir un infinito silencioso.
No hay ninguno (`min|A_p| = 324.9 a₀³` en R'=111).

## 3. Barrido fino alrededor de la resonancia

`A_p` con ambos métodos:

| R [a₀] | A_p lineal | A_p inverso | factor |
|---|---|---|---|
| 552 | +141 778 | +140 493 | 0.99 |
| 556 | +353 527 | +353 499 | 1.00 |
| 558 | +365 160 | +440 679 | 1.21 |
| 560 | +101 585 | +761 466 | 7.50 |
| **562** | −162 644 | **+2 746 566** | 16.89 |
| **564** | −427 529 | **−1 729 079** | 4.04 |
| 566 | −675 747 | −675 143 | 1.00 |
| 570 | −373 025 | −353 549 | 0.95 |
| 574 | −133 010 | −132 677 | 1.00 |

**Fuera del hueco los dos métodos coinciden al 1 %**; dentro, el método inverso
reproduce la divergencia (+2.7×10⁶ justo antes del polo, −1.7×10⁶ justo
después) mientras que el lineal la aplana absurdamente.

Energía del autovalor más bajo del bloque M_J=0:

| R [a₀] | E lineal [GHz] | E inverso [GHz] | ¿dato tabulado real? |
|---|---|---|---|
| 554 | −67.7 | −67.7 | sí |
| 556 | −76.9 | −59.5 | sí |
| 558 | −70.3 | −51.8 | **no (hueco)** |
| 560 | −82.0 | −563.5 | **no (hueco)** |
| 562 | −1914.1 | −209.8 | **no (hueco)** |
| 564 | −4980.7 | **−20 149.4** | **no (hueco)** |
| **566** | **−7790.1** | **−7783.1** | **sí** |
| 568 | −6393.5 | −6299.5 | sí |
| 570 | −4208.9 | −3989.0 | sí |
| 572 | −2236.8 | −2109.3 | sí |

### Respuesta a la pregunta del punto 3

**El pozo no se atenúa: se profundiza y se vuelve formalmente no acotado.**

- El valor anterior de −6472 GHz (en R=565, dentro del hueco) pasa a −20 149 GHz
  en R=564 con la interpolación correcta. La profundidad exacta depende de lo
  cerca que caiga la malla del polo en R=562.77 a₀, y **diverge al acercarse**.
- Pero, y esto es lo decisivo, **en R = 566, 568, 570, 572 a₀ —todos fuera del
  hueco, con `A_p` leído de la tabla— los dos métodos coinciden y dan
  −7790, −6300, −3989, −2109 GHz.** El pozo profundo está en los datos, no en
  la interpolación.
- La anchura tampoco cambia sustancialmente: la estructura sigue confinada a
  ~15–20 a₀, consistente con la FWHM de 13.2 a₀ ya medida.

## 4. Consecuencia: qué estaba mal en la conclusión anterior

`analysis_resonancia_onda_p.md` §6.2 dice: «la profundidad de −6472 GHz no es un
resultado, es un artefacto de interpolación». **Es falso.** Sólo era artefacto
el valor concreto en ese punto de malla. La existencia y el orden de magnitud
del pozo (miles de GHz) se reproducen con nodos tabulados genuinos.

La explicación correcta es la que se planteó como hipótesis de partida y
descarté demasiado rápido: **el pseudopotencial de Fermi de rango cero diverge
en la resonancia de forma**. `V ∝ A_p → ∞` no es una patología numérica, es la
señal de que la aproximación de rango cero deja de valer justo ahí. Es
exactamente el régimen para el que existen las correcciones de rango efectivo
(Omont) que hasta ahora hemos evitado implementar.

Que la Fig. 2 del paper no muestre un pozo de THz sigue sin explicarse. Ahora
las opciones son: (a) el paper usa correcciones de rango efectivo, (b) su
ventana de energía excluye las curvas butterfly, (c) su malla en R no resuelve
la región. Ya **no** vale la opción «es un artefacto nuestro de interpolación».

## 5. Estado de las tareas

| punto | estado |
|---|---|
| 1. justificación con fórmula | hecho (§1), con verificación empírica |
| 2. implementación con opción de vuelta atrás | hecho (§2) |
| 3. barrido fino y comparación | hecho (§3) |
| 4. curva BOP completa regenerada | **HECHO** (2026-08-19, ronda siguiente): 223.4 s, misma malla en R. Fuera de la resonancia los dos métodos coinciden a mediana 3.7×10⁻⁶ GHz y máximo 0.194 GHz. Ver `analysis_ventana_exclusion_resonancia.md` §1 |
| 5. tests | hecho, sin regresiones |

## 6. Archivos tocados

- `src/trimero/hamiltonians/fermi_krb.py` — `ScatteringLengths(p_interpolation=…)`,
  tablas en ε, comprobación de ceros de `A_p`.
- `plots/bop_curve_MJ0_linear.npz` — copia de la curva con el método antiguo.
- `docs/analysis_interpolacion_polo_Ap.md` — este documento.

## 7. Pendiente

1. ~~Regenerar la curva BOP completa con `p_interpolation="inverse"` y comparar
   con `plots/bop_curve_MJ0_linear.npz`~~ → **HECHO**, ver
   `analysis_ventana_exclusion_resonancia.md` §1.
2. ~~**Valorar en serio la corrección de rango efectivo (Omont)**~~ →
   **CERRADO POR DECISIÓN DEL USUARIO (2026-08-19), NO por implementación.**
   Ver la nota final de este documento.
3. Conseguir `A_p(k)` tabulado de una fuente independiente para acotar la
   magnitud cerca del polo. **Sigue abierto**, pero ya no bloquea nada: la
   región donde importaría está excluida de las comparaciones.

## 8. Nota final: cómo se cerró esto (2026-08-19)

**Decisión explícita del usuario: opción (b).** Se acepta la divergencia del
pseudopotencial de rango cero en la resonancia de forma p como **limitación
física conocida y documentada** del modelo, y **no se implementa la corrección
de rango efectivo de Omont (1977)**. La comparación cuantitativa con la Fig. 2
de González-Férez 2015 se acota a la región fuera de la resonancia, mediante una
ventana de exclusión explícita `R = 562.771 ± 2×13.187 a₀` calculada del propio
polo y de la FWHM de |A_p| (`ScatteringLengths.p_resonance_window()`).

Es decir: el §4 de este documento —«es exactamente el régimen para el que
existen las correcciones de rango efectivo (Omont)»— **sigue siendo cierto y
sigue sin implementarse**, ahora a propósito y por escrito, no por omisión.

El detalle del criterio, la comparación numérica, la tabla de características
fuera de la ventana y la nota de limitación redactada para un manuscrito están
en **`docs/analysis_ventana_exclusion_resonancia.md`**.

## Referencias

- Expansión de rango efectivo de onda p: `k³ cot δ_p = -1/a_p + ½ r_p k²`.
- Omont, A., *J. Physique* **38**, 1343 (1977) — correcciones más allá del rango cero.
- González-Férez, Sadeghpour & Schmelcher, *New J. Phys.* **17**, 013021 (2015).
- `docs/analysis_resonancia_onda_p.md` (⚠️ su §6.2 queda corregida por este documento),
  `docs/analysis_curva_bop_MJ0.md`.
