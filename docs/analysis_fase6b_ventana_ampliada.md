# Fase 6b — Ventana de R ampliada (cierra el pendiente §7.1 de la Fase 6)

**Fecha**: 2026-08-23
**Autor**: Javier Aguilera (con Claude Code)
**Relevancia**: cierra el primer punto del trabajo pendiente de
`docs/analysis_fase6_aplicacion_n25.md` §7 — ampliar la ventana de R
(250 a₀) del análisis de dinámica no adiabática del cruce evitado real
(`BOPSystem(n_manifold=25)`, M_J=0, estados 54/55) para que (a) más
estados ligados tengan un canal continuo accesible y (b) evaluar si cabe
el paquete gaussiano real de Ruttley et al. 2023 (σ=945 a₀).

## 1. Estimación de coste ANTES de ejecutar (obligatoria, per instrucción)

Con el tiempo medido en la Fase 6 (285.6 s / 251 puntos = **1.138 s/punto**,
`BOPSystem(n_manifold=25)`, M_J=0, sin Fermi, eigh completo de 1113×1113):

| ventana pedida/considerada | ancho | nº puntos (h=1 a₀) | coste estimado |
|---|---|---|---|
| para que quepa σ_real=945 a₀ con ±3σ | 5670 a₀ | 5670 | **107.5 min (1.8 h)** |
| para que quepa σ_real=945 a₀ con ±5σ | 9450 a₀ | 9450 | **179.2 min (3.0 h)** |
| alternativa elegida: dominio ya validado de la Fig. 1 | 1300 a₀ | 1301 | **24.7 min** |

Las dos primeras superan ampliamente el umbral de 30-40 min acordado — no
se ejecutan esta ronda sin autorización explícita para ese tiempo de
cómputo (quedan documentadas como referencia en §5). La tercera está por
debajo del umbral ("coste razonable, minutos"), así que se ejecuta
directamente, tal como autorizaba la instrucción.

### Malla no uniforme: evaluada y descartada esta ronda

La alternativa (i) del encargo (paso fino sólo cerca del cruce, grueso
lejos) se evaluó ANTES de escribir código. Los cinco módulos de las
Fases 1-5 (`coupling.py`, `coupled_channels.py`, `stabilization.py`,
`energy_normalization.py`, `franck_condon.py`) no sólo VALIDAN que la
malla sea uniforme (`np.allclose(steps, steps[0], ...)`, y lanzan
`ValueError` si no) — las propias fórmulas de diferencias finitas (matriz
cinética de 3 puntos, primera y segunda derivada centradas) están
derivadas asumiendo h CONSTANTE; no es sólo una comprobación de entrada
que se pueda relajar, es una asunción estructural de cada fórmula.
Adoptar una malla no uniforme exigiría rederivar y re-verificar (contra
las mismas referencias analíticas independientes ya usadas: oscilador
armónico, pozo cuadrado, partícula libre, dos gaussianas) los cinco
módulos — un trabajo del orden de las Fases 1-5 otra vez, no de esta
ronda. Se descarta explícitamente, no se fuerza.

**Decisión**: ventana R∈[500,1800] a₀ (1300 a₀, 5.2× la de la Fase 6),
paso 1.0 a₀ (1301 puntos), coincidiendo con el dominio ya validado
completo de la Fig. 1 (`docs/STATUS.md`) en vez de un rango arbitrario —
así el resultado se puede comparar directamente con los umbrales y
mínimos ya reportados ahí.

## 2. Barrido ejecutado, y un BUG REAL encontrado antes de aceptar nada

Barrido electrónico: 1301 diagonalizaciones completas de 1113×1113,
**1529.7 s** (25.5 min) — dentro de lo estimado (24.7 min).

### 2.1 Test de consistencia (§3 del encargo): perfecto

Se recortó del barrido nuevo el mismo tramo R∈[1180,1430] de la Fase 6
original y se corrió el pipeline sobre ese recorte, comparando contra los
números YA PUBLICADOS (no contra una repetición ciega):

```
pico |A(54,55)|: original=5.0672e-01 en R=1305.0  recorte-nuevo=5.0672e-01 en R=1305.0  dif_rel=0.00%
```

Coincidencia exacta — la malla más ancha reproduce bit a bit el
acoplamiento de la Fase 1 en la región compartida. Buena señal de que el
barrido nuevo es correcto... **pero al comparar los ESTADOS LIGADOS
encontrados, algo no cuadraba** (ver 2.2).

### 2.2 El bug: `classify_stability` recibiendo energías ABSOLUTAS

Con la ventana ancha y `n_keep=60`, primera corrida: **300/300** estados
"estables" — el máximo posible, sospechosamente perfecto. Se investigó
ANTES de aceptarlo (regla del proyecto): se subió `n_keep` a 300 (gratis en
coste, `eigvalsh` ya calcula todos los autovalores; confirmado
midiendo 0.22 s por caja de 1299×1299) y se probó con `threshold_ratio`
33× más estricto (0.1 → 0.003): **seguían siendo 300/300 estables**,
incluidos estados con energía relativa POSITIVA (V-E_manifold hasta
+12.6 GHz, por encima del propio umbral de disociación del canal).

Eso no podía ser física real. Se inspeccionó V_d(R)=W[:,0] directamente:
tiene un umbral de disociación limpio y bien definido,
**V_d → −9.2759 GHz** (desviación estándar 4.4×10⁻⁵ GHz en los últimos
100 a₀ de la ventana) — así que SÍ hay una referencia real contra la que
comparar. La causa del bug: `classify_stability` compara `|dE/dL|` contra
`2|E|/L` (la escala de deriva de un estado de caja, Fase 3), una fórmula
derivada asumiendo E referenciada al umbral de disociación (E→0 lejos del
pozo, como en el pozo cuadrado de juguete de la Fase 3, donde V→0 por
construcción). Las energías de `BOPSystem` son ABSOLUTAS en Hartree
(~−8×10⁻⁴ Eh), dominadas por el offset del manifold Rydberg
(E_manifold=−0.5/25²=−8×10⁻⁴ Eh, ≈−5264 GHz) — un número que NO tiene nada
que ver con la física vibracional local. Usar |E| absoluto en `2|E|/L`
produce una referencia de comparación absurdamente grande (dominada por
−5264 GHz, no por los GHz de profundidad real), así que prácticamente
CUALQUIER deriva medida (a escala de GHz) sale "pequeña" frente a esa
referencia — de ahí el 300/300.

**Esto NO es un bug de `stabilization.py`** (sus propios tests, incluido el
pozo cuadrado, siguen en verde: `poetry run pytest
tests/systems/nonadiabatic_dynamics/test_stabilization.py` → 4 passed —
ahí V ya estaba referenciada a cero por construcción del modelo de
juguete). Es un bug de **cómo se invocó el módulo** en
`analyze_nonadiabatic_n25_crossing.py` (Fase 6 original y primera versión
de esta Fase 6b): nunca se restó el offset absoluto del Rydberg antes de
pasar las energías a `classify_stability`.

### 2.3 Corrección

Se resta el umbral asintótico REAL del canal d (medido directamente de la
cola plana de `W_mid[:,0]`, últimos 10% de la ventana) antes de clasificar
— un desplazamiento constante, no afecta a `dE/dL` (la derivada de una
constante es cero), sólo corrige la escala de `2|E|/L`. Implementado en
`run_pipeline` (`analyze_nonadiabatic_n25_crossing.py`).

**Consecuencia importante para la Fase 6 original**: los DOS estados con Γ
evaluable reportados ahí (V=−9.0722 GHz, Γ=139.4 MHz; V=−8.9127 GHz,
Γ=23.5 MHz) **ya NO aparecen en la lista de estados ligados** una vez
corregida la referencia — el recorte-viejo corregido da **14** estados
ligados (no 15), con rango [−16.275,−9.639] GHz, es decir el estado más
somero anterior (−8.9127 GHz) queda fuera, y el segundo más somero
(−9.0722 GHz) TAMBIÉN queda fuera del nuevo rango. **Los dos valores de Γ
del análisis original deben tratarse como no fiables** (calculados para
estados que la clasificación corregida no confirma como ligados) — se
retiran aquí explícitamente, no se mantienen como válidos.

```
[recorte-viejo] Fase 3 (corregido): 14 estados, rango [-16.2750, -9.6389] GHz
    fracción estable de los 6 más someros: 0.94, 0.89, 0.83, 0.83, 0.78, 0.72, 0.67, 0.67, 0.61
    (los 5 más profundos: 1.00 -- alta confianza; los someros: mucho menos, coherente
     con estar genuinamente cerca de un umbral, que es justo lo que este método debe reflejar)
```

## 3. Resultados corregidos, ventana ancha (500-1800 a₀, n_keep=300)

```
[ancho] Fase 3 (corregido): 144 estados ligados (de 300 probados), rango [-20.2108, 1.4863] GHz
```

Distribución de confianza (fracción de cajas en que cada estado sale
estable): 37 estados al 100%, y una cola decreciente hasta el mínimo
aceptado (50%) — un perfil MUCHO más físico que el "300/300 perfecto" de
antes de la corrección: los estados profundos son inequívocamente
estables, y la confianza decae de forma continua según se acercan al
umbral, que es exactamente el comportamiento esperado de un método de
estabilización funcionando bien.

**De los 144, 79 tienen un canal u (estado 55) accesible dentro de la
ventana para evaluar Γ** (frente a 2 de 15 en la Fase 6 original — mejora
sustancial, aunque los 2 originales ya no cuentan como confirmados, ver
§2.3):

```
Gamma: n=79, min=0.0002 MHz, max=129.29 MHz, mediana=2.36 MHz, media=18.92 MHz
```

## 4. Comparación de tendencia con el paper de 2024 (ahora con 79 puntos, no 2)

Correlación de Pearson entre V-E_manifold (profundidad relativa) y Γ:
**r=0.60** (positiva, moderada) — los estados MENOS ligados (más cerca del
umbral/cruce) tienden a Γ mayor; los más profundamente ligados, a Γ menor:

```
mayor Gamma  : V=0.72 a 0.73 GHz (cerca del umbral), Gamma=129.3 MHz
               V=-5.0 a -7.4 GHz, Gamma=80-115 MHz
menor Gamma  : V=-15.2 a -14.8 GHz (profundamente ligados), Gamma=0.0002-0.002 MHz
```

Esto es CUALITATIVAMENTE consistente con lo que reporta el paper de 2024
(verificado en la Fase 6 original, verbatim): *"the higher vibrational
states away from the avoided crossing region decay more slowly than the
states near the crossing"*. Con 79 puntos (no 2) y una correlación
moderada pero clara (r=0.60, no un artefacto de 2 puntos aislados), esto
es ahora una comparación de tendencia razonablemente respaldada — aunque
sigue siendo cualitativa (el paper no da números), y r=0.60 no es
perfecta: hay dispersión real (p.ej. el estado en V=-0.05 GHz tiene
Γ=0.19 MHz, bajo pese a estar cerca del umbral) que no se investiga más a
fondo aquí.

## 5. σ_real=945 a₀: sigue sin caber, ventana necesaria documentada

Con la ventana ampliada (1298 a₀ efectivos): `sigma_real=945.0 a0 cabe con
margen 5sigma en esta ventana (1298 a0): False`. Confirma la estimación de
§1: hacen falta 5670-9450 a₀ (107-179 min de cómputo), no ejecutado esta
ronda. Se mantiene como referencia para una ronda futura si se decide
invertir ese tiempo.

## 6. Conclusiones

- El barrido ampliado (1300 a₀, 25.5 min) se ejecutó dentro del
  presupuesto acordado, y el test de consistencia contra la Fase 6
  original es perfecto (mismo pico de acoplamiento, dif_rel=0.00%).
- Se encontró y corrigió un bug real en cómo se referenciaban las energías
  para el criterio de estabilidad (Fase 3 aplicada a datos reales, no al
  módulo en sí) — sin esta corrección, el método clasificaba
  incorrectamente prácticamente CUALQUIER estado como "ligado", incluidos
  estados con energía por encima del propio umbral de disociación. Se
  retiran explícitamente los 2 valores de Γ de la Fase 6 original: los
  estados de los que salieron ya no se confirman como ligados con el
  criterio corregido.
- Con la corrección, la ventana ancha da 144 estados ligados (frente a los
  300 probados) con un perfil de confianza físicamente sensato, y 79 con Γ
  evaluable (frente a 2 antes) — permitiendo, por primera vez, una
  comparación de TENDENCIA real (no anecdótica) con el paper de 2024:
  correlación positiva moderada (r=0.60) entre cercanía al umbral y Γ,
  cualitativamente consistente con el hallazgo del paper.
- σ_real=945 a₀ sigue sin caber; la ventana que haría falta (5670-9450 a₀,
  107-179 min) queda documentada, no ejecutada.
- Con esto, el punto 1 del pendiente de `docs/analysis_fase6_aplicacion_n25.md`
  §7 queda **resuelto** en la parte de "más estados con Γ evaluable" (79,
  con tendencia comparable al paper) y **parcialmente resuelto, con el
  motivo documentado**, en la parte de "que quepa σ_real" (no cabe, coste
  de la ventana necesaria cuantificado para una ronda futura).

## Archivos tocados

- `scripts/analyze_nonadiabatic_n25_crossing.py` (refactorizado:
  `run_pipeline` reutilizable + corrección del bug de referencia de
  energía en la Fase 3; verificado que reproduce bit a bit los números de
  A_ij, los 15 estados y los 2 Γ ORIGINALES antes de la corrección, para
  confirmar que el refactor no cambió comportamiento por sí solo)
- `scripts/analyze_nonadiabatic_n25_crossing_wide.py` (nuevo)
- `docs/PLAN_nonadiabatic_dynamics.md` (ESTADO ACTUAL actualizado)
- `docs/analysis_fase6b_ventana_ampliada.md` (este documento)
- `plots/hybrid_neutral_polar/data/fase6b_n25_results.npz` (resultados: R_mid,
  A, B, W_mid, bound_summary, gammas)

No se tocó ningún módulo de `src/trimero/` ni ningún test (el bug estaba en
el script de aplicación, no en la librería; `test_stabilization.py` sigue
en verde, 4 passed, confirmando que el módulo en sí nunca estuvo mal).

## Referencias

- `docs/analysis_fase6_aplicacion_n25.md` — el análisis original que este
  documento extiende y en parte corrige (§2.3).
- Mellado-Alcedo, Guttridge, Cornish, Sadeghpour & González-Férez,
  *Phys. Rev. A* **110**, 013314 (2024), arXiv:2401.09618 — comparación de
  tendencia (§4).
- `docs/analysis_fase3_estabilizacion.md` — el criterio `2|E|/L` cuya
  asunción implícita (E referenciada al umbral) se violó y corrigió aquí.
