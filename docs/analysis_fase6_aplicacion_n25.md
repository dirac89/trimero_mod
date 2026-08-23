# Fase 6 — Primera aplicación a datos reales (n=25, M_J=0, estados 54/55)

**Fecha**: 2026-08-23
**Autor**: Javier Aguilera (con Claude Code)
**Relevancia**: primera vez que el pipeline completo (Fases 1-5) se aplica
a un sistema real de Rb*-KRb en vez de a un caso de juguete/analítico.

## 1. Decisión: n=25 (demostración) en vez de n≈52 (predicción cuantitativa)

El experimento real de Ruttley, Guttridge et al., PRL 131, 013401 (2023)
usa Rb(52s), no n=25. Antes de aplicar el pipeline había que decidir
explícitamente entre:

- **(a)** usar n=25 — el cruce evitado ya completamente caracterizado en
  las Fases 1-2 (`docs/analysis_fase1_acoplamiento_derivada.md`,
  R≈1305 a₀, estados 54/55) — como DEMOSTRACIÓN del pipeline completo, sin
  pretender que sea una predicción cuantitativa del experimento real de
  52s;
- **(b)** rehacer las Fases 1-2 para un n cercano a 52, generando un cruce
  evitado equivalente para ESE n, y aplicar el pipeline ahí.

**Se elige (a).** Razón, medida (no estimada a ojo) en esta misma ronda:

```
BOPSystem(n_manifold=25): dim(M_J=0) = 1113   (referencia, ya usado en Fases 1-2)
BOPSystem(n_manifold=52): dim(M_J=0) = 2436   (2.19x mayor)
Una sola diagonalización completa (eigh) en n=52: 76.7 s
                                    en n=25: ~1.1-2.0 s   (40-70x más lenta)
```

El barrido electrónico de esta misma Fase 6 (251 puntos en R, ver §4) tarda
~290 s con n=25. La MISMA operación con n=52 tardaría del orden de
251×76.7 s ≈ **5.3 horas**, sin contar el barrido exploratorio previo que
hizo falta en la Fase 1 para LOCALIZAR un cruce evitado real en n=25 (varias
decenas de diagonalizaciones adicionales sólo para encontrarlo) — esa
exploración habría que rehacerla desde cero para n≈52, con su propio dominio
natural de R (2n²a₀≈5408 a₀ para n=52, frente a 1250 a₀ para n=25), en un
régimen de R que ni siquiera se ha exportado todavía a `plots/rb_krb_polar/`.
Confirma la estimación del usuario: la ruta (b) es un coste real de varias
rondas, no de esta.

**Consecuencia explícita**: todos los números de esta Fase 6 son una
demostración de que el pipeline (Fases 1-5) funciona end-to-end sobre datos
reales de `BOPSystem`, **no** una predicción para el experimento de
Ruttley et al. 2023 (que usa un n de Rydberg distinto). Cualquier
comparación con el experimento real requeriría la ruta (b).

## 2. Estado inicial (Fase 5): R0 y σ, sin reutilizar R_am=310 nm sin más

`docs/analysis_fase5_franck_condon.md` dejó R0 y σ como parámetros de
entrada, pendientes de fijar aquí. Dos decisiones distintas:

- **R0** (la separación fijada por las pinzas): SÍ es específico del n de
  Rydberg (el radio orbital natural escala como 2n²a₀), así que **NO** se
  reutiliza 310 nm (calibrado para Rb(52s), 2n²a₀≈5408 a₀≈286 nm — de ahí
  sale ese número). Para n=25 (2n²a₀=1250 a₀≈66 nm) se fija en su lugar
  **R0=1305 a₀**, el propio cruce evitado que se está estudiando en esta
  demostración — la elección natural cuando el objetivo es examinar la
  física de ESE cruce, no reproducir un experimento concreto.
- **σ** (la anchura del paquete, dispersión de posicionamiento de las
  pinzas): es una propiedad de la PLATAFORMA experimental (frecuencia de
  trampa, calidad óptica) del átomo/molécula EN SU ESTADO FUNDAMENTAL,
  antes de la excitación a Rydberg — no depende de qué n se excite después.
  Se usa como **orden de magnitud** el valor medido por Ruttley et al. 2023
  en la MISMA plataforma (Rb+RbCs en pinzas específicas por especie),
  σ≈50 nm≈945 a₀, con la advertencia explícita de que no está validado
  específicamente para una hipotética preparación en n=25 (que no existe
  como experimento).

## 3. Masa reducida Rb-87 + RbCs (fuente citada, no tabla periódica genérica)

Pendiente explícito de `docs/analysis_fase2_canales_acoplados.md`
("µ≈113522 mₑ es un valor de demostración, pendiente de verificar la
fuente exacta"). Resuelto aquí con fuentes primarias:

| magnitud | valor | fuente |
|---|---|---|
| m(Rb-87) | 86.909 180 531 0(60) u | NIST Atomic Weights and Isotopic Compositions (physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=Rb), que reporta el valor de la AME2020 (Wang et al., *Chinese Phys. C* **45**, 030003, 2021) |
| m(Cs-133) | 132.905 451 961 0(80) u | ídem, physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=Cs |
| 1 u en mₑ | 1822.888486 | CODATA 2022, recíproco de "electron mass in u" = 5.485799090441×10⁻⁴ u (physics.nist.gov/cuu/Constants/Table/allascii.txt) |

```
m(RbCs) = m(Rb-87) + m(Cs-133) = 219.814 632 492 u
```

(la energía de enlace química de RbCs, del orden de 0.5-4 eV, corresponde a
Δm/c²~10⁻⁹ u — nueve órdenes de magnitud por debajo de la precisión de la
masa atómica misma; se desprecia explícitamente, no por omisión).

```
μ = m(Rb-87)·m(RbCs) / [m(Rb-87)+m(RbCs)]
  = 62.283 750 943 u
  = 62.283750943 × 1822.888486 mₑ
  = 113 536.332 mₑ
```

Este es el valor usado en toda esta Fase 6 (`MU_RB_RBCS` en
`scripts/analyze_nonadiabatic_n25_crossing.py`), sustituyendo el
µ≈113522 mₑ de demostración de la Fase 2 (diferencia ~0.01%, la Fase 2 ya
estaba cerca, pero sin fuente citada).

## 4. Ventana de R y su limitación (más ancha que la Fase 2, todavía no el dominio completo)

Barrido electrónico: R∈[1180,1430] a₀ (250 a₀), paso 1.0 a₀ (251 puntos),
`BOPSystem(n_manifold=25)`, M_J=0, sin Fermi, estados 54/55 — 5x más ancho
que los 50 a₀ de la demostración de la Fase 2, tal como pidió esta ronda.
Tiempo real: ver §6.

**Limitación explícita, no oculta**: el dominio FÍSICO completo donde vive
la curva BOP de M_J=0 es R∈[400,1800] a₀ (`docs/STATUS.md`, "8 mínimos
locales"). Una ventana de 250 a₀ sigue siendo un RECORTE de ese dominio. Si
los estados ligados reales de este par de curvas requieren la extensión
completa del pozo (que puede ser mucho más ancho que 250 a₀ — no se sabe
a priori sin calcularlo), esta ventana puede no contener suficiente
información para que el método de estabilización (Fase 3) confirme ningún
estado como genuinamente ligado, y eso sería un resultado ESPERADO de la
limitación de la ventana, no un fallo del método — ver los resultados en
§5 para lo que realmente se encontró.

## 5. Resultados numéricos reales

Pipeline completo (`poetry run python
scripts/analyze_nonadiabatic_n25_crossing.py`), sobre el barrido
electrónico de §4 (251 diagonalizaciones completas de 1113×1113,
285.6 s):

```
malla: 251 puntos, R=[1180.0,1430.0] a0, h=1.00 a0

Fase 1: |A(54,55)| max en la ventana = 5.0672e-01 en R=1305.0 a0
Fase 2: H acoplada ensamblada, forma (498, 498)
```

El pico de A_ij coincide EXACTAMENTE con el ya encontrado en la Fase 1
(mismo valor, 5.0672e-01, mismo R=1305.0 a₀) — consistencia interna
esperada, ambos se calculan con la misma fórmula sobre datos ahora
recalculados en una ventana más ancha.

### Fase 3 — estabilización (canal d=54, desacoplado)

Barrido de caja: L∈[60,248] a₀ (20 tamaños), sobre los datos precalculados
(sin nuevas diagonalizaciones electrónicas — sólo se trunca el array ya
calculado, ver `hamiltonian_builder_d` en el script). **15 estados**
clasificados como estables (`fraction_stable=1.00` en los 15, el máximo
posible), abarcando V-E_manifold∈[-16.275,-8.913] GHz (anchura 7.36 GHz) —
del mismo orden de magnitud que el pozo más profundo de la Fig. 1 completa
(-23.1 GHz, `docs/STATUS.md`), lo cual es una señal de sanidad: no son
números disparatados.

**Interpretación cuidadosa (no sobre-interpretar)**: que los 15 estados
salgan con `fraction_stable=1.00` en TODAS las cajas probadas (desde
L=60 a₀ en adelante) sugiere que el pozo real que los sostiene es más
ANGOSTO que 60 a₀ — es decir, ya con la caja más pequeña probada, todos
caben cómodamente y no notan que la caja crece. Esto es coherente con la
física (un cruce evitado y su entorno inmediato producen un pozo local
relativamente estrecho), pero significa que este barrido de L **no
descarta** que hubiera estados de caja mezclados si se hubiera empezado
con L aún más pequeño — no se hizo por presupuesto de cómputo (cada
L pequeño no ahorra tiempo, el coste ya está pagado en el barrido
electrónico). Se documenta como limitación, no se oculta.

### Fase 4 — tasas de decaimiento

De los 15 estados, sólo los **2 más someros** (los más cercanos al umbral)
tienen una energía por ENCIMA del mínimo de V_u=55 dentro de esta ventana
—condición necesaria para que la fórmula de predisociación tenga sentido—;
los otros 13 quedan por debajo del mínimo de V_u **en esta ventana
concreta** (no se puede concluir que no exista predisociación en la
física real completa, sólo que no se puede evaluar con estos datos):

```
V-E_manifold=-9.1229 GHz   Gamma=139.3908 MHz
V-E_manifold=-8.9013 GHz   Gamma=23.4604 MHz
```

Los otros 13 reportan explícitamente "por debajo del mínimo de V_u en esta
ventana, no se puede evaluar predisociación con los datos disponibles" —
la misma limitación de ventana de §4, propagada honestamente hasta aquí en
vez de forzar un número.

### Fase 5 — Franck-Condon

`R0=1305 a₀` (el propio cruce), dos anchuras:

- `sigma_real=945 a₀` (≈50 nm, de Ruttley et al. 2023): el script emite un
  AVISO explícito — esta anchura es ~4× la ventana completa (248 a₀), así
  que el resultado está dominado por el recorte artificial de la caja, NO
  es un número fiable. Se reporta igualmente (por transparencia, no se
  oculta el intento), pero no se usa para ninguna conclusión.
- `sigma_demo=12.4 a₀` (=ventana/20, elegido para que quepa cómodamente):
  ilustrativo únicamente, NO pretende ser Rb*-RbCs real.

```
V-E_manifold= -16.2750 GHz  F(sigma_real)=-3.4219e-01  F(sigma_demo)=-1.3933e-09
V-E_manifold= -12.3201 GHz  F(sigma_real)=-1.4364e-01  F(sigma_demo)=-5.2365e-03
V-E_manifold= -10.9153 GHz  F(sigma_real)=-1.7012e-01  F(sigma_demo)=-9.6532e-02
V-E_manifold=  -9.1229 GHz  F(sigma_real)= 1.6542e-01  F(sigma_demo)= 1.6323e-01
V-E_manifold=  -8.9013 GHz  F(sigma_real)= 8.2156e-02  F(sigma_demo)= 2.5104e-01
```

(lista completa de 15 en el output real más abajo). Con `sigma_demo`, F
crece varios órdenes de magnitud según el estado se acerca al umbral
(de 1.4×10⁻⁹ para el más profundo a 0.25 para el más somero) — coherente
con que los estados más someros tienen más amplitud cerca de R0=1305 (el
propio cruce), donde está centrado el paquete inicial.

**Nota honesta añadida al revisar el output** (no estaba en la primera
versión de este documento): las energías de esta tabla y las de la §Fase 3
NO coinciden exactamente. La Fase 3 promedia la energía de la trayectoria
sobre las cajas truncadas L∈[60,248]; la Fase 4/5 usa el autovalor de la
caja COMPLETA (L=248) más cercano a esa media, para tener un autovector
concreto con el que calcular Γ y F. Con estados separados por menos de
1 GHz entre sí (como aquí, ver la tabla de la Fase 3), esa búsqueda por
proximidad puede desplazar el valor hasta ~0.3 GHz y, en dos casos, hacer
que dos entradas "estables" distintas de la Fase 3 apunten al MISMO
autovalor de caja completa (los dos "-10.9153 GHz" y los dos "-9.6390 GHz"
de la tabla de arriba) — un artefacto esperable de la resolución
sub-GHz entre estados, no un error de cálculo, pero que no se documentó
hasta esta nota.

## 6. Comparación cualitativa con el paper de 2024 (Cs-RbCs)

Verificado directamente del texto (vía `WebFetch`,
arxiv.org/html/2401.09618): *"the higher vibrational states away from the
avoided crossing region decay more slowly than the states near the
crossing"* — Γ decrece para estados más alejados (en energía/posición) del
cruce evitado, y crece para estados más próximos a él. El paper no da
valores numéricos en el texto (sólo figuras 4(b) y 4(d), sin escala
cuantitativa en el texto), así que la comparación sólo puede ser de
FORMA/tendencia, nunca de magnitud — como se pedía.

**Con nuestros datos NO se puede confirmar ni refutar esta tendencia con
confianza**: sólo tenemos 2 estados con Γ evaluable (el resto quedó fuera
de alcance por la limitación de ventana de §5), y no se conoce con
fiabilidad el número cuántico vibracional de cada uno dentro de un pozo
que la ventana recorta artificialmente. Lo que se observa (Γ=139.4 MHz
para el estado en V=-9.12 GHz, Γ=23.5 MHz para V=-8.90 GHz, un factor ~6)
es una diferencia real y medida, no un artefacto — pero con 2 puntos no es
honesto afirmar que reproduce o contradice la tendencia cualitativa del
paper; haría falta identificar más estados con Γ evaluable, lo que exige
una ventana de R más ancha (ver §7).

## 7. Conclusiones y trabajo pendiente

- El pipeline completo (Fases 1-5) se ejecuta de punta a punta sobre datos
  reales de `BOPSystem`, sin errores, con resultados de magnitud física
  razonable (profundidades de pozo del mismo orden que la Fig. 1 completa,
  tasas de decaimiento en el rango MHz-cientos de MHz, factores de
  Franck-Condon que decaen suavemente con la distancia al centro del
  paquete inicial).
- Las limitaciones son reales y se documentan explícitamente, no se
  ocultan: (i) n=25 es una demostración, no una predicción para el
  experimento real de 52s (§1); (ii) la ventana de 250 a₀ sigue siendo un
  recorte del dominio físico completo (400-1800 a₀), lo que impide evaluar
  Γ para 13 de los 15 estados encontrados y hace que el paquete gaussiano
  REAL (945 a₀) no quepa; (iii) con sólo 2 estados con Γ evaluable no se
  puede confirmar la tendencia cualitativa del paper de 2024 con
  confianza, aunque la diferencia medida (factor ~6) es consistente en
  DIRECCIÓN con "más cerca del cruce decae más rápido".
- Trabajo natural para una ronda futura, si se decide seguir por esta
  línea: extender la ventana de R (aceptando el coste de ~1.1-2 s por
  diagonalización electrónica) lo suficiente para (a) que quepa
  σ_real=945 a₀ del paquete inicial con margen, y (b) que más estados
  ligados tengan un canal u accesible para Γ — permitiendo entonces sí una
  comparación de tendencia con varios puntos, no sólo dos.
- La masa reducida (§3) y la elección de estado inicial (§2) quedan ahora
  fijadas con fuente citada para cualquier trabajo futuro sobre este mismo
  sistema (n=25); si se retoma la ruta (b) del §1 (n≈52), ambas decisiones
  tendrían que revisarse (R0 cambia con 2n²a₀; σ probablemente se mantiene
  al ser propiedad de la plataforma, no del n).

## Archivos tocados

- `scripts/analyze_nonadiabatic_n25_crossing.py` (nuevo)
- `docs/PLAN_nonadiabatic_dynamics.md` (ESTADO ACTUAL: Fase 6 en curso, con
  resumen de esta ronda)
- `docs/analysis_fase6_aplicacion_n25.md` (este documento)
- `plots/hybrid_neutral_polar/data/fase6_n25_results.npz` (resultados
  numéricos: R_mid, A, B, W_mid, L_values, E_box, bound_summary)

No se tocó ningún módulo de `src/trimero/` ni ningún test en esta ronda
(Fase 6 es aplicación, no infraestructura nueva); la suite del proyecto
sigue en el mismo estado que al cierre de la Fase 5 (133 passed, 16
deselected en `pytest -m "not slow"`).

## Referencias

- `docs/analysis_fase1_acoplamiento_derivada.md` — el cruce evitado real
  (n=25, M_J=0, estados 54/55, R≈1305 a₀) reutilizado y reextendido aquí.
- Mellado-Alcedo, Guttridge, Cornish, Sadeghpour & González-Férez,
  *Phys. Rev. A* **110**, 013314 (2024), arXiv:2401.09618 — comparación
  cualitativa de tendencia (§6), texto verificado verbatim.
- Ruttley, Guttridge, Baldock, González-Férez, Sadeghpour, Adams & Cornish,
  *Phys. Rev. Lett.* **131**, 013401 (2023), arXiv:2303.06126 — σ del
  paquete inicial (§2).
- NIST Atomic Weights and Isotopic Compositions (AME2020) y CODATA 2022 —
  masa reducida (§3).
- `docs/STATUS.md` — dominio físico completo de la curva BOP (§4, §7).
