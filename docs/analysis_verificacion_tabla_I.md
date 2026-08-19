# Verificación contra la Tabla I de González-Férez 2015: defectos cuánticos de Rb

**Fecha**: 2026-08-19
**Autor**: Javier Aguilera
**Relevancia**: Primera comparación directa de nuestros números con valores
publicados. No requiere construir base ni diagonalizar: sólo `Atom.E_Rb()`.
Aísla un sistemático real en la serie s antes de que contamine cualquier
comparación posterior de curvas BOP.

## Resumen

De los 8 niveles de la Tabla I, **6 concuerdan a ±0.004 GHz (≤0.002 %)**. Los
dos de la serie s fallan por −0.299 GHz (27s) y +0.265 GHz (28s). Se demuestra
que ambos desajustes piden **el mismo desplazamiento constante en δ₀(ns)**
(+6.19×10⁻⁴), lo que descarta δ₂, la dependencia en n y la masa reducida.
**Es una discrepancia real, no redondeo.**

## 1. Método

Conversión CODATA 2018: `1 E_h = 6.579683920502×10⁶ GHz`.

⚠️ **No** se usa `EhtoGHz = 6.579683920729e9` del código legado
(`hamiltonians/trimer.py`): ese factor es en realidad Hartree→**MHz** mal
etiquetado, hallazgo previo de la sesión de refactor.

Referencia: `E(24, l=3) = -8.680555555556×10⁻⁴ E_h` (manifold cuasi-degenerado,
hidrogenoide, defecto cuántico despreciable para l≥3).

## 2. Tabla I — comparación

`ΔE_nl = |E_nl − E_24,3|` en GHz.

| n | l | E_nl [E_h] | n* efectivo | nuestro [GHz] | paper [GHz] | dif [GHz] | rel |
|---|---|---|---|---|---|---|---|
| 25 | 3 | −8.000000000000000e−04 | 25.0000000 | 447.78404 | 447.78 | +0.004 | 0.001 % |
| 28 | 0 | −8.084804394978180e−04 | 24.8685374 | 391.98543 | 391.72 | **+0.265** | **0.068 %** |
| 26 | 2 | −8.226318540067425e−04 | 24.6537076 | 298.87360 | 298.87 | +0.004 | 0.001 % |
| 27 | 1 | −8.431955228646972e−04 | 24.3512274 | 163.57116 | 163.57 | +0.001 | 0.001 % |
| 24 | 3 | −8.680555555555555e−04 | 24.0000000 | 0.00000 | 0.00 | +0.000 | — |
| 27 | 0 | −8.776457784040770e−04 | 23.8685133 | 63.10064 | 63.40 | **−0.299** | **0.472 %** |
| 25 | 2 | −8.936519942688221e−04 | 23.6537928 | 168.41648 | 168.42 | −0.004 | 0.002 % |
| 26 | 1 | −9.169637842916073e−04 | 23.3511842 | 321.80069 | 321.80 | +0.001 | 0.000 % |

**Agrupado por serie**:

| serie | diferencias [GHz] |
|---|---|
| l=3 (hidrogenoide) | +0.004 |
| l=1 (p) | +0.001, +0.001 |
| l=2 (d) | +0.004, −0.004 |
| **l=0 (s)** | **+0.265, −0.299** |

El fallo está **confinado a l=0**.

## 3. Procedencia de los defectos cuánticos de `systems/atom.py`

Los coeficientes son la **media aritmética de estructura fina** de valores
j-resueltos de la literatura estándar de Rb. Comprobado aritméticamente:

| serie | δ₀ en `atom.py` | media j-resuelta | dif δ₀ | δ₂ código | media δ₂ | dif δ₂ |
|---|---|---|---|---|---|---|
| s | 3.1311804 | 3.1311804 | +0.0000000 | 0.1745312 | 0.1784 | **−0.00387** |
| p | 2.6482793 | 2.6482793 | +0.0000000 | 0.2925324 | 0.2925 | +0.00003 |
| d | 1.3472787 | 1.3472785 | +0.0000002 | −0.5994376 | −0.5997 | +0.00026 |

Valores j-resueltos usados:

| estado | δ₀ | δ₂ | fuente |
|---|---|---|---|
| ns₁/₂ | 3.1311804 | 0.1784 | Li, Mourachko, Noel & Gallagher, PRA **67**, 052502 (2003) |
| np₁/₂ | 2.6548849 | 0.2900 | Han, Jamil, Norum, Tanner & Gallagher, PRA **74**, 054502 (2006) |
| np₃/₂ | 2.6416737 | 0.2950 | ídem |
| nd₃/₂ | 1.3480948 | −0.6054 | Li+ 2003 |
| nd₅/₂ | 1.3464622 | −0.5940 | Li+ 2003 |

- p y d coinciden con la media a **2×10⁻⁷**: no es casualidad, es la procedencia.
- `δ₀(ns) = 3.1311804` coincide **exacto** con Li+ 2003 (l=0 no tiene
  estructura fina que promediar).
- `δ₂(ns) = 0.1745312` **no** coincide con Li+ 2003 (0.1784), −2.2 %.

## 4. El desajuste NO es δ₂

Sensibilidad de la serie s a δ₂:

```
  n=27: δ₂=0.1745312 ->  63.10064 GHz ;  δ₂=0.1784 ->  63.10392 GHz ;  cambio = +0.0033 GHz  (falta +0.299)
  n=28: δ₂=0.1745312 -> 391.98543 GHz ;  δ₂=0.1784 -> 391.98276 GHz ;  cambio = -0.0027 GHz  (falta +0.265)
```

δ₂ mueve **0.003 GHz**: dos órdenes de magnitud por debajo del desajuste.
Además el signo del cambio en n=28 va en dirección contraria a la necesaria.
Corregir δ₂ al valor de Li+ 2003 sería correcto por consistencia bibliográfica,
pero **no resuelve nada aquí**.

## 5. El desajuste ES δ₀ — evidencia

Ajustando independientemente el δ₀(ns) que reproduciría cada nivel del paper:

| n | δ₀ requerido | desplazamiento vs 3.1311804 |
|---|---|---|
| 27 | 3.1317990 | +0.0006186 |
| 28 | 3.1318008 | +0.0006204 |

**Los dos niveles piden el mismo desplazamiento constante**, coincidentes entre
sí a 1.8×10⁻⁶. Un desplazamiento constante e **independiente de n** en δ₀ es
exactamente la firma de un valor de δ₀ distinto; un error en δ₂ o en la
expansión de Rydberg-Ritz daría un desplazamiento que varía con n.

El paper corresponde a **δ₀(ns) ≈ 3.13180** frente a nuestro 3.1311804.

### Controles descartados

| Hipótesis | Descartada porque |
|---|---|
| Redondeo de la tabla publicada | 0.299 GHz sobre valores dados a 2 decimales; y sería aleatorio, no sistemático en una sola serie |
| δ₂(ns) | mueve 0.003 GHz (§4) |
| Dependencia en n | el desplazamiento requerido es constante (§5) |
| Masa reducida (⁸⁷Rb vs Rydberg infinito) | `R_M/R_∞ = 0.9999936879` reescala **todas** las energías por igual: mueve 0.0004 GHz y afectaría también a p, d y f, que sí concuerdan |
| Conversión de unidades | l=1,2,3 usan la misma conversión y concuerdan a 0.001 % |

### Hipótesis sobre el origen (SIN VERIFICAR)

⚠️ No se ha podido consultar la referencia [20] del paper. Marinescu, Sadeghpour
& Dalgarno, PRA **49**, 982 (1994) es un trabajo de **potencial modelo**, no un
ajuste Rydberg-Ritz. Los autovalores de un potencial modelo difieren de un
ajuste Rydberg-Ritz sobre todo en los canales más penetrantes en el core, y
`l=0` es el más penetrante de todos, mientras que `l≥1` penetran mucho menos.
Que el desajuste esté **confinado exactamente a l=0** encaja con esa
explicación. Queda como hipótesis a comprobar contra la fuente.

## 6. Impacto sobre nuestros resultados

- Magnitud: 0.299 GHz = 0.47 % del espaciado 27s–manifold, 5.2×10⁻⁵ de la
  energía absoluta del 27s (−5775.6 GHz).
- Efecto: desplaza la posición del 27s respecto del manifold n=24, y con ella
  los cruces evitados 27s–manifold, en esa cantidad.
- Frente a los desplazamientos de decenas–centenares de GHz calculados en
  `analysis_pseudopotencial_fermi_krb.md`, es un sistemático **pequeño pero
  registrable**. No invalida nada; hay que tenerlo en cuenta al comparar
  posiciones de resonancia.

## 7. Límites libres de campo con excitación rotacional

Verificación trivial hecha explícita: `N(N+1)` para N=5 es `5·6 = 30` ✓ y para
N=6 es `6·7 = 42` ✓.

Con `B(KRb) = 1.114 GHz`:

| límite | nuestro [GHz] | con el valor del paper [GHz] |
|---|---|---|
| ΔE₂₇ₛ + 30·B | 63.10064 + 33.420 = **96.5206** | 63.40 + 33.420 = 96.8200 |
| ΔE₂₇ₛ + 42·B | 63.10064 + 46.788 = **109.8886** | 63.40 + 46.788 = 110.1880 |

La diferencia es la misma −0.299 GHz de §2, arrastrada íntegra: el término
`N(N+1)·B` no añade error propio y `B` es un input común a ambos cálculos.

## 8. Qué más sería verificable barato

⚠️ Esta sección clasifica **candidatos propuestos por el usuario**, no un
inventario del paper: no se dispone de su texto completo, sólo de los ocho
valores de la Tabla I.

| Dato | ¿Verificable ya? | Qué haría falta |
|---|---|---|
| `B(KRb)=1.114 GHz`, `d(KRb)=0.566 D` | **No** | son *inputs* nuestros, no un chequeo independiente |
| Dimensiones de base (568, 49, 27832, 1016, 59 bloques) | **Sí, ya hecho** | verificado en rondas anteriores |
| Convergencia en N: `E(N≤6) − E(N≤8) < 2×10⁻⁶` rel. | **Sí, barato** | `CoupledBasis(N_max=8)` + rediagonalizar M_J=0, ~1 min. Única afirmación cuantitativa del paper ya recogida en el doc de diseño y aún sin comprobar |
| ⟨cosθ⟩ máximo (orientación) | **Sí, barato** | `cos_theta_element` ya verificado; ⟨ψ\|cosθ_d\|ψ⟩ sobre un autovector son ~10 líneas. ⟨cos²θ⟩ (alineamiento) también, vía `gaunt` con Y₂₀. Requiere saber a qué estado y a qué R se refiere el valor publicado |
| Cruce evitado ~1.8 GHz | **Moderado** | barrido en R (~50–100 puntos × ~1 s) para localizar el cruce y medir el gap |
| Espaciados vibracionales | **No** | requiere curvas BOP en malla fina de R **y** resolver el movimiento nuclear en R |
| Constantes rotacionales B₀, B₇ del trímero | **No** | ídem: son constantes de los niveles vibracionales del ULRM |

Restricción que afecta a los tres últimos: nuestras curvas sólo existen para
`R ≲ 1152 a₀`, límite de dominio del remapeo k(R) (punto de retorno clásico
externo de n=24). Si el paper cita estructura vibracional a R mayores, no se
alcanza con las tablas actuales.

## 9. APLICADO: override local δ₀(ns) = 3.13180 (decisión del usuario, opción 1)

`systems/atom.py` **no se ha tocado**. Se añadió `systems/rb_defects.py` con el
override, y un parámetro opcional `delta0_ns=None` en
`hamiltonians/charge_dipole.rydberg_diagonal`, `hamiltonians/fermi_krb.n_star_of_l`
y `FermiPseudopotential.__init__`. Con `None` el comportamiento es idéntico al
anterior; `energy_rb(n, l, None)` delega literalmente en `Atom.E_Rb()`.

**Verificación de no-regresión**: `pytest -m "not slow"` → 37 passed;
`pytest -m slow` (golden files G1-G4) → 2 passed en 331.48 s. Ningún resultado
cambia.

### Tabla I recalculada

| n | l | E_nl [E_h] | n* | nuevo [GHz] | antes [GHz] | paper | dif |
|---|---|---|---|---|---|---|---|
| 25 | 3 | −8.000000000000000e−04 | 25.0000000 | 447.78404 | 447.78404 | 447.78 | +0.004 |
| 28 | 0 | −8.085207285236110e−04 | 24.8679178 | 391.72034 | 391.98543 | 391.72 | **+0.000** |
| 26 | 2 | −8.226318540067425e−04 | 24.6537076 | 298.87360 | 298.87360 | 298.87 | +0.004 |
| 27 | 1 | −8.431955228646972e−04 | 24.3512274 | 163.57116 | 163.57116 | 163.57 | +0.001 |
| 24 | 3 | −8.680555555555555e−04 | 24.0000000 | 0.00000 | 0.00000 | 0.00 | +0.000 |
| 27 | 0 | −8.776913467604867e−04 | 23.8678936 | 63.40046 | 63.10064 | 63.40 | **+0.000** |
| 25 | 2 | −8.936519942688221e−04 | 23.6537928 | 168.41648 | 168.41648 | 168.42 | −0.004 |
| 26 | 1 | −9.169637842916073e−04 | 23.3511842 | 321.80069 | 321.80069 | 321.80 | +0.001 |

Por serie: s `+0.0003, +0.0005` · p `+0.0012, +0.0007` · d `+0.0036, −0.0035` ·
f `+0.0040`. **Peor |dif| sobre las 8 = 0.0040 GHz.** La serie s pasa a ser la
que mejor concuerda; el residuo de d y f es compatible con el redondeo a dos
decimales de la tabla publicada. `n*(27s)`: 23.8685133 → 23.8678936.

## 10. Convergencia en N (N_max=6 vs N_max=8), con δ₀ corregido

`dim(M_J=0)`: 1016 (N≤6) → 1640 (N≤8). Diez autovalores más bajos de
`H_a + H_mol + V_Fermi`:

| R [a₀] | peor \|ΔE\| [MHz] | peor \|ΔE\|/\|E\| | peor \|ΔE\|/\|E−E_man\| | ¿< 2×10⁻⁶ en rel(E)? |
|---|---|---|---|---|
| 600 | 18.972 | **3.072×10⁻⁶** | 4.081×10⁻⁵ | **NO** |
| 900 | 0.771 | 1.342×10⁻⁷ | 2.547×10⁻⁵ | sí |

- A `R=900` se cumple con holgado margen (1.3×10⁻⁷).
- A `R=600` **no se cumple**: 3.07×10⁻⁶, un factor 1.5 por encima del criterio.
  Lo aporta un solo estado (k=9, 18.97 MHz); los otros nueve están entre
  10⁻¹¹ y 10⁻¹³ relativo. Es sensibilidad de un estado casi degenerado, no un
  fallo global de la base.
- ⚠️ **La afirmación del paper es ambigua en el denominador.** Referido a la
  energía absoluta, se cumple en R=900 y falla por poco en R=600. Referido a la
  energía de ligadura (`E − E_manifold`, que es la magnitud físicamente
  relevante), es 2-4×10⁻⁵ en ambos casos, más de diez veces el criterio.

## 11. Orientación ⟨cos θ_d⟩ (bloque M_J=0)

Operador `⟨i|cosθ_d|j⟩ = δ_{ll'}δ_{mm'}·cos_theta_element(N,M_N,N',M_N')`,
construido con la función ya verificada; `‖C‖_F = 17.1957`, simétrico.

Dominio: el par **27s–manifold** (energía media) tiene punto de retorno a
**1145.6 a₀**, más corto que el del manifold solo (1151.5 a₀). Ese es el límite
real del escaneo — un bug de la primera versión del script, que usaba el
manifold y reventaba en R=1150.

Rama suave (dos criterios de seguimiento coinciden), monótona decreciente:

| R [a₀] | 110 | 200 | 290 | 380 | 400 | 480 | 560 |
|---|---|---|---|---|---|---|---|
| ⟨cos θ_d⟩ | +0.930 | +0.885 | +0.843 | +0.807 | +0.770 | +0.713 | +0.668 |

- **Cruza 0.78 entre R = 380 y 400 a₀.**
- **Máximo en el dominio = 0.930 en R = 110 a₀** — pero es el **borde del
  escaneo**, no un máximo interior: la curva decrece monótonamente, así que el
  "máximo" lo fija dónde se corta, y R≈110 es justo el extremo interno del
  dominio del remapeo (105.6 a₀), donde la aproximación semiclásica y el
  pseudopotencial puntual son menos fiables.
- Para `R ≳ 600` **el seguimiento del estado deja de ser fiable**: el criterio
  "más bajo con carácter de manifold" salta entre estados casi degenerados
  (⟨cos⟩ oscila +0.13, +0.69, +0.24, −0.11, +0.57…) y el seguimiento adiabático
  con pasos de 40 a₀ tampoco resuelve un espectro tan denso. No se cita ningún
  valor de esa región.

**No se puede confirmar ni refutar el 0.78 del paper con este cálculo**: nuestra
curva no tiene máximo interior en el rango accesible, y el 0.78 publicado
corresponde presumiblemente a un estado y un R concretos que no se pueden
identificar sin más contexto del artículo. Lo que sí se puede decir es que el
valor 0.78 cae dentro del rango que produce nuestra curva, en R ≈ 390 a₀.

## 12. Acciones pendientes

1. ~~Decisión del usuario sobre δ₀(ns)~~ — **RESUELTO**: opción 1, override
   local en `systems/rb_defects.py` (§9).
2. `δ₂(ns) = 0.1745312` no coincide con Li+ 2003 (0.1784). Corregirlo sería
   consistente bibliográficamente aunque sólo mueva 0.003 GHz. **No se ha
   tocado**: cambiar `atom.py` altera `trimer.py`, congelado por golden files.
3. Verificar la hipótesis del potencial modelo contra Marinescu+ 1994.

## Referencias

- González-Férez, Sadeghpour & Schmelcher, *New J. Phys.* **17**, 013021 (2015), Tabla I. DOI: 10.1088/1367-2630/17/1/013021 · arXiv:1406.6549
- Li, Mourachko, Noel & Gallagher, *Phys. Rev. A* **67**, 052502 (2003) — δ(ns), δ(nd)
- Han, Jamil, Norum, Tanner & Gallagher, *Phys. Rev. A* **74**, 054502 (2006) — δ(np)
- Marinescu, Sadeghpour & Dalgarno, *Phys. Rev. A* **49**, 982 (1994) — potencial modelo, ref. [20] del paper
- CODATA 2018 (NIST): 1 E_h = 6.579683920502×10¹⁵ Hz
