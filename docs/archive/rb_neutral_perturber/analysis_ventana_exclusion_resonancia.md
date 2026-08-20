# Ventana de exclusión de la resonancia de forma p: criterio, comparación y cierre

> ⚠️ **AVISO DE PREMISA (2026-08-20).** Este documento trata el
> **pseudopotencial de Fermi**, que **NO forma parte del Hamiltoniano de
> Aguilera-Fernández et al. 2015 / González-Férez et al. 2015** para Rb*-KRb.
> Su Ec. 1 es `H_ad = H_A + H_mol`, con KRb como **dipolo puntual**, no como
> centro de dispersión de contacto. Su contenido técnico **sigue siendo
> válido** —y sigue aplicando a la línea de perturbador NEUTRO
> (Aguilera-Fernández 2016)— pero **su uso para comparar con esos dos papers
> partía de una premisa equivocada**. Ver
> `docs/analysis_fig1_carga_dipolo_sin_fermi.md` §1.

> ⚠️ **Nota añadida 2026-08-19.** Todo este documento usa la base electrónica
> **incompleta** (manifold n=24 + 27s, dim 1064→1016), anterior a leer
> arXiv:1507.07972. La **ventana de exclusión en sí no cambia** —depende sólo del
> manifold y de `rvsAP.dat`, no de los vecinos: R_polo = 562.771 a₀,
> FWHM = 13.187 a₀, ventana [536.40, 589.15] a₀, verificado en
> `tests/simulation/test_bop_system.py::test_s1_n24_base_correcta`—. Lo que sí
> cambia es el **dominio**: con la base correcta llega a 1090.16 a₀, no a
> 1138.92, porque el par más ligado pasa a ser 26p-26p. Y la identificación de la
> curva por índice fijo k=6 deja de ser válida con la base completa; ver
> `docs/analysis_base_correcta_3_vecinos.md`. Las curvas y la figura de este
> documento están, por tanto, **pendientes de rehacer**.

**Fecha**: 2026-08-19
**Autor**: Javier Aguilera
**Relevancia**: Cierra el bloque abierto en `analysis_curva_bop_MJ0.md` →
`analysis_resonancia_onda_p.md` → `analysis_interpolacion_polo_Ap.md`. Fija por
decisión explícita que la divergencia del pseudopotencial de rango cero en la
resonancia de forma p se **acepta como limitación conocida** y se **excluye** de
toda comparación cuantitativa, sin implementar la corrección de rango efectivo.
**Tipo**: analysis

## Resumen

Tres resultados.

1. **El cambio de interpolación no tocó nada fuera de la resonancia.** Con la
   curva BOP de M_J=0 regenerada con `p_interpolation="inverse"` y comparada
   punto a punto con la versión antigua (`"linear"`), la discrepancia máxima
   fuera de la ventana de resonancia es **0.194 GHz** sobre las 60 curvas
   adiabáticas, con **mediana 3.7×10⁻⁶ GHz**; dentro de la ventana es
   **4.6×10³ GHz**. El cambio está confinado donde debía estarlo.
2. **La ventana de exclusión queda definida por un criterio calculado, no a
   ojo**: `R ∈ [536.40, 589.15] a₀`, es decir `R_polo ± 2×FWHM` con
   `R_polo = 562.771 a₀` (cero del interpolante de `1/A_p`) y
   `FWHM = 13.187 a₀`. Está implementada en
   `ScatteringLengths.p_resonance_window()` y cubierta por un test.
3. **Pero ±2×FWHM no basta para la comparación con la Fig. 2**, y hay que
   decirlo: la ventana elimina la *divergencia*, no su *cola*. Justo fuera del
   borde derecho la adiabática k=6 aún vale **−659.6 GHz**, y no vuelve a caber
   en la ventana de energía de la figura (|E| ≤ 120 GHz) hasta **R = 680 a₀**.
   Ese tramo se marca aparte como *banda de cola butterfly* y se excluye también
   de la estadística.

## Palabras clave

- resonancia de forma p, ³P^o de e⁻–Rb(5S)
- pseudopotencial de Fermi de rango cero, divergencia
- estados butterfly
- ventana de exclusión, FWHM
- Omont 1977 (corrección no implementada)

## 1. Punto 4 pendiente: curva BOP regenerada con `inverse`

`poetry run python scripts/run_bop_curve.py` con el valor por defecto
`p_interpolation="inverse"`. 211 puntos, 220 diagonalizaciones, **223.4 s**.
La malla en R sale **idéntica** a la del barrido antiguo (mismos 11
refinamientos, mismos R), de modo que la comparación es punto a punto y no
requiere reinterpolar nada.

Comparación contra `plots/bop_curve_MJ0_linear.npz`, con la ventana de exclusión
de la §2 (201 puntos fuera, 10 dentro):

| curva | max‖Δ‖ **fuera** | mediana‖Δ‖ **fuera** | max‖Δ‖ **dentro** |
|---|---|---|---|
| adiabática k=6 (carácter manifold) | 1.939×10⁻¹ GHz | 3.655×10⁻⁶ GHz | 4.577×10³ GHz |
| diabática (seguida por solapamiento) | 4.922×10⁻⁴ GHz | 5.299×10⁻⁶ GHz | 1.803×10⁻¹ GHz |
| **las 60 adiabáticas a la vez** | **1.939×10⁻¹ GHz** | — | **4.577×10³ GHz** |

Dónde está ese máximo de 0.194 GHz y qué lo produce:

| R [a₀] | max‖Δ‖ sobre las 60 curvas [GHz] |
|---|---|
| 590.0 | 1.939×10⁻¹ |
| 600.0 | 1.306×10⁻¹ |
| 605.0 | 1.176×10⁻¹ |
| 610.0 | 8.005×10⁻² |
| 595.0 | 5.541×10⁻² |

Los cinco puntos con más discrepancia son **los cinco inmediatamente pegados al
borde derecho de la ventana**: es el flanco de la resonancia, donde `A_p` todavía
vale ~6×10⁴ a₀³. Alejándose, la diferencia cae cuatro órdenes de magnitud. En
todo el resto del dominio (R < 536 a₀ y R > 620 a₀) los dos métodos coinciden a
nivel de **10⁻⁶–10⁻² GHz**, muy por debajo de cualquier cifra que se compare con
el paper.

**Conclusión del punto 4: confirmado, el cambio de interpolación no alteró nada
fuera de la región problemática.**

## 2. La ventana de exclusión, con el criterio explícito

Implementada en `ScatteringLengths.p_resonance_window(n_star, margin_factor)`.
Los tres ingredientes se calculan de `rvsAP.dat`, ninguno se escribe a mano:

### 2.1 Centro = el polo que ve el propio modelo

No el punto medio del hueco de malla (que daba 561.8 a₀ en
`analysis_resonancia_onda_p.md`, un número que dependía de dónde acabara la
malla), sino el **cero del interpolante lineal de `1/A_p` en ε = k²/2 entre los
dos nodos donde `A_p` cambia de signo**. Es literalmente el polo que produce el
modo `p_interpolation="inverse"` en producción: si el modelo diverge en algún
sitio, diverge exactamente ahí.

```
eps_polo = 0.00090887 E_h = 24.732 meV
R_polo   = 562.771 a0
```

El ajuste lineal global sobre 20 nodos (§1 de `analysis_interpolacion_polo_Ap.md`)
da 24.711 meV → R = 563.01 a₀. Los dos valores distan **0.24 a₀**, un 2 % de la
FWHM: irrelevante para dónde poner la ventana. Se adopta el local por ser el que
usa el código.

### 2.2 Anchura = FWHM de |A_p|

Extensión en R de los nodos de la tabla con `|A_p| ≥ max|A_p|/2`:

```
R in [555.994, 569.181] a0   ->   FWHM = 13.187 a0
```

Reproduce los 13.2 a₀ ya medidos. **Salvedad honesta**: `max|A_p|` es finito sólo
porque la malla no cae sobre el polo, así que esta FWHM es una **escala operativa
de anchura**, no la Γ de la resonancia en sentido estricto. Sirve para
dimensionar una ventana, no para citarla como anchura física.

### 2.3 Margen y ventana

```
margen  = ± 2 × FWHM = ± 26.375 a0
VENTANA EXCLUIDA:  R ∈ [536.40, 589.15] a0     (anchura 52.75 a0, 10 de 211 puntos)
```

Es el **5.1 % del dominio** [105.6, 1138.9] a₀.

### 2.4 Diagnóstico: ¿por qué 2 y no otro factor?

| factor | ventana [a₀] | n puntos fuera | min E(k=6) fuera [GHz] | max‖Δ‖ fuera [GHz] |
|---|---|---|---|---|
| 1 | [549.58, 575.96] | 205 | −1049.86 | 2.411 |
| 1.5 | [542.99, 582.55] | 203 | −816.17 | 0.581 |
| **2** | **[536.40, 589.15]** | **201** | **−659.61** | **0.194** |
| 3 | [523.21, 602.33] | 195 | −410.09 | 0.118 |
| 4 | [510.02, 615.52] | 190 | −285.43 | 0.029 |
| 6 | [483.65, 641.89] | 179 | −180.40 | 0.010 |
| 8 | [457.27, 668.27] | 169 | −134.06 | 0.003 |

Dos lecturas de esta tabla, y conviene no mezclarlas:

- **Para la reproducibilidad del cambio de interpolación**, ±2×FWHM ya es
  suficiente: 0.194 GHz de discrepancia máxima, y sólo en el punto pegado al
  borde.
- **Para la comparación de energías con la Fig. 2**, ningún factor razonable
  basta: la columna `min E(k=6)` baja de −1050 a −134 GHz entre f=1 y f=8 sin
  estabilizarse. La influencia energética de la resonancia decae mucho más
  despacio que |A_p|, porque a partir del polo la curva k=6 **ya no es la del
  manifold: es el estado butterfly** que la resonancia genera.

Se adopta **f = 2** como ventana declarada —el criterio es sobre `A_p`, que es
donde está la divergencia— y la cola se trata por separado en la §3.

## 3. La banda de cola butterfly

Criterio, también explícito: `R > R_hi` y `|E(k=6) − E_man| > 120 GHz`, siendo
120 GHz el borde inferior de la ventana de energía comparable a la Fig. 2 que ya
usábamos en `scripts/run_bop_curve.py`.

```
BANDA DE COLA:  R ∈ (589.15, 680.00] a0   —   19 puntos
```

En ella la adiabática k=6 pasa de −659.6 GHz a −120.6 GHz. No se excluye por
divergencia (ahí `A_p` es finito y está tabulado), sino porque **la etiqueta
«curva adiabática más baja con carácter de manifold» apunta en ese tramo al
estado butterfly, no al objeto que dibuja la Fig. 2**. Se dibuja rayada y se
excluye también de la estadística del criterio de 20 GHz.

## 4. Características de la curva, excluyendo la ventana

Adiabática k=6, con la ventana fuera. Los extremos locales se buscan **por
separado en cada segmento contiguo**, para que el hueco de la exclusión no cree
mínimos ni máximos falsos:

| # | segmento | R_min [a₀] | E−E_man [GHz] | prof. vs máx. derecho [GHz] | ΔR al anterior [a₀] |
|---|---|---|---|---|---|
| 0 | izq (R < 536.4) | 235.00 | −24.9547 | 10.6575 | — |
| 1 | izq | 400.00 | −15.3433 | 0.0238 | 165.00 |
| 2 | izq | 455.00 | −15.9004 | 0.1332 | 55.00 |
| 3 | der (R > 589.1) | 905.00 | −40.0682 | 9.8539 | 450.00 \* |
| 4 | der | 1083.28 | −30.2291 | — | 178.28 |

\* salto que atraviesa la ventana excluida; **no cuenta** como separación física.

- Ningún mínimo local cae dentro de la banda de cola (589–680 a₀): ahí la curva
  es monótona. La tabla de características **no cambia** por incluir o no esa
  banda; sólo cambia la estadística de la §4.2.
- Separación entre mínimos **contiguos** (descontando el salto marcado):
  media **132.76 a₀**, σ **55.25 a₀**, CV = **0.416** → sigue siendo **NO
  regular**. Con sólo tres separaciones válidas, esta cifra es indicativa y no
  soporta una afirmación fuerte sobre regularidad del patrón.
- El mínimo global de la curva **fuera de la ventana** es −659.6 GHz en
  R = 590 a₀, y es cola butterfly; **fuera de ventana y de cola** es
  −114.0 GHz en R = 685 a₀ (todavía flanco de la misma estructura: el criterio
  de corte es energético, así que este número queda por construcción pegado al
  umbral de 120 GHz). Dentro de la ventana la curva llega a **−11 049 GHz**,
  valor sin significado físico.

### 4.2 «shifted more than 20 GHz for R ≲ 1200 a₀»

Nuestro dominio llega sólo a R = 1138.9 a₀, así que la comprobación cubre casi
todo el rango que enuncia el paper.

| curva | sin excluir nada | excluida la ventana | excluidas ventana + cola |
|---|---|---|---|
| adiabática k=6 | 72.99 % | 73.63 % | **70.88 %** |
| diabática (solapamiento) | 15.64 % | 16.42 % | **18.13 %** |

Los 53 puntos de la adiabática k=6 que quedan **por debajo** de 20 GHz están
confinados a `R ∈ [275, 535] a₀`.

**La afirmación del paper se sostiene con la adiabática y no con la diabática, y
esa conclusión es robusta frente a la exclusión**: las tres columnas de la fila
k=6 caen dentro de un margen de 2.8 puntos porcentuales. Es decir, la resonancia
no es la que hace que el 73 % del dominio esté desplazado más de 20 GHz — se
puede quitar entera y la cifra apenas se mueve. Ése es el resultado que
buscábamos aislar.

## 5. Figura final

`plots/bop_curve_MJ0_excluded.png`, generada por
`scripts/analyze_resonance_window.py`.

- **Panel superior (rango completo)**: se **muestra** la divergencia, para que se
  vea qué se está excluyendo y por qué. Banda gris = ventana excluida; banda
  rayada = cola butterfly; línea de puntos = polo de `A_p` en 562.8 a₀.
- **Panel inferior (ventana de energía de la Fig. 2)**: las curvas se **cortan**
  dentro de la ventana, pero la región queda sombreada y rotulada con el motivo
  y la referencia al documento, no simplemente recortada sin explicación.
- El rótulo dentro de la propia figura dice: *«REGIÓN EXCLUIDA — divergencia
  conocida del pseudopotencial de rango cero en la resonancia de forma p
  (rayado: cola butterfly, tampoco comparable), ver
  docs/analysis_interpolacion_polo_Ap.md»*.

## 6. Nota de limitación, redactada para un manuscrito

> El pseudopotencial de Fermi de rango cero empleado aquí diverge en la
> resonancia de forma ³P^o de e⁻–Rb(5S), que nuestras tablas de dispersión
> sitúan en E_r = 24.7 meV (frente a ≈ 23 meV en la literatura) y que en la
> coordenada internuclear corresponde a R ≈ 562.8 a₀. Se trata de una limitación
> conocida de la aproximación de rango cero, subsanable mediante el tratamiento
> dependiente de la energía y de rango efectivo de Omont [Omont 1977], que no se
> ha implementado en este trabajo; en consecuencia excluimos de toda comparación
> cuantitativa la ventana R = 562.8 ± 26.4 a₀, definida como el polo de A_p más
> menos dos veces la anchura a media altura de |A_p| (13.2 a₀). Esta salvedad
> afecta a la profundidad que nuestro modelo asigna a los estados butterfly en
> esa región, no a su existencia: son un fenómeno bien establecido, predicho por
> Hamilton, Greene y Sadeghpour [Hamilton 2002] y observado experimentalmente por
> Niederprüm y colaboradores [Niederprüm 2016].

## Conclusiones y aplicación al proyecto

1. El cambio a `p_interpolation="inverse"` está **verificado como inocuo** fuera
   de la resonancia (mediana 3.7×10⁻⁶ GHz de diferencia). Queda cerrado el punto
   4 que había quedado sin ejecutar.
2. La región problemática está ahora **acotada con un criterio reproducible y
   testado**, no con límites de eje ad hoc.
3. La comparación con la Fig. 2 de González-Férez 2015 puede hacerse sobre el
   94.9 % del dominio, y la afirmación de los 20 GHz **se sostiene con la curva
   adiabática independientemente de la exclusión**.
4. **No se ha implementado ninguna corrección de rango efectivo**, por decisión
   explícita. La divergencia se documenta y se cita, no se corrige.

## Estado de los pendientes heredados

| documento | pendiente | estado |
|---|---|---|
| `analysis_interpolacion_polo_Ap.md` | 1. regenerar curva BOP con `inverse` y comparar | **HECHO** (§1) |
| `analysis_interpolacion_polo_Ap.md` | 2. valorar corrección de rango efectivo (Omont) | **CERRADO por decisión**: no se implementa; se documenta y se excluye la región |
| `analysis_interpolacion_polo_Ap.md` | 3. `A_p(k)` de fuente independiente | sigue abierto (no bloquea nada ahora) |
| `analysis_resonancia_onda_p.md` | 2. interpolación consciente del polo | **HECHO** en la ronda anterior |
| `analysis_resonancia_onda_p.md` | 3. valorar rango efectivo | **CERRADO por decisión** |
| `analysis_curva_bop_MJ0.md` | 1. ¿el butterfly es físico o artefacto? | **RESUELTO**: es la divergencia del rango cero, limitación conocida del modelo; región excluida |
| todos | M_J = 1 | **sigue abierto**, sin abordar |

## Archivos tocados

- `src/trimero/hamiltonians/fermi_krb.py` — `ScatteringLengths.p_resonance_window()` (nuevo).
- `scripts/analyze_resonance_window.py` — nuevo: comparación, ventana, características y figura.
- `tests/hamiltonians/test_fermi_krb.py` — `test_f4_p_resonance_window` (nuevo).
- `plots/bop_curve_MJ0.npz`, `plots/bop_curve_MJ0.png` — regenerados con `inverse`.
- `plots/bop_curve_MJ0_excluded.png` — figura final con la región marcada.

`poetry run pytest -m "not slow"` → **38 passed**, sin regresiones.

## Referencias

- Omont, A., *J. Physique* **38**, 1343 (1977) — correcciones más allá del rango
  cero (dependencia en energía / rango efectivo). **Vía de corrección NO
  implementada**, citada como tal.
- Hamilton, E. L., Greene, C. H. & Sadeghpour, H. R., *J. Phys. B* **35**, L199
  (2002) — estados butterfly inducidos por la resonancia de forma p.
  https://lweb.cfa.harvard.edu/~hrs/itamp/JPBL2002.pdf
- Niederprüm, T. *et al.*, *Nat. Commun.* **7**, 12820 (2016) — observación de
  moléculas de Rydberg butterfly pendulares.
  https://www.nature.com/articles/ncomms12820
- González-Férez, R., Sadeghpour, H. R. & Schmelcher, P., *New J. Phys.* **17**,
  013021 (2015), Fig. 2 — objeto de la comparación.
- `docs/analysis_interpolacion_polo_Ap.md`, `docs/analysis_resonancia_onda_p.md`,
  `docs/analysis_curva_bop_MJ0.md`.
