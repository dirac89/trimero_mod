# Curva BOP de M_J=0: seguimiento, resultado y una señal de alarma

> ⚠️ **AVISO DE PREMISA (2026-08-20).** Las curvas de este documento incluyen el
> **pseudopotencial de Fermi**, que NO forma parte del Hamiltoniano de
> Aguilera-Fernández et al. 2015 para Rb*-KRb. El §2 —seguimiento adiabático vs
> diabático, y que el k-ésimo autovalor ordenado es la k-ésima curva
> adiabática— sigue siendo metodológicamente válido, pero **los números y el
> «butterfly de −6.5 THz» son de un modelo que no es el de ese paper**. Ver
> `docs/analysis_fig1_carga_dipolo_sin_fermi.md`.

> ⚠️ **Nota añadida 2026-08-19.** Base electrónica **incompleta** (manifold n=24
> + 27s). Con la base del paper —que añade 25d y 26p— el §2 sigue siendo cierto
> («el k-ésimo autovalor ordenado ES la k-ésima curva adiabática») pero deja de
> bastar: los estados de 25d y 26p caen en el rango de energía relevante y la
> curva de índice fijo **cambia de carácter** al atravesar sus cruces evitados.
> Hay que identificar la curva por CARÁCTER (peso de manifold > 50 % en cada R).
> Con la base incompleta las dos coincidían, que es por lo que aquí funcionó.
> Ver `docs/analysis_base_correcta_3_vecinos.md` §5.

**Fecha**: 2026-08-19
**Autor**: Javier Aguilera
**Relevancia**: Primera curva de potencial Born-Oppenheimer propia, comparable
con la Fig. 2 de González-Férez 2015. Destapa que el pseudopotencial de onda p
domina el espectro de una forma que no encaja con la figura publicada.

## Resumen

Se implementó `simulation/bop_tracking.trace_curve` (seguimiento por solapamiento
máximo con refinamiento adaptativo del paso) y se trazó la curva de M_J=0 sobre
todo el dominio válido. **El seguimiento funciona como se especificó pero produce
una curva DIABÁTICA**, no la adiabática que muestra una figura BOP. La curva
adiabática correcta se obtiene gratis del espectro ordenado. El resultado
destapa un estado *butterfly* de onda p a −6.5 THz que no encaja con la figura
publicada.

## 1. Dominio: hay que calcularlo, no elegirlo

Dos bugs sucesivos en el guardado del dominio, ambos por elegir a mano el par de
estados que fija el límite:

1. Usar `n*=24` (manifold): falla en R=1150 por el par mixto 27s-manifold.
2. Usar la media 27s-manifold: falla en R=1143 por el par **27s-27s**, el más
   ligado de toda la base.

Corregido con `FermiPseudopotential.domain_bounds()`, que recorre **todos** los
pares (l₁,l₂) y devuelve el intervalo más restrictivo, con margen de 1e-9
relativo para el redondeo del borde:

```
R ∈ [105.61, 1138.92] a₀     (el par crítico es 27s-27s)
```

## 2. El seguimiento por solapamiento da una curva diabática

`trace_curve` con paso base 5 a₀, umbral de solapamiento 0.7 y hasta 6
bisecciones:

- 211 puntos, 220 diagonalizaciones, 164.6 s
- solapamiento |⟨v|v_prev⟩|²: **mínimo 0.7099, media 0.9587**
- **11 refinamientos** por solapamiento bajo, concentrados en R ≈ 920 y
  R ≈ 1002–1018 a₀. En R=1018.1 agotó las 6 bisecciones (paso 0.156 a₀) sin
  superar 0.7: ahí hay un racimo de estados casi degenerados, no un cruce
  evitado resoluble.

⚠️ **El índice ordenado del estado seguido cambia 22 veces** (recorre 16 índices
distintos, de k=2 a k=25). Dentro de un bloque M_J las curvas adiabáticas no se
cruzan, así que 22 cambios no son cruces σ_v genuinos: **el seguimiento por
solapamiento está saltando entre curvas adiabáticas, es decir, sigue el carácter
diabático del estado**. Es un resultado correcto pero NO es lo que dibuja una
figura de curvas BOP.

**Consecuencia práctica**: para curvas adiabáticas no hace falta seguimiento
ninguno. El k-ésimo autovalor ordenado ES la k-ésima curva adiabática, por la
regla de no cruce. El seguimiento por solapamiento sirve para lo contrario:
seguir un estado diabático a través de cruces evitados.

## 3. Validación fuerte: los umbrales asintóticos

En el borde del dominio (R=1138.9 a₀) las curvas deben tender a
27s+KRb(N) y manifold+KRb(N). Asignación por peso del autovector, inequívoca
(>99 %):

| k | E−E_man [GHz] | peso manifold | asignación | umbral esperado [GHz] |
|---|---|---|---|---|
| 0 | −63.589 | 0.0019 | 27s, N=0 | −63.40 |
| 1 | −61.373 | 0.0024 | 27s, N=1 | −61.17 |
| 2 | −56.914 | 0.0024 | 27s, N=2 | −56.72 |
| 3 | −50.231 | 0.0024 | 27s, N=3 | −50.03 |
| 4 | −41.320 | 0.0025 | 27s, N=4 | −41.12 |
| 5 | −30.181 | 0.0026 | 27s, N=5 | −29.98 |
| **6** | **−25.718** | **0.9949** | **manifold** | (aún desplazada) |
| 7 | −16.785 | 0.0027 | 27s, N=6 | −16.60 |

Los siete umbrales del 27s se reproducen con un residuo constante de ≈ 0.19 GHz
— la curva todavía no ha llegado a su asíntota a R=1138.9 a₀, que es el borde
del dominio, no el infinito.

**La curva análoga a la roja gruesa de la Fig. 2 es la adiabática k=6**: la más
baja con carácter de manifold en el borde.

## 4. Características de la curva k=6

| R [a₀] | E−E_man [GHz] | profundidad vs máximo derecho [GHz] |
|---|---|---|
| 235.0 | −24.955 | 10.657 |
| 400.0 | −15.343 | 0.024 |
| 455.0 | −15.900 | 0.133 |
| **565.0** | **−6472.043** | **6432.572** |
| 905.0 | −40.068 | 9.854 |
| 1083.3 | −30.229 | — |

Separación entre mínimos: media 169.7 a₀, σ 95.7 a₀, **CV = 0.564 → el patrón
NO es regular** en nuestro cálculo.

### Contraste con el texto del paper ("shifted more than 20 GHz for R ≲ 1200 a₀")

| curva | fracción del dominio con \|E−E_man\| > 20 GHz |
|---|---|
| adiabática k=6 | **73.46 %** (los puntos por debajo están confinados a R ∈ [275, 550] a₀) |
| diabática seguida por solapamiento | 15.64 % |

La adiabática es ampliamente consistente con la afirmación; la diabática no, lo
que confirma que la diabática no es el objeto que describe el paper.

## 5. ⚠️ SEÑAL DE ALARMA: el butterfly de onda p a −6.5 THz

El rasgo dominante de nuestro espectro es un pozo de **−6472 GHz en R = 565 a₀**,
que arrastra hacia abajo las ~12 curvas más bajas simultáneamente. Es un estado
*butterfly*: lo genera el término de onda p, cuyo volumen de dispersión alcanza
`A_p ≈ −6.8×10⁵ a₀³` en la tabla.

Los estados butterfly son físicamente reales en las ULRM, pero **−6.5 THz es dos
órdenes de magnitud más profundo que cualquier cosa que la Fig. 2 pueda mostrar**
(su ventana es de decenas de GHz). Dos lecturas posibles, sin resolver:

1. El paper restringe la ventana de energía y simplemente no dibuja las curvas
   butterfly.
2. Nuestro tratamiento de `A_p` las está exagerando. R=565 a₀ remapea a
   R' = 763 a₀ en la tabla de n=35, justo donde `A_p` atraviesa la resonancia de
   forma p. **Es exactamente el punto donde la aproximación del remapeo k(R)
   (§3.2 de `analysis_pseudopotencial_fermi_krb.md`) es más frágil**: un error
   pequeño en `k` se traduce en un error enorme en `A_p` cerca de la resonancia.

Esto **no se puede decidir sin la figura o sin regenerar las tablas para n=24**
(opción (ii) del bloqueante original). Queda como el punto abierto más
importante antes de dar por buena ninguna comparación de curvas.

## 6. Archivos

- `src/trimero/simulation/bop_tracking.py` — `trace_curve` con refinamiento adaptativo.
- `src/trimero/hamiltonians/fermi_krb.py` — `FermiPseudopotential.domain_bounds()` (nuevo).
- `scripts/run_bop_curve.py` — barrido, características y figura.
- `plots/bop_curve_MJ0.png`, `plots/bop_curve_MJ0.npz`.

## 7. Pendiente

> **Actualizado 2026-08-19.** El punto 1 está resuelto.

1. ~~Decidir si el butterfly de −6.5 THz es físico o artefacto del remapeo
   (§5).~~ → **RESUELTO**: no es artefacto del remapeo ni de la interpolación.
   Es la **divergencia del pseudopotencial de Fermi de rango cero en la
   resonancia de forma p**, una limitación conocida del modelo. Con la
   interpolación correcta de `1/A_p` el pozo se profundiza hasta −11 049 GHz y
   es formalmente no acotado. Por decisión del usuario **no se implementa la
   corrección de rango efectivo (Omont 1977)**: la región se excluye con una
   ventana explícita `R ∈ [536.40, 589.15] a₀`. Ver
   `docs/analysis_ventana_exclusion_resonancia.md`.
   ⚠️ Con ello, la cifra de **−6472 GHz en R = 565 a₀** de §4 y §5 queda
   **obsoleta**: era el valor con interpolación lineal a través del polo.
2. M_J=1 (no abordado, a la espera de validar M_J=0).
3. La curva k=6 aún no ha alcanzado su asíntota en el borde del dominio: no
   podemos ver el régimen R > 1139 a₀ con las tablas actuales.
4. La estadística de regularidad de §4 (CV = 0.564) se recalculó excluyendo la
   ventana de resonancia: **CV = 0.416 sobre tres separaciones válidas**, sigue
   siendo «no regular» pero con muy poca base estadística.

## Referencias

- González-Férez, Sadeghpour & Schmelcher, *New J. Phys.* **17**, 013021 (2015), Fig. 2.
- `docs/analysis_pseudopotencial_fermi_krb.md` — aproximaciones del remapeo k(R).
