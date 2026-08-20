# La resonancia de forma de onda p: ¿física o artefacto de la tabla?

> ⚠️ **AVISO DE PREMISA (2026-08-20).** Este documento trata el
> **pseudopotencial de Fermi**, que **NO forma parte del Hamiltoniano de
> Aguilera-Fernández et al. 2015 / González-Férez et al. 2015** para Rb*-KRb.
> Su Ec. 1 es `H_ad = H_A + H_mol`, con KRb como **dipolo puntual**, no como
> centro de dispersión de contacto. Su contenido técnico **sigue siendo
> válido** —y sigue aplicando a la línea de perturbador NEUTRO
> (Aguilera-Fernández 2016)— pero **su uso para comparar con esos dos papers
> partía de una premisa equivocada**. Ver
> `docs/analysis_fig1_carga_dipolo_sin_fermi.md` §1.

> ⚠️ **CORREGIDO por `analysis_interpolacion_polo_Ap.md` (2026-08-19).** La
> conclusión 2 del resumen y la §6.2 —«la profundidad de −6472 GHz es un
> artefacto de interpolación»— son **incorrectas**. Sólo lo era el valor en ese
> punto de malla concreto; nodos tabulados genuinos (R=566, 568, 570 a₀) dan
> pozos igual de profundos. La causa real es la divergencia del pseudopotencial
> de rango cero en la resonancia.

**Fecha**: 2026-08-19
**Autor**: Javier Aguilera
**Relevancia**: La curva BOP de M_J=0 tiene un pozo de −6.5 THz en R≈565 a₀, dos
órdenes de magnitud más profundo que la ventana de la Fig. 2 de González-Férez
2015. Antes de tocar nada del Hamiltoniano hay que saber si ese pozo es física o
un defecto de nuestros datos de entrada.

## Resumen

**Las dos cosas a la vez, y se pueden separar limpiamente:**

1. **La resonancia es física.** El polo de `A_p` en `rvsAP.dat` está en
   **E = 24.82 meV**, frente a los **23 meV** que la literatura da para la
   resonancia de forma ³P^o de e⁻–Rb(5S). Un 7.9 % de diferencia en la posición
   de una resonancia es buen acuerdo: la tabla es consistente con un cálculo
   genuino de dispersión e⁻–Rb, no con un artefacto.

2. **El valor concreto de −6472 GHz NO es fiable.** `rvsAP.dat` tiene un
   **hueco de 15 a₀ en su malla** (frente a 3 a₀ de mediana) **exactamente sobre
   el polo**, entre R'=750 (A_p=+410652) y R'=765 (A_p=−680412). `np.interp`
   une esos dos puntos con una recta, es decir **interpola linealmente a través
   de una divergencia**. El mínimo de nuestra curva, en R=565 a₀, cae dentro de
   ese hueco: su `A_p = −560219` no es una lectura de la tabla, es el valor de
   una recta trazada entre dos lados opuestos de un polo.

No hizo falta invocar la limitación del pseudopotencial de rango cero para
explicar el número: falla antes, en la interpolación de los datos de entrada.
Eso no quita que la limitación de rango cero exista y siga siendo relevante para
la magnitud real.

## 1. Búsqueda de una fuente independiente de A_p(k)

**No se consiguió una tabla o parametrización independiente de A_p(k).** El
artículo de Bahrim & Thumm en Phys. Rev. A **61**, 022722 (2000) devuelve HTTP
403 (de pago), y el PDF de Hamilton, Greene & Sadeghpour, *J. Phys. B* **35**,
L199 (2002) —el trabajo que introduce los estados *butterfly* inducidos por
esta resonancia— se descargó pero su texto no es extraíble.

Lo que **sí** se obtuvo, de fuentes secundarias, es la **posición** de la
resonancia, que es lo que permite el contraste pedido:

> «Calculations by I. I. Fabrikant revealed the presence of a low-energy
> electron-Rb(5S) p-wave shape resonance for an energy of **23 meV**, which was
> predicted by A. Johnston and P. Burrow.»

y, en otra fuente, «a metastable ³P negative ion exists at fairly low energies
(e.g. **~0.03 eV** in rubidium)». Ambas cifras acotan la resonancia en el rango
23–30 meV.

## 2. Posición de la resonancia en NUESTRA tabla

Método: localizar el cambio de signo de `A_p` en `rvsAP.dat` (firma del polo,
donde δ_p(k)=π/2) y convertir esa R' a k con la relación semiclásica de n=35 que
ya usamos en el remapeo, `k² = 2(E₃₅ + 1/R')`.

```
  |A_p| máximo: -680411.6125 a0^3 en R' = 765.00 a0  ->  k = 0.042403 a.u.
  cambio de signo entre R'=750.00 (A_p=+410651.6) y R'=765.00 (A_p=-680411.6)
  polo en R' = 757.5  ->  k = 0.042708 a.u.  ->  E = k²/2 = 24.82 meV
  cota por el hueco:  E in [24.46, 25.18] meV
```

| | E de la resonancia |
|---|---|
| nuestra tabla `rvsAP.dat` | **24.82 meV** (acotado a [24.46, 25.18]) |
| literatura (Fabrikant; Johnston & Burrow) | **23 meV** |
| diferencia | **+1.82 meV (+7.9 %)** |

**Veredicto: la posición coincide razonablemente.** Es evidencia de que la
resonancia en nuestra tabla es el fenómeno físico real, no un defecto de la
tabla ni del remapeo k(R). El 7.9 % es del orden de lo que separa distintos
cálculos de estructura del ion negativo Rb⁻.

Esto además **refuerza indirectamente la procedencia** de `rvsAP.dat`, que
seguía siendo un punto abierto (`analysis_procedencia_rvsAS_rvsAP.md`): una
tabla que sitúa la resonancia ³P^o a 24.8 meV difícilmente puede ser otra cosa
que un cálculo de dispersión e⁻–Rb.

## 3. El hueco de malla sobre el polo

Al inspeccionar la malla apareció algo que no esperábamos:

```
  paso mediano de la malla: 3.0 a0
  paso máximo:             15.0 a0, entre R'=750.0 y R'=765.0
     R'=750.0  ->  A_p = +410651.6
     R'=765.0  ->  A_p = -680411.6
```

El único hueco anómalo de toda la tabla está **exactamente sobre el polo**. Es
comprensible que el generador saltase esos puntos —ahí `A_p` diverge— pero tiene
una consecuencia grave para el consumidor: `np.interp` no sabe que hay una
divergencia y traza una recta de +410652 a −680412.

### Anchura, en nuestra coordenada R (n=24)

Mapeando con `R = 1/(1/R' + E₂₄ − E₃₅)`:

| magnitud | en R' (tabla n=35) | en R (nuestro, n=24) |
|---|---|---|
| polo | 757.5 a₀ | **561.8 a₀** |
| FWHM de `A_p` | [750, 774] a₀ | [557.7, 570.8] a₀ → **13.2 a₀** |
| hueco de malla | [750, 765] a₀ | [557.7, 565.9] a₀ → **8.3 a₀** |

La FWHM ocupa el **1.27 %** de nuestro dominio [105.6, 1138.9] a₀; el hueco
interpolado, el **0.80 %**. Es una estructura **muy estrecha**.

### Qué puntos de nuestro barrido están contaminados

```
  3 de 211 puntos del barrido tienen A_p interpolado a través del polo:
    R = 555.0, 560.0, 565.0

  mínimo de la curva adiabática k=6: -6472.0 GHz en R = 565.0  ->  R' = 763.3
  ¿está dentro del hueco?  SÍ
  A_p interpolado ahí = -560218.8  (recta entre +410652 y -680412)
```

**El mínimo global de la curva BOP cae dentro del hueco.** Su profundidad de
−6472 GHz está construida sobre un valor de `A_p` que la tabla no contiene.

## 4. Respuesta a la pregunta del punto 3 (sensibilidad)

La resonancia es estrecha: **13.2 a₀ de FWHM** sobre un dominio de 1033 a₀. Con
el paso de 5 a₀ del barrido caen dentro 2–3 puntos. Las consecuencias:

- Un cálculo directo para n=24, sin pasar por el remapeo de una tabla de n=35,
  la ubicaría en un R ligeramente distinto: el desplazamiento del 7.9 % en la
  posición en energía se traduce en varios a₀ en R, comparable a la propia
  anchura.
- Una malla de R más gruesa que ~10 a₀ **puede saltarse el pico por completo** o
  registrarlo con una profundidad arbitraria, según dónde caigan los nodos.
- Es por tanto perfectamente plausible que la Fig. 2 del paper no muestre este
  pozo: o su malla no lo resuelve, o su ventana de energía lo deja fuera, o su
  `A_p` está tabulado de forma que no requiere interpolar sobre el polo.

## 5. Ventana de energía (punto 4), formalizada

`scripts/run_bop_curve.py` pasa a tener opciones en vez de límites de ejes ad
hoc:

```
--reuse                reutiliza plots/bop_curve_MJ0.npz sin rebarrer R
--step S               paso base en R [a0]
--ymin / --ymax        ventana de energía [GHz]  (por defecto -120, 25)
--ap-threshold X       |A_p| [a0^3] por encima del cual el punto se marca
--exclude-resonance    elimina de la curva los puntos con A_p interpolado
                       a través del polo
--out RUTA             fichero PNG de salida
```

Y en la librería, dos métodos nuevos de `FermiPseudopotential`:

- `large_gaps(factor=3.0)` — devuelve los intervalos de R' donde el paso de la
  malla supera `factor` veces la mediana. En `rvsAP.dat`: `[(750.0, 765.0)]`.
- `bridges_gap(R, factor=3.0)` — True si el R' remapeado de **algún** par
  (l₁,l₂) cae dentro de un hueco. (Primera versión sólo probaba un par y
  detectaba 1 punto en vez de 3; corregido.)

Salida con los valores por defecto:

```
  ventana: [-120.0, 25.0] GHz
  umbral |A_p| = 1.0e+05 a0^3 -> 6 de 211 puntos marcados
  huecos de malla en la tabla: [(750.0, 765.0)]
  puntos cuyo A_p se interpola A TRAVÉS del polo: 3  -> R = 555.0, 560.0, 565.0
  curva adiabática más baja de carácter manifold en R_max: k = 6
```

El PNG marca la zona contaminada con una banda amarilla y los puntos de `|A_p|`
alto con círculos naranjas, de modo que la parte «normal» de la curva se puede
comparar con la Fig. 2 sin que el butterfly la domine visualmente.

## 6. Conclusiones honestas

1. **La resonancia de forma p es física** (posición a 7.9 % de la literatura) y
   los estados butterfly que genera son un fenómeno real de las ULRM.
2. **La profundidad de −6472 GHz no es un resultado, es un artefacto de
   interpolación**: procede de una recta trazada a través de un polo por un
   hueco de la malla de entrada. La profundidad verdadera en R≈562 a₀ queda
   **indeterminada** con estas tablas.
3. **No se pudo obtener A_p(k) de una fuente independiente** — sólo la posición
   de la resonancia. Para acotar la magnitud haría falta el dato tabulado.
4. **No se ha implementado ninguna corrección de rango efectivo**, según lo
   pedido. Sigue siendo una vía abierta, pero ahora sabemos que **no es el
   primer problema a resolver**: antes está la interpolación de los datos.
5. La estrechez de la resonancia (13.2 a₀) explica por qué el resultado es tan
   sensible y por qué la Fig. 2 puede no mostrarla.

## 7. Pendiente

> **Actualizado 2026-08-19.** Los puntos 2 y 3 están cerrados; ver
> `docs/analysis_interpolacion_polo_Ap.md` y
> `docs/analysis_ventana_exclusion_resonancia.md`.

1. Conseguir `A_p(k)` tabulado de Bahrim & Thumm o Fabrikant (fuente de pago o
   acceso institucional) para acotar la magnitud, no sólo la posición.
   **Sigue abierto**; ya no bloquea nada, porque la región donde importaría está
   excluida de toda comparación cuantitativa.
2. ~~Sustituir `np.interp` por una interpolación consciente del polo~~ →
   **HECHO**: `p_interpolation="inverse"` interpola `1/A_p` en ε=k²/2. El hueco
   de malla dejó de importar, pero **el pozo butterfly no era artefacto suyo**
   (ver el aviso al principio de este documento).
3. ~~Sólo después: valorar si hace falta corrección de rango efectivo.~~ →
   **CERRADO POR DECISIÓN DEL USUARIO (2026-08-19): no se implementa.** La
   divergencia se acepta como limitación conocida del pseudopotencial de rango
   cero, se documenta citando Omont (1977) como la vía no tomada, y la región
   afectada se excluye con una ventana explícita
   `R ∈ [536.40, 589.15] a₀` = `R_polo ± 2×FWHM`.
4. M_J=1, aún sin abordar.

### Nota sobre la anchura citada en §3

La «FWHM de 13.2 a₀» de este documento es la que se ha reutilizado para
dimensionar la ventana de exclusión, recalculada con
`ScatteringLengths.p_resonance_window()` → **13.187 a₀**. Es una **escala
operativa**: `max|A_p|` es finito sólo porque la malla no cae sobre el polo, así
que no debe citarse como la anchura Γ de la resonancia. El centro sí se corrigió:
**562.771 a₀** (cero del interpolante de `1/A_p`) en vez de los 561.8 a₀ del
punto medio del hueco.

## 8. Archivos tocados

- `src/trimero/hamiltonians/fermi_krb.py` — `large_gaps()`, `bridges_gap()` (nuevos).
- `scripts/run_bop_curve.py` — `argparse` con ventana de energía, umbral de `A_p`,
  exclusión de puntos contaminados y `--reuse`; el plot marca ambas cosas.
- `plots/bop_curve_MJ0.png` — regenerado con el marcado.
- `docs/analysis_resonancia_onda_p.md` — este documento.

`pytest -m "not slow"` → 37 passed, sin regresiones.

## Referencias

- González-Férez, Sadeghpour & Schmelcher, *New J. Phys.* **17**, 013021 (2015).
- Bahrim & Thumm, *Phys. Rev. A* **61**, 022722 (2000) — R-matrix de Dirac para
  e⁻–Rb/Cs/Fr. https://journals.aps.org/pra/abstract/10.1103/PhysRevA.61.022722 (403, de pago)
- Hamilton, Greene & Sadeghpour, *J. Phys. B* **35**, L199 (2002) — estados
  butterfly inducidos por la resonancia de forma. https://lweb.cfa.harvard.edu/~hrs/itamp/JPBL2002.pdf (texto no extraíble)
- Engel, F., *Ultracold chemistry of a Rydberg atom in a rubidium-87 BEC*, Masterarbeit,
  Univ. Stuttgart (2016) — cita la resonancia a 23 meV (Fabrikant; Johnston & Burrow).
  https://www.pi5.uni-stuttgart.de/documents/abgeschlossene-arbeiten/2016-Engel-Felix-Ultracold-chemistry-of-a-Rydberg-atom-in-a-rubidium-87-BEC-MSC.pdf
- Schlagmüller et al., *Probing a scattering resonance in Rydberg molecules with a BEC*,
  arXiv:1510.07003. https://arxiv.org/abs/1510.07003
- Niederprüm et al., *Observation of pendular butterfly Rydberg molecules*,
  Nat. Commun. **7**, 12820 (2016). https://www.nature.com/articles/ncomms12820
- `docs/analysis_curva_bop_MJ0.md`, `docs/analysis_pseudopotencial_fermi_krb.md`.
