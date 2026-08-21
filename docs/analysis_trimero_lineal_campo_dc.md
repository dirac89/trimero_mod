# Trímero Rydberg lineal simétrico en campo DC — validación contra Aguilera-Fernández et al. (2016)

**Fecha**: 2026-08-21
**Autor**: Javier Aguilera
**Relevancia**: Reactiva la línea del **perturbador neutro** aplicándola al paper
para el que se escribió el código legado. Es la primera vez que el
pseudopotencial de Fermi de este repositorio se compara contra *su* sistema
físico, y no contra Rb\*-KRb.
**Tipo**: analysis

---

## Resumen

El pseudopotencial de Fermi de `systems/rb_neutral_perturber/` se escribió para
Aguilera-Fernández, Schmelcher & González-Férez, *Ultralong-range triatomic
Rydberg molecules in an electric field*, **J. Phys. B 49, 124002 (2016)**
([arXiv:1601.05049](https://arxiv.org/abs/1601.05049)). Tras varias rondas en
las que se le intentó aplicar al sistema polar Rb\*-KRb —donde **no
interviene**— esta ronda lo lleva a su sitio: Rb\*(n=35, l≥3) con **dos átomos
neutros de Rb** en θ=0 y θ=π, más un campo eléctrico DC.

La ventaja decisiva es que **n=35 es el manifold nativo de las tablas**
`rvsAS.dat` / `rvsAP.dat`: su R máximo, 2448 a₀, coincide con el punto de
retorno clásico 2n² = 2450 a₀. No hace falta remapeo k(R), ni ventana de
exclusión de la resonancia p, ni límite de dominio: la tabla se lee tal cual y
además se evalúa **en sus propios nodos**, así que tampoco se interpola.

El resultado reproduce el paper con precisión inesperada. Seis anclas
cuantitativas del texto caen dentro del 3 %, incluidas dos que coinciden
prácticamente al alfiler: el cruce Π sin campo en **1060.5 a₀** (paper: ≈1060) y
el mínimo de la Π más baja en **1116 a₀** (paper: ≈1115).

## Palabras clave

- ULRM triatómica, trímero Rydberg
- pseudopotencial de Fermi, ondas s y p
- efecto Stark DC, límite de Inglis-Teller
- simetría gerade / ungerade, Σ y Π
- Rb(n=35, l≥3), 38s / 37p / 36d

---

## 1. Obtención de la fuente

`WebFetch` sobre `arxiv.org/abs/1601.05049` sólo devolvió el resumen, y sobre
`arxiv.org/pdf/1601.05049` devolvió el PDF binario. El PDF quedó guardado en
disco y se extrajo su texto con `pypdf` (8 páginas, 35 054 caracteres). **Se
dispone del texto completo**, no sólo del resumen. No hay tabla numérica en el
paper: todo lo comparable son valores sueltos citados en §III.A y las escalas
de los ejes de las Figs. 3-5.

### Erratum detectado en el propio paper

El texto y los pies de las Figs. 4-7 dicen **V/m**; las etiquetas dentro de los
paneles dicen **V/cm**. Son incompatibles y la correcta es **V/m**:

| criterio | V/m | V/cm |
|---|---|---|
| ensanchamiento Stark del manifold n=35 a F=500 | ±11.4 GHz | ±1142 GHz |
| eje de las Figs. 4-5 | −20 … 0 GHz ✓ | ✗ fuera de escala ×50 |
| límite de Inglis-Teller 1/(3n⁵) = 6.35·10⁻⁹ u.a. ≈ 3264 V/m | F=500 = 0.15·F_IT, subcrítico: abanicos Stark resueltos como en las Figs. 4-5 ✓ | ✗ F=500 V/cm = 15·F_IT, manifolds adyacentes solapados: esa estructura no sobreviviría |

Los cálculos de este documento usan **V/m**.

> **Corrección (2026-08-21)**: esta fila decía antes «1/(3n⁵) = 326 V/m …
> F=500 justo por encima». El número llevaba un deslizamiento de factor 10:
> con n=35 y 1 u.a. = 5.14220674763·10¹¹ V/m, F_IT = 6.3466·10⁻⁹ u.a. =
> **3264 V/m (32.6 V/cm)**; con la convención de semianchura igual a la
> separación inter-manifold (~147 GHz) el umbral sube a ≈6400 V/m. Consecuencia:
> **F=500 V/m está claramente POR DEBAJO del límite (×6.5), no justo por
> encima**, y bajo lectura V/cm estaría ×15 por encima. La elección V/m se
> mantiene por las dos primeras filas, que son mediciones directas. Los
> «several avoided crossings» del paper son cruces en R entre APC individuales
> y niveles vestidos (§III.A), no un efecto de solape de manifolds: este
> cálculo los reproduce a 0.15·F_IT sin contradicción. El paper no menciona
> Inglis-Teller en ningún punto (verificado contra el texto completo, ar5iv);
> nótese además que su propia §III.2 rotula «300 V/cm y 500 V/cm», mientras
> III.1 y los pies de figura dicen V/m.

## 2. El Hamiltoniano del paper y su correspondencia con el código

El paper (Ecs. 1-4):

```math
H = H_0 + \mathbf{F}\cdot\mathbf{r} + V(\mathbf{r},\mathbf{R}_1) + V(\mathbf{r},\mathbf{R}_2)

V(\mathbf{r},\mathbf{R}_i) = 2\pi A_s[k(R_i)]\,\delta^3(\mathbf{r}-\mathbf{R}_i)
   + 6\pi A_p^3[k(R_i)]\,\overleftarrow{\nabla}\delta^3(\mathbf{r}-\mathbf{R}_i)\overrightarrow{\nabla}

A_s(k) = -\tan\delta_0(k)/k,\qquad A_p^3(k) = -\tan\delta_1(k)/k^3

E_\mathrm{kin} = k^2/2 = 1/R_i - 1/2n^2,\qquad \mathbf{F} = F\hat{Z}
```

| término del paper | dónde está | estado |
|---|---|---|
| `H₀` (energías Rydberg con defecto cuántico) | `systems/rb_atom.py::Atom.E_Rb` | ya servía |
| `F·r` con **F ∥ Z** | `systems/rb_atom.py::Atom.Vfield` | **ya servía tal cual**, ver §2.1 |
| `V(r,R₁)` con R₁ en θ=0 | `rb_neutral_perturber/fermi_krb.py` | ya servía |
| `V(r,R₂)` con R₂ en θ=π | **nuevo**, `linear_trimer.py` | ver §3 |
| base n=35 + 38s/37p/36d | `basis/quantum.py::CoupledBasis(N_max=0)` | reutilizada |

### 2.1 `Vfield` no necesitaba generalizarse

`Atom.Vfield(li, lj, mi, mj, radial, strength)` calcula
`F · ⟨n l|r|n' l'⟩ · ⟨l m|cosθ|l' m'⟩`, que es exactamente `F·r = F·z` con
**F ∥ Z**, el único caso que trata el paper. El factor angular
`Angular_dc_field` se comprobó término a término contra
`⟨l,m|cosθ|l±1,m⟩ = √((l_>² − m²)/((2l_>−1)(2l_>+1)))`, y el test **T8** lo
valida de forma independiente: con un manifold hidrogenoide puro reproduce la
escalera Stark analítica `(3/2)·n·q·F`, `q = −(n−1), …, +(n−1)` de dos en dos,
a `rtol=1e-9`.

**No hubo que tocarlo.**

### 2.2 Elementos dipolares radiales

`data/Wavefunction/exp_val_r.txt` contiene ⟨l|r|l+1⟩ para l = 0..33, y su
estructura estaba sin documentar. Se ha determinado:

| índice | es | valor |
|---|---|---|
| 0 | ⟨38s\|r\|37p⟩ | +1278.215 |
| 1 | ⟨37p\|r\|36d⟩ | −1626.510 |
| 2 | ⟨36d\|r\|35f⟩ | +1614.194 |
| l ≥ 3 | ⟨35,l\|r\|35,l+1⟩ | −(3/2)·35·√(35²−(l+1)²) |

Los valores con l ≥ 3 coinciden con la forma cerrada hidrogenoide **incluido el
signo** (que es el convenio de `fermi_krb.hydrogenic_R`, luego el mismo que usa
V_Fermi), pero sólo a 3.3·10⁻⁴ relativo en el peor caso (l=7): residuo de
integración numérica del generador original. El módulo usa por eso la **forma
cerrada** para l ≥ 3 y reserva el fichero para los tres cruzados, que no tienen
forma cerrada porque sus n\* no son enteros.

## 3. Simetría real del sistema: **m_l**, no M_J

Ésta era la pregunta del punto 2 del encargo, y la respuesta condicionaba todo
lo demás. **El buen número cuántico es m_l.** Razones:

1. **No hay rotor.** M_J = m_l + M_N era el buen número en Rb\*-KRb porque el
   dímero polar tenía momento angular rotacional N. Aquí los perturbadores son
   átomos neutros sin estructura interna: M_N ≡ 0 y M_J colapsa a m_l.
2. **El eje Z es único para todo.** Es a la vez el eje de los dos
   perturbadores (θ=0, π) y la dirección del campo. Los tres términos de H
   conservan m_l por separado:
   - V_Fermi con el perturbador **sobre el eje** sólo conecta m₁ = m₂: en θ=0/π
     los armónicos esféricos con |m|≥2 y sus gradientes se anulan;
   - `F·r = F·z` tiene Δm = 0, Δl = ±1.

Luego **H es bloque-diagonal en m_l incluso con campo**, y la nomenclatura
molecular del paper es literalmente m_l:

| símbolo | m_l | contenido |
|---|---|---|
| Σ | 0 | 35 estados (l = 0..34); s + p |
| Π | ±1 | 34 estados (l = 1..34); **sólo onda p** |
| Δ, Φ, … | \|m_l\| ≥ 2 | V_Fermi ≡ 0: sólo el manifold vestido por el campo |

Se implementa reutilizando `CoupledBasis(N_max=0)`, que degenera a estados
(l, m_l, 0, 0) y a bloques etiquetados por M_J = m_l. No es un apaño: es la
afirmación de que este sistema **es** el límite «sin rotor» de la base
compartida. Verificado en el test **T1**.

### 3.1 Los dos perturbadores se reducen a un factor de paridad

Como R⃗₂ = −R⃗₁ y ψ_{lm}(−r⃗) = (−1)^l ψ_{lm}(r⃗), también
(∇ψ)(−r⃗) = (−1)^{l+1}(∇ψ)(r⃗). Los dos términos del pseudopotencial son
bilineales en ψ o en ∇ψ, luego los dos signos extra del gradiente se cancelan:

```math
\langle l_1 m|V(\mathbf{R}_2)|l_2 m\rangle = (-1)^{l_1+l_2}\,
\langle l_1 m|V(\mathbf{R}_1)|l_2 m\rangle
\;\Longrightarrow\;
V_\text{total} = \left[1 + (-1)^{l_1+l_2}\right] V(\theta=0)
```

que vale **2·V(θ=0) si l₁+l₂ es par y CERO si es impar**. Ésa es exactamente la
separación gerade/ungerade que el paper cita de su Ref. [15]: sin campo el
bloque m_l se parte en l pares y l impares sin acoplamiento entre ellos, y el
término F·r (Δl = ±1) es el **único** que los conecta. De ahí que «the electric
field couples the adiabatic electronic states with gerade and ungerade
symmetry» (§III.A).

El test **T3** no da esto por bueno: evalúa el pseudopotencial por fuerza bruta
sobre ψ y ∇ψ **numéricos** en (0,0,±R), sin reutilizar ninguna forma cerrada, y
comprueba el factor a `rel=1e-6`. De paso valida las fórmulas de `fermi_krb`
contra la evaluación numérica directa (`rel=2e-4`).

## 4. Sin remapeo: n=35 lee la tabla nativa

Se comprobó lo que pedía el encargo. `ScatteringLengths` **sí** forzaba el
remapeo aunque no hiciera falta, por dos vías:

1. `remap_R_from_energy` reconstruía `R' = 1/(E + 1/R − E_table)`. Con
   `E == E_table` eso es la identidad algebraica `1/(1/R) = R`, pero la ida y
   vuelta en coma flotante movía el extremo R=2448 fuera de la tabla por 1 ulp
   y `scattering_from_energy` lo rechazaba. **Cortocircuitado**: si
   `E == E_table` devuelve `R` exactamente. Es evitar reconstruir un número que
   ya se tiene, no un cambio de física.
2. `FermiPseudopotential` usaba `scattering_pair`, es decir la **energía media
   del par (l₁,l₂)**, lo que metía los n\* de los vecinos (34.87, 34.35, 34.65)
   y volvía a activar el remapeo. El paper usa **un único k(R) por punto**, el
   del número cuántico principal Rydberg. Añadido `uniform_n_star`.

Con `uniform_n_star = n_star_table = 35` el remapeo es la identidad y, además,
las curvas se evalúan **en los nodos de la propia tabla** (`table_R()`), igual
que hace el camino legado (`R = As[row,0]`). Resultado: **cero interpolación de
A_s y A_p** en producción. Test **T9**.

## 5. Resultados

Comando:

```bash
poetry run python scripts/compute_trimer_curves.py --symmetry Sigma
poetry run python scripts/compute_trimer_curves.py --symmetry Pi
poetry run python scripts/compute_trimer_curves.py --symmetry Sigma --fields 0 --s-wave-only --rmin 500
```

483 nodos nativos, R = 1002 … 2448 a₀, F = 0, 100, 300, 500 V/m
(0, 1.945e-10, 5.834e-10, 9.724e-10 u.a.). Salidas en
`plots/trimer_lineal_{Sigma,Pi}_n35.{npz,png}`.

### 5.1 Comparación CUANTITATIVA con el texto

Todo lo que el paper da como número, y lo que sale:

| ancla (cita del paper) | paper | calculado | desv. |
|---|---|---|---|
| Rb(38s) por debajo del manifold (Fig. 3(b), eje −20…0) | ≈ −20 GHz | **−20.267 GHz** | — |
| «the crossing of the field-free APCs at R ≈ 1060 a₀» (Π) | 1060 a₀ | **1060.5 a₀** | 0.05 % |
| «the minimum appearing at R ≈ 1115 a₀ for the lowest lying Π-APC» | 1115 a₀ | **1116 a₀** | 0.09 % |
| ese mínimo «is shifted 0.3 GHz … for F = 500 V/m» | 0.3 GHz | **−0.244 GHz** | 19 % |
| «several avoided crossings close to R ≈ 1500 a₀» (Σ) | ≈1500 a₀ | **1503 / 1524 / 1557 a₀** | ≤ 4 % |
| Fig. 3(a) sólo onda s, eje −10…0 GHz | fondo −10 | **−9.871 GHz** | apura el eje |
| Fig. 3(c) Π s+p, eje −35…0 GHz | fondo −35 | **−34.267 GHz** | apura el eje |
| Stark lineal del manifold a R grande, ±(3/2)n(n−1−\|m\|)F | ±11.08 GHz (Π, 500 V/m) | **[−11.02, +10.97]** | ≤ 2.2 % |
| «quadratic Stark shift of the 38s» | ∝ F² | **1 : 9.0 : 25 (0.011/0.098/0.272 GHz)** | < 2 % |
| «resonance of the p-wave scattering length at R ≈ 780 a₀» | 780 a₀ | **759.3 a₀** (polo de la tabla) | **2.7 %** |

Observaciones sobre las dos filas que no cierran del todo:

* **Desplazamiento Stark del mínimo Π (0.244 vs 0.3 GHz).** El paper da una sola
  cifra significativa y no da el signo. El cálculo lo baja (más ligado), que es
  lo que hace la repulsión de niveles sobre un estado situado por debajo del
  manifold. Compatible.
* **Resonancia de onda p (759 vs 780 a₀).** Es una propiedad de **la tabla**,
  no del Hamiltoniano: el polo de `rvsAP.dat` está en 759.3 a₀, con FWHM
  operativa [747, 771] a₀. 780 queda fuera de esa ventana. O el paper redondea
  a la baja de forma generosa, o las tablas del repositorio no son exactamente
  las que produjeron las figuras. **No se ha resuelto**; queda registrado en
  `test_resonancia_de_onda_p_de_la_tabla`.

### 5.2 Comparación CUALITATIVA con las Figs. 3-5

Cada afirmación del paper y lo que hace el cálculo:

| afirmación (§III.A) | ¿se reproduce? |
|---|---|
| sólo onda s: dos APC (gerade/ungerade) que oscilan alrededor de la del dímero y convergen a ella | **sí**, Fig. 3(a) |
| con onda p: **dos APC Σ adicionales** se separan del manifold (4 en total) | **sí**, exactamente 4 |
| «their slope becomes pronounced for R ≲ 1200 a₀» por la resonancia p | **sí**: de −16.9 GHz en 1200 a₀ a −41.6 GHz en 1002 a₀ |
| dos APC sufren cruces evitados con el estado Rb(5s)Rb(38s)Rb(5s) | **sí**, en R ≈ 1.1–1.2 ·10³ a₀ |
| la APC del 38s «remains approximately constant for larger values of R» | **sí**, plana en −20.28 GHz |
| a R grande los pares s-dominado y p-dominado **degeneran** y convergen al dímero | **sí**: en 2400 a₀ los gaps son 0.074 y 0.001 GHz |
| Π: **la onda p es la responsable** de esas APC | **sí**: con `--s-wave-only` el bloque Π queda vacío |
| Π: dos APC con pozos que admiten varios estados vibracionales | **sí**, pozo de −34.3 GHz en 1116 a₀ |
| el campo **rompe la degeneración** gerade/ungerade a R grande | **sí**, test T6 |
| el campo saca **APC adicionales** del manifold (acopla l y l±1) | **sí**, el abanico Stark |
| F=100: una sube y otra baja; F≥300 **todas** bajan | **sí** |
| F=500: dos APC se **funden** con el abanico vestido, quedan sólo dos bien separadas | **sí** |
| el cruce evitado con el 38s **se ensancha** al crecer F | **sí** |
| las Π son «weakly affected» para 1000 ≲ R ≲ 1250, con Stark cuadrático apenas visible | **sí**: los dos pozos profundos son indistinguibles entre los cuatro paneles |
| a R grande **todas** las APC del manifold se desplazan **linealmente** en F | **sí**, §5.1 |

No se ha calculado la configuración **asimétrica** (θ₁=θ₂=π, R₁\*=1200 a₀ fijo)
ni la **planar** (θ₂ = π−θ₁, φ=0), que son §III.B y §IV del paper. La geometría
está encapsulada en `parity_factor`: extenderla es sustituir ese factor por los
armónicos esféricos evaluados en cada θᵢ, con A_s/A_p distintos para cada R_i.

## 6. Aproximaciones explícitas y su sensibilidad

| # | aproximación | efecto medido |
|---|---|---|
| 1 | radiales de 38s/37p/36d: se usan las **tabuladas** (spline cúbico) en vez de hidrogenoides de n entero | cambiar a hidrogenoides mueve Σ hasta **2.6 GHz** en el peor punto a 500 V/m (mediana 0.005 GHz); Π es insensible (0.002 GHz) |
| 2 | el **signo relativo** entre las radiales tabuladas y los tres dipolos cruzados de `exp_val_r.txt` se asume consistente | sin campo es un **gauge exacto**: 0.000 GHz. Con campo, hasta **1.16 GHz** en el peor punto de Σ a 500 V/m (mediana 0.006 GHz); Π < 0.01 GHz |
| 3 | rango R ≤ 2448 a₀ | las figuras del paper llegan a 2600 a₀; ese tramo exigiría extrapolar A_s, A_p al límite k→0 y **no se ha hecho** |
| 4 | defecto cuántico del 35f despreciado | es la misma decisión del paper («We neglect through the quantum defect of the 35f Rydberg state») |

Las aproximaciones 1 y 2 sólo muerden en los cruces evitados de Σ con campo, no
en las anclas de §5.1 (todas Π o de campo nulo). Cuantificadas en
`test_el_convenio_de_signo_de_los_vecinos_es_gauge_sin_campo`.

## 7. Bugs encontrados en los datos y el código del legado

Se documentan aquí y **no se propagan** a la capa moderna. No se corrigen en
`trimer.py`/`fermi_potentials.py`, que están congelados por los golden files.

1. **`rvsDR38s.dat` vale 10⁻⁶ veces la derivada real de `rvsR38s.dat`.**
   Comprobado contra `np.gradient`: la razón tabla/numérica es 1.0011e−6, frente
   a 1.0010 para `rvsDR37p` y `rvsDR36d`, que sí son correctas. En `trimer.py`
   la contribución de onda p del 38s entra al cuadrado, luego está suprimida por
   10⁻¹²: es efectivamente **cero**.
2. **Desalineación de mallas.** `rvsR38s.dat` y `rvsDR38s.dat` tienen 780 filas
   (malla uniforme 111..2448, paso 3), mientras `rvsAS`, `rvsAP`, `rvsR37p` y
   `rvsR36d` tienen 776: les faltan las 4 filas del hueco 750→765 a₀, sobre la
   resonancia p. `trimer.py` las indexa **por número de fila**, así que a partir
   de la fila 214 lee el 38s con un desfase de 12 a₀ (en la fila 297, donde
   arranca su bucle: R=1014 para todo y R=1002 para el 38s).
3. **`EhtoGHz = 6.579683920729e9` en `trimer.py` es 10³ demasiado grande.**
   1 E_h = 6.5797e15 Hz = 6.5797e6 GHz. Ya estaba documentado como bug congelado
   en `.claude/CLAUDE.md`; aquí queda el valor correcto,
   `HARTREE_TO_GHZ = 6.579683920502e6`.

La capa moderna evita los tres: interpola por valor de R (no por índice), saca
la derivada del spline de la función de onda, y usa la constante CODATA 2018.

---

## Conclusiones y aplicación al proyecto

1. **La línea del perturbador neutro vuelve a estar vigente.** El
   pseudopotencial de Fermi de este repositorio reproduce el paper para el que
   se escribió, con seis anclas cuantitativas dentro del 3 % y catorce
   afirmaciones cualitativas confirmadas. Deja de ser código «verde y
   congelado».
2. **`Vfield` de `rb_atom.py` servía tal cual**: es F·z con F ∥ Z, el único
   caso del paper, y queda validado de forma independiente por la escalera
   Stark analítica.
3. **La simetría de este sistema es m_l**, no M_J. Σ ≡ m_l=0, Π ≡ |m_l|=1, y
   |m_l| ≥ 2 no siente el perturbador en absoluto.
4. **Dos perturbadores simétricos = un factor de paridad.** Reduce el sistema a
   un solo elemento de matriz por par (l₁,l₂), y explica de dónde sale la
   estructura gerade/ungerade.
5. **n=35 no necesita remapeo.** Se ha implementado el atajo (identidad exacta)
   y el `uniform_n_star` del convenio del paper. En producción no se interpola
   ni A_s ni A_p: se evalúa en los nodos de la tabla.
6. **Siguiente paso natural**: las geometrías asimétrica y planar (§III.B y §IV
   del paper), sustituyendo `parity_factor` por los armónicos esféricos
   evaluados en cada θᵢ. También, extrapolar A_s/A_p al límite k→0 para llegar
   a los 2600 a₀ de las figuras.

### Qué se añadió

| fichero | qué |
|---|---|
| `src/trimero/systems/rb_neutral_perturber/linear_trimer.py` | **nuevo**: `SymmetricLinearTrimer`, `NeighborRadials`, `field_au` |
| `src/trimero/systems/rb_neutral_perturber/fermi_krb.py` | **aditivo**: `uniform_n_star`, `radial_fn`/`dradial_fn`, atajo identidad en `remap_R_from_energy` |
| `scripts/compute_trimer_curves.py` | **nuevo**: script de producción |
| `tests/systems/rb_neutral_perturber/test_linear_trimer.py` | **nuevo**: 18 tests analíticos (T1-T9) |
| `tests/systems/rb_neutral_perturber/test_regression_trimer_2016.py` | **nuevo**: 11 anclas del paper, marcadas `slow` |

Los 48 tests previos siguen pasando, **golden files incluidos**
(`66 passed in 472.72s` antes de añadir los de regresión). El sistema polar
`rb_krb_polar/` no se ha tocado.

## Referencias

- [arXiv:1601.05049](https://arxiv.org/abs/1601.05049) — Aguilera-Fernández,
  Schmelcher & González-Férez, *Ultralong-range triatomic Rydberg molecules in
  an electric field*, J. Phys. B **49**, 124002 (2016). Texto completo obtenido
  y extraído con `pypdf`.
- Ref. [15] del paper — orbitales adaptados por simetría, origen de la
  descomposición gerade/ungerade en sumas sobre l par / impar.
- Ref. [16] del paper — el segundo átomo hace posible un trímero ligado.
- [`docs/archive/rb_neutral_perturber/analysis_pseudopotencial_fermi_krb.md`](archive/rb_neutral_perturber/analysis_pseudopotencial_fermi_krb.md)
  — la derivación de las formas cerradas de `fermi_krb`.
- [`docs/archive/rb_neutral_perturber/analysis_procedencia_rvsAS_rvsAP.md`](archive/rb_neutral_perturber/analysis_procedencia_rvsAS_rvsAP.md)
  — la procedencia de las tablas sigue sin verificarse; §5.1 añade el dato de
  que su polo p está en 759 a₀ y no en los 780 a₀ que cita el paper.
- [`docs/archive/rb_neutral_perturber/analysis_interpolacion_polo_Ap.md`](archive/rb_neutral_perturber/analysis_interpolacion_polo_Ap.md)
  — interpolación de 1/A_p; **aquí no se usa**, porque se evalúa en los nodos.

## Notas adicionales

- Las figuras van a `plots/`, que **sí** se commitea (ver `.claude/CLAUDE.md`).
- El paper trae un erratum V/m ↔ V/cm entre texto y etiquetas de figura; §1.
- El barrido completo de Σ o Π (483 nodos × 4 campos) cuesta ~40 s.
