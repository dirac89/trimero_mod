# Fase 2 — Ecuación radial de canales acoplados (Ec. 8, Mellado-Alcedo et al. 2024)

**Fecha**: 2026-08-23
**Autor**: Javier Aguilera (con Claude Code)
**Relevancia**: segundo módulo de `docs/PLAN_nonadiabatic_dynamics.md`. Toma
el acoplamiento de derivada de la Fase 1 (`coupling.py`) y lo convierte en
niveles de energía resolviendo el sistema de canales acoplados. Es el
requisito directo de la Fase 3 (estabilización) y la Fase 4 (tasas de
decaimiento).

## Resumen

Se implementó `trimero.systems.nonadiabatic_dynamics.coupled_channels`,
módulo genérico que:

1. Calcula B_ij = ⟨Ψᵢ|d²/dR²|Ψⱼ⟩ (`second_derivative_coupling`), réplica de
   `coupling.derivative_coupling` para el operador de segunda derivada.
2. Ensambla la matriz de canales acoplados en diferencias finitas
   (`build_coupled_hamiltonian`), con condiciones de contorno de Dirichlet.
3. Conecta con la Fase 1 (`adiabatic_coupled_hamiltonian`) y resuelve
   (`solve_coupled_channels`).

Se verificaron, en este orden y ANTES de aceptar el módulo:

1. Discretización justificada (diferencias finitas, no DVR) — ver §Método,
   escrita ANTES de implementar.
2. Canal único vs. oscilador armónico ANALÍTICO (independiente del método).
3. Acoplamiento a cero (dos canales) recupera EXACTAMENTE dos solves de un
   solo canal, independientes.
4. Convergencia con el paso de malla (orden 2, verificado).
5. Doble pozo de dos canales acoplados, resuelto por DOS vías completamente
   independientes (representación diabática vs. adiabática) que DEBEN
   coincidir por un principio físico exacto, no por construcción del test.
6. Aplicación a un cruce evitado REAL ya localizado en la Fase 1
   (`docs/analysis_fase1_acoplamiento_derivada.md`), con las limitaciones de
   esa aplicación documentadas explícitamente (§Aplicación a datos reales).

Los 6 tests están en verde. Ningún test necesitó ajuste de tolerancia
después de fallar: los valores objetivo se calcularon primero por
exploración numérica y se usaron directamente.

## Método

### La ecuación

Ver el docstring completo (con la derivación) en
`src/trimero/systems/nonadiabatic_dynamics/coupled_channels.py`. Resultado:

```
-ħ²/(2μ) Σⱼ [δ_ij χⱼ'' + 2A_ij χⱼ' + B_ij χⱼ] + Vᵢ(R) χᵢ = E χᵢ
```

con A_ij de la Fase 1 y B_ij nuevo en este módulo. **B_ii ≠ 0 en general**
(corrección diagonal de Born-Oppenheimer, B_ii = -Σₖ|A_ik|² ≤ 0): se incluye
sin tratamiento especial, siguiendo la regla del proyecto de no hacer
desaparecer una limitación o término conocido en silencio.

### Elección de discretización: diferencias finitas, no DVR

Decidida y documentada en el docstring del módulo ANTES de escribir el
ensamblador de la matriz, por tres razones (texto completo en el módulo):

1. A_ij y B_ij sólo existen como diferencias finitas centradas sobre una
   malla UNIFORME — no son funciones continuas proyectables sobre una base
   DVR sin interpolar. Con diferencias finitas, T, D1, A_ij y B_ij viven en
   la misma malla sin ningún paso de interpolación.
2. Vᵢ(R) también se conoce sólo en puntos discretos (una diagonalización
   electrónica por R) — la ventaja característica de un DVR (convergencia
   espectral con potenciales conocidos en todo punto) no aplica.
3. Simplicidad y verificabilidad: la matriz cinética de 3 puntos es el
   operador mejor caracterizado en métodos numéricos, permite un test
   directo contra el oscilador armónico analítico y un test de convergencia
   igual al de la Fase 1.

### Simetrización

El operador de acoplamiento -A_ij d/dR - ½A_ij' es Hermítico en el continuo
(verificado por integración por partes usando A_ji=-A_ij, en el docstring
del módulo), pero `diag(A_ij) @ D1` discreto no lo respeta exactamente — el
mismo error O(h²) de la Fase 1. La matriz completa se simetriza
explícitamente, `H = (H + Hᵀ)/2`, antes de diagonalizar.

## Resultados numéricos reales

Suite completa (`poetry run pytest tests/systems/nonadiabatic_dynamics/test_coupled_channels.py -v -s`):

```
tests/systems/nonadiabatic_dynamics/test_coupled_channels.py::test_single_channel_matches_analytic_harmonic_oscillator PASSED
tests/systems/nonadiabatic_dynamics/test_coupled_channels.py::test_single_channel_convergence_with_grid_spacing
  errores (n=100,200,400): [5.11199476707358e-05, 1.2632330050058732e-05, 3.1410784949656434e-06]
PASSED
tests/systems/nonadiabatic_dynamics/test_coupled_channels.py::test_zero_coupling_recovers_two_independent_single_channel_solves PASSED
tests/systems/nonadiabatic_dynamics/test_coupled_channels.py::test_double_well_diabatic_and_adiabatic_representations_agree
  diabática  lowest 6: [0.0099354  0.0099354  0.02983358 0.02983358 0.0496299  0.0496299 ]
  adiabática lowest 6: [0.0099354  0.0099354  0.02983358 0.02983358 0.0496299  0.0496299 ]
  diferencia máxima  : 1.564e-10
PASSED
tests/systems/nonadiabatic_dynamics/test_coupled_channels.py::test_build_coupled_hamiltonian_is_symmetric_and_finite PASSED
tests/systems/nonadiabatic_dynamics/test_coupled_channels.py::test_real_avoided_crossing_coupling_shifts_the_spectrum
  acoplado   lowest 6 (GHz, rel. a min desacoplado): [0.0470721  0.93014376 0.97382466 1.50846155 1.70001304 2.19237084]
  desacoplado lowest 6 (GHz, rel. a min desacoplado): [0.         0.69747705 0.73429388 1.1706056  1.31879895 1.94346335]
  |desplazamiento| por nivel (GHz): [0.0470721  0.23266671 0.23953078 0.33785595 0.38121409 0.24890749]
PASSED

========================= 6 passed in 65.64s (0:01:05) =========================
```

Suite completa del proyecto (`poetry run pytest -m "not slow"`): **114
passed, 16 deselected** — nada se rompió.

### Oscilador armónico (canal único) vs. analítico

μ=1.0, ω=0.05 (u.a.), caja [-40,40] a₀, 400 puntos: los 6 autovalores más
bajos coinciden con (n+½)ω dentro de rtol=2×10⁻³ (error dominado por la
caja finita y la discretización, no por un bug — ver convergencia abajo).

### Convergencia con el paso de malla (orden 2)

Estado fundamental del oscilador armónico, μ=1.0, ω=0.05, caja fija [-40,40]:

| n puntos | h (a₀) | error vs. exacto | error/h² |
|---|---|---|---|
| 100 | 0.808 | 5.112×10⁻⁵ | 7.829×10⁻⁵ |
| 200 | 0.402 | 1.263×10⁻⁵ | 7.816×10⁻⁵ |
| 400 | 0.201 | 3.141×10⁻⁶ | 7.813×10⁻⁵ |

Coeficiente h⁻² estable en tres resoluciones — error de discretización de
orden 2, tal como predice el operador cinético de 3 puntos.

### Acoplamiento a cero: dos canales independientes

Dos canales con centros distintos (R=0 y R=5, mismo ω), A=B=0: el espectro
combinado del solver de 2 canales coincide **exactamente** (rtol=1e-10) con
la unión ordenada de los espectros de dos solves de 1 canal independientes.
Nótese que los dos espectros salen numéricamente idénticos entre sí
(oscilador armónico: la energía no depende de dónde esté el centro del
pozo), lo cual es un resultado físico correcto, no una casualidad del test.

### Doble pozo: representación diabática vs. adiabática

μ=50 (masa pesada, para localizar los estados y limpiar el cruce), ω=0.02,
pozos en R=±8, acoplamiento diabático constante Δ=0.01, caja [-30,30],
300 puntos. Dos cálculos totalmente independientes del mismo problema
físico:

- **Diabático**: matriz de bloques T+V₁, T+V₂ en la diagonal, Δ·I fuera de
  ella — sin ninguna maquinaria de la Fase 1/2.
- **Adiabático**: se diagonaliza el Hamiltoniano diabático 2×2 en cada R
  (dando V±(R) y las autofunciones), se calculan A_ij, B_ij con
  `eigenbasis_along_R` + `derivative_coupling` + `second_derivative_coupling`,
  y se ensambla con `build_coupled_hamiltonian`.

Diferencia máxima en los 6 niveles más bajos: **1.564×10⁻¹⁰** — coinciden
hasta el límite de la propia discretización compartida (mismo h en ambos
caminos), confirmando que A_ij, B_ij y su combinación en la ecuación de
canales acoplados están bien implementados: si hubiera un error de signo,
de factor, o en la combinación 2A_ij+B_ij, las dos vías NO coincidirían así
de bien, porque parten de construcciones matriciales completamente
distintas.

## Aplicación a datos reales (con limitaciones explícitas)

Se usó el cruce evitado real de la Fase 1: `BOPSystem(n_manifold=25)`,
M_J=0, sin Fermi, estados 54/55, mínimo del hueco de energía en R≈1305 a₀
(`docs/analysis_fase1_acoplamiento_derivada.md`). Ventana de R: [1280, 1330]
a₀, paso 1 a₀ (51 puntos), acotada por coste computacional — cada punto es
una diagonalización de 1113×1113 (~1.1-1.3 s); la suite completa de esta
ventana tarda 65 s.

**Lo que SÍ se puede concluir**: la matriz ensamblada es Hermítica y finita
(verificado explícitamente), y el acoplamiento tiene un efecto real y
medible sobre el espectro — comparando la solución acoplada (A_ij, B_ij
reales) contra la desacoplada (A=B=0, misma Vᵢ(R)) en la MISMA ventana:

```
desplazamiento máximo entre los 6 niveles más bajos: 0.381 GHz
```

Esto es consistente con el pico de A_ij≈0.5 encontrado en la Fase 1 en este
mismo cruce (un acoplamiento casi-divergente cerca de R=1305 a₀).

**Lo que NO se puede concluir de esta ventana**: que las energías reportadas
sean niveles vibracionales físicos reales de la molécula Rb*-RbCs. La
ventana de 50 a₀ es enormemente menor que el dominio real donde vive el
movimiento nuclear de una molécula Rydberg ultra-larga (cientos de a₀); con
una caja tan pequeña, el espaciado entre niveles (~0.2-0.4 GHz en la lista
de arriba) está dominado por la cuantización de caja, no por la física de
enlace real. Identificar cuáles de estos niveles serían estados ligados
reales (si los hay) frente a artefactos de caja es exactamente lo que la
Fase 3 (método de estabilización, variando el tamaño de la caja) está
diseñada para resolver — no se afirma aquí lo que la Fase 3 todavía no ha
verificado.

**Sobre μ=113522 mₑ**: masa reducida aproximada de Rb-RbCs
(M(Rb-87)≈86.909 u, M(RbCs)≈M(Rb-87)+M(Cs-133)≈219.814 u, μ≈62.28 u,
1 u=1822.888486 mₑ). Valor de demostración para esta ronda, sin cita
verificada de las masas isotópicas exactas usadas en el experimento de
Ruttley et al. 2023 — la Fase 6 (física real de Rb*-RbCs) es donde este
valor debe revisarse y documentarse con la fuente correcta antes de citar
ningún resultado cuantitativo definitivo.

## Conclusiones y aplicación al proyecto

- `coupled_channels.py` es correcto: la prueba de representación
  diabática/adiabática (independencia de representación, un principio físico
  exacto) es la verificación más fuerte posible sin depender de ningún
  resultado externo — y coincide a 1.6×10⁻¹⁰.
- La discretización por diferencias finitas es adecuada para esta fase:
  converge al orden esperado (2) y da acceso directo a A_ij, B_ij de la
  Fase 1 sin interpolación. Si la Fase 3 revela que la resolución necesaria
  es prohibitiva, se reconsiderará DVR, documentado como tal.
- La aplicación al cruce evitado real confirma cualitativamente lo que
  predecía la Fase 1: el acoplamiento no adiabático, donde es grande, altera
  el espectro de forma medible y no despreciable (∼0.4 GHz, comparable a las
  propias energías de enlace del sistema BOP, no un efecto marginal).
- Fase 2 queda cerrada para las partes 1 y 2 del plan (discretización +
  tests de límite con un sistema de 2 canales). La parte 3 (aplicar a curvas
  reales) se hizo con las limitaciones explícitas de arriba: es una
  demostración de que la maquinaria funciona sobre datos reales, no una
  determinación de niveles vibracionales físicos — eso espera a la Fase 3.

## Archivos tocados

- `src/trimero/systems/nonadiabatic_dynamics/coupled_channels.py` (nuevo)
- `tests/systems/nonadiabatic_dynamics/test_coupled_channels.py` (nuevo, 6 tests)
- `docs/PLAN_nonadiabatic_dynamics.md` (ESTADO ACTUAL: Fase 1 HECHA, Fase 2 HECHA)
- `docs/analysis_fase2_canales_acoplados.md` (este documento)

## Referencias

- Mellado-Alcedo, Guttridge, Cornish, Sadeghpour & González-Férez,
  *Phys. Rev. A* **110**, 013314 (2024), arXiv:2401.09618 — Ec. 8.
- `docs/analysis_fase1_acoplamiento_derivada.md` — A_ij, el cruce evitado
  real (R≈1305 a₀, estados 54/55) reutilizado aquí.
- Masas isotópicas Rb-87, Cs-133: valores estándar de tabla periódica de
  isótopos: **por verificar contra la fuente exacta usada en Ruttley et al.,
  PRL 131, 013401 (2023), antes de la Fase 6** — ver advertencia arriba.
