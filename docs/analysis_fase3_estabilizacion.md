# Fase 3 — Método de estabilización (Hazi & Taylor, 1970)

**Fecha**: 2026-08-23
**Autor**: Javier Aguilera (con Claude Code)
**Relevancia**: tercer módulo de `docs/PLAN_nonadiabatic_dynamics.md`. Toma
la ecuación de canales acoplados de la Fase 2 y decide, de forma
automática, qué autovalores son estados ligados reales frente a artefactos
de la caja finita — el requisito directo de la Fase 4 (tasas de
decaimiento no adiabático), que sólo tiene sentido aplicarse a estados
ligados de verdad.

## Resumen

Se implementó `trimero.systems.nonadiabatic_dynamics.stabilization`, módulo
genérico (no conoce la ecuación de canales acoplados directamente, sólo un
`hamiltonian_builder(L) -> (R_grid, H)`):

1. `stabilization_scan`: diagonaliza en una malla de tamaños de caja L.
2. `track_trajectories`: empareja autovalores de cajas vecinas por energía
   más cercana (algoritmo húngaro), construyendo trayectorias continuas.
3. `classify_stability`: criterio AUTOMÁTICO (no visual) — compara la
   deriva medida dE/dL contra la deriva esperada de un estado de caja en esa
   misma energía y longitud, 2|E|/L (fórmula del pozo infinito).
4. `stable_state_summary`: decide por trayectoria si es un estado ligado
   real.

Se verificó, en este orden y ANTES de aceptar el módulo:

1. La propia referencia independiente (pozo cuadrado finito, ecuación
   trascendente) se contrastó contra la fórmula analítica de conteo
   ⌊k₀·2a/π⌋+1 — y esa verificación **atrapó un bug real en la primera
   versión de la referencia** (ver §Un bug encontrado antes de aceptar nada).
2. `track_trajectories` mantiene la identidad de cada rama a través de un
   cruce evitado sintético (nunca degenerado, hueco mínimo conocido).
3. Con la referencia corregida (3 estados ligados exactos, energías
   conocidas), el clasificador automático encuentra EXACTAMENTE 3 estados
   ligados, con energías que coinciden dentro de la discretización, y
   clasifica correctamente TODOS los estados de energía positiva (por
   encima del umbral de disociación) como estados de caja.

No se aplicó a datos reales de Rb*-RbCs en esta ronda — instrucción
explícita: primero hay que confiar en el método contra un caso conocido.

## Un bug encontrado antes de aceptar nada

La primera versión de `bound_states_square_well` (la referencia
independiente) escaneaba la energía en una malla y buscaba cambios de signo
de la ecuación trascendente `κ - k·tan(ka)` (par) / `κ + k·cot(ka)` (impar).
Con V0=0.1, a=8.0, μ=1.0 "encontró" **5** raíces:

```
E=-8.833047e-02 (par)
E=-8.072343e-02 (par)     <- espuria
E=-5.482889e-02 (impar)
E=-2.289372e-02 (impar)   <- espuria
E=-8.088245e-03 (par)
```

La fórmula analítica de conteo del pozo finito, N = ⌊k₀·2a/π⌋+1 con
k₀=√(2μV0), da **3** para estos parámetros (k₀=0.4472, k₀·2a/π=2.278,
⌊2.278⌋+1=3) — no 5. La discrepancia delató el bug: un polo de `tan`/`cot`
también produce un cambio de signo aparente en un escaneo ingenuo de
`np.diff(np.sign(...))`, así que dos de las cinco "raíces" eran cruces de
polo, no ceros reales de la ecuación. Corregido buscando cada raíz dentro de
una única rama monótona de tan/cot (entre polos consecutivos, acotando el
intervalo de búsqueda de `brentq` a esa rama), la referencia da exactamente
las 3 raíces reales, coincidiendo con la fórmula de conteo:

```
E = -8.833046956086295e-02  (par)
E = -5.482888535439445e-02  (impar)
E = -8.088245288672372e-03  (par)
```

Esto es exactamente la disciplina que pide el documento maestro: no se
aceptó el resultado del módulo de estabilización contra una referencia sin
verificar esa referencia de forma independiente primero (aquí, contra la
fórmula de conteo analítica).

## Método

### Criterio automático de estabilidad

Para un pozo infinito 1D de longitud L, Eₙ(L) ≈ n²π²ħ²/(2μL²), así que
dEₙ/dL = -2Eₙ/L: la escala de variación TÍPICA de un estado de caja en esa
energía y longitud. Un punto de trayectoria se clasifica ESTABLE si

```
|dE/dL| < umbral · (2|E|/L),   umbral = 0.1 por defecto
```

y una trayectoria completa se acepta como ligada si es estable en al menos
la mitad de los L escaneados Y en el L más grande. Es un criterio numérico,
reproducible, sin inspección visual de ningún diagrama.

### Emparejamiento de trayectorias

A diferencia de la Fase 1 (donde se sigue la fase de autovectores de
dimensión constante), aquí la dimensión de H crece con L, así que se seguye
por PROXIMIDAD EN ENERGÍA entre cajas vecinas (algoritmo húngaro,
`scipy.optimize.linear_sum_assignment`), válido por la regla de no-cruce:
autovalores de la misma simetría, dependientes de un único parámetro
continuo, no se cruzan exactamente — con un paso de L fino, la trayectoria
más cercana en energía es la física correcta.

## Resultados numéricos reales

Suite completa (`poetry run pytest tests/systems/nonadiabatic_dynamics/test_stabilization.py -v -s`):

```
tests/systems/nonadiabatic_dynamics/test_stabilization.py::test_square_well_reference_matches_analytic_count_formula PASSED
tests/systems/nonadiabatic_dynamics/test_stabilization.py::test_track_trajectories_keeps_identity_through_avoided_crossing PASSED
tests/systems/nonadiabatic_dynamics/test_stabilization.py::test_classify_stability_rejects_nonuniform_l_grid PASSED
tests/systems/nonadiabatic_dynamics/test_stabilization.py::test_stabilization_identifies_exact_number_and_energy_of_bound_states
  estados ligados encontrados: 3 (referencia: 3)
    numérico=-8.805280e-02  exacto=-8.833047e-02  dif=2.777e-04  fracción_estable=1.000
    numérico=-5.382147e-02  exacto=-5.482889e-02  dif=1.007e-03  fracción_estable=1.000
    numérico=-6.691792e-03  exacto=-8.088245e-03  dif=1.396e-03  fracción_estable=0.831
PASSED

============================== 4 passed in 1.16s ===============================
```

Suite completa del proyecto (`poetry run pytest -m "not slow"`): **118
passed, 16 deselected** — nada se rompió (114 de la Fase 1+2, +4 de esta
ronda).

### El caso de prueba: pozo cuadrado finito

μ=1.0, V0=0.1 Eh, semi-anchura a=8.0 a₀ (3 estados ligados exactos). Barrido
de caja: L∈[30,150] a₀, paso 2.0 a₀ (61 tamaños de caja), malla espacial de
paso fijo h=0.25 a₀ (creciendo en número de puntos con L), `n_keep=12`
autovalores más bajos guardados en cada L.

Diferencias numérico vs. exacto: 2.8×10⁻⁴, 1.0×10⁻³, 1.4×10⁻³ — creciendo
para el estado más somero, coherente con que su función de onda es más
extendida y más sensible a la discretización espacial h=0.25 (no se hizo
un test de convergencia en h en esta ronda; queda como trabajo natural si
se necesita más precisión en una fase posterior).

Todas las trayectorias de energía POSITIVA (por encima del umbral de
disociación E=0; hay 9 en la lista de `n_keep=12`, ya que 3 son las
ligadas) se clasificaron correctamente como NO ligadas
(`fraction_stable=0.0` en todas, verificado explícitamente en el test).

## Conclusiones y aplicación al proyecto

- El método de estabilización, con un criterio de clasificación
  completamente automático (sin inspección visual), reproduce EXACTAMENTE
  el número de estados ligados de un caso de referencia independiente, con
  energías que coinciden dentro del error esperado de la discretización
  espacial.
- El emparejamiento por energía más cercana (húngaro) es robusto a través
  de cruces evitados (verificado con un modelo analítico donde el hueco
  mínimo es conocido): no confunde identidades entre ramas que se acercan
  mucho sin fundirse.
- Encontrar y corregir el bug de la propia referencia (§arriba) antes de
  confiar en cualquier resultado del módulo es, en sí mismo, el resultado
  más importante de esta ronda desde el punto de vista de la disciplina del
  proyecto: confirma que "verificar la referencia independiente contra algo
  TODAVÍA más independiente" (aquí, la fórmula analítica de conteo) no es
  un paso opcional.
- Fase 3 queda cerrada para lo que pedía esta ronda: criterio automático +
  verificación contra un caso completamente conocido. Deliberadamente NO se
  aplicó a los datos reales de Rb*-RbCs — eso es la Fase 4 en adelante,
  cuando además haga falta decidir qué estados ligados alimentan el cálculo
  de tasas de decaimiento.

## Archivos tocados

- `src/trimero/systems/nonadiabatic_dynamics/stabilization.py` (nuevo)
- `tests/systems/nonadiabatic_dynamics/test_stabilization.py` (nuevo, 4 tests)
- `docs/PLAN_nonadiabatic_dynamics.md` (ESTADO ACTUAL: Fase 1-3 HECHAS)
- `docs/analysis_fase3_estabilizacion.md` (este documento)

## Referencias

- Hazi, A. U. & Taylor, H. S., *Phys. Rev. A* **1**, 1109 (1970) — método de
  estabilización.
- Griffiths, D. J., *Introduction to Quantum Mechanics* — pozo cuadrado
  finito, ecuación trascendente y fórmula de conteo de estados ligados.
- `docs/analysis_fase2_canales_acoplados.md` — `build_coupled_hamiltonian`,
  reutilizado aquí sin modificar para ensamblar cada caja del barrido.
