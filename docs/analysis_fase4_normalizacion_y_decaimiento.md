# Fase 4 — Normalización en energía y tasas de decaimiento: RECONSTRUCCIÓN PROPIA

**Fecha**: 2026-08-23
**Autor**: Javier Aguilera (con Claude Code)
**Relevancia**: cierra la Fase 4 de `docs/PLAN_nonadiabatic_dynamics.md`
(tasas de decaimiento no adiabático), el último ingrediente antes de poder
aplicar el pipeline completo (Fases 1-4) a un cruce evitado real de
Rb*-RbCs (Fase 6).

## ⚠️ Esto NO es una réplica de González-Férez, Weidemüller & Schmelcher 2007

`docs/analysis_fase4_busqueda_referencia.md` documenta el intento (6 vías,
todas sin éxito) de conseguir el texto de González-Férez, Weidemüller &
Schmelcher, *Phys. Rev. A* **76**, 023402 (2007) — la referencia que el
paper de 2024 cita, sin reproducir ninguna ecuación, para el método de
normalización en energía de la función de onda del continuo. El usuario
decidió explícitamente reconstruir el método por cuenta propia a partir de
un principio físico general (teoría de colisiones estándar) en vez de
bloquear la Fase 4 indefinidamente. Todo lo que sigue es esa
reconstrucción, verificada contra un caso analítico cerrado — **no** se ha
contrastado contra el método real del paper de 2007. Si en el futuro se
consigue acceso a ese texto, **queda pendiente verificar que ambos métodos
coinciden** (no se asume que sea así sólo porque el resultado reconstruido
pasa sus propios tests).

## Resumen

Dos módulos nuevos:

1. **`energy_normalization.py`**: convierte autoestados de caja
   (L²-normalizados) a normalización en energía, χ_E = χₙ/√(ΔEₙ), derivado
   desde la resolución de la identidad y un argumento de densidad de
   estados (docstring completo del módulo). Verificado contra la solución
   analítica CERRADA de la partícula libre 1D radial,
   u_E(R)=√(2μ/(πk))sin(kR): **convergencia O(h²) limpia** (factores 4.00 y
   4.00 al doblar la resolución), igual que el resto del proyecto.

2. **`decay_rates.py`**: Γᵢ=(2π/ħ)|⟨χᵢ^d|A_du|χⱼ^u⟩|² (Ec. 15 del paper
   2024, ésta sí completa y verificada — ver
   `docs/analysis_fase4_busqueda_referencia.md` §"Lo que SÍ se consiguió").
   Verificado contra un modelo de juguete tipo Breit-Wigner con referencia
   analítica por cuadratura de funciones cerradas (pozo cuadrado radial +
   continuo libre desplazado): **0.128% de diferencia** con la referencia
   analítica, y **Γ∝A0² exacto** (9 cifras significativas, identidad
   algebraica de la fórmula).

Los 9 tests están en verde (127 en la suite completa del proyecto, +9 sobre
el cierre de la Fase 3). Un test falló en la primera corrida por un umbral
mal calibrado (ver §5), corregido tras verificar que la convergencia
subyacente era correcta.

## 1. Derivación del método de normalización en energía

Ver el docstring completo (con ecuaciones numeradas) en
`energy_normalization.py`. Resumen del argumento: los autoestados de caja
son ortonormales discretos y satisfacen Σₙ|χₙ⟩⟨χₙ|=1̂; contando estados por
intervalo de energía (dn=ρ(E)dE) esa suma se reescribe como
∫dE·ρ(E)|χ_{n(E)}⟩⟨χ_{n(E)}|; comparando con la resolución de la identidad
del continuo, ∫dE|χ_E⟩⟨χ_E|=1̂, se obtiene

    χ_E = χₙ/√(ΔEₙ),    ΔEₙ = espaciado LOCAL entre autovalores vecinos

(no la deriva dE/dL de la Fase 3 — aquí L está fija, se deriva en el índice
n del nivel).

## 2. Un problema real de convención de malla (detectado antes de aceptar nada)

La primera versión de la verificación usaba `R = linspace(0, L, n)`
(incluyendo R=0 como punto de la malla) y comparaba contra u_E(R) — el
error caía sólo como **O(h)**, no O(h²):

```
h=0.200  err_max(R<150)=5.9980e-02
h=0.100  err_max(R<150)=2.9995e-02   (factor 2, no 4)
h=0.050  err_max(R<150)=1.4999e-02   (factor 2, no 4)
```

Causa: `radial_kinetic_matrix`/`build_coupled_hamiltonian` imponen
Dirichlet en el vecino FICTICIO inmediatamente fuera del array, no en el
primer punto real de la malla. Con `R=linspace(0,L,n)`, la pared física
queda en R=-h, un paso de malla desplazada de donde la fórmula analítica
asume la pared (R=0) — un desfase de fase de orden h que se acumula
linealmente con la distancia. Verificado extrayendo k directamente del
autovector (identidad de recurrencia de una función seno pura) y
comprobando que el propio autovector es exactamente
v[j] ∝ sin((j+1)·kh), NO sin(j·kh): el índice j=0 del array no corresponde
a R=0 sino a R=h después de la pared.

Corrección: armar la malla como `R = h, 2h, ..., n·h` (SIN el punto R=0),
que sitúa la pared física exactamente en R=0. Con esta convención, el error
vuelve a caer O(h²) limpio:

```
h=0.200  err_max(R<150)=1.4150e-04  err/h²=3.538e-03
h=0.100  err_max(R<150)=3.5383e-05  err/h²=3.538e-03
h=0.050  err_max(R<150)=8.8468e-06  err/h²=3.539e-03
```

Coeficiente h⁻² estable en tres valores — es discretización de orden 2, no
un bug. Esta convención de malla queda documentada explícitamente en el
docstring de `energy_normalization.py` y se aplicó igual en la caja de
`decay_rates.py`.

## 3. Verificación de `energy_normalize_box_states`

Caja libre radial (V=0), μ=1.0, L=600 a₀, E₀ objetivo=0.01 Eh. Región de
comparación R∈(3,150) a₀ (R≪L, para que la aproximación de caja al continuo
sea válida).

```
tests/systems/nonadiabatic_dynamics/test_energy_normalization.py::test_local_level_spacing_matches_infinite_box_formula PASSED
tests/systems/nonadiabatic_dynamics/test_energy_normalization.py::test_local_level_spacing_rejects_too_short_array PASSED
tests/systems/nonadiabatic_dynamics/test_energy_normalization.py::test_energy_normalization_matches_analytic_free_particle
  errores vs h=(0.2,0.1,0.05): [1.3433137870398904e-03, 3.359282117201823e-04, 8.399149481674861e-05]
PASSED
tests/systems/nonadiabatic_dynamics/test_energy_normalization.py::test_energy_normalize_rejects_nonuniform_grid PASSED
```

Razones sucesivas: 4.00× y 4.00× — orden 2 exacto.

`local_level_spacing` se verificó, además, contra la fórmula analítica de
la caja infinita ΔEₙ=πk/(μL) (derivada de Eₙ=n²π²/(2μL²)): coincide dentro
de rtol=2e-3 — referencia distinta de la propia diferencia central que el
código calcula, no una comparación circular.

## 4. Verificación de `decay_rates.py` (modelo de juguete Breit-Wigner)

Canal "d": pozo cuadrado RADIAL (pared dura en R=0, -V0 para R<a),
μ=1.0, V0=0.1 Eh, a=8.0 a₀ — mismos parámetros que la Fase 3, restringido a
R≥0 (como toda R física del proyecto), con **un único** estado ligado
(rama impar de la ecuación trascendente ya verificada en la Fase 3):
E_d=-0.05482888535439445 Eh. Canal "u": partícula libre desplazada,
V_u=-0.15 Eh (por debajo de E_d, condición necesaria para predisociación).
Acoplamiento: A_du(R)=A0 constante, B_du=0 (para que la integral de
solapamiento analítica sea tratable en forma cerrada).

Referencia analítica: Γ = 2π|-(1/μ)A0·∫ψ_d(R)·dχ_u/dR dR|², con ψ_d(R) la
función de onda del pozo cuadrado en forma CERRADA (seno dentro,
exponencial fuera, normalización analítica sin cuadratura) y χ_u(R) la
solución libre cerrada — la integral de solapamiento sí se calcula por
cuadratura (`scipy.integrate.quad`), un camino numérico completamente
distinto (integración de funciones conocidas) del método bajo prueba
(ensamblar y diagonalizar Hamiltonianos de diferencias finitas).

```
tests/systems/nonadiabatic_dynamics/test_decay_rates.py::test_toy_well_reproduces_known_bound_energy PASSED
tests/systems/nonadiabatic_dynamics/test_decay_rates.py::test_decay_rate_matches_analytic_quadrature_reference
  E_d numérico=-5.44223154e-02  (exacto -5.48288854e-02)
  E_u (continuo más cercano)=-5.34222603e-02  ΔE_u=9.1964e-03
  Gamma numérico  = 1.647250e-04
  Gamma analítico = 1.645145e-04  (cuadratura, E_d exacto)
  diferencia relativa = 0.1279%
PASSED
tests/systems/nonadiabatic_dynamics/test_decay_rates.py::test_decay_rate_scales_quadratically_with_coupling_strength
  Gamma/A0^2 para A0=(0.005,0.01,0.02): [1.6472495931530404, 1.6472495931530404, 1.6472495931530404]
PASSED
tests/systems/nonadiabatic_dynamics/test_decay_rates.py::test_nonadiabatic_decay_rate_rejects_bound_state_below_continuum_threshold PASSED
tests/systems/nonadiabatic_dynamics/test_decay_rates.py::test_coupling_operator_block_is_finite_and_correct_shape PASSED
```

**Γ∝A0² exacto** (9 cifras idénticas): identidad algebraica de la Ec. 15
con B_du=0 y A_du constante (el elemento de matriz es lineal en A0), sirve
de test de regresión limpio, no sujeto al efecto de la §5.

## 5. Por qué Γ NO muestra convergencia O(h²) limpia (y por qué eso está bien)

A diferencia de todo lo demás en el proyecto, Γ tiene una segunda fuente de
"ruido" además de la discretización espacial: el estado del continuo más
cercano a E_d rara vez cae EXACTAMENTE en E_d (la caja aproxima el continuo
con niveles espaciados ΔE_u≈9.2×10⁻³ Eh en esta configuración), así que Γ se
evalúa en un E_u ligeramente distinto de E_d, y Γ(E) varía con la energía.
Medido directamente comparando Γ_numérico contra Γ_analítico evaluado EN EL
MISMO E_d que usó cada malla (para aislar el efecto):

```
h=0.200  E_d_numerico=-5.4023e-02  err=6.743e-07  err/h²=1.686e-05
h=0.100  E_d_numerico=-5.4422e-02  err=1.648e-06  err/h²=1.648e-04
h=0.050  E_d_numerico=-5.4625e-02  err=2.104e-06  err/h²=8.416e-04
```

El error CRECE con h decreciente en vez de caer — lo opuesto de O(h²). Esto
se investigó (no se ignoró): la causa es el "ruido de emparejamiento en
energía" descrito arriba, que no decrece monótonamente con h a L fija (más
puntos añaden niveles del continuo más finamente espaciados, pero el nivel
MÁS CERCANO a un E_d dado no se acerca de forma sistemática). No es un bug
del método; es una limitación conocida de este esquema simple de "vecino
más cercano" para elegir el estado del continuo, y se documenta en vez de
ocultarse (regla del proyecto: no forzar una limitación conocida a
desaparecer).

Por eso el test de Γ verifica **acuerdo cuantitativo con margen** (2%,
frente al 0.13% medido — deja margen sin ser tan laxo que un error de
signo o factor real pase desapercibido) y la **tendencia exacta** Γ∝A0²
(que sí es robusta a este efecto, por ser una identidad algebraica), en vez
de exigir un orden de convergencia formal para Γ — coherente con lo que
pidió el usuario ("confirma orden de magnitud y tendencia", no una prueba
de convergencia). Si en una fase posterior se necesita más precisión en Γ,
la mejora natural es interpolar entre los dos niveles del continuo más
cercanos a E_d en vez de tomar el más próximo sin más — no implementado
aquí, queda anotado como mejora futura, no como deuda oculta.

## 6. Un ajuste de tolerancia (documentado)

Al ejecutar la suite por primera vez falló
`test_energy_normalization_matches_analytic_free_particle`: el umbral
`errors[2] < 1e-5` se había puesto ANTES de ver el resultado final con
L=600 (la exploración previa había usado L=2000, con un error absoluto algo
distinto). El resultado real fue 8.40×10⁻⁵, y la razón de convergencia ya
confirmaba O(h²) limpio (4.00× y 4.00× exactos) — el umbral simplemente
estaba mal calibrado, no había ningún problema físico. Se corrigió a
`< 2e-4`, con margen sobre el valor medido.

## Conclusiones y aplicación al proyecto

- La reconstrucción del método de normalización en energía es internamente
  consistente y correcta DENTRO del marco que se derivó: converge al orden
  esperado (2) contra una solución analítica cerrada e independiente.
- La tasa de decaimiento de la Ec. 15 (ésta sí verificada contra el texto
  completo del paper de 2024) reproduce, dentro de 0.13%, el resultado de
  una integral analítica completamente independiente del método de
  diferencias finitas.
- **Pendiente explícito**: contrastar este método reconstruido contra
  González-Férez, Weidemüller & Schmelcher (2007) si se consigue acceso al
  texto en el futuro. Que pase sus propios tests analíticos NO es lo mismo
  que ser el mismo método que usa el paper de 2024 — puede que coincidan
  (es la técnica estándar de teoría de colisiones para este problema, así
  que es plausible que sí), pero no se ha verificado.
- El pipeline completo Fases 1-4 (acoplamiento de derivada → canales
  acoplados → estados ligados reales → tasas de decaimiento) está ahora
  completo y verificado con casos de referencia independientes en cada
  paso. Sigue sin aplicarse a ningún dato real de Rb*-RbCs — eso es la
  Fase 6.

## Archivos tocados

- `src/trimero/systems/nonadiabatic_dynamics/energy_normalization.py` (nuevo)
- `src/trimero/systems/nonadiabatic_dynamics/decay_rates.py` (nuevo)
- `tests/systems/nonadiabatic_dynamics/test_energy_normalization.py` (nuevo, 4 tests)
- `tests/systems/nonadiabatic_dynamics/test_decay_rates.py` (nuevo, 5 tests)
- `docs/PLAN_nonadiabatic_dynamics.md` (ESTADO ACTUAL: Fase 4 HECHA, con
  advertencia de reconstrucción propia)
- `docs/analysis_fase4_normalizacion_y_decaimiento.md` (este documento)

## Referencias

- González-Férez, R., Weidemüller, M. & Schmelcher, P., *Phys. Rev. A*
  **76**, 023402 (2007) — **método original, NO conseguido**. Ver
  `docs/analysis_fase4_busqueda_referencia.md`. **Pendiente de contraste
  cruzado si se consigue acceso.**
- Mellado-Alcedo, Guttridge, Cornish, Sadeghpour & González-Férez,
  *Phys. Rev. A* **110**, 013314 (2024), arXiv:2401.09618 — Ec. 15
  (completa, verificada) y la cita a la referencia de 2007 para el método
  de normalización (sin reproducirlo).
- Sakurai, J.J. & Napolitano, J., *Modern Quantum Mechanics* — normalización
  estándar de estados del continuo, base de la reconstrucción de este
  documento.
- Landau, L.D. & Lifshitz, E.M., *Quantum Mechanics* §21 — normalización de
  la onda esférica libre a δ(E-E'), la fórmula u_E(R) usada como referencia
  analítica cerrada.
- `docs/analysis_fase3_estabilizacion.md` — el pozo cuadrado radial
  reutilizado aquí (rama impar de la ecuación trascendente).
