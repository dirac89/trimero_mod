# Fase 5 — Factores de Franck-Condon

**Fecha**: 2026-08-23
**Autor**: Javier Aguilera (con Claude Code)
**Relevancia**: quinto módulo de `docs/PLAN_nonadiabatic_dynamics.md`.
Estima la probabilidad de formar la molécula de Rydberg por excitación de
dos fotones, conectando el estado inicial (átomo + molécula en sus
trampas) con los estados finales ya calculados en las Fases 2-3 (canales
acoplados, estados ligados reales).

## Corrección de numeración respecto al plan original

El documento maestro llamaba a esto "Ec. A9 del paper 2024" (nombre
provisional, antes de poder leer el texto). Verificado contra
`arxiv.org/html/2401.09618`: la ecuación real es la **Ec. (17)**, en el
Apéndice A (que va de la Ec. 9 a la 17). Se corrige aquí para que quede
trazable.

## 1. Qué estado inicial se usa para Rb*-RbCs (decisión explícita, no copiada)

El paper de 2024 define ψ_scat para SU protocolo (Cs-RbCs), con dos
escenarios: (a) pinzas ópticas FUSIONADAS (Cs y RbCs comparten una trampa
combinada, Cs(42s), separación <200 nm) — el estado inicial natural es el
fundamental del movimiento relativo en la trampa fusionada (oscilador
armónico centrado en R≈0); (b) pinzas SEPARADAS (Cs(74s), ~500 nm) — el
paper menciona el escenario sin dar su ψ_scat.

**Para Rb*-RbCs, ninguno de los dos aplica directamente.** El sistema real
YA DEMOSTRADO es un tercer caso: Ruttley, Guttridge, Baldock,
González-Férez, Sadeghpour, Adams & Cornish, "Observation of Rydberg
blockade due to the charge-dipole interaction between an atom and a polar
molecule", *Phys. Rev. Lett.* **131**, 013401 (2023), arXiv:2303.06126 (en
abierto). Hechos verificados directamente del texto (vía `WebFetch`, no
supuestos):

- "**Species-specific tweezers are used to control the separation between
  the atom and molecule**" — Rb y RbCs están en pinzas ÓPTICAS
  INDEPENDIENTES, específicas por especie, **no fusionadas**.
- La separación átomo-molécula demostrada para el bloqueo Rydberg
  (transición Rb→52s) fue **R_am=310(40) nm**.
- Dispersión experimental de la alineación relativa entre pinzas, tiro a
  tiro: "**an estimated standard deviation of 50 nm in each coordinate**".
- Esquema de excitación: dos fotones, |g⟩→|6p₃/₂⟩→|r⟩, con |r⟩=|52s⟩.

Esto es físicamente distinto del escenario de pinzas fusionadas: no hay una
trampa combinada con un oscilador armónico de movimiento relativo centrado
en R=0. El átomo y la molécula ocupan trampas INDEPENDIENTES separadas por
una distancia R0 FIJADA experimentalmente (el parámetro de control del
experimento), con una incertidumbre de posicionamiento relativo σ medida
directamente. El modelo adaptado a ESTE sistema (no al del paper) es un
**paquete de ondas gaussiano centrado en la separación fijada R0**, de
anchura σ dada por la incertidumbre de posicionamiento medida — no un
oscilador armónico centrado en R=0:

```
u_scat(R) = (πσ²)^(-1/4) · exp[-(R-R0)²/(2σ²)]
```

`gaussian_wavepacket(R_grid, R0, sigma)` la implementa. **R0 y σ son
parámetros de entrada de este módulo, NO fijados a Rb*-RbCs en esta
ronda**: fijar R0 al cruce evitado real (Fase 1, R≈1305 a₀) y σ al valor de
Ruttley et al. (50 nm ≈ 945 a₀) es tarea de la Fase 6, que además exige
antes decidir a qué n de Rydberg corresponde el análisis (52s del paper de
2023 no es necesariamente el n relevante para el cruce ya estudiado en las
Fases 1-2, que usa n=25) — no se mezclan ambas decisiones aquí.

## 2. Convención: por qué "R dR" y no "dR"

Ver la derivación completa en el docstring de `franck_condon.py`. Resumen:
la Ec. (17) tiene un factor extra R en la medida de integración porque
ψ_scat(R) está escrita en la forma NO reducida (3D, ∫|ψ_scat|²R²dR=1),
mientras que χᵢᵏ(R) usa la forma REDUCIDA (u=R·R_l, ∫|χ|²dR=1) que emplea
el resto del proyecto (Fases 1-4). Con u_scat(R)≡R·ψ_scat(R), la Ec. (17)
por canal se reescribe F = Ω_ns∫χᵢᵏ(R)*·C(R)·u_scat(R)dR — la misma
convención reducida de siempre, así que `gaussian_wavepacket` devuelve
directamente u_scat (no ψ_scat), normalizado vía ∫u²dR=1.

## 3. Verificación: solapamiento de dos gaussianas (referencia analítica cerrada)

El estado fundamental de un oscilador armónico (ya verificado contra la
energía analítica en la Fase 2, E₀=½ω) es exactamente gaussiano, con
anchura σ_HO=1/√(μω). Solapado contra `gaussian_wavepacket` (el modelo de
estado inicial de esta fase) centrado en un punto distinto, el resultado
debe coincidir con la fórmula CERRADA del solapamiento de dos gaussianas:

```
⟨φₐ|φᵦ⟩ = √(2σₐσᵦ/(σₐ²+σᵦ²)) · exp[-(Rₐ-Rᵦ)²/(2(σₐ²+σᵦ²))]
```

resultado de libro de texto, independiente de `franck_condon_factor`
(camino de cuadratura discreta) salvo por compartir el mismo estado final
calculado por diferencias finitas.

Parámetros: μ=1.0, ω=0.05, oscilador centrado en R_center=100 a₀ (lejos de
la pared R=0, σ_HO=4.472 a₀), paquete inicial centrado en R0=110 a₀,
σ_scat=6.0 a₀ (elección de demostración de la máquina, no de Rb*-RbCs real
— eso es Fase 6).

## Resultados numéricos reales

Suite completa (`poetry run pytest tests/systems/nonadiabatic_dynamics/test_franck_condon.py -v -s`):

```
tests/systems/nonadiabatic_dynamics/test_franck_condon.py::test_gaussian_wavepacket_is_normalized PASSED
tests/systems/nonadiabatic_dynamics/test_franck_condon.py::test_gaussian_wavepacket_rejects_nonuniform_grid PASSED
tests/systems/nonadiabatic_dynamics/test_franck_condon.py::test_franck_condon_matches_two_gaussian_overlap_closed_form
  F_exacto = 4.00858575e-01
  errores vs h=(0.2,0.1,0.05): [3.0647522500149016e-05, 7.660697576350461e-06, 1.915100405702752e-06]
PASSED
tests/systems/nonadiabatic_dynamics/test_franck_condon.py::test_franck_condon_coupling_weight_scales_result PASSED
tests/systems/nonadiabatic_dynamics/test_franck_condon.py::test_franck_condon_rejects_nonuniform_grid PASSED
tests/systems/nonadiabatic_dynamics/test_franck_condon.py::test_franck_condon_vanishes_for_far_separated_states
  F_near=4.0085e-01  F_far=2.1595e-82
PASSED

============================== 6 passed in 37.73s ==============================
```

Suite completa del proyecto (`poetry run pytest -m "not slow"`): **133
passed, 16 deselected** — nada se rompió (127 al cierre de la Fase 4, +6 de
esta ronda). Los 6 tests pasaron en la primera ejecución: los parámetros
(σ_scat, R0, resolución de malla) se calibraron por exploración numérica
antes de fijar el test, siguiendo la misma disciplina de las rondas
anteriores.

### Convergencia (orden 2)

Razones sucesivas de error: 4.00× y 4.00× exactos — coincide con el orden 2
esperado de la diferencia centrada usada en la matriz cinética
(`build_coupled_hamiltonian`, Fase 2) que genera el estado final.

### Sanidad de magnitud (decaimiento con la distancia)

Con el paquete inicial desplazado ~45σ_HO del centro del oscilador
(R0=300 en vez de 110), F cae de 0.401 a 2.16×10⁻⁸² — la cola gaussiana
del solapamiento se comporta exactamente como se espera, confirmando que
no hay ninguna fuga numérica ni normalización rota que sostenga un
solapamiento artificialmente grande a distancia.

## Conclusiones y aplicación al proyecto

- El estado inicial para Rb*-RbCs se definió explícitamente a partir del
  protocolo REAL ya demostrado (Ruttley et al. 2023, no una copia del
  escenario de pinzas fusionadas del paper de 2024): paquete de ondas
  gaussiano en la separación fijada por las pinzas, con anchura de la
  incertidumbre de posicionamiento medida — una adaptación física genuina,
  no una traducción mecánica de la fórmula del paper.
- La maquinaria de solapamiento (`franck_condon_factor`) es correcta:
  coincide con una referencia analítica cerrada e independiente, converge
  al orden esperado, y se comporta razonablemente en los límites (peso
  lineal, decaimiento con la distancia).
- Con esto, las Fases 1-5 del plan quedan completas y verificadas cada una
  con un caso de referencia independiente. Sigue sin aplicarse nada a
  datos reales de Rb*-RbCs — eso es la Fase 6, que además tendrá que
  resolver una decisión pendiente de esta fase: a qué n de Rydberg (y por
  tanto qué R0) corresponde el análisis, antes de fijar R0 y σ a valores
  reales.

## Archivos tocados

- `src/trimero/systems/nonadiabatic_dynamics/franck_condon.py` (nuevo)
- `tests/systems/nonadiabatic_dynamics/test_franck_condon.py` (nuevo, 6 tests)
- `docs/PLAN_nonadiabatic_dynamics.md` (ESTADO ACTUAL: Fase 5 HECHA)
- `docs/analysis_fase5_franck_condon.md` (este documento)

## Referencias

- Mellado-Alcedo, Guttridge, Cornish, Sadeghpour & González-Férez,
  *Phys. Rev. A* **110**, 013314 (2024), arXiv:2401.09618 — Ec. (17)
  (verificada verbatim contra el texto).
- Ruttley, Guttridge, Baldock, González-Férez, Sadeghpour, Adams & Cornish,
  *Phys. Rev. Lett.* **131**, 013401 (2023), arXiv:2303.06126 — el
  protocolo REAL de Rb*-RbCs (pinzas específicas por especie, separación
  controlada R_am=310(40) nm, dispersión de alineación 50 nm) que motiva
  el modelo de estado inicial de esta fase.
