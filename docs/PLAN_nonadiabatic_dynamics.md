# PLAN — Dinámica no adiabática Rb*-RbCs (réplica de Mellado-Alcedo et al. 2024)

CONTEXTO MAESTRO — LÉELO ENTERO ANTES DE HACER NADA. Este documento es un plan
de trabajo persistente para una investigación de varias rondas. Es posible que
la sesión se interrumpa (fin de tokens) y se retome más tarde, posiblemente
por otra instancia sin memoria de esta conversación — por eso este documento
tiene que ser autosuficiente. Se mantiene actualizado al final de cada ronda
de trabajo con una sección "ESTADO ACTUAL" que dice exactamente qué fase está
hecha, cuál está en curso, y qué falta — así, si se reactiva sin contexto de
esta conversación, lo primero que hay que hacer es LEER ese fichero antes de
preguntar nada.

═══════════════════════════════════════════════════════════════════════════
CONTEXTO DEL PROYECTO (para que se entienda dónde encaja esto)
═══════════════════════════════════════════════════════════════════════════

Repositorio: dirac89/trimero_mod, rama migrate-python. Paquete src/trimero/,
organizado en systems/rb_krb_polar/ (Rydberg + molécula polar, carga-dipolo,
validado contra Aguilera-Fernández et al. 2015 y González-Férez et al. 2015),
systems/rb_neutral_perturber/ (Rydberg + átomos neutros, pseudopotencial de
Fermi s+p, validado contra Aguilera-Fernández et al. 2016), y
systems/hybrid_neutral_polar/ (combinación de ambos, primer sistema tetra-
atómico del proyecto — resultado: el perturbador neutro NO orienta la
molécula polar como lo haría un campo externo, efecto débil y no monótono).

Todo el trabajo de esta sesión sigue una disciplina estricta que hay que
mantener:
  - Cada ronda de trabajo se documenta en un fichero docs/analysis_*.md nuevo
    (contexto, método, resultados numéricos REALES —nunca solo un resumen
    narrativo—, conclusiones honestas, archivos tocados, referencias).
  - Los tests de límite y de sanidad se escriben y verifican ANTES de aceptar
    ningún resultado físico nuevo. Si un test falla, SE PARA Y SE REPORTA —
    nunca se ajusta el test ni el criterio para que pase.
  - Cuando un cálculo mejora en precisión (y se puede verificar
    independientemente, p.ej. contra aritmética exacta), se actualiza la
    referencia/golden file y se documenta por qué cambió — no se relaja la
    tolerancia para acomodar el número antiguo.
  - Las curvas de potencial se identifican por CARÁCTER electrónico (peso en
    el manifold relevante), nunca por "el autovalor más bajo sin más" — ya
    hemos tenido varias rondas donde ese atajo llevó a interpretar mal la
    física (ver docs/analysis_hibrido_caracter_E0.md como caso de estudio).
  - La convergencia numérica (N_max del rotor, paso de malla, tamaño de caja)
    se verifica explícitamente antes de citar cualquier cifra como definitiva.
  - No se fuerza ninguna limitación conocida a desaparecer: se documenta y se
    decide explícitamente si se acepta o se corrige (ver el manejo de la
    resonancia de forma p / pseudopotencial de rango cero como precedente).

═══════════════════════════════════════════════════════════════════════════
OBJETIVO DE ESTA INVESTIGACIÓN
═══════════════════════════════════════════════════════════════════════════

Replicar, para el sistema Rb*-RbCs (ya completamente validado en
systems/rb_krb_polar/), la maquinaria de dinámica no adiabática de:

  Mellado-Alcedo, Guttridge, Cornish, Sadeghpour & González-Férez,
  "Ultralong-range Cs-RbCs Rydberg molecules: non-adiabaticity of dipole
  moments", Phys. Rev. A 110, 013314 (2024), arXiv:2401.09618

Ese paper hace el análisis para Cs*-RbCs. Nadie ha hecho el análisis
equivalente para Rb*-RbCs — que es la combinación de especies que SÍ se ha
demostrado experimentalmente en el bloqueo de Rydberg de Ruttley, Guttridge
et al., PRL 131, 013401 (2023). Ese es el hueco que llenamos.

Hamiltoniano de referencia del paper 2024 (Ec. 1-4), para Rb*-RbCs:

    H = H_A + H_mol + H_int
    H_int = -d·F_ryd(R,r) + 2π a_S(k) δ³(r-R)

El segundo término (contacto de onda s) usa una aproximación de umbral MUY
simple, sin tablas: a_S = 1/√(2·E_A), con E_A la afinidad electrónica de
RbCs (E_A = 0.478 ± 0.020 eV, Eaton et al., Chem. Phys. Lett. 193, 141,
1992) — un valor CONSTANTE (no depende de k), y es una propiedad de RbCs,
no del átomo de Rydberg, así que sirve igual para Rb* que para Cs*. NO
requiere ninguna tabla de dispersión nueva.

═══════════════════════════════════════════════════════════════════════════
PLAN DE FASES
═══════════════════════════════════════════════════════════════════════════

--- FASE 0: Este documento ---
Escribir este plan en docs/PLAN_nonadiabatic_dynamics.md. Actualizar su
sección "ESTADO ACTUAL" (al final de este documento) después de cada ronda.

--- FASE 1: Acoplamiento no adiabático de derivada ⟨Ψᵢ|d/dR|Ψⱼ⟩ ---
Módulo nuevo, GENÉRICO (no depende de qué sistema genera las curvas):
src/trimero/systems/nonadiabatic_dynamics/coupling.py

1. Función que, dado un conjunto de R y, para cada R, la matriz de
   autovectores de un Hamiltoniano ya diagonalizado (mismo orden de base en
   todos los R), calcule ⟨Ψᵢ(R)|d/dR|Ψⱼ(R)⟩ por diferencias finitas
   centradas. Debe aceptar tanto una malla precalculada como una función
   solve(R)->(w,V) para generarla internamente.
2. CRÍTICO: fijar el signo de cada autovector en R_k maximizando el
   solapamiento con R_{k-1} antes de derivar (numpy.linalg.eigh no garantiza
   continuidad de fase entre puntos vecinos). Documentar esto como el paso
   más delicado.
3. Tests, en este orden, ANTES de aceptar nada:
   a. Antisimetría exacta ⟨i|d/dR|j⟩ = -⟨j|d/dR|i⟩ (identidad matemática).
   b. Diagonal exactamente nula ⟨i|d/dR|i⟩ = 0.
   c. Acoplamiento pequeño y suave lejos de cruces evitados (usar curvas ya
      calculadas de rb_krb_polar, tramo sin cruces cercanos).
   d. Pico de acoplamiento en un cruce evitado real ya documentado en el
      proyecto (por ejemplo, alguno visto al seguir la curva de carácter del
      híbrido o al cruzar un umbral de vecino en el sistema polar).
   e. Convergencia con el paso de diferencias finitas (2-3 pasos distintos).
4. NO implementar todavía la ecuación radial acoplada ni nada posterior.

--- FASE 2: Ecuación radial de canales acoplados (Ec. 8 del paper) ---
src/trimero/systems/nonadiabatic_dynamics/coupled_channels.py

1. Discretizar y resolver el sistema de N canales acoplados (empezar por 2):
   T + V_i + A_ii en la diagonal, A_ij fuera de diagonal (T = energía
   cinética radial reducida, A_ij del módulo de la Fase 1 más el término de
   energía cinética cruzado ⟨Ψᵢ|T|Ψⱼ⟩). Elegir método de discretización
   (diferencias finitas o DVR) y JUSTIFICAR la elección.
2. Tests de límite ANTES de física real:
   a. Acoplamiento puesto a cero -> debe recuperar exactamente los estados
      ligados de un solo canal, resueltos independientemente con un solver
      de referencia simple (Numerov o diferencias finitas estándar de un
      pozo 1D) — no comparar contra sí mismo, usar una referencia distinta.
   b. Verificar con un potencial de juguete (doble pozo simple con solución
      semi-analítica conocida o muy bien caracterizada en la literatura de
      métodos numéricos) antes de aplicarlo a ninguna curva real de RbCs.
3. Aplicar a un par de curvas REALES de Rb*-RbCs con un cruce evitado ya
   identificado (de la Fase 1d). Reportar los niveles vibracionales
   resultantes.

--- FASE 3: Método de estabilización (Hazi-Taylor 1970) ---
src/trimero/systems/nonadiabatic_dynamics/stabilization.py

1. Repetir la resolución de la Fase 2 variando el límite inferior de la
   malla radial R_min (o el tamaño de caja), y trazar el diagrama de
   estabilización.
2. Criterio automático (no solo visual) para identificar qué niveles son
   estados ligados reales (energía estable frente a R_min) frente a estados
   de caja (energía que se mueve sistemáticamente).
3. Test: verificar en el caso de acoplamiento cero que TODOS los estados
   ligados conocidos del canal inferior se identifican correctamente como
   estables, y que aparecen estados de caja claramente inestables por
   encima del umbral de disociación.

--- FASE 4: Tasas de decaimiento no adiabático ---
src/trimero/systems/nonadiabatic_dynamics/decay_rates.py

1. Implementar la normalización en energía de la función de onda del
   continuo del canal inferior (referencia: González-Férez, Weidemüller &
   Schmelcher, PRA 76, 023402, 2007 — CONSEGUIR y leer ese paper primero,
   no improvisar el método).
2. Calcular Γᵢ = (2π/ℏ)|⟨χᵢ^d|A_du|χⱼ^u⟩|² para los estados ligados de la
   Fase 3.
3. Test de sanidad: verificar que Γ crece cuando el acoplamiento de la
   Fase 1 es mayor (cerca de un cruce evitado más estrecho/fuerte) y es
   pequeño lejos de cruces — comparación cualitativa mínima antes de
   cualquier cifra cuantitativa.

--- FASE 5: Factores de Franck-Condon ---
src/trimero/systems/nonadiabatic_dynamics/franck_condon.py

Réplica de la Ec. A9 del paper 2024, adaptada al estado inicial relevante
para Rb*-RbCs (definir explícitamente qué estado inicial se usa, análogo al
modelo de excitación de dos fotones del paper, antes de calcular nada).

--- FASE 6: Física real de Rb*-RbCs ---
1. Añadir el término de contacto de onda s a rb_krb_polar (constante
   a_S=1/√(2·0.478 eV), SIN tablas), como término OPCIONAL (no romper nada
   existente). Verificar que su efecto es pequeño (el paper reporta ~48.6
   MHz de diferencia máxima) — si sale muy distinto, parar y reportar.
2. Identificar cruces evitados reales en el sistema Rb*-RbCs ya calculado
   (candidatos: los ya vistos en el trabajo de validación de la Fig. 1 y en
   el híbrido) y aplicar el pipeline completo (Fases 1-5) a al menos dos
   casos, en analogía directa con las Figs. 2-7 del paper 2024.
3. Comparar cualitativamente contra los hallazgos del paper para Cs-RbCs
   (estructura de niveles ligados, dependencia de las tasas de decaimiento
   con el número cuántico vibracional, hibridación de ondas parciales) —
   NO se espera coincidencia cuantitativa (especies distintas), sí
   consistencia cualitativa.

--- FASE 7: Escritura ---
Actualizar (o crear, si se decide que merece manuscrito propio en vez de
extender el existente) la página de Notion correspondiente con Resultados y
Discusión reales, siguiendo el mismo formato que
"DRAFT Manuscript — Ultralong-Range Rb*-Rb-RbCs Hybrid Rydberg Molecule" ya
existente. Esto se coordina en el chat, no es tarea de código.

═══════════════════════════════════════════════════════════════════════════
CÓMO RETOMAR SI SE INTERRUMPE LA SESIÓN
═══════════════════════════════════════════════════════════════════════════

Si se reactiva sin el contexto de esta conversación:
1. Lee docs/PLAN_nonadiabatic_dynamics.md completo (este documento).
2. Lee su sección "ESTADO ACTUAL" para saber qué fase toca.
3. Lee los docs/analysis_*.md más recientes relacionados con
   nonadiabatic_dynamics para ver el último resultado numérico real.
4. Corre pytest completo para confirmar que el estado del repo coincide con
   lo que dice la documentación antes de seguir avanzando.
5. Continúa desde la fase indicada, actualizando este mismo plan al terminar
   cada ronda.

═══════════════════════════════════════════════════════════════════════════
ESTADO ACTUAL (actualizar esta sección al final de cada ronda)
═══════════════════════════════════════════════════════════════════════════

Fase 0: HECHA (este documento).

Fase 1: HECHA (docs/analysis_fase1_acoplamiento_derivada.md). Módulo
`src/trimero/systems/nonadiabatic_dynamics/coupling.py`
(`fix_eigenvector_signs`, `eigenbasis_along_R`, `derivative_coupling`),
11 tests en verde en `tests/systems/nonadiabatic_dynamics/test_coupling.py`
(incluye 1 marcado `slow`, ~12 s la suite completa). Resultado numérico real
verificado: pico de acoplamiento |A(54,55)|=0.5067 en el cruce evitado real
de `BOPSystem(n_manifold=25)` M_J=0 en R≈1305 a₀ (hueco de energía mínimo
≈4.98e-9 Eh), frente a 5.98e-3 lejos de él (razón 84.75×); convergencia de
orden 2 verificada explícitamente con el paso de diferencias finitas.
Cerrada por consistencia con los datos de esa ronda; el pico en el cruce
evitado real se acepta por consistencia con datos ya validados en rondas
anteriores, sin reconstruir el sistema completo esta vez.

Fase 2: HECHA (docs/analysis_fase2_canales_acoplados.md). Módulo
`src/trimero/systems/nonadiabatic_dynamics/coupled_channels.py`
(`second_derivative_coupling`, `radial_kinetic_matrix`,
`first_derivative_matrix`, `build_coupled_hamiltonian`,
`adiabatic_coupled_hamiltonian`, `solve_coupled_channels`), 6 tests en verde
en `tests/systems/nonadiabatic_dynamics/test_coupled_channels.py` (1 marcado
`slow`, ~66 s). Discretización: diferencias finitas (justificado en el
docstring del módulo antes de implementar), simetrizada explícitamente.
Verificación más fuerte: doble pozo de 2 canales resuelto por dos vías
independientes (representación diabática vs. adiabática) coincide a
1.564e-10 (principio físico exacto, no ajustado). Acoplamiento a cero
recupera exactamente (rtol=1e-10) dos solves de un solo canal
independientes, distintos del propio método. Aplicado también al cruce
evitado real de la Fase 1 (R≈1305 a₀, estados 54/55): desplazamiento máximo
medido 0.381 GHz frente a la solución desacoplada, en una ventana de R
limitada (50 a₀) cuyas limitaciones quedan documentadas explícitamente — NO
se afirman niveles vibracionales físicos reales, eso es tarea de la Fase 3.
Masa reducida μ≈113522 mₑ (Rb-RbCs) es un valor de demostración, pendiente
de verificar la fuente exacta en la Fase 6.

Fase 3: HECHA (docs/analysis_fase3_estabilizacion.md). Módulo
`src/trimero/systems/nonadiabatic_dynamics/stabilization.py`
(`stabilization_scan`, `track_trajectories`, `classify_stability`,
`stable_state_summary`), 4 tests en verde en
`tests/systems/nonadiabatic_dynamics/test_stabilization.py` (~1.2 s, nada
`slow`). Criterio de estado ligado automático (no visual): |dE/dL| <
umbral·(2|E|/L), la escala de deriva de un estado de caja en esa energía y
longitud. Verificado contra un pozo cuadrado finito (μ=1.0, V0=0.1 Eh,
a=8.0 a₀, 3 estados ligados exactos por ecuación trascendente): el
clasificador encuentra EXACTAMENTE 3 estados ligados (ni de más ni de
menos), energías dentro de 1.4e-3 del valor exacto, y clasifica
correctamente todas las trayectorias de energía positiva como no ligadas.
Nota importante: la primera versión de la referencia independiente (pozo
cuadrado) tenía un bug (2 raíces espurias de un cruce de polo de tan/cot,
5 en vez de 3 estados) que se detectó ANTES de aceptar nada, contrastándola
contra la fórmula analítica de conteo ⌊k₀·2a/π⌋+1 — no se aplicó el método
a datos reales de Rb*-RbCs en esta ronda (instrucción explícita: primero
validar contra un caso completamente conocido).

Fase 4: HECHA — con advertencia (docs/analysis_fase4_normalizacion_y_decaimiento.md).
Tras no conseguir González-Férez, Weidemüller & Schmelcher PRA 76, 023402
(2007) (docs/analysis_fase4_busqueda_referencia.md, 6 vías agotadas), el
usuario decidió explícitamente RECONSTRUIR el método de normalización en
energía desde primeros principios (teoría de colisiones estándar) en vez
de bloquear la fase. ⚠️ NO es una réplica del paper de 2007 — es una
reconstrucción propia, verificada sólo contra un caso analítico cerrado
(partícula libre), no contra el método real de esa referencia. Pendiente
explícito: contrastar contra González-Férez et al. 2007 si se consigue
acceso en el futuro.

Módulos: `src/trimero/systems/nonadiabatic_dynamics/energy_normalization.py`
(`local_level_spacing`, `energy_normalize_box_states`) y
`decay_rates.py` (`coupling_operator_block`, `coupling_matrix_element`,
`decay_rate`, `nonadiabatic_decay_rate`). 9 tests en verde
(`test_energy_normalization.py` + `test_decay_rates.py`, ~26 s, nada
`slow`). Resultados clave: normalización en energía converge O(h²) limpio
(factores 4.00×, 4.00×) contra la solución analítica cerrada de la
partícula libre radial u_E(R)=√(2μ/(πk))sin(kR); Γ del modelo de juguete
Breit-Wigner (pozo cuadrado radial + continuo libre desplazado, A_du
constante) coincide con la referencia analítica por cuadratura dentro de
0.128%, y Γ∝A0² se cumple exactamente (9 cifras). Nota importante
documentada: Γ NO converge O(h²) limpio como el resto del proyecto —hay un
segundo efecto (el continuo más cercano a E_d rara vez cae exactamente en
E_d, "ruido de emparejamiento en energía")—, investigado y explicado en el
documento, no ignorado. Detectado y corregido en esta ronda: la malla debe
empezar en R=h (no en R=0) para que la pared de Dirichlet caiga en el punto
físico correcto — con R=linspace(0,L,n) el error caía sólo O(h), no O(h²).

Pipeline Fases 1-4 (acoplamiento de derivada → canales acoplados → estados
ligados reales → tasas de decaimiento) completo y verificado con
referencias independientes en cada paso. Todavía NO aplicado a ningún dato
real de Rb*-RbCs — eso es la Fase 6.

Fase 5: HECHA (docs/analysis_fase5_franck_condon.md). Módulo
`src/trimero/systems/nonadiabatic_dynamics/franck_condon.py`
(`gaussian_wavepacket`, `franck_condon_factor`), 6 tests en verde en
`tests/systems/nonadiabatic_dynamics/test_franck_condon.py` (~38 s, nada
`slow`). Corrección de numeración: la fórmula del plan ("Ec. A9") es en
realidad la Ec. (17) del paper 2024, verificada verbatim
(arxiv.org/html/2401.09618). Estado inicial para Rb*-RbCs definido
explícitamente (NO copiado del escenario de pinzas fusionadas del paper):
paquete de ondas gaussiano centrado en la separación FIJADA por pinzas
específicas por especie, INDEPENDIENTES (no fusionadas) — el protocolo REAL
ya demostrado en Ruttley, Guttridge et al., PRL 131, 013401 (2023),
arXiv:2303.06126 (R_am=310(40) nm, dispersión de alineación 50 nm, ambos
verificados directamente del texto). Verificado contra la fórmula cerrada
del solapamiento de dos gaussianas (estado fundamental de oscilador
armónico de la Fase 2 vs. el paquete gaussiano de esta fase): convergencia
O(h²) limpia (4.00×, 4.00×). R0 y σ quedan como parámetros de entrada, NO
fijados a valores reales de Rb*-RbCs en esta ronda — eso, y decidir a qué n
de Rydberg corresponde el análisis, es tarea de la Fase 6.

Pipeline Fases 1-5 completo y verificado con referencias independientes en
cada paso.

Fase 6: HECHA (docs/analysis_fase6_aplicacion_n25.md).
Primera aplicación del pipeline completo a datos reales:
`BOPSystem(n_manifold=25)`, M_J=0, estados 54/55 (el cruce de la Fase 1),
ventana ampliada a R∈[1180,1430] a₀ (250 a₀, 5× la de la Fase 2).
Decisión explícita: n=25 se usa como DEMOSTRACIÓN del pipeline, no como
predicción cuantitativa del experimento real de Ruttley et al. 2023 (que
usa Rb(52s)) — medido en esta ronda: una sola diagonalización en n=52 tarda
76.7 s (dim=2436) frente a ~1.1-2 s en n=25 (dim=1113), así que repetir
las Fases 1-2 para n≈52 sería obra de varias rondas, no de ésta. Masa
reducida resuelta con fuente citada: μ=113536.332 mₑ (NIST/AME2020 para
m(Rb-87) y m(Cs-133), CODATA 2022 para u→mₑ), sustituyendo el valor de
demostración de la Fase 2. Estado inicial (Fase 5): paquete gaussiano en
R0=1305 a₀ (el propio cruce, escalado a 2n²a₀ de n=25, NO se reutiliza
R_am=310nm que es específico de n=52), σ=945 a₀≈50 nm de Ruttley et al.
2023 (propiedad de la plataforma, no del n, con la advertencia explícita
de que no está validado para n=25).

Resultados: 15 estados clasificados ligados por estabilización (V-E_manifold
∈[-16.28,-8.91] GHz, mismo orden que el pozo más profundo de la Fig. 1
completa, −23.1 GHz); de ellos sólo 2 tienen un canal continuo accesible
dentro de la ventana para evaluar Γ (139.4 MHz y 23.5 MHz) — los otros 13
quedan fuera de alcance por el tamaño de la ventana, limitación explícita
y documentada, no oculta. El paquete gaussiano REAL (σ=945 a₀) no cabe en
la ventana de 248 a₀ (se reporta con aviso explícito; se añadió un σ_demo
ilustrativo que sí cabe). Comparación cualitativa con el paper de 2024
(Cs-RbCs): el paper reporta que los estados más cercanos al cruce decaen
más rápido — nuestros 2 puntos (factor ~6 entre ellos) son consistentes en
DIRECCIÓN pero insuficientes en número para confirmar la tendencia con
confianza. No se tocó ningún módulo de `src/trimero/` ni test en esta
ronda; suite del proyecto sin cambios (133 passed, 16 deselected).

Fase 6b: HECHA (docs/analysis_fase6b_ventana_ampliada.md). Cierra el punto
1 del pendiente de la Fase 6 (§7): ventana ampliada a R∈[500,1800] a₀
(1300 a₀, dominio completo de la Fig. 1), 1301 diagonalizaciones, 25.5 min
(estimado 24.7 min antes de ejecutar, dentro del umbral de 30-40 min
acordado; la alternativa de σ_real=945 a₀ con margen 5σ habría costado
107-179 min, no ejecutada, documentada como referencia). Test de
consistencia contra la Fase 6 original: pico de A(54,55) idéntico
(dif_rel=0.00%) en el tramo compartido.

**Bug real encontrado y corregido**: `classify_stability` (Fase 3) recibía
energías ABSOLUTAS de Hartree (dominadas por el offset del manifold
Rydberg, ~-5264 GHz) en vez de energías referenciadas al umbral de
disociación real del canal (medido: V_d→-9.2759 GHz, plano con desviación
4.4e-5 GHz) — sin corregir, el criterio `2|E|/L` quedaba con una escala de
comparación absurda y clasificaba CUALQUIER estado como "ligado" (300/300,
incluso con energía positiva, ni con threshold_ratio 33x más estricto
cambiaba). Es un bug de cómo se invocó el módulo en el script de Fase 6,
NO del módulo `stabilization.py` (sus tests, incluido el pozo cuadrado,
siguieron en verde durante todo el proceso). **Consecuencia: los 2
valores de Γ de la Fase 6 original (139.4 y 23.5 MHz) se retiran
explícitamente** — los estados de los que salieron ya no se confirman
como ligados con el criterio corregido (14 estados en el recorte-viejo,
no 15; rango [-16.275,-9.639] GHz).

Resultados corregidos: ventana ancha (n_keep=300) → 144 estados ligados
(de 300 probados) con perfil de confianza gradual y físicamente sensato
(37 al 100%, cola decreciente hasta el mínimo 50%) en vez del "300/300
perfecto" pre-corrección; de ellos, **79 tienen Γ evaluable** (frente a 2
antes), min=0.0002 MHz, max=129.29 MHz, mediana=2.36 MHz. Con 79 puntos,
correlación de Pearson r=0.60 entre profundidad relativa y Γ (menos
ligado → Γ mayor), cualitativamente consistente con el paper de 2024
("estados cerca del cruce decaen más rápido") — ahora una comparación de
tendencia real, no anecdótica de 2 puntos. σ_real=945 a₀ sigue sin caber
en esta ventana; se mantiene la estimación de 5670-9450 a₀ (107-179 min)
para una ronda futura. No se tocó ningún módulo de `src/trimero/`; suite
del proyecto sin cambios.

Fase 7: EN CURSO. Escritura — actualizar (o crear) la página de Notion
correspondiente con los Resultados y la Discusión reales de las Fases 1-6,
siguiendo el formato del "DRAFT Manuscript — Ultralong-Range Rb*-Rb-RbCs
Hybrid Rydberg Molecule" ya existente. Por diseño del propio plan maestro
("Esto se coordina en el chat, no es tarea de código"), esta fase no tiene
entregable de código ni de test — necesita que el usuario aporte o
autorice el acceso a la página de Notion (no hay integración de Notion
disponible en esta sesión) y decida, en el chat, si el contenido de las
Fases 1-6 se integra en ese manuscrito existente o merece uno propio.
