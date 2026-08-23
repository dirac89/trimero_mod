# Plan maestro — figuras de publicación (Rb*-KRb y el híbrido Rb*-Rb-RbCs)

**Fecha de creación**: 2026-08-23
**Autor**: Javier Aguilera (via Claude Code)
**Relevancia**: plan persistente de varias rondas; puede interrumpirse y
retomarse por otra instancia sin memoria de esta conversación.

> Este documento es el plan de trabajo tal cual se acordó en el chat. Su
> sección **ESTADO ACTUAL** (al final) es lo primero que hay que leer si te
> reactivan sin contexto.

> ⚠️ **Corrección 2026-08-23, tras cerrar la Fase B**: el texto original de
> este plan (abajo, dictado tal cual por el usuario) llama "Rb*-RbCs" al
> sistema `rb_krb_polar` y le atribuye las constantes moleculares del RbCs
> (d=1.225 D, B=490.17 MHz). Es un error: `rb_krb_polar` corre con
> constantes de **KRb** (`B_KRB_GHZ=1.114`, `D_KRB_DEBYE=0.566`,
> Ni et al. 2008/2009, ver `src/trimero/systems/rb_krb_polar/
> charge_dipole.py`), consistente con las referencias que lo validan
> (Aguilera-Fernández 2015, González-Férez 2015 — ambos papers de KRb). Las
> constantes de RbCs (d=1.225 D, B=490.17 MHz) son del sistema híbrido
> `hybrid_neutral_polar` (`D_RBCS_DEBYE`/`B_RBCS_MHZ` en `hybrid_system.py`),
> un módulo distinto. El sistema es **Rb\*-KRb**, no Rb\*-RbCs; se corrigen
> las menciones sueltas de "RbCs" referidas a `rb_krb_polar` en el resto de
> este documento y en `docs/analysis_faseA_*.md` / `docs/analysis_faseB_*.md`
> y sus figuras. No cambia ningún número: es una corrección de etiqueta, no
> de física — los cálculos de las Fases A y B ya usaban las constantes de
> KRb correctamente, sólo el nombre en los títulos estaba mal.

═══════════════════════════════════════════════════════════════════════════
CONTEXTO DEL PROYECTO
═══════════════════════════════════════════════════════════════════════════

Repositorio: dirac89/trimero_mod, rama migrate-python, paquete src/trimero/.
Sistemas existentes y ya validados:
  - systems/rb_krb_polar/: H_A+H_mol (carga-dipolo), validado contra
    Aguilera-Fernández 2015 (Fig. 1, n=25) y González-Férez 2015 (Tabla I).
    Base: manifold(n,l≥3)+(n+1)d+(n+2)p+(n+3)s. [Corregido 2026-08-23:
    constantes de KRb, no RbCs — d=0.566 D, B=1.114 GHz, Ni et al.
    2008/2009. Ver nota de corrección al inicio de este documento.]
  - systems/rb_neutral_perturber/: pseudopotencial de Fermi s+p, validado
    contra Aguilera-Fernández 2016 (geometría lineal simétrica, campo DC ya
    implementado ahí, Vfield en rb_atom.py).
  - systems/hybrid_neutral_polar/: Rb*(n=35)-Rb(neutro,θ=π)-RbCs(θ=0).
    Resultado ya documentado: el perturbador neutro NO orienta RbCs como un
    campo externo (efecto débil, ≲5%, no monótono en R1).
  - systems/nonadiabatic_dynamics/: pipeline completo (acoplamiento de
    derivada, canales acoplados, estabilización, tasas de decaimiento,
    Franck-Condon) verificado en 6 fases (+ Fase 6b, ventana ampliada,
    corrección importante de un bug de referencia de energía absoluta —
    ver docs/analysis_fase6b_ventana_ampliada.md), aplicado como
    demostración a n=25.
  - systems/visualization/geometry_diagram.py: esquemas de geometría
    molecular reutilizables, ya construido.

Disciplina a mantener (igual que toda la sesión): cada figura/resultado
nuevo se verifica antes de aceptarse (comparación con casos ya validados,
límites físicos conocidos, órdenes de magnitud razonables); toda limitación
se documenta explícitamente, nunca se oculta ni se fuerza un resultado;
cada ronda se documenta en docs/analysis_*.md con resultados numéricos
REALES, no solo narrativa; si un resultado sospechosamente "perfecto"
aparece (ver el precedente de la Fase 6b: 300/300 estados "estables"),
se investiga antes de aceptarlo, no se reporta sin más.

═══════════════════════════════════════════════════════════════════════════
OBJETIVO
═══════════════════════════════════════════════════════════════════════════

Producir un conjunto de figuras, en el estilo de Aguilera-Fernández et al.
2015/2016 (IOP JPCS 635, 012023 y JPB 49, 124002), aplicadas a los sistemas
Rb*-KRb (rb_krb_polar) y al híbrido Rb*-Rb-RbCs, con el fin de establecer
una conclusión cualitativa/comparativa: **¿qué probabilidad relativa hay de
que estas configuraciones formen realmente una molécula de Rydberg ligada?**
(vía factores de Franck-Condon, ya implementados en
systems/nonadiabatic_dynamics/franck_condon.py — NO se busca una cifra
absoluta cuantitativa con n≈52 en esta ronda, sino una comparación entre
configuraciones: distintos n, presencia/ausencia de campo eléctrico,
sistema polar puro vs. híbrido).

Alcance decidido con el usuario:
  - Sistemas: rb_krb_polar (Rb*-KRb) Y el híbrido (Rb*-Rb-RbCs). NO el
    perturbador neutro puro en esta ronda.
  - Campo ELÉCTRICO: sí, se añade a rb_krb_polar (no existe ahí todavía).
    Campo MAGNÉTICO: diferido, no se implementa en este plan.
  - Conclusión final: cualitativa/comparativa vía Franck-Condon relativos,
    NO la corrida cuantitativa cara de n≈52.

═══════════════════════════════════════════════════════════════════════════
PLAN DE FASES
═══════════════════════════════════════════════════════════════════════════

--- FASE 0: Este documento ---
Escribir este plan en docs/PLAN_figuras_publicacion.md. Actualizar su
sección ESTADO ACTUAL después de cada ronda.

--- FASE A: Curvas BOP comparadas entre varios n (rb_krb_polar) ---
Ya existe la maquinaria (scripts/compute_bop_curve.py, BOPSystem) y datos
para n=25 (Fig. 1 ya reproducida). Objetivo: generalizar a un conjunto de n
comparables (candidatos: n=24,25,26,27 — confirma con domain_bounds() que
todos son computacionalmente razonables antes de comprometerte a la lista
completa; si alguno es prohibitivo, dilo y ajusta el conjunto, no fuerces).
Un plot con las curvas de M_J=0 superpuestas o en paneles, mismo estilo que
la Fig. 1 ya validada (energía vs R, cero en el manifold libre). Verifica
ANTES de nada las tendencias ya conocidas de la sesión (pozos se desplazan
a mayor R con n creciente, profundidad decrece con n) como test de sanidad
del propio barrido nuevo, no solo como resultado a reportar.

--- FASE B: Orientación ⟨cos θ_d⟩ vs R para varios n ---
Ya se calculó ⟨cosθ_d⟩ para n=24 (validación contra González-Férez 2015,
cruce en R≈390 a₀) en rondas anteriores de esta sesión — recupera ese
código/resultado en vez de rehacerlo desde cero. Generalízalo al mismo
conjunto de n de la Fase A. Un plot comparativo (⟨cosθ_d⟩ vs R, una curva
por n). Ten en cuenta la limitación ya documentada entonces: el seguimiento
del estado se vuelve poco fiable en regiones de alta densidad de estados —
no fuerces una curva limpia donde no la hay, documenta dónde se vuelve
ambiguo el seguimiento.

--- FASE C: Campo eléctrico DC en rb_krb_polar ---
Añadir H_ext = e·r·F_ext - d·F_ext (Ec. 5 de González-Férez et al. 2015,
YA VERIFICADA en el texto del paper en una ronda muy temprana de esta
sesión — no la rederives, cítala) como término OPCIONAL nuevo (no rompas
nada existente; por defecto F_ext=0 debe reproducir exactamente los
resultados ya validados, bit a bit, como test de no-regresión obligatorio
antes de seguir). Dos piezas:
  1. e·r·F_ext: acoplamiento Stark del electrón Rydberg — verifica primero
     si esto ya está parcialmente cubierto por algo existente en
     rb_atom.py (Vfield, que sí existe para el sistema neutro) o si hay
     que generalizarlo/reimplementarlo para la base de rb_krb_polar.
  2. -d·F_ext: acoplamiento Stark del dipolo permanente [corregido: de KRb,
     no RbCs — ver nota de corrección al inicio] — más
     simple, un término directo sobre el rotor.
Tests de límite ANTES de física real: F_ext=0 reproduce exactamente lo ya
validado (test de no-regresión); a F_ext grande, verifica que domina sobre
el resto del Hamiltoniano de forma razonable (comparación de escala, no
solo "no explota").
Con el término verificado: curvas BOP a varios F_ext (candidatos de partida
100, 300, 500 V/m, mismo orden que usó Aguilera-Fernández 2016 para el
sistema neutro — convierte a u.a. con el factor CODATA ya usado en esa
ronda, 5.14220674763e11 V/m por u.a.), mismo estilo de figura que la de
la ronda del trímero lineal con campo DC.

--- FASE D: Curvas BOP y ⟨cosθ⟩ del híbrido (Rb*-Rb-RbCs) ---
Ya existe la maquinaria (systems/hybrid_neutral_polar/) y el resultado ya
documentado (efecto débil del perturbador neutro). Objetivo aquí es
puramente GRÁFICO: producir las figuras de curvas BOP (ya semi-existen, ver
plots/hybrid_neutral_polar/) en el estilo consistente con las Fases A-C
(mismos ejes, misma paleta, para que sean comparables visualmente), y una
figura nueva de ⟨cosθ⟩ del híbrido vs R2 (a distintos R1 fijos) si la
maquinaria de la Fase B se puede reutilizar aquí sin gran esfuerzo — si
requiere trabajo sustancial nuevo, dilo y decide si vale la pena en esta
ronda o se pospone.

--- FASE E: Conclusión comparativa — probabilidad relativa de formación ---
Usa franck_condon.py (ya construido y verificado) para comparar, de forma
CUALITATIVA/COMPARATIVA (no cuantitativa absoluta con n≈52), los factores
de Franck-Condon relativos entre:
  - distintos n (Fase A/B)
  - con y sin campo eléctrico (Fase C)
  - sistema polar puro vs. híbrido (Fase D)
Con esto, sintetiza una tabla o figura resumen que responda: ¿qué
configuración(es) tienen mayor probabilidad relativa de formar una
molécula de Rydberg ligada, dentro de lo que este modelo puede decir sin
la corrida cuantitativa cara? Sé explícito sobre qué es una comparación
relativa válida y qué NO se puede concluir sin más trabajo (p.ej. no se
puede dar una probabilidad absoluta sin n≈52 real).

--- FASE F: Consolidación ---
Reunir las figuras finales (formato consistente, listas para manuscrito) y
actualizar la página de Notion del manuscrito con una sección de
Resultados gráficos. Esto se coordina en el chat, no es tarea de código.

═══════════════════════════════════════════════════════════════════════════
CÓMO RETOMAR SI SE INTERRUMPE LA SESIÓN
═══════════════════════════════════════════════════════════════════════════

1. Lee docs/PLAN_figuras_publicacion.md completo.
2. Lee su sección ESTADO ACTUAL.
3. Lee los docs/analysis_*.md más recientes relacionados.
4. Corre pytest completo para confirmar que el repo coincide con lo
   documentado.
5. Continúa desde la fase indicada.

═══════════════════════════════════════════════════════════════════════════
ESTADO ACTUAL
═══════════════════════════════════════════════════════════════════════════

Fase 0: **COMPLETA** (este documento).

Fase A: **COMPLETA** (2026-08-23). Curvas BOP M_J=0 calculadas y verificadas
para n=24,25,26,27 (n=25 reutilizado del cálculo previo, sin recalcular).
Script nuevo `scripts/compare_bop_curves_n.py` (no toca `compute_bop_curve.py`,
que ya aceptaba `--n-manifold` arbitrario). Figura:
`plots/rb_krb_polar/figures/fig_compare_n_MJ0.png`. Análisis completo con números
reales en `docs/analysis_faseA_curvas_bop_varios_n.md`. Tendencia de escala
confirmada (con un matiz importante documentado ahí: el "pozo más profundo"
reportado por `compute_bop_curve.py` está anclado en el borde R=400 a₀ de la
ventana fija para los cuatro n — la tendencia física real de desplazamiento
a mayor R se ve en el panel reducido R/2n² de la figura, donde las cuatro
curvas colapsan). `pytest -m "not slow"` en verde (133 passed) antes de
empezar; `pytest` completo también en verde al cierre (**149 passed en
644.69s**, ningún golden se movió) — Fase A cerrada sin deuda pendiente.

Fase B: **cálculo y análisis completos, cierre formal pendiente** (2026-08-23).
⟨cosθ_d⟩ calculado para n=24,25,26,27 (M_J=0, fermi=False, R∈[100,1800] a₀).
Código recuperado de la ronda anterior y reproducido bit a bit para n=24
(confirmado); se detectó y corrigió una dependencia desactualizada
(`fermi=True` por defecto — la premisa equivocada que ya se había corregido
para las curvas BOP pero no para orientación). Acotación [-1,1] verificada
explícitamente en los 4 n. Diagnóstico cuantitativo de ambigüedad por n
hecho: escala con 2n² (razón ≈0.53-0.56), igual que el codo de la Fase A.
Cruce de 0.78 (González-Férez 2015) en R≈365-395 a₀ para los 4 n. Figura
`plots/rb_krb_polar/figures/fig_orientation_compare_n_MJ0.png` generada. Análisis
completo en `docs/analysis_faseB_orientacion_varios_n.md`.

**Pendiente antes de dar la Fase B por cerrada** (checkpoint por límite de
uso, no por bloqueo): correr `poetry run pytest` completo (como se hizo al
cierre de la Fase A) y confirmar 149/149 en verde — los scripts nuevos
(`compute_orientation_curve.py`, `compare_orientation_n.py`) no tocan
`src/`, así que no se espera regresión, pero no se ha confirmado todavía.
Si retomas esto: corre el pytest completo, pega el resultado real en el
análisis y en el chat, y entonces sí marca Fase B como COMPLETA.

Fase C-F: NO EMPEZADAS.
