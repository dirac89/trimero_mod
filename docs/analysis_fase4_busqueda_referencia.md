# Fase 4 — Búsqueda de la referencia (González-Férez, Weidemüller & Schmelcher, PRA 76, 023402, 2007): bloqueada

**Fecha**: 2026-08-23
**Autor**: Javier Aguilera (con Claude Code)
**Relevancia**: la Fase 4 de `docs/PLAN_nonadiabatic_dynamics.md` (tasas de
decaimiento no adiabático) exige, por instrucción explícita del propio plan
maestro y reforzada por el usuario en esta ronda, conseguir y leer
González-Férez, Weidemüller & Schmelcher, *Phys. Rev. A* **76**, 023402
(2007) para el método de normalización en energía de la función de onda del
continuo ANTES de escribir ningún código — "no improvisar el método". Este
documento registra el intento (búsquedas web reales, no un resumen
narrativo) y por qué se para aquí sin implementar nada.

## Resumen

No se consiguió el texto completo del paper. Se intentó por seis vías
independientes, documentadas abajo con el resultado real de cada una. Según
la instrucción explícita del plan maestro ("Si no consigues el texto
completo, PARA y repórtalo explícitamente en vez de aproximar el método por
tu cuenta"), la Fase 4 queda BLOQUEADA: no se ha escrito ningún módulo
`decay_rates.py`, ni se ha implementado ninguna aproximación del método de
normalización en energía.

## Intentos realizados (orden cronológico, con resultado real)

### 1. Búsqueda web del título y confirmación de la referencia

`WebSearch`: "González-Férez Weidemüller Schmelcher Phys Rev A 76 023402
2007". Resultado: se confirma el título exacto — **"Photoassociation of
cold heteronuclear dimers in static electric fields"** — y que la
referencia (revista, volumen, página, año) es correcta. No se obtuvo texto.

### 2. DOI directo

`WebFetch` a `https://doi.org/10.1103/PhysRevA.76.023402` → redirige a
`https://link.aps.org/doi/10.1103/PhysRevA.76.023402` → **HTTP 403
Forbidden**. El sitio de APS bloquea el acceso automático (paywall, como es
de esperar sin suscripción institucional).

### 3. arXiv, búsqueda por título exacto

`WebFetch` a la API de arXiv
(`export.arxiv.org/api/query?search_query=ti:"Photoassociation of cold
heteronuclear dimers"`) → **`<opensearch:totalResults>0</opensearch:totalResults>`**.
Sin resultados.

### 4. arXiv, búsqueda por los dos autores juntos (sin restringir título)

`WebFetch` a `export.arxiv.org/api/query?search_query=au:Gonzalez-Ferez
AND au:Weidemuller` → **0 resultados**, en cualquier tema. Esto es más
fuerte que el punto 3: confirma que González-Férez y Weidemüller no tienen
NINGÚN preprint conjunto en arXiv — el paper de 2007 nunca se depositó ahí
(común en algunos grupos/épocas para papers sólo de revista).

### 5. Página del grupo de teoría en Heidelberg (donde trabajaba Schmelcher
en 2007)

Los resultados de `WebSearch` apuntaban a
`physi.uni-heidelberg.de/Forschung/ka/Theorie/Publications/` como posible
listado con PDFs. `WebFetch` → **HTTP 404 Not Found**, en esa URL y en
variantes sin la barra final y en el directorio padre. El sitio se ha
reestructurado desde 2007 (Schmelcher se trasladó a Hamburgo hace años); la
ruta indexada por el buscador ya no existe.

### 6. Repositorio institucional de la Universidad de Granada (digibug)

González-Férez está afiliada actualmente a la UGR, y su paper de 2024 (el
que se está replicando en este plan) SÍ está en abierto ahí:
`digibug.ugr.es/bitstream/handle/10481/93714/PhysRevA.110.013314.pdf`. Se
buscó el paper de 2007 por el mismo camino (`WebSearch` con `digibug.ugr.es`
+ términos del título, y `WebFetch` directo al buscador simple de digibug)
→ sin resultado para este artículo (el de 2024 sí aparece, el de 2007 no;
probable que el mandato de depósito en abierto de la UGR no cubra
retroactivamente un paper de una época en la que la autora ni siquiera
estaba afiliada allí — publicado durante o después de su etapa en
Heidelberg, antes de su incorporación a Granada).

### 7. Semantic Scholar (API)

`WebFetch` a la API de búsqueda de Semantic Scholar → **HTTP 429 Too Many
Requests**, en dos intentos con una espera entre medias. No se pudo
consultar `openAccessPdf`.

## Lo que SÍ se consiguió (colateral, útil de todos modos)

Al buscar el paper de 2007 se recuperó el HTML completo del paper de 2024
(`arxiv.org/html/2401.09618`, el que da origen a todo este plan), y de ahí
tres datos que sirven para la Fase 4 en cuanto se desbloquee:

1. **La cita exacta** tal como aparece en el paper de 2024: "R.
   González-Férez, M. Weidemüller, and P. Schmelcher, Phys. Rev. A 76,
   023402 (2007)" — sin número de arXiv (confirma el punto 4: no hay
   preprint).
2. **Confirmación independiente de la Fase 2**: la Ec. 8 del paper 2024 y la
   definición del operador de acoplamiento que da justo después,

   ```
   A_ij = ⟨Ψᵢ|T|Ψⱼ⟩ - (ħ²/m)⟨Ψᵢ|d/dR|Ψⱼ⟩ d/dR
   ```

   coinciden EXACTAMENTE con la ecuación derivada de forma independiente en
   `docs/analysis_fase2_canales_acoplados.md` (el término ⟨Ψᵢ|T|Ψⱼ⟩ es
   -ħ²/(2m)·B_ij, y el segundo término es -ħ²/m·A_ij·d/dR, con la misma A_ij
   de la Fase 1) — una confirmación externa, tardía pero bienvenida, de que
   la Fase 2 no tiene un error de signo o de factor.
3. La única mención del método de normalización en el paper de 2024 es esta
   frase (Apéndice A, cerca de la Ec. 15): *"Note that the scattering state
   χ^d_i is energy normalized, whereas computationally it is L² normalized.
   A numerical way to obtain these energy normalized wave functions is
   described in Ref. [González-Férez et al. 2007]."* — es decir, el propio
   paper de 2024 remite el método COMPLETO a la referencia de 2007 sin
   reproducir ni una ecuación. No hay atajo: hace falta el texto de 2007
   para tener el método real, no una aproximación de lo que probablemente
   dice.

## Por qué no se improvisa un sustituto

Existen métodos ESTÁNDAR de teoría de la dispersión para convertir una
función de onda del continuo normalizada en caja (L²) a normalización en
energía (por ejemplo, vía el factor de densidad de estados dk/dE en el
límite asintótico, o el emparejamiento de fase con una onda plana/esférica
libre) que aparecen en textos de mecánica cuántica estándar. Pero el plan
maestro y el usuario piden específicamente EL método de esta referencia, no
"un método estándar plausible" — y la razón para esa exigencia (evitar que
una ronda futura repita el patrón de "premisa equivocada" que ya afectó a
este proyecto más de una vez, según el propio `CLAUDE.md`) se aplica igual
aquí: sin el texto, no se puede saber si el método de 2007 tiene un detalle
concreto (convención de normalización, tratamiento del borde de la caja,
corrección por el rango finito del potencial, etc.) distinto del genérico
de libro de texto. Aproximarlo sería exactamente el tipo de atajo que este
documento maestro prohíbe.

## Opciones para desbloquear (decisión del usuario, no de esta ronda)

1. **El usuario consigue el PDF** (acceso institucional, biblioteca, o
   contacto directo con los autores) y lo pega o lo deja en el repositorio;
   se retoma la Fase 4 leyéndolo primero, como pide el plan.
2. **El usuario autoriza explícitamente** usar un método estándar de
   normalización en energía de un texto de mecánica cuántica/teoría de la
   dispersión distinto (a elegir y citar con el mismo rigor que las fases
   anteriores), documentando que es una sustitución deliberada del método
   de la referencia original, no una réplica de ella — cambiaría el
   alcance de "réplica de Mellado-Alcedo et al. 2024" que motiva todo el
   plan, así que es una decisión que corresponde al usuario, no a esta
   sesión.
3. **Se aparca la Fase 4** y se continúa con otro trabajo mientras se
   consigue la referencia por otra vía.

Esta sesión no elige entre estas opciones: se para y se reporta, tal como
pide el documento maestro.

## Archivos tocados

- `docs/PLAN_nonadiabatic_dynamics.md` (ESTADO ACTUAL: Fase 4 BLOQUEADA)
- `docs/analysis_fase4_busqueda_referencia.md` (este documento)

Ningún archivo de `src/` ni `tests/` se tocó en esta ronda: no hay
`decay_rates.py`, no hay tests nuevos. La suite del proyecto sigue en el
mismo estado que al cierre de la Fase 3 (118 passed, 16 deselected en
`pytest -m "not slow"`).

## Referencias

- González-Férez, R., Weidemüller, M. & Schmelcher, P., *Phys. Rev. A*
  **76**, 023402 (2007), DOI: 10.1103/PhysRevA.76.023402 — **texto
  completo NO conseguido en esta ronda**, ver intentos arriba.
- Mellado-Alcedo, Guttridge, Cornish, Sadeghpour & González-Férez,
  *Phys. Rev. A* **110**, 013314 (2024), arXiv:2401.09618 — cita la
  referencia de 2007 para el método, sin reproducirlo; confirma
  independientemente la Ec. 8/operador A_ij ya usados en
  `docs/analysis_fase2_canales_acoplados.md`.
