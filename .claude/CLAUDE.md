# Trimero Atómico - Instrucciones de Proyecto

## Descripción General

Estructura electrónica de **moléculas Rydberg de largo alcance**: un átomo de Rb
en estado Rydberg perturbado por uno o varios compañeros, resuelta por
diagonalización del Hamiltoniano en una base acoplada.

⚠️ **El repositorio cubre SIETE paquetes de sistemas físicos.** Es lo primero
que hay que saber. Confundirlos ya costó rondas enteras de trabajo: hay dos
documentadas en `docs/PLAN_figuras_publicacion.md` (KRb calculado creyéndose
RbCs, y después una "corrección" que cambió la etiqueta en vez del cálculo).

| paquete (`src/trimero/systems/`) | Hamiltoniano | script de producción | estado |
|---|---|---|---|
| `polar_rydberg/` | `H_A + B N² − d·(F_ion + F_elec)` | — (motor genérico) | vigente |
| `rb_krb_polar/` | idem, molécula KRb | `compute_bop_curve.py --molecule krb` | vigente |
| `rb_rbcs_polar/` | idem, molécula RbCs | `compute_bop_curve.py --molecule rbcs`, `compute_field_curves.py` | vigente |
| `double_polar_rydberg/` | `H_A + Σᵢ[B Nᵢ² − dᵢ·F_ryd(Rᵢ)] + V_dd + H_F` | `compute_double_rbcs_curves.py` | vigente |
| `rb_neutral_perturber/` | `H_A + F_ext·r + V_Fermi` | `compute_trimer_curves.py` | vigente |
| `hybrid_neutral_polar/` | `H_A + H_mol(R₂) + V_Fermi^π(R₁)` | `compute_hybrid_curves.py` | vigente |
| `nonadiabatic_dynamics/` | ecuación nuclear sobre BOP y acoplamientos derivados | `analyze_nonadiabatic_n25_crossing[_wide].py` | vigente |

Dentro de `rb_neutral_perturber/`, `trimer.py` y `fermi_potentials.py` son el
**camino legado congelado**, protegido bit a bit por golden files. No se
extiende: se replica en la capa moderna (`linear_trimer.py`, `fermi_krb.py`).

El grafo de dependencias entre paquetes, los invariantes y la deuda técnica
están en [`ARCHITECTURE.md`](ARCHITECTURE.md); es la fuente de verdad del
código y este documento no debe contradecirlo.

**El pseudopotencial de Fermi NO interviene en los sistemas polares puros.** Si
te ves añadiendo `V_Fermi` a una curva de Rb*-KRb o Rb*-RbCs, para y lee
`docs/STATUS.md`.

Y al revés: en el sistema neutro el buen número cuántico es **`m_l`**, no `M_J`
(no hay rotor). Σ ≡ m_l=0, Π ≡ |m_l|=1.

**Status**: migración C++ → Python completada. Rama de trabajo: `migrate-python`.
El estado científico vigente, con sus anclas numéricas verificadas, está en
[`docs/STATUS.md`](../docs/STATUS.md).

## Documentos de entrada

| documento | para qué |
|---|---|
| [`docs/STATUS.md`](../docs/STATUS.md) | **empieza aquí**: la física vigente en una página |
| [`.claude/ARCHITECTURE.md`](ARCHITECTURE.md) | el código: capas, dependencias, invariantes, deuda |
| [`docs/INDEX.md`](../docs/INDEX.md) | catálogo completo de documentación |
| [`plots/README.md`](../plots/README.md) | dónde vive cada resultado |

## Estructura del Proyecto

```
src/trimero/
├── mathlib/          angular.py (3j/Gaunt), special.py, laplacian.py
├── basis/            quantum.py (CoupledBasis), radial.py (RadialBasis)
├── simulation/       bop_tracking.py (trace_curve)
├── visualization/    geometry_diagram.py — no importa sistemas
└── systems/          rb_atom.py y polar_molecule.py (compartidos)
                      + los siete paquetes de la tabla de arriba

scripts/              ~14 scripts de producción y análisis (ver §Cómo Ejecutar)
scripts/archive/      los 12 scripts de las rondas de exploración
tests/                goldens del legado en
                      systems/rb_neutral_perturber/characterization/
data/Wavefunction/    .dat de entrada (sólo los sistemas con Fermi los usan)
docs/                 referencia activa + STATUS.md + INDEX.md
plots/<sistema>/      data/ (.npz) y figures/ (.png); histórico en plots/archive/
```

Para el número exacto de tests **no lo copies aquí**: pregúntaselo a pytest.

```bash
poetry run pytest --collect-only -q | tail -1
```

## Cómo Ejecutar

```bash
poetry install

# BOP polar — la molécula es OBLIGATORIA, no hay default
poetry run python scripts/compute_bop_curve.py --molecule rbcs --n-manifold 25 --mj 0 1
poetry run python scripts/compute_orientation_curve.py --molecule rbcs --n-manifold 25

# Campo DC (polar) y dos moléculas polares
poetry run python scripts/compute_field_curves.py --molecule rbcs
poetry run python scripts/compute_double_rbcs_curves.py --geometry both

# Perturbador neutro: trímero lineal simétrico en campo DC
poetry run python scripts/compute_trimer_curves.py --symmetry Sigma

# Híbrido y dinámica no adiabática
poetry run python scripts/compute_hybrid_curves.py
poetry run python scripts/analyze_nonadiabatic_n25_crossing.py

# Convergencia y comparativas (baratas; no recalculan barridos)
poetry run python scripts/check_polar_convergence.py --molecule rbcs
poetry run python scripts/compare_bop_curves_n.py --molecule rbcs

# Tests: rápidos mientras iteras, completo antes de commitear
poetry run pytest -m "not slow"    # ~25 s
poetry run pytest                  # ~11 min, incluye los goldens del legado
```

Cada script escribe en `plots/<sistema>/data/*.npz` y `plots/<sistema>/figures/*.png`.
El camino legado (`Trimer_energies_field`) produce `Trimer_R_sp_wave_*.dat`.

## Reglas de Desarrollo

### 1. Estilo y Formato
- **Python**: 3.13+, sigue PEP 8
- **Dependencias**: Gestiona con Poetry (`pyproject.toml`)
- **Imports**: Ordena alfabéticamente, sin imports circulares
- **Tipos**: Usa type hints en funciones públicas

### 2. Cambios de Código
- **Respeta la frontera entre sistemas**: un paquete de `systems/` no importa de
  otro paquete hermano. Hoy hay excepciones documentadas como deuda en
  `ARCHITECTURE.md` (`BOPSystem` → `fermi_krb`; `rb_defects` usado desde el lado
  neutro); no añadas la siguiente.
- **`mathlib/`, `basis/` y `visualization/` son compartidos**: no deben conocer
  un sistema concreto.
- **No refactorices sin necesidad**: el código está estructurado, mantén el diseño
- **Mantén la trazabilidad**: los `print()` de depuración están permitidos si
  marcan puntos clave
- **Datos de entrada**: valida que los archivos en `data/Wavefunction/` existan
  antes de usarlos

### 3. Sistema y molécula, siempre explícitos
Antes de lanzar cualquier cálculo, **declara y contrasta con `docs/STATUS.md`**:
sistema, molécula, `n`, `M_J` (o simetría) y campo. No hay atajo: la tabla de
arriba tiene siete filas y varias comparten motor.

- **Prohibidos los defaults silenciosos en parámetros físicos.** Un default no
  aparece en la línea de órdenes que se copia al documento, así que nada delata
  el error. Fue exactamente la causa raíz de la ronda perdida de KRb/RbCs
  (`--molecule` con `default="krb"`). Los cuatro scripts polares llevan hoy
  `--molecule ... required=True`; los scripts nuevos hacen lo mismo.
- Los `.npz` polares graban la especie y sus constantes (`molecule`,
  `B_hz`, `d_debye`, `schema_version`). Un `.npz` sin ese campo **se rechaza**,
  no se supone la molécula.

### 4. El coste se anuncia antes de calcular
Un barrido cuesta minutos de diagonalizaciones (~1.14 s/punto para n=25; 281
puntos ≈ 5 min; la ventana ancha de la Fase 6b fueron ~25 min).

- **Estima y comunica los minutos ANTES de ejecutar**, no después.
- Por encima de ~10 min, **confirma con el usuario** antes de lanzar.
- Los barridos largos van en background, no bloqueando la sesión.
- Valida primero con un rango corto (`--rmin 400 --rmax 500 --step 50`).
- Al terminar, **audita el `.npz`** (`/dataset-check`): autovalores reales, sin
  NaN, peso de manifold ~1 en el borde superior, `overlap` por encima del
  umbral. Si el peso o el solapamiento caen, la curva cambió de carácter y el
  resultado no es lo que crees — es lo que provocó los commits de recuperación
  `a0730b1` y `b85295d`.

### 5. Testing y Validación
- `poetry run pytest -m "not slow"` mientras iteras; **completo antes de commitear**
- **Los golden files son ley.** `tests/systems/rb_neutral_perturber/characterization/`
  fija el camino legado bit a bit (`rtol=1e-12`). Si un cambio los mueve, es un
  **cambio de física**: para y repórtalo, no ajustes el golden para que encaje.
  Conservan a propósito dos bugs del legado (unidades en `EhtoGHz`, asimetría de
  matriz) para que el golden siga siendo fiel al C++ original.
  `.claude/settings.json` los protege además con una regla `deny`.
- Tocar `mathlib/` o `basis/` es tocar la capa de la que dependen los goldens de
  todos los sistemas: suite **completa** antes de commitear. Precedente: el
  `wigner_3j` canónico (`586776d`, `642fabd`) obligó a regenerar goldens y a
  rehacer las fases 1-2 del híbrido.
- `tests/systems/rb_krb_polar/test_regression_fig1.py` ancla los números
  verificados de la Fig. 1: −23.100 GHz, E(1800 a₀) = −0.338 GHz, 8 mínimos.

### 6. Documentación
- Comenta el **POR QUÉ**, no el **QUÉ** (el código es autodocumentado)
- Documenta cambios en física o matemática que no sean obvios
- Si cambias la estructura de paquetes, actualiza `.claude/ARCHITECTURE.md`
- Si cambias la física vigente, actualiza `docs/STATUS.md`
- Mantén el README.md actualizado si cambias la interfaz pública

### 7. Rectificar sin borrar
Cuando un resultado nuevo revoca un documento anterior:

- **marca el viejo como revocado**, con fecha y con qué lo sustituye; no lo
  reescribas en silencio ni lo borres;
- **documenta el mecanismo del error**, no sólo el error. Lo reutilizable es
  «un default silencioso no aparece en la orden que copias», no «era KRb».
- Registra el documento nuevo en `docs/INDEX.md`.

Es el patrón que ya usan `docs/PLAN_figuras_publicacion.md` (rectificación del
2026-08-24) y `docs/analysis_faseA_curvas_bop_varios_n.md` §0.

### 8. Documentación de Investigación
**REGLA IMPORTANTE**: Toda información de investigación debe documentarse en
formato Markdown y guardarse en `docs/` (usa `/research-doc`).

#### Qué Documentar
- **Papers, artículos científicos**: Resumen, palabras clave, fórmulas relevantes, enlace/DOI
- **Búsquedas web**: Temas investigados, resultados útiles, fuentes confiables
- **Análisis de datos**: Interpretación de resultados, patrones observados
- **Notas teóricas**: Derivaciones, explicaciones de conceptos física
- **Metodología**: Decisiones de diseño y su justificación científica
- **Benchmarks y comparaciones**: Resultados con código C++ original, optimizaciones probadas

#### Formato de Archivos
Guardar en `docs/` con nomenclatura clara:
- `docs/research_<tema>.md` — Investigación general sobre un tema
- `docs/paper_<titulo_corto>.md` — Resumen de un paper
- `docs/analysis_<tipo>.md` — Análisis de datos
- `docs/theory_<concepto>.md` — Notas teóricas

#### Estructura Mínima de Cada Documento
```markdown
# Título de la Investigación

**Fecha**: YYYY-MM-DD
**Autor**: Nombre
**Relevancia**: Por qué es importante para el proyecto

## Resumen
[1-2 párrafos resumiendo lo investigado]

## Contenido Principal
[Secciones según el tipo de investigación]

## Conclusiones y Aplicación al Proyecto
[Cómo se aplica a trimero_mod]

## Referencias
- [Fuente 1]: enlace/DOI
```

### 9. Rama y Commits
- **Rama principal**: `master`
- **rama de trabajo**: `migrate-python` (activa)
- **Commits**: Mensajes claros en inglés o español, referencian la física si es relevante

### 10. Qué NO se commitea
- **`graphify-out/` NO se commitea.** Está en `.gitignore` y fue purgado del
  histórico el 2026-08-21. Es un artefacto **derivado y regenerable**: sale de
  `/graphify` a partir del propio repositorio, pesa ~2.8 MB en 100 ficheros y
  se reescribe entero en cada reconstrucción, así que versionarlo sólo añadía
  ruido y peso sin aportar nada que no se pueda regenerar.
  - Se mantiene **en local**: no lo borres del disco, sólo no entra en git.
  - Para regenerarlo tras un cambio grande de estructura: `/graphify`
    (reconstrucción completa) o `/graphify . --update` (incremental).
  - Un hook de `.claude/settings.json` bloquea el `git add` que lo arrastre.
- **`plots/` SÍ se commitea**, y nunca va al `.gitignore`: las figuras y sus
  `.npz` son resultados verificados, no artefactos regenerables baratos
  (un barrido cuesta ~5 min de diagonalizaciones).

## Dependencias Principales

| Librería | Versión | Propósito |
|----------|---------|-----------|
| numpy | ^2.3.1 | Matrices, operaciones numéricas |
| scipy | ^1.16.0 | Álgebra lineal (diagonalización) |
| matplotlib | ^3.10.3 | Visualización de resultados |

## Extensión del Proyecto

### Agregar Nueva Física
**Primero decide de qué paquete es.** Ese es el punto de todo el árbol.

- **Polar (KRb, RbCs o una molécula nueva)**: añade la especie al catálogo
  `systems/polar_molecule.py` y usa `polar_rydberg/PolarBOPSystem`, el motor
  genérico. Los términos del Hamiltoniano viven en
  `rb_krb_polar/charge_dipole.py` (nombre histórico, física genérica: ver
  Deuda §1 de `ARCHITECTURE.md`). Otro manifold es sólo `--n-manifold`.
- **Dos moléculas polares**: `double_polar_rydberg/` (base de dos rotores,
  Hamiltoniano disperso, geometrías simétrica y unilateral).
- **Perturbador neutro**: `rb_neutral_perturber/linear_trimer.py`
  (`SymmetricLinearTrimer`) sobre `fermi_krb.py`. Toda la geometría está en
  `parity_factor`: las configuraciones asimétrica y planar del paper se añaden
  ahí. El camino legado está **congelado**: replica en la capa moderna.
- **Híbrido neutro+polar**: `hybrid_neutral_polar/`.
- **Dinámica nuclear sobre las BOP**: `nonadiabatic_dynamics/`.
- **Matemáticas nuevas**: `mathlib/angular.py` (álgebra angular) o
  `mathlib/special.py` (funciones especiales). No metas contexto físico ahí.

### Paralelización
Cada punto de R es independiente. `compute_double_rbcs_curves.py` ya usa
`ProcessPoolExecutor`; el resto se paraleliza igual, manteniendo la salida
ordenada por R. Es donde está el 95 % del tiempo (~1.1 s/punto).
Ojo: el seguimiento por solapamiento **sí** consume los resultados en orden.

### Análisis de Resultados
- Los `.npz` polares guardan `R`, `E` (GHz, relativa al cero de energía), `K`
  (índice de la curva), `W` (peso de manifold), `spectrum` y los metadatos de
  molécula/base/criterio. Los del sistema doble añaden `overlap`,
  `warning_overlap` y `COS1`/`COS2`.
- Los `.dat` del camino legado contienen `R, eigenvalue_1, eigenvalue_2, ...`
- Auditoría estándar tras cada barrido: `/dataset-check`.

## Skills y Agentes Disponibles

| comando | para qué |
|---|---|
| `/sweep` | protocolo del barrido caro: declarar sistema → estimar coste → ejecutar → auditar |
| `/dataset-check` | audita un `.npz` de resultados (carácter, solapamiento, NaN, metadatos) |
| `/quick-test` | validación rápida (suite `not slow` + barrido corto) |
| `/physics-review` | revisa un cambio de física antes de gastar CPU en él |
| `/research-doc` | documenta investigación en `docs/` y la registra en `INDEX.md` |
| `/graphify` | reconstruye el grafo de conocimiento (`graphify-out/`, local) |

Subagentes en `.claude/agents/`: `physics-reviewer` (revisión de física, sólo
lectura), `dataset-auditor` (auditoría de `.npz`), `docs-curator` (cierre
documental de una ronda).

## Preguntas Frecuentes (en código)

**¿Qué pasa si falta un archivo `.dat`?**
- La simulación fallará con `FileNotFoundError`. Revisa `data/Wavefunction/` y
  nombres exactos. Sólo los sistemas con pseudopotencial de Fermi los necesitan;
  los polares puros no leen longitudes de dispersión.

**¿Cómo cambio parámetros físicos?**
- Por línea de órdenes del script de producción del sistema (`--molecule`,
  `--n-manifold`, `--mj`, `--rmin/--rmax/--step`, `--weight`, `--fields`…), o
  construyendo el sistema directamente (`PolarBOPSystem`, `RbTwoRbCsSystem`,
  `SymmetricLinearTrimer`, `HybridNeutralPolar`).
- Camino legado: `Trimer_energies_field(n1, dc_field_au)`.

**¿Por qué `hamiltonian()` tiene un flag `fermi`?**
- Es deuda técnica, no diseño. `fermi=False` no significa «apagar un término de
  este sistema»: significa que el pseudopotencial **no forma parte** del modelo
  polar. Los sistemas polares nuevos usan `PolarBOPSystem`, donde el flag no
  existe. Ver `.claude/ARCHITECTURE.md` §Deuda técnica.

**Un golden file ha cambiado, ¿lo regenero?**
- **No.** Es un cambio de física disfrazado. Para y repórtalo.

**La curva me sale discontinua o el peso de manifold cae, ¿subo el umbral?**
- No sin entender por qué. El índice espectral no identifica un objeto físico a
  través de cruces evitados; para eso está el seguimiento por solapamiento con
  subpaso adaptativo (`simulation/bop_tracking.py`,
  `--catastrophic-overlap`/`--min-substep`). Audita con `/dataset-check` antes
  de tocar umbrales.

**¿Puedo manejar valores complejos en la matriz?**
- Sí, scipy soporta matrices complejas. Asegúrate que los autovalores sean
  reales si es esperado físicamente.

---

**Autor Original**: Migración a Python por Javier Aguilera
**Email**: aguilerajavier58@gmail.com
