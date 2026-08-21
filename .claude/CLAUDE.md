# Trimero Atómico - Instrucciones de Proyecto

## Descripción General

Estructura electrónica de **moléculas Rydberg de largo alcance**: un átomo de Rb
en estado Rydberg perturbado por un compañero, resuelta por diagonalización del
Hamiltoniano en una base acoplada.

⚠️ **El repositorio cubre DOS sistemas físicos distintos.** Es lo primero que
hay que saber; mezclarlos ya produjo varias rondas con premisa equivocada.

| | `rb_krb_polar` | `rb_neutral_perturber` |
|---|---|---|
| perturbador | KRb, **polar** | átomo/molécula **neutra** |
| interacción | carga-dipolo, `−d·F_ryd` | pseudopotencial de Fermi (contacto) |
| Hamiltoniano | `H_ad = H_A + H_mol` | `H = H_a + V_Fermi` |
| estado | **vigente** | verde y congelado, línea en pausa |

**El pseudopotencial de Fermi NO interviene en el sistema polar.** Si te ves
añadiendo `V_Fermi` a una curva de Rb*-KRb, para y lee `docs/STATUS.md`.

**Status**: migración C++ → Python completada; reorganizado por sistema físico
el 2026-08-20. Rama: `migrate-python`.

## Documentos de entrada

| documento | para qué |
|---|---|
| [`docs/STATUS.md`](../docs/STATUS.md) | **empieza aquí**: la física vigente en una página |
| [`.claude/ARCHITECTURE.md`](ARCHITECTURE.md) | el código: capas, dependencias, invariantes, deuda |
| [`docs/INDEX.md`](../docs/INDEX.md) | catálogo completo de documentación |

## Estructura del Proyecto

```
src/trimero/
├── mathlib/          angular.py (3j/Gaunt), special.py, laplacian.py
├── basis/            quantum.py (CoupledBasis), radial.py (RadialBasis)
├── simulation/       bop_tracking.py (trace_curve)
└── systems/
    ├── rb_atom.py                 Atom.E_Rb() — defectos cuánticos, COMPARTIDO
    ├── rb_krb_polar/              charge_dipole.py, bop_system.py, rb_defects.py
    └── rb_neutral_perturber/      fermi_krb.py, fermi_potentials.py, trimer.py

scripts/compute_bop_curve.py       único script de producción
scripts/archive/                   los 12 scripts de exploración
tests/{basis,systems}/             48 tests; goldens del legado en
                                   systems/rb_neutral_perturber/characterization/
data/Wavefunction/                 .dat de entrada
docs/                              referencia activa + STATUS.md
docs/archive/rb_neutral_perturber/ la saga del pseudopotencial (otro sistema)
plots/                             sólo lo vigente; el resto en plots/archive/
```

## Cómo Ejecutar

```bash
poetry install

# Curvas BOP de Rb*-KRb (el camino vigente)
poetry run python scripts/compute_bop_curve.py --n-manifold 25 --mj 0 1

# Tests: rápidos mientras iteras, completo antes de commitear
poetry run pytest -m "not slow"    # ~25 s
poetry run pytest                  # ~6 min, incluye los goldens del legado
```

`compute_bop_curve.py` escribe `plots/fig1_ad_MJ<..>_n<n>.npz` y su PNG.
El camino legado (`Trimer_energies_field`) produce `Trimer_R_sp_wave_*.dat`.

## Reglas de Desarrollo

### 1. Estilo y Formato
- **Python**: 3.13+, sigue PEP 8
- **Dependencias**: Gestiona con Poetry (`pyproject.toml`)
- **Imports**: Ordena alfabéticamente, sin imports circulares
- **Tipos**: Usa type hints en funciones públicas

### 2. Cambios de Código
- **Respeta la frontera entre sistemas**: nada de `rb_krb_polar/` debería
  importar de `rb_neutral_perturber/` ni al revés. Hoy hay una excepción
  (`BOPSystem` → `fermi_krb`), documentada como deuda en `ARCHITECTURE.md`;
  no añadas la segunda.
- **`mathlib/` y `basis/` son compartidos**: no deben conocer un sistema concreto.
- **No refactorices sin necesidad**: el código está estructurado, mantén el diseño
- **Mantén la trazabilidad**: los `print()` de depuración están permitidos si
  marcan puntos clave
- **Datos de entrada**: valida que los archivos en `data/Wavefunction/` existan
  antes de usarlos

### 3. Testing y Validación
- `poetry run pytest -m "not slow"` mientras iteras; **completo antes de commitear**
- **Los golden files son ley.** `tests/systems/rb_neutral_perturber/characterization/`
  fija el camino legado bit a bit (`rtol=1e-12`). Si un cambio los mueve, es un
  **cambio de física**: para y repórtalo, no ajustes el golden para que encaje.
  Conservan a propósito dos bugs del legado (unidades en `EhtoGHz`, asimetría de
  matriz) para que el golden siga siendo fiel al C++ original.
- `tests/systems/rb_krb_polar/test_regression_fig1.py` ancla los números
  verificados de la Fig. 1: −23.100 GHz, E(1800 a₀) = −0.338 GHz, 8 mínimos.
- Valida con parámetros pequeños antes de lanzar barridos largos
  (un barrido de 281 puntos son ~5 min)

### 4. Documentación
- Comenta el **POR QUÉ**, no el **QUÉ** (el código es autodocumentado)
- Documenta cambios en física o matemática que no sean obvios
- Si cambias la estructura de paquetes, actualiza `.claude/ARCHITECTURE.md`
- Si cambias la física vigente, actualiza `docs/STATUS.md`
- Mantén el README.md actualizado si cambias la interfaz pública

### 5. Documentación de Investigación
**REGLA IMPORTANTE**: Toda información de investigación debe documentarse en formato Markdown y guardarse en `docs/`

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
- [Fuente 2]: enlace/DOI
```

#### Ejemplos
- `docs/research_fermi_potentials.md` — Investigación sobre potenciales Fermi
- `docs/paper_diagonalization_methods.md` — Resumen de métodos numéricos para diagonalización
- `docs/theory_spherical_harmonics.md` — Notas sobre armónicos esféricos
- `docs/analysis_eigenvalue_convergence.md` — Análisis de convergencia numérica

### 6. Rama y Commits
- **Rama principal**: `master`
- **rama de trabajo**: `migrate-python` (activa, para completar migración)
- **Commits**: Mensajes claros en inglés o español, referencian la física si es relevante

### 7. Qué NO se commitea
- **`graphify-out/` NO se commitea.** Está en `.gitignore` y fue purgado del
  histórico el 2026-08-21. Es un artefacto **derivado y regenerable**: sale de
  `/graphify` a partir del propio repositorio, pesa ~2.8 MB en 100 ficheros y
  se reescribe entero en cada reconstrucción, así que versionarlo sólo añadía
  ruido y peso sin aportar nada que no se pueda regenerar.
  - Se mantiene **en local**: no lo borres del disco, sólo no entra en git.
  - Para regenerarlo tras un cambio grande de estructura: `/graphify`
    (reconstrucción completa) o `/graphify . --update` (incremental).
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
**Primero decide de qué sistema es.** Ese es el punto de todo el árbol.

- **Sistema polar (vigente)**: nuevos términos de `H_mol` en
  `systems/rb_krb_polar/charge_dipole.py`; el montaje de la base y el cero de
  energía en `bop_system.py`. Otro manifold es sólo `--n-manifold`: no hay nada
  cableado a n=24 ni n=25.
- **Perturbador neutro**: `systems/rb_neutral_perturber/fermi_krb.py` (versión
  vectorizada). El camino legado (`trimer.py`, `fermi_potentials.py`) está
  **congelado**: no lo extiendas, replica en la capa moderna.
- **Matemáticas nuevas**: `mathlib/angular.py` (álgebra angular) o
  `mathlib/special.py` (funciones especiales). No metas contexto físico ahí.

### Paralelización
Cada punto de R es independiente: el barrido de `compute_bop_curve.py` se
paraleliza con `multiprocessing` o `joblib` sin más cuidado que mantener la
salida ordenada por R. Es donde está el 95 % del tiempo (~1.1 s/punto).

### Análisis de Resultados
- `compute_bop_curve.py` guarda `.npz` con `R`, `E` (GHz, relativa al cero de
  energía), `K` (índice de la curva), `W` (peso de manifold) y `spectrum`.
- Los `.dat` del camino legado contienen `R, eigenvalue_1, eigenvalue_2, ...`
- Valida que los autovalores sean reales y que el peso de manifold sea ~1 en el
  borde superior: si no, la curva ha cambiado de carácter y el resultado no es
  lo que crees.

## Skills Disponibles

- `/run-simulation`: Ejecuta la simulación principal con parámetros configurables
- `/quick-test`: Prueba rápida con parámetros pequeños
- `/physics-review`: Revisa cambios de física antes de validarlos

## Preguntas Frecuentes (en código)

**¿Qué pasa si falta un archivo `.dat`?**
- La simulación fallará con `FileNotFoundError`. Revisa `data/Wavefunction/` y nombres exactos.

**¿Cómo cambio parámetros físicos?**
- Camino vigente: opciones de `scripts/compute_bop_curve.py` (`--n-manifold`,
  `--mj`, `--rmin/--rmax/--step`, `--weight`), o construye un `BOPSystem`
  directamente.
- Camino legado: llama `Trimer_energies_field(n1, dc_field_au)`.

**¿Por qué `hamiltonian()` tiene un flag `fermi`?**
- Es deuda técnica, no diseño. `fermi=False` no significa «apagar un término de
  este sistema»: significa que el pseudopotencial **no forma parte** del modelo
  polar. Ver `.claude/ARCHITECTURE.md` §Deuda técnica.

**Un golden file ha cambiado, ¿lo regenero?**
- **No.** Es un cambio de física disfrazado. Para y repórtalo.

**¿Puedo manejar valores complejos en la matriz?**
- Sí, scipy soporta matrices complejas. Asegúrate que los autovalores sean reales si es esperado físicamente.

---

**Autor Original**: Migración a Python por Javier Aguilera  
**Email**: aguilerajavier58@gmail.com
