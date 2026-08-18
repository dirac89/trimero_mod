# Trimero Atómico - Instrucciones de Proyecto

## Descripción General

Este proyecto simula la física de un **trimero atómico** (tres átomos interactuantes) usando:
- Potenciales de Fermi para la interacción atómica
- Diagonalización de matrices Hamiltonianas
- Análisis de estados propios en función de campos eléctricos

**Status**: Migración completada de C++ a Python (rama: `migrate-python`)

## Estructura del Proyecto

```
src/
├── main.py              # Punto de entrada, llama simulación principal
├── trimer.py            # Lógica principal: Trimer_energies_field()
├── atom.py              # Clase Atom: encapsula física atómica
├── fermi_potentials.py  # Clase FermiPotentials: cálculo de potenciales
├── math_aux.py          # Funciones matemáticas especiales (armónicos esféricos, etc.)
└── laplacian.py         # Interfaz de funciones matemáticas

data/
└── Wavefunction/        # Archivos de datos (*.dat, *.txt)
```

## Cómo Ejecutar

```bash
poetry install
poetry run python src/main.py
```

Salida: Archivos `Trimer_R_sp_wave_*.dat` con autovalores por radio.

## Reglas de Desarrollo

### 1. Estilo y Formato
- **Python**: 3.13+, sigue PEP 8
- **Dependencias**: Gestiona con Poetry (`pyproject.toml`)
- **Imports**: Ordena alfabéticamente, sin imports circulares
- **Tipos**: Usa type hints en funciones públicas

### 2. Cambios de Código
- **No refactorices sin necesidad**: El código está estructurado, mantén el diseño
- **Mantén la trazabilidad**: Los `print()` de depuración están permitidos si marcan puntos clave
- **Datos de entrada**: Siempre valida que los archivos en `data/Wavefunction/` existan antes de usar
- **Salida**: Guarda resultados en archivos `.dat` con formato consistente (R, autovalores)

### 3. Testing y Validación
- Usa la función `test_trimer_energies_field()` en `src/main.py` para pruebas rápidas
- Valida con parámetros pequeños (ej: `n1=5`, `dc_field_au=0.1`)
- Prueba cambios en física antes de ejecutar simulaciones grandes

### 4. Documentación
- Comenta el **POR QUÉ**, no el **QUÉ** (el código es autodocumentado)
- Documenta cambios en física o matemática que no sean obvios
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

## Dependencias Principales

| Librería | Versión | Propósito |
|----------|---------|-----------|
| numpy | ^2.3.1 | Matrices, operaciones numéricas |
| scipy | ^1.16.0 | Álgebra lineal (diagonalización) |
| matplotlib | ^3.10.3 | Visualización de resultados |

## Extensión del Proyecto

### Agregar Nueva Física
1. Extiende `Atom` en `src/atom.py` o `FermiPotentials` en `src/fermi_potentials.py`
2. Usa `math_aux.py` para funciones matemáticas complejas
3. Actualiza `Trimer_energies_field()` si cambias la matriz Hamiltoniana

### Paralización
El bucle principal de `Trimer_energies_field()` puede paralelizarse con `multiprocessing` o `joblib`, pero mantén la salida de datos ordenada por radio.

### Análisis de Resultados
- Los archivos `.dat` generados contienen: `R, eigenvalue_1, eigenvalue_2, ...`
- Usa matplotlib para graficar energías vs radio
- Valida que los autovalores sean reales y ordenados ascendentemente

## Skills Disponibles

- `/run-simulation`: Ejecuta la simulación principal con parámetros configurables
- `/quick-test`: Prueba rápida con parámetros pequeños
- `/physics-review`: Revisa cambios de física antes de validarlos

## Preguntas Frecuentes (en código)

**¿Qué pasa si falta un archivo `.dat`?**
- La simulación fallará con `FileNotFoundError`. Revisa `data/Wavefunction/` y nombres exactos.

**¿Cómo cambio parámetros físicos?**
- Edita `src/main.py` o llama `Trimer_energies_field()` con otros valores en `n1`, `dc_field_au`, etc.

**¿Puedo manejar valores complejos en la matriz?**
- Sí, scipy soporta matrices complejas. Asegúrate que los autovalores sean reales si es esperado físicamente.

---

**Autor Original**: Migración a Python por Javier Aguilera  
**Email**: aguilerajavier58@gmail.com
