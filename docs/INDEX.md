# Índice de Documentación de Investigación

Este directorio contiene toda la documentación de investigación, análisis teórico y notas científicas relacionadas con el proyecto Trimero.

**Última actualización**: 2026-08-18

## Documentación Disponible

### Investigación General (research_*)
- Investigación sobre temas específicos relevantes para el proyecto

### Papers y Referencias (paper_*)
- Resúmenes de artículos científicos
- Métodos numéricos
- Publicaciones relevantes

### Análisis y Resultados (analysis_*)
- Análisis de datos de simulaciones
- Comparaciones de métodos
- Validación de resultados

### Teoría y Derivaciones (theory_*)
- Notas sobre conceptos teóricos
- Derivaciones matemáticas
- Explicaciones de fórmulas

### Comparaciones (comparison_*)
- Comparación entre implementaciones
- Benchmarks de rendimiento
- Validación C++ vs Python

---

## Documentos Existentes

| Documento | Contenido |
|---|---|
| `AUDIT_MIGRACIÓN_PYTHON.md` | Auditoría de la migración C++ → Python |
| `DISEÑO_HAMILTONIANO_KRBC.md` / `_v2.md` | Diseño del Hamiltoniano Rb*-KRb; base acoplada y bloqueo por M_J |
| `analysis_procedencia_rvsAS_rvsAP.md` | ⚠️ Punto abierto: procedencia de los `.dat` de entrada |
| `analysis_validacion_carga_dipolo.md` | `B·N²` + campo del ion Rb⁺: derivación, 4 tests analíticos, escalado 1/R⁴ |
| `analysis_interpolacion_polo_Ap.md` | Interpolación de 1/A_p a través del polo; ⚠️ CORRIGE la conclusión de `analysis_resonancia_onda_p.md`: el pozo butterfly no era artefacto |
| `analysis_resonancia_onda_p.md` | Resonancia de forma p: posición confirmada (24.8 vs 23 meV) pero el pozo de −6.5 THz es artefacto de interpolar sobre un hueco de malla |
| `analysis_curva_bop_MJ0.md` | ⚠️ Curva BOP M_J=0: seguimiento adiabático vs diabático, y un butterfly de onda p a −6.5 THz sin explicar |
| `analysis_verificacion_tabla_I.md` | ⚠️ Tabla I del paper: 6/8 niveles concuerdan a ±0.004 GHz; sistemático real de 0.3 GHz en la serie s, localizado en δ₀(ns) |
| `analysis_pseudopotencial_fermi_krb.md` | Pseudopotencial de Fermi s+p, remapeo k(R), dominancia frente al carga-dipolo |
| `analysis_campo_electron_rydberg.md` | Campo del electrón Rydberg (Ec. A.6-A.10): expansión multipolar, 9 tests, validación contra cuadratura 2D, dominancia frente al ion |
| `superpowers/specs/2026-08-18-refactor-estructura-design.md` | Diseño del refactor a paquete `trimero`: capas, ABC `Hamiltonian`, goldens y plan de 11 pasos. Incluye dos hallazgos medidos: bug de unidades ×1000 en `EhtoGHz` y coste de 34 días para `n1=35` |

---

## Cómo Agregar Nueva Documentación

Usa el skill `/research-doc`:

```
/research-doc

Usuario: Acabo de leer un paper sobre optimización de matrices...
```

El skill te guiará para:
1. Elegir el tipo de documento
2. Crear la estructura Markdown
3. Guardar en `docs/` con nombre apropiado
4. Actualizar este índice

## Estructura Recomendada

Cada documento debe tener:
- **Encabezado principal** (# Título)
- **Metadatos** (Fecha, Autor, Relevancia)
- **Resumen** (1-2 párrafos)
- **Contenido** (secciones según el tema)
- **Conclusiones** (aplicación a trimero)
- **Referencias** (fuentes con DOI/URL)

## Ejemplo de Nombres

```
docs/
├── research_fermi_potentials.md
├── paper_scipy_eigensolver.md
├── theory_spherical_harmonics.md
├── analysis_convergence_n35.md
├── comparison_cpp_python_perf.md
└── INDEX.md (este archivo)
```

## Búsqueda Rápida

**¿Cómo buscar documentación existente?**
```bash
# Buscar por tipo
ls docs/research_*.md
ls docs/paper_*.md
ls docs/theory_*.md

# Buscar por tema
grep -l "fermi\|potential" docs/*.md
grep -l "diagonalización\|eigenvalue" docs/*.md
```

---

**Política de Documentación**: Toda investigación, análisis y hallazgos importantes deben quedar documentados aquí para referencia futura.
