---
name: research-doc
description: Crear y guardar documentación de investigación en docs/ siguiendo el formato estándar
---

# Skill: research-doc

## Propósito
Documentar cualquier investigación (papers, análisis, teoría, web) en formato Markdown estructurado dentro de `docs/`.

## Trigger
- Usuario dice: "documenta la investigación", "crea research doc", "guarda lo investigado"
- O: `/research-doc`

## Flujo

1. **Identificar tipo de investigación**:
   - ¿Paper/artículo científico?
   - ¿Análisis de datos?
   - ¿Notas teóricas?
   - ¿Búsqueda web?
   - ¿Comparación de métodos?

2. **Recopilar información**:
   - Resumen del tema
   - Conceptos clave
   - Fórmulas/ecuaciones relevantes
   - Fuentes (DOI, URL, referencias)

3. **Crear documento Markdown**:
   - Ubicación: `docs/`
   - Nombre: `{tipo}_{tema_corto}.md`
   - Estructura: encabezado, resumen, contenido, conclusiones, referencias

4. **Guardar y registrar**:
   - Archivo creado en `docs/`
   - Commit con mensaje: "docs: add research on [tema]"
   - Actualizar índice si existe `docs/INDEX.md`

## Plantilla Estándar

```markdown
# {Título descriptivo}

**Fecha**: {YYYY-MM-DD}  
**Autor**: {Nombre}  
**Relevancia**: {Por qué es importante para trimero}  
**Tipo**: {paper|analysis|theory|research|comparison}  

## Resumen
{1-2 párrafos resumiendo lo investigado}

## Palabras Clave
- keyword1
- keyword2
- keyword3

## Contenido Principal

### Sección 1
{Contenido según el tema}

### Sección 2
{Desarrollo}

## Fórmulas Clave
```math
E = \hbar \omega
H = H_0 + V_{int}
```

## Conclusiones y Aplicación al Proyecto
{Cómo se conecta con trimero_mod, qué cambios o mejoras habilita}

## Referencias
- [Nombre Corto](URL/DOI) — Descripción
- [Paper XYZ](doi.org/...) — Año, autor, resumen

## Notas Adicionales
{Links, Follow-ups, temas relacionados}
```

## Ejemplos de Nombres

| Tema | Nombre del Archivo |
|------|------------------|
| Potenciales de Fermi | `docs/research_fermi_potentials.md` |
| Paper sobre diagonalización | `docs/paper_scipy_eigensolver.md` |
| Análisis de convergencia | `docs/analysis_eigenvalue_convergence.md` |
| Notas sobre armónicos esféricos | `docs/theory_spherical_harmonics.md` |
| Comparación C++ vs Python | `docs/comparison_cpp_python_performance.md` |
| Búsqueda sobre campos eléctricos | `docs/research_electric_field_interaction.md` |

## Ejemplo Completo

```markdown
# Potenciales de Fermi en Sistemas Trimoleculares

**Fecha**: 2026-08-18  
**Autor**: Javier Aguilera  
**Relevancia**: Base teórica de FermiPotentials() en src/fermi_potentials.py  
**Tipo**: research  

## Resumen
Los potenciales de Fermi modelan la interacción entre átomos fermionicos 
mediante un decaimiento exponencial. Este documento resume la formulación 
matemática y parámetros típicos.

## Palabras Clave
- Fermi gas
- Scattering length
- Contacto interaction
- Ultracold atoms

## Contenido Principal

### Formulación Matemática
V(r) = -V₀ * exp(-r/a₀)

### Parámetros Típicos
- a₀ (longitud de dispersión): 50-500 Bohr
- V₀ (profundidad): 0.1-1.0 Hartree

## Aplicación al Proyecto
FermiPotentials usa esta fórmula en get_potential(r_au).
Parámetros se ajustan según el tipo de átomo (Rb, K, Li).

## Referencias
- [Quantum Mechanics of Atomic Systems](doi.org/...) — Cohen-Tannoudji et al., 2019
- [Fermi Contact in Ultracold Gases](URL) — Review, 2020
```

## Cuándo Usar

✓ **Usar este skill cuando**:
- Investigues un concepto nuevo para el proyecto
- Leas un paper relevante
- Analices resultados de simulaciones
- Explores optimizaciones o métodos alternativos

✗ **No necesario para**:
- Cambios de código simples (usa git commit)
- Conversaciones rápidas sobre bugs
- Decisiones de triviales

## Ventajas

- **Trazabilidad**: Queda registrado QUÉ se investigó y CUÁNDO
- **Reutilización**: Futuros trabajos usan el mismo conocimiento
- **Contexto**: Próximas personas entienden las decisiones de diseño
- **Reproducibilidad**: Referencias completamente documentadas

## Notas

- Los documentos se versionan con git
- Usa Markdown para todo: tablas, código, LaTeX para fórmulas
- Enlaza con otros docs usando `[referencia](path_archivo.md)`
- Revisa `docs/` regularmente para temas relacionados que expandir
