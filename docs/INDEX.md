# Índice de Documentación de Investigación

Este directorio contiene toda la documentación de investigación, análisis teórico y notas científicas del proyecto. El repositorio cubre **dos sistemas físicos distintos** —Rb*-KRb
polar y perturbador neutro—; ver `STATUS.md` para la separación.

**Última actualización**: 2026-08-21

> **Empieza por [`STATUS.md`](STATUS.md)** — estado vigente del proyecto en una
> página. Este índice es el catálogo completo; `STATUS.md` dice qué es actual.

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

### Vigentes

| Documento | Contenido |
|---|---|
| `STATUS.md` | 🟢 **PUNTO DE ENTRADA**: modelo físico vigente de Rb*-KRb, base, criterio de carácter, δ₀(ns), deuda técnica y enlaces a la referencia activa |
| `AUDIT_MIGRACIÓN_PYTHON.md` | Auditoría de la migración C++ → Python |
| `DISEÑO_HAMILTONIANO_KRBC.md` / `_v2.md` | Diseño del Hamiltoniano Rb*-KRb; base acoplada y bloqueo por M_J |
| `analysis_validacion_carga_dipolo.md` | `B·N²` + campo del ion Rb⁺: derivación, 4 tests analíticos, escalado 1/R⁴ |
| `analysis_fig1_carga_dipolo_sin_fermi.md` | ✅ **EL CÁLCULO BUENO para Rb*-KRb**: `H_ad = H_A + H_mol`, sin pseudopotencial de Fermi. Curvas M_J=0 y M_J=1 en R∈[400,1800] a₀ completas, sin remapeo ni ventanas. Explica qué documentos anteriores quedan con premisa equivocada |
| `analysis_base_correcta_3_vecinos.md` | ⚠️ **CORRIGE la base electrónica de toda la sesión**: el paper usa manifold + (n+1)d + (n+2)p + (n+3)s, no sólo (n+3)s. dim(M_J=0) 1016→1064; el dominio del remapeo se acorta ~49 a₀; convergencia y orientación apenas cambian |
| `analysis_verificacion_tabla_I.md` | ⚠️ Tabla I del paper: 6/8 niveles concuerdan a ±0.004 GHz; sistemático real de 0.3 GHz en la serie s, localizado en δ₀(ns) |
| `analysis_campo_electron_rydberg.md` | Campo del electrón Rydberg (Ec. A.6-A.10): expansión multipolar, 9 tests, validación contra cuadratura 2D, dominancia frente al ion |
| `analysis_wigner3j_orden_canonico.md` | 🔧 **Arreglo de `wigner_3j` (orden canónico)**: causa raíz de la asimetría de hermiticidad a n=35 (2.8e−13 → 3.0e−19 rel). Simetría bit a bit bajo permutaciones y volteo de m, validado contra Racah racional exacto |
| `analysis_hibrido_fases12_golden_regen.md` | ✅ **Cierre de la ronda wigner + Fases 1-2 del híbrido**: golden `fig1_ad_MJ0_n25.npz` regenerado (Opción A; mejora de precisión 1.12×, no corrección de física; tabla old→new), números de §5 verificados, NaN de `dg_integrals` diagnosticado como uso fuera de dominio (0 en producción), y primeras curvas E0(R2) del sistema Rb*-Rb-RbCs a R1∈{600,900,1100} a₀ (`scripts/compute_hybrid_curves.py`). Goldens MJ0 y MJ1 regenerados (K idéntico bit a bit en ambos) |
| `superpowers/specs/2026-08-18-refactor-estructura-design.md` | Diseño del refactor a paquete `trimero`: capas, ABC `Hamiltonian`, goldens y plan de 11 pasos. Incluye dos hallazgos medidos: bug de unidades ×1000 en `EhtoGHz` y coste de 34 días para `n1=35` |

### Vigentes — sistema de perturbador NEUTRO

| Documento | Contenido |
|---|---|
| `analysis_trimero_lineal_campo_dc.md` | 🟢 **EL CÁLCULO BUENO para el perturbador neutro**: trímero lineal simétrico Rb(5s)Rb(35,l≥3)Rb(5s) en campo DC, validado contra Aguilera-Fernández 2016 con 6 anclas cuantitativas dentro del 3 %. Establece que la simetría es `m_l` (no M_J), que los dos perturbadores son un factor de paridad, y que n=35 lee las tablas sin remapeo. Documenta 3 bugs de datos del legado |

### Archivados — sistema de perturbador NEUTRO

Escritos cuando se creía que el pseudopotencial aplicaba a Rb*-KRb. Técnicamente
correctos pero con esa premisa; **no aplicables al sistema polar**. Siguen
siendo la referencia de la derivación del pseudopotencial y de la procedencia de
los datos, ahora usada por `analysis_trimero_lineal_campo_dc.md`.

| Documento | Contenido |
|---|---|
| `archive/rb_neutral_perturber/analysis_procedencia_rvsAS_rvsAP.md` | ⚠️ Punto abierto: procedencia de los `.dat` de entrada |
| `archive/rb_neutral_perturber/analysis_extension_dominio_fig1.md` | ❌ **Extensión del dominio DESCARTADA con medida**: poner V_Fermi=0 más allá del remapeo mete un salto de 12-14 GHz (41-43 % de la ligadura). Incluye las curvas M_J=0 de n=24 y n=25 con base completa y las tres tendencias de la Fig. 1(a) |
| `archive/rb_neutral_perturber/analysis_ventana_exclusion_resonancia.md` | ✅ **Cierre del bloque de la resonancia p**: ventana de exclusión `R ∈ [536.4, 589.1] a₀` con criterio calculado, verificación de que `inverse` no alteró nada fuera, tabla de características re-hecha y nota de limitación para manuscrito. Omont 1977 **no** implementado, por decisión |
| `archive/rb_neutral_perturber/analysis_interpolacion_polo_Ap.md` | Interpolación de 1/A_p a través del polo; ⚠️ CORRIGE la conclusión de `analysis_resonancia_onda_p.md`: el pozo butterfly no era artefacto |
| `archive/rb_neutral_perturber/analysis_resonancia_onda_p.md` | Resonancia de forma p: posición confirmada (24.8 vs 23 meV); ⚠️ su §6.2 está corregida — el pozo profundo NO era artefacto de interpolación |
| `archive/rb_neutral_perturber/analysis_curva_bop_MJ0.md` | Curva BOP M_J=0: seguimiento adiabático vs diabático, y el butterfly de onda p (explicado ya como divergencia del rango cero) |
| `archive/rb_neutral_perturber/analysis_pseudopotencial_fermi_krb.md` | Pseudopotencial de Fermi s+p, remapeo k(R), dominancia frente al carga-dipolo |

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
