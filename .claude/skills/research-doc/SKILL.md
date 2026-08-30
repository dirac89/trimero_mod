---
name: research-doc
description: Documenta una investigación o una ronda de cálculo en docs/, la registra en docs/INDEX.md y marca como revocado lo que el resultado nuevo invalide. Usa esto al cerrar cualquier ronda con conclusiones.
---

# Skill: research-doc

## Propósito

Que quien llegue sin contexto dentro de un mes sepa **qué es cierto hoy** sin
reconstruir la historia. `docs/` tiene ~35 documentos y varias rectificaciones
encadenadas: un documento sin registrar en el índice, o una conclusión revocada
sin marcar, es una trampa para la siguiente ronda.

## Trigger

`/research-doc`, o: «documenta esto», «guarda lo investigado», y al cerrar
cualquier ronda de cálculo con conclusiones.

## Flujo

### 1. Redactar
Ubicación `docs/`, nomenclatura por tipo:

| tipo | nombre |
|---|---|
| investigación general | `research_<tema>.md` |
| resumen de paper | `paper_<titulo_corto>.md` |
| análisis de datos/ronda | `analysis_<tema>.md` |
| notas teóricas | `theory_<concepto>.md` |
| plan de varias rondas | `PLAN_<tema>.md` |

Si la ronda fue un cálculo, el documento **debe** incluir:

- **la identidad completa**: sistema, molécula, `n`, `N_max`, `M_J`/simetría, campo;
- **la orden exacta** que lo produjo, con todos los flags explícitos;
- **la ruta del `.npz`** resultante;
- **las anclas numéricas** verificadas, con su valor y contra qué se comparan;
- lo que el auditor (`/dataset-check`) haya señalado, incluidos los tramos
  dudosos. No se maquilla un resultado con reservas.

### 2. Registrar en `docs/INDEX.md`
**Una entrada nueva no existe hasta que está en el índice.** Respeta la sección
y el estilo de las entradas ya presentes.

### 3. Rectificar sin borrar
Si el resultado revoca algo anterior (regla §7 de `.claude/CLAUDE.md`):

- localiza el documento afectado (`grep` por el número, la figura o la afirmación);
- añade **arriba** un aviso fechado: qué queda revocado, por qué, y qué lo sustituye;
- **no borres ni reescribas el texto revocado**;
- documenta **el mecanismo**, no sólo el error. Lo reutilizable es «un default
  silencioso no aparece en la orden que copias al documento», no «era KRb».
- si lo revocado tocaba `STATUS.md`, actualízalo también.

Formato a imitar: la rectificación del 2026-08-24 en
`docs/PLAN_figuras_publicacion.md` y `docs/analysis_faseA_curvas_bop_varios_n.md` §0.

### 4. Verificar `docs/STATUS.md`
Es la página de entrada: la tabla de sistemas, las anclas y la lista de «lo que
no está hecho» deben seguir siendo verdad tras esta ronda.

Para los pasos 2-4 puedes delegar en el subagente **`docs-curator`**, que edita
sólo dentro de `docs/`.

## Plantilla

```markdown
# {Título descriptivo}

**Fecha**: {YYYY-MM-DD}
**Autor**: {Nombre}
**Relevancia**: {por qué importa para trimero_mod}
**Tipo**: {paper|analysis|theory|research|comparison}

## Resumen
{1-2 párrafos}

## Identidad del cálculo        ← si aplica
sistema, molécula, n, N_max, M_J/simetría, campo, rango de R
orden exacta ejecutada · fichero .npz producido

## Contenido Principal
{secciones según el tema}

## Anclas numéricas
| magnitud | referencia | calculado |

## Conclusiones y Aplicación al Proyecto
{qué habilita, qué revoca}

## Referencias
- [Nombre](DOI/URL) — autor, año
```

## Cuándo Usar

✓ Al cerrar una ronda de cálculo con conclusiones; al leer un paper relevante;
al investigar un concepto nuevo; al tomar una decisión de diseño con
justificación científica.

✗ No hace falta para un bugfix trivial o una conversación rápida.
