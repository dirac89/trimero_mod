# .claude — Configuración de Proyecto Trimero

Directorio de configuración específica del proyecto para Claude Code.

## Contenido

### Documentación Principal

- **[CLAUDE.md](CLAUDE.md)** — Instrucciones, reglas y guía de desarrollo
  - La tabla de los **siete paquetes de sistemas** y su script de producción
  - Reglas de estilo, de coste y de rectificación documental
  - Preguntas frecuentes

- **[ARCHITECTURE.md](ARCHITECTURE.md)** — Documentación de arquitectura
  - **Fuente de verdad del código**: capas, grafo real de dependencias
  - Invariantes, deuda técnica conocida y puntos de extensión

- **[settings.json](settings.json)** — Configuración de Claude Code
  - Modelo, permisos (`allow`/`deny`) y hooks del proyecto

### Skills Personalizadas

Ubicadas en `skills/`:

| Skill | Comando | Propósito |
|-------|---------|-----------|
| [sweep](skills/sweep/SKILL.md) | `/sweep` | Barrido caro: declarar sistema → estimar coste → ejecutar → auditar |
| [dataset-check](skills/dataset-check/SKILL.md) | `/dataset-check` | Audita un `.npz` de resultados |
| [quick-test](skills/quick-test/SKILL.md) | `/quick-test` | Validación rápida |
| [physics-review](skills/physics-review/SKILL.md) | `/physics-review` | Revisa cambios de física |
| [research-doc](skills/research-doc/SKILL.md) | `/research-doc` | Documenta investigación en `docs/` |

### Subagentes

Ubicados en `agents/`:

| Agente | Propósito |
|--------|-----------|
| [physics-reviewer](agents/physics-reviewer.md) | Revisión de física sobre un diff; sólo lectura, no edita |
| [dataset-auditor](agents/dataset-auditor.md) | Auditoría numérica de un `.npz` de resultados |
| [docs-curator](agents/docs-curator.md) | Cierre documental de una ronda: `analysis_*`, `INDEX.md`, revocaciones |

### Hooks

Ubicados en `hooks/`, referenciados desde `settings.json`:

| Hook | Cuándo | Qué hace |
|------|--------|----------|
| `block-graphify-out.sh` | antes de un `Bash` | **Bloquea** un `git add`/`git commit` que arrastre `graphify-out/` |
| `shared-layer-warning.sh` | tras editar un fichero | Avisa si se tocó `mathlib/` o `basis/`: de ahí dependen los goldens |
| `remind-full-suite.sh` | antes de `git commit` | Recuerda la suite completa (~11 min), no sólo `-m "not slow"` |

## Cómo Usar

1. **Para entender el proyecto** (en este orden):
   - [`docs/STATUS.md`](../docs/STATUS.md) — la física vigente. **Empieza aquí.**
   - [ARCHITECTURE.md](ARCHITECTURE.md) — el código y sus fronteras
   - [CLAUDE.md](CLAUDE.md) — reglas de desarrollo

   ⚠️ El repositorio cubre **siete paquetes de sistemas físicos** con motores
   parcialmente compartidos. Confundirlos ya costó dos rondas de trabajo,
   documentadas en `docs/PLAN_figuras_publicacion.md`.

2. **Para ejecutar simulaciones**: `/sweep` (declara sistema y molécula, estima
   el coste antes de gastar CPU) o `/quick-test` para validación rápida.

3. **Antes de hacer cambios de física**: `/physics-review`.

4. **Al cerrar una ronda**: `/research-doc`, y registra en `docs/INDEX.md`.

## Estructura Visual

```
trimero_mod/
├── .claude/                          # ← TÚ ESTÁS AQUÍ
│   ├── CLAUDE.md                     # Reglas y guía
│   ├── ARCHITECTURE.md               # Fuente de verdad del código
│   ├── settings.json                 # Modelo, permisos, hooks
│   ├── README.md                     # Este archivo
│   ├── agents/                       # physics-reviewer, dataset-auditor, docs-curator
│   ├── hooks/                        # los tres scripts de la tabla de arriba
│   └── skills/                       # sweep, dataset-check, quick-test, …
├── src/trimero/
│   ├── mathlib/  basis/              # compartidos: de aquí dependen los goldens
│   ├── simulation/  visualization/
│   └── systems/                      # los siete paquetes + rb_atom, polar_molecule
├── scripts/                          # producción y análisis (~14)
├── scripts/archive/                  # las rondas de exploración
├── tests/                            # goldens del legado en …/characterization/
├── data/Wavefunction/                # entrada (sólo sistemas con Fermi)
├── docs/STATUS.md                    # ← la física vigente, en una página
├── plots/<sistema>/{data,figures}/   # resultados verificados; se commitean
└── pyproject.toml
```

## Próximos Pasos

✓ Configuración completada
→ Lee [`docs/STATUS.md`](../docs/STATUS.md) para saber qué es lo vigente
→ `poetry run pytest -m "not slow"` para validar el proyecto (~25 s)
→ `/sweep` para tu primer cálculo

---

**Última actualización**: 2026-08-30
