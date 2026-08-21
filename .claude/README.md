# .claude — Configuración de Proyecto Trimero

Directorio de configuración específica del proyecto para Claude Code.

## Contenido

### Documentación Principal

- **[CLAUDE.md](CLAUDE.md)** — Instrucciones, reglas y guía de desarrollo del proyecto
  - Cómo ejecutar y extender la simulación
  - Reglas de estilo y contribución
  - Preguntas frecuentes

- **[ARCHITECTURE.md](ARCHITECTURE.md)** — Documentación de arquitectura
  - **Los dos sistemas físicos y por qué están separados**
  - Capas, grafo real de dependencias, invariantes
  - Deuda técnica conocida y puntos de extensión

- **[settings.json](settings.json)** — Configuración de Claude Code
  - Modelo de IA preferido
  - Permisos de herramientas
  - Idioma (español)

### Skills Personalizadas

Ubicadas en `skills/`:

| Skill | Comando | Propósito |
|-------|---------|-----------|
| [run-simulation](skills/run-simulation/SKILL.md) | `/run-simulation` | Ejecuta la simulación con parámetros |
| [quick-test](skills/quick-test/SKILL.md) | `/quick-test` | Prueba rápida de validación |
| [physics-review](skills/physics-review/SKILL.md) | `/physics-review` | Revisa cambios de física |
| [research-doc](skills/research-doc/SKILL.md) | `/research-doc` | Documenta investigación en `docs/` |

## Cómo Usar

1. **Para entender el proyecto** (en este orden):
   - [`docs/STATUS.md`](../docs/STATUS.md) — la física vigente. **Empieza aquí.**
   - [ARCHITECTURE.md](ARCHITECTURE.md) — el código y sus fronteras
   - [CLAUDE.md](CLAUDE.md) — reglas de desarrollo

   ⚠️ El repositorio cubre **dos sistemas físicos distintos** (Rb*-KRb polar y
   perturbador neutro). Confundirlos ya costó varias rondas de trabajo con
   premisa equivocada.

2. **Para ejecutar simulaciones**:
   - Usa `/run-simulation` con parámetros específicos
   - O usa `/quick-test` para validación rápida

3. **Antes de hacer cambios de física**:
   - Invoca `/physics-review` para revisar la matemática
   - Asegúrate de coherencia en unidades y fórmulas

## Estructura Visual

```
trimero_mod/
├── .claude/                          # ← TÚ ESTÁS AQUÍ
│   ├── CLAUDE.md                     # Reglas y guía
│   ├── ARCHITECTURE.md               # Documentación técnica
│   ├── settings.json                 # Configuración
│   ├── README.md                     # Este archivo
│   └── skills/
│       ├── run-simulation/
│       ├── quick-test/
│       └── physics-review/
├── src/trimero/
│   ├── mathlib/                      # primitivas matemáticas
│   ├── basis/                        # CoupledBasis, RadialBasis (compartido)
│   ├── simulation/                   # trace_curve
│   └── systems/
│       ├── rb_atom.py                # defectos cuánticos de Rb (compartido)
│       ├── rb_krb_polar/             # ← sistema VIGENTE
│       └── rb_neutral_perturber/     # ← el otro sistema, congelado
├── scripts/compute_bop_curve.py      # único script de producción
├── scripts/archive/                  # los 12 de exploración
├── tests/{basis,systems}/            # 48 tests
├── data/Wavefunction/                # Archivos de entrada
├── docs/STATUS.md                    # ← la física vigente, en una página
├── docs/archive/                     # material del otro sistema
├── plots/{rb_krb_polar,rb_neutral_perturber}/  # espejo de systems/, lo vigente
├── pyproject.toml                    # Dependencias (Poetry)
└── ...
```

## Próximos Pasos

✓ Configuración completada  
→ Lee [`docs/STATUS.md`](../docs/STATUS.md) para saber qué es lo vigente  
→ `poetry run pytest -m "not slow"` para validar el proyecto (~25 s)  
→ `poetry run python scripts/compute_bop_curve.py --n-manifold 25 --mj 0` para tu primer cálculo  

---

**Última actualización**: 2026-08-20
