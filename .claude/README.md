# .claude — Configuración de Proyecto Trimero

Directorio de configuración específica del proyecto para Claude Code.

## Contenido

### Documentación Principal

- **[CLAUDE.md](CLAUDE.md)** — Instrucciones, reglas y guía de desarrollo del proyecto
  - Cómo ejecutar y extender la simulación
  - Reglas de estilo y contribución
  - Preguntas frecuentes

- **[ARCHITECTURE.md](ARCHITECTURE.md)** — Documentación de arquitectura
  - Capas del sistema: Datos, Física, Simulación, I/O, Control
  - Flujo de ejecución y dependencias
  - Puntos de extensión y optimización

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

1. **Para entender el proyecto**:
   - Lee [CLAUDE.md](CLAUDE.md) (reglas y flujo)
   - Lee [ARCHITECTURE.md](ARCHITECTURE.md) (diseño técnico)

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
├── src/
│   ├── main.py
│   ├── trimer.py
│   ├── atom.py
│   ├── fermi_potentials.py
│   ├── math_aux.py
│   └── laplacian.py
├── data/
│   └── Wavefunction/                 # Archivos de entrada
├── README.md                         # Documentación general del proyecto
├── pyproject.toml                    # Dependencias (Poetry)
└── ...
```

## Próximos Pasos

✓ Configuración completada  
→ Lee [CLAUDE.md](CLAUDE.md) para entender las reglas  
→ Usa `/quick-test` para validar el proyecto  
→ Usa `/run-simulation` para tu primer cálculo  

---

**Última actualización**: 2026-08-18
