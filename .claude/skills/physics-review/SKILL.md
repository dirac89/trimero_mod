---
name: physics-review
description: Revisa un cambio de física (Hamiltoniano, potenciales, base, unidades, álgebra angular) antes de gastar CPU en él. No ejecuta barridos — analiza y emite veredicto.
---

# Skill: physics-review

## Propósito

Detectar inconsistencias en un cambio de física **antes** de ejecutar
simulaciones costosas o de commitear algo que mueva un golden.

## Trigger

`/physics-review`, o: «revisa la física», «valida el cambio de potencial»,
«¿esto está bien planteado?».

## Flujo

### Paso 0 — ¿de qué sistema es el cambio? (obligatorio, primero)

Siete paquetes, motores parcialmente compartidos. Este paso es el que más rondas
ha salvado; ver la tabla completa en `.claude/CLAUDE.md` y el grafo en
`.claude/ARCHITECTURE.md`.

| paquete | Hamiltoniano |
|---|---|
| `polar_rydberg/` (genérico), `rb_krb_polar/`, `rb_rbcs_polar/` | `H_A + B N² − d·(F_ion + F_elec)` |
| `double_polar_rydberg/` | `H_A + Σᵢ[B Nᵢ² − dᵢ·F_ryd(Rᵢ)] + V_dd + H_F` |
| `rb_neutral_perturber/` | `H_A + F_ext·r + V_Fermi` (legado congelado: `trimer.py`, `fermi_potentials.py`) |
| `hybrid_neutral_polar/` | `H_A + H_mol(R₂) + V_Fermi^π(R₁)` |
| `nonadiabatic_dynamics/` | ecuación nuclear sobre las BOP |
| compartido | `mathlib/`, `basis/`, `rb_atom.py`, `polar_molecule.py` |

Tres banderas rojas inmediatas:

- **`V_Fermi` en un sistema polar puro** → ✗. El pseudopotencial de contacto
  modela un perturbador **neutro**.
- **`M_J` usado en el sistema neutro** → ✗. Ahí no hay rotor: el buen número
  cuántico es `m_l` (Σ ≡ m_l=0, Π ≡ |m_l|=1).
- **Molécula ambigua** → para. KRb: B=1.114 GHz, d=0.566 D. RbCs: B=490.17 MHz,
  d=1.225 D. Confundirlas costó una ronda entera
  (`docs/PLAN_figuras_publicacion.md`).

### Paso 1 — Delegar en el subagente

Para cualquier revisión que no sea trivial, delega en **`physics-reviewer`**
(`.claude/agents/physics-reviewer.md`): leer Hamiltonianos consume mucho
contexto, y el subagente devuelve el veredicto sin arrastrarlo a la sesión.
Pásale el diff o los ficheros tocados.

### Paso 2 — Checklist

1. **Hamiltoniano**: hermítico, diagonal real, dimensión coherente, bloqueado
   por el buen número cuántico correcto.
2. **Base**: polar → manifold `(n,l≥3)` + (n+1)d + (n+2)p + (n+3)s (**tres**
   vecinos; única definición en `rb_defects.neighbor_levels()`). Neutro →
   n=35 (l≥3) + 38s + 37p + 36d.
3. **Potenciales**: Fermi real, → 0 en r→∞, sin divergencia en r→0, evaluado en
   los nodos de tabla. Carga-dipolo con el escalado 1/R⁴ del ion.
4. **Álgebra angular**: convención canónica de `wigner_3j`/Gaunt. Cambiarla
   mueve goldens de todos los sistemas (precedente: `586776d`).
5. **Unidades**: Bohr y Hartree dentro; GHz **sólo** en la salida; cero de
   energía declarado.
6. **Estabilidad**: sin divergencias, autovalores reales donde la física lo pide.
7. **Goldens**: ⚠️ si el cambio puede mover
   `tests/systems/rb_neutral_perturber/characterization/`, **para y repórtalo**.
   No se regenera un golden para que encaje.

### Paso 3 — Salida

```
SISTEMA: <paquete> · molécula <X> · buen nº cuántico <M_J|m_l>
CAPA COMPARTIDA: sí/no

✓ Hamiltoniano: hermítico, dim 1113 para M_J=0 — correcto
⚠ Base: (n+2)p reescrito a mano en scripts/foo.py:42; usa neighbor_levels()
✗ Unidades: la conversión a GHz entra en el bucle (charge_dipole.py:118)

GOLDENS EN RIESGO: no
VEREDICTO: válido con reservas
SIGUIENTE PASO: /quick-test, y barrido corto antes del completo
```

## Notas

- **No ejecuta barridos.** Comprobaciones numéricas de segundos, sí.
- Si el veredicto es válido: `/quick-test` y luego `/sweep`.
