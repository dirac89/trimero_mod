---
name: quick-test
description: Validación rápida tras un cambio — suite sin los tests lentos y, si el cambio toca un sistema, un barrido corto de pocos puntos. Segundos, no minutos.
---

# Skill: quick-test

## Propósito

Comprobar en menos de un minuto que un cambio no ha roto nada, **antes** de
gastar minutos en un barrido completo o en la suite entera.

## Trigger

`/quick-test`, o: «prueba rápida», «valida el flujo», «¿sigue pasando todo?».

## Flujo

### 1. La suite rápida — siempre lo primero

```bash
poetry run pytest -m "not slow"    # ~25 s
```

Deselecciona los tests de caracterización end-to-end (~3.5 min cada uno). Para
saber cuántos son ahora mismo, sin copiar cifras a mano:

```bash
poetry run pytest --collect-only -q | tail -1              # total
poetry run pytest -m "not slow" --collect-only -q | tail -1  # rápidos / total
```

### 2. Barrido corto — sólo si el cambio toca un sistema

Pocos puntos, sin figura. **La molécula es obligatoria**, no hay default:

```bash
# polar (KRb o RbCs)
poetry run python scripts/compute_bop_curve.py \
    --molecule rbcs --n-manifold 25 --mj 0 \
    --rmin 400 --rmax 500 --step 50 --no-plot

# perturbador neutro
poetry run python scripts/compute_trimer_curves.py --symmetry Sigma --rmax 600
```

Escribe en `plots/<sistema>/data/`: si no quieres tocar los datasets buenos, usa
`--npz-dir` (polar) o un `--out` temporal.

### 3. Validación

- ✓ La suite rápida pasa entera
- ✓ La matriz se construye con la dimensión esperada
- ✓ Autovalores reales y ordenados; sin NaN ni Inf
- ✓ El peso de manifold `W` está por encima del umbral en los puntos calculados

### 4. Reporte

Tiempo, tests pasados/fallados, dimensión del bloque, primeros autovalores y
`W`. Si algo falla, **el fallo tal cual** — no lo resumas a «un test rojo».

## Qué NO es esto

- **No sustituye a la suite completa antes de commitear** (`poetry run pytest`,
  ~11 min, incluye los goldens del legado). Y si el cambio toca `mathlib/` o
  `basis/`, la completa es obligatoria: de esa capa dependen los goldens de
  todos los sistemas.
- **No es la vía para validar el camino legado.**
  `Trimer_energies_field(n1, dc_field_au)` es código **congelado**; sus goldens
  ya lo cubren en la suite `slow`. No lo llames como validación de rutina.
- Un golden que se mueve no se regenera: se reporta. Ver `.claude/CLAUDE.md` §5.

## Cuándo Usar

Después de un cambio de física o matemáticas, al empezar a trabajar en el
repositorio, y antes de lanzar cualquier barrido largo.
