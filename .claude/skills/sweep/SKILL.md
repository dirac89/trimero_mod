---
name: sweep
description: Protocolo del barrido caro — declara sistema y molécula, estima y anuncia el coste, ejecuta en background y audita el .npz. Usa esto para cualquier cálculo de curvas (BOP, orientación, campo, trímero, híbrido), no lances los scripts a mano.
---

# Skill: sweep

## Propósito

Un barrido son minutos de diagonalizaciones (~1.14 s/punto para n=25). Los dos
modos de tirarlos a la basura están documentados en este repositorio:
calcularlos **del sistema equivocado** (`docs/PLAN_figuras_publicacion.md`) y
darlos por buenos cuando **la rama se perdió** (commits `a0730b1`, `b85295d`).
Este protocolo cierra los dos.

## Trigger

`/sweep`, o: «calcula la curva…», «lanza el barrido», «saca la BOP de…».

## Paso 1 — Declarar la identidad del cálculo (antes de nada)

Escribe explícitamente, y **contrasta con `docs/STATUS.md`**:

- **sistema** (paquete de `src/trimero/systems/`)
- **molécula** — `krb` o `rbcs`; nunca implícita
- **`n`**, **`N_max`**, **`M_J`** o simetría (`Sigma`/`Pi`), **campo**
- **rango**: `rmin`, `rmax`, `step`

Si alguno no está claro en la petición, **pregunta**; no lo supongas. Un
parámetro físico supuesto en silencio es exactamente lo que costó la ronda
KRb/RbCs.

| sistema | script |
|---|---|
| polar (BOP) | `compute_bop_curve.py --molecule <krb\|rbcs>` |
| polar (orientación) | `compute_orientation_curve.py --molecule <..>` |
| polar en campo DC | `compute_field_curves.py --molecule <..>` |
| dos moléculas polares | `compute_double_rbcs_curves.py --geometry <symmetric\|unilateral\|both>` |
| perturbador neutro | `compute_trimer_curves.py --symmetry <Sigma\|Pi>` |
| híbrido | `compute_hybrid_curves.py` |
| dinámica no adiabática | `analyze_nonadiabatic_n25_crossing[_wide].py` |

Comparativas (`compare_*.py`, `plot_orientation_alignment.py`) **no recalculan**:
leen `.npz` existentes. No pasan por este protocolo.

## Paso 2 — Estimar y ANUNCIAR el coste

```
n_puntos = (rmax − rmin)/step + 1        (× nº de bloques M_J × nº de campos)
coste    ≈ n_puntos × 1.14 s             (n=25, N_max=6; escala con la dimensión)
```

Los defaults de `compute_bop_curve.py` (400→1800 a₀, paso 5, `--mj 0 1`) son
**281 puntos × 2 bloques ≈ 11 min**, no un cálculo instantáneo.

- Di los minutos **antes** de ejecutar, con la cuenta a la vista.
- **> ~10 min: pide confirmación** al usuario antes de lanzar.
- ¿Dudas del setup? Barrido corto primero:
  `--rmin 400 --rmax 500 --step 50 --no-plot` (~30 s).
- ¿Existe ya el `.npz`? `--reuse` en vez de rebarrer.

## Paso 3 — Ejecutar

- Lánzalo **en background** si pasa de un par de minutos, para no bloquear la
  sesión, y repite en la orden **todos los flags explícitos** (esa línea se
  copia luego al documento de la ronda).
- Vigila que no aparezcan avisos de solapamiento o de peso durante la corrida:
  son la señal temprana de rama perdida.

## Paso 4 — Auditar antes de dar nada por bueno

Delega en el subagente **`dataset-auditor`** (o `/dataset-check`) sobre el `.npz`
recién escrito. No presentes resultados ni los commitees antes de eso.

Si el auditor devuelve tramos con `overlap` bajo o `W` por debajo del umbral:
recalcula ese tramo con `--min-substep` menor **antes** de aceptar la curva; no
subas el umbral para que encaje.

## Paso 5 — Cerrar

- Resume: identidad del cálculo, coste real, fichero(s) escrito(s), anclas
  numéricas y lo que dijo el auditor.
- Si la ronda produce conclusiones físicas, `/research-doc` (o el subagente
  `docs-curator`) para dejarlas en `docs/` y en `docs/INDEX.md`.

## Ejemplo

```
Usuario: saca la BOP de RbCs para n=27

Sistema: rb_rbcs_polar (polar puro, sin Fermi) · molécula rbcs (B=490.17 MHz,
d=1.225 D) · n=27, N_max=6, M_J = 0 y 1 · R ∈ [400,1800] a0, paso 5.
Coste: 281 puntos × 2 bloques × ~1.3 s/punto ≈ 12 min. ¿Lanzo?

  poetry run python scripts/compute_bop_curve.py \
      --molecule rbcs --n-manifold 27 --n-max 6 --mj 0 1 \
      --rmin 400 --rmax 1800 --step 5

→ plots/rb_rbcs_polar/data/fig1_ad_MJ0_n27.npz  (+ MJ1, + PNG)
→ auditoría: apto; W ∈ [0.71, 0.99], sin NaN, K constante salvo en R=1520 a0.
```
