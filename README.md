# trimero_mod

Estructura electrónica de **moléculas Rydberg de largo alcance**: un átomo de Rb
en estado Rydberg perturbado por un compañero, resuelta por diagonalización del
Hamiltoniano en una base acoplada. Originalmente C++, migrado a Python.

> **Empieza por [`docs/STATUS.md`](docs/STATUS.md)** — la física vigente en una
> página. Para el código, [`.claude/ARCHITECTURE.md`](.claude/ARCHITECTURE.md).

## ⚠️ Este repositorio cubre DOS sistemas físicos distintos

Es lo primero que hay que saber. Estuvieron mezclados en el mismo espacio de
nombres hasta la reorganización del 2026-08-20, y esa mezcla produjo varias
rondas de trabajo con premisa equivocada.

| | `rb_krb_polar` | `rb_neutral_perturber` |
|---|---|---|
| perturbador | KRb, **polar** | átomo/molécula **neutra** |
| interacción | carga-dipolo, `−d·F_ryd` | dispersión de contacto (pseudopotencial de Fermi) |
| Hamiltoniano | `H_ad = H_A + H_mol` | `H = H_a + V_Fermi` |
| referencia | Aguilera-Fernández 2015 / González-Férez 2015 | Aguilera-Fernández 2016 |
| estado | **vigente** | verde y congelado, línea en pausa |

**El pseudopotencial de Fermi no interviene en el sistema polar.** Si te ves
añadiendo `V_Fermi` a una curva de Rb\*-KRb, para y lee `docs/STATUS.md`.

## Estructura

```
src/trimero/
├── mathlib/          angular.py (3j/Gaunt), special.py, laplacian.py
├── basis/            quantum.py (CoupledBasis), radial.py (RadialBasis)
├── simulation/       bop_tracking.py (trace_curve)
└── systems/
    ├── rb_atom.py                 Atom.E_Rb() — defectos cuánticos, COMPARTIDO
    ├── rb_krb_polar/              charge_dipole.py, bop_system.py, rb_defects.py
    └── rb_neutral_perturber/      fermi_krb.py, fermi_potentials.py, trimer.py

scripts/compute_bop_curve.py       único script de producción
scripts/archive/                   los 12 scripts de exploración
tests/{basis,systems}/             48 tests
data/Wavefunction/                 .dat de entrada (sólo los usa el perturbador neutro)
docs/                              referencia activa + STATUS.md
docs/archive/                      material del otro sistema
plots/                             sólo lo vigente; el resto en plots/archive/
graphify-out/                      grafo de conocimiento del repositorio
```

## Instalación y uso

```sh
poetry install
```

**Curvas BOP de Rb\*-KRb** — el camino vigente:

```sh
poetry run python scripts/compute_bop_curve.py --n-manifold 25 --mj 0 1
```

Opciones: `--n-manifold`, `--mj`, `--rmin/--rmax/--step`, `--weight`,
`--no-plot`, `--reuse`. Escribe `plots/fig1_ad_MJ<..>_n<n>.npz` y su PNG.
Un barrido completo (281 puntos) tarda ~5 min: cada punto es una
diagonalización de 1113×1113.

**Camino legado** (perturbador neutro, congelado):

```python
from trimero.systems.rb_neutral_perturber.trimer import Trimer_energies_field
Trimer_energies_field(n1=5, dc_field_au=0.1)
```

Necesita los `.dat` de `data/Wavefunction/` (`rvsAS.dat`, `rvsAP.dat`,
`rvsR38s.dat`, `rvsR36d.dat`, `rvsR37p.dat`, sus derivadas `rvsDR*.dat` y
`exp_val_r.txt`). Produce `Trimer_R_sp_wave_*.dat` con `R` y los autovalores.

## Tests

```sh
poetry run pytest -m "not slow"    # 45 tests, ~25 s — mientras iteras
poetry run pytest                  # 48 tests, ~6 min — antes de commitear
```

Los `slow` son los **golden files** que congelan el camino legado bit a bit
(`rtol=1e-12`). Conservan a propósito dos bugs del C++ original (un factor de
unidades en `EhtoGHz` y una asimetría de matriz) para que el golden siga siendo
fiel. **Si un cambio mueve un golden, es un cambio de física: para y repórtalo,
no regeneres el golden.**

`tests/systems/rb_krb_polar/test_regression_fig1.py` ancla los números
verificados de la Fig. 1 (n=25, M_J=0): profundidad −23.100 GHz,
E(1800 a₀) = −0.338 GHz, 8 mínimos locales.

## Dependencias

Python ≥ 3.13 · numpy ^2.3 · scipy ^1.16 · matplotlib ^3.10 · Poetry · pytest

## Extensión

**Primero decide de qué sistema es** — ese es el punto de toda la estructura.

- **Polar (vigente)**: términos de `H_mol` en `rb_krb_polar/charge_dipole.py`;
  montaje de base y cero de energía en `bop_system.py`. Otro manifold es sólo
  `--n-manifold`: no hay nada cableado a n=24 ni n=25.
- **Perturbador neutro**: `rb_neutral_perturber/fermi_krb.py`. El camino legado
  (`trimer.py`, `fermi_potentials.py`) está congelado — no lo extiendas.
- **Matemáticas**: `mathlib/`. No metas contexto físico ahí.

Cada punto de R es independiente, así que el barrido se paraleliza con
`multiprocessing` o `joblib` sin más cuidado que mantener la salida ordenada.

## Documentación

Toda la investigación se documenta en `docs/` (ver [`docs/INDEX.md`](docs/INDEX.md)).
`docs/archive/rb_neutral_perturber/` conserva la saga del pseudopotencial:
técnicamente correcta, pero de **otro sistema físico**.

## Contacto

Migración a Python: Javier Aguilera — aguilerajavier58@gmail.com
