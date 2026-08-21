# trimero_mod

Estructura electrónica de **moléculas Rydberg de largo alcance**: un átomo de Rb
en estado Rydberg perturbado por un compañero, resuelta por diagonalización del
Hamiltoniano en una base acoplada. Originalmente C++, migrado a Python.

> **Empieza por [`docs/STATUS.md`](docs/STATUS.md)** — la física vigente en una
> página. Para el diseño del código, [`.claude/ARCHITECTURE.md`](.claude/ARCHITECTURE.md).

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

---

## Estructura de directorios

```
trimero_mod/
│
├── src/trimero/                  el paquete (20 ficheros)
│   ├── mathlib/                  primitivas matemáticas, sin contexto físico
│   │   ├── angular.py                wigner_3j, gaunt        → lo usa el lado polar
│   │   ├── special.py                Spherical, DRnl, DOlm, DPhilm, hydrogenicR
│   │   └── laplacian.py              reexporta special.py
│   ├── basis/                    enumeración de la base — COMPARTIDO
│   │   ├── quantum.py                CoupledBasis, QuantumBasisBlock
│   │   │                             |l m_l⟩ ⊗ |N M_N⟩, bloqueada por M_J
│   │   └── radial.py                 RadialBasis — tablas radiales hidrogenoides
│   ├── simulation/
│   │   └── bop_tracking.py           trace_curve — agnóstico del sistema
│   └── systems/
│       ├── rb_atom.py            Atom.E_Rb() — defectos cuánticos de Rb.
│       │                         COMPARTIDO: fuente de verdad de los dos sistemas
│       ├── rb_krb_polar/         ◄── SISTEMA VIGENTE
│       │   ├── charge_dipole.py      RydbergElectronField, ChargeDipoleHamiltonian,
│       │   │                         rydberg_diagonal, constantes de KRb
│       │   ├── bop_system.py         BOPSystem — monta base + H para un manifold
│       │   └── rb_defects.py         DELTA0_NS_PAPER, neighbor_levels(), n_star_nl()
│       └── rb_neutral_perturber/ ◄── el otro sistema, congelado
│           ├── fermi_krb.py          ScatteringLengths, FermiPseudopotential
│           ├── fermi_potentials.py   FermiPotentials       ─┐ legado C++,
│           └── trimer.py             Trimer_energies_field ─┘ golden files
│
├── scripts/
│   ├── compute_bop_curve.py      ÚNICO script de producción
│   └── archive/                  los 12 scripts de exploración (imports al día,
│                                 pero su física no es necesariamente la vigente)
│
├── tests/                        48 tests (17 ficheros)
│   ├── conftest.py                   pone src/ en el path
│   ├── basis/                        enumeración de la base (compartida)
│   └── systems/
│       ├── rb_krb_polar/             charge_dipole, rydberg_field, bop_system,
│       │                             regression_fig1
│       └── rb_neutral_perturber/     fermi_krb
│           └── characterization/     golden files G1→G4 del legado
│
├── data/Wavefunction/            22 ficheros .dat/.txt de entrada.
│                                 SÓLO los usa el perturbador neutro
├── docs/                         5 documentos de referencia activa + STATUS.md
│   └── archive/                      material del otro sistema físico
├── plots/                        sólo lo vigente: 1 PNG + los 2 .npz que lo respaldan
│   └── archive/                      11 figuras y datos de rondas superadas
│
├── graphify-out/                 grafo de conocimiento. LOCAL, NO SE VERSIONA
│                                 (regenerable con /graphify; ver CLAUDE.md §7)
└── Trimer_R_sp_wave_N35_R_300_{GHz,au}.dat
                                  salida del camino legado, versionada como muestra
```

**Qué se versiona y qué no**: `plots/` **sí** — sus figuras y `.npz` son
resultados verificados que cuestan minutos de diagonalizaciones.
`graphify-out/` **no** — se rehace entero con un comando.

---

## Cómo sacar resultados

### 0. Instalación

```sh
poetry install
```

Python ≥ 3.13 · numpy ^2.3 · scipy ^1.16 · matplotlib ^3.10 · pytest (dev).

### 1. Comprueba que el repo está sano

```sh
poetry run pytest -m "not slow"     # 45 tests, ~25 s
```

Si esto falla, para: no tiene sentido calcular nada todavía.

### 2. Calcula una curva BOP de Rb\*-KRb

Todo pasa por un único comando:

```sh
poetry run python scripts/compute_bop_curve.py --n-manifold 25 --mj 0 1
```

| opción | por defecto | qué hace |
|---|---|---|
| `--n-manifold` | `25` | n del manifold cuasi-degenerado (l ≥ 3). No hay nada cableado a un n concreto |
| `--mj` | `0 1` | bloques M_J a calcular (acepta varios) |
| `--rmin` / `--rmax` / `--step` | `400` / `1800` / `5` | malla en R, en a₀ |
| `--weight` | `0.5` | peso mínimo de manifold para aceptar la curva |
| `--npz-dir` | `plots` | dónde se guardan/leen los `.npz` |
| `--reuse` | — | reutiliza el `.npz` si existe en vez de rebarrer R |
| `--no-plot` / `--out` | — | omitir la figura / cambiar su ruta |
| `--ymin` / `--ymax` | `-25` / `1` | límites del eje de energía en la figura |
| `--system` | `polar` | `neutral` está reservado y aún lanza `NotImplementedError` |

**Coste**: ~1,1 s por punto de R. El barrido por defecto son 281 puntos
(diagonalizaciones de 1113×1113) ≈ **5 min por bloque M_J**. Para tantear
primero, acorta la malla:

```sh
poetry run python scripts/compute_bop_curve.py \
    --n-manifold 25 --mj 0 --rmin 400 --rmax 600 --step 50 --no-plot
```

### 3. Qué sale

Por pantalla, la composición de la base, el cero de energía, los umbrales
asintóticos y la forma de la curva:

```
  manifold n=25 (l=3..24) + 26d + 27p + 28s
  delta0_ns = 3.1318   cero de energía: E(n=25, l>=3) + KRb(N=0) = -8.000000000000e-04 E_h

  umbrales asintóticos 28s + KRb(N):
    ΔE(28s) + 30B (N=5) = -22.6437 GHz
    ΔE(28s) + 42B (N=6) =  -9.2757 GHz

    M_J = 0
      pozo MÁS PROFUNDO   :  -23.1002 GHz en R = 400.0 a0
      E(R = 1800 a0)      :   -0.3376 GHz   (k = 55, peso = 1.0000)
      mínimos locales     : 8
```

En disco, `plots/fig1_ad_MJ<M_J>_n<n>.npz` con cinco arrays:

| clave | forma | qué es |
|---|---|---|
| `R` | (n_R,) | malla radial en a₀ |
| `E` | (n_R,) | la curva BOP en GHz, relativa a `E_manifold + KRb(N=0)` |
| `K` | (n_R,) | índice del autovalor que la curva ocupa en cada R |
| `W` | (n_R,) | peso de carácter de manifold del estado elegido |
| `spectrum` | (n_R, 250) | los 250 autovalores más bajos, para el fondo de la figura |

Más el PNG con los dos paneles, el resto del bloque en gris, la curva de
carácter en negro y los umbrales rotacionales.

### 4. Lee los datos

```python
import numpy as np
d = np.load("plots/fig1_ad_MJ0_n25.npz")
R, E, W = d["R"], d["E"], d["W"]

print(f"pozo más profundo: {E.min():.3f} GHz en R = {R[E.argmin()]:.0f} a0")
print(f"E(R_max) = {E[-1]:.3f} GHz")
```

**Comprobación obligatoria**: `W` debe valer ~1 en el borde superior. Si no, la
curva cambió de carácter por el camino y lo que estás mirando no es lo que
crees. Es el error que motivó el criterio de carácter en vez de índice fijo.

### 5. Verifica contra los números de referencia

Para n=25, M_J=0, R ∈ [400, 1800] a₀ el resultado está verificado de forma
independiente y anclado por un test:

| magnitud | valor |
|---|---|
| pozo más profundo | **−23.100 GHz** (R = 400 a₀) |
| E(R = 1800 a₀) | **−0.338 GHz** |
| mínimos locales | **8** (R = 450, 515, 590, 670, 765, 875, 1015, 1155 a₀) |

```sh
poetry run pytest tests/systems/rb_krb_polar/test_regression_fig1.py
```

Si tu barrido no reproduce esos números, **algo se rompió**: no ajustes la
referencia para que encaje.

### 6. Camino legado (perturbador neutro)

Congelado. Se ejecuta desde Python, escribe en el directorio actual y necesita
los `.dat` de `data/Wavefunction/`:

```python
from trimero.systems.rb_neutral_perturber.trimer import Trimer_energies_field
Trimer_energies_field(n1=5, dc_field_au=0.1)
```

Produce `Trimer_R_sp_wave_*.dat` con `R` seguido de los autovalores por fila.
⚠️ El fichero rotulado `GHz` contiene en realidad MHz: es un bug del C++
original que **se conserva a propósito** para que los golden files sigan siendo
fieles.

---

## Tests

```sh
poetry run pytest -m "not slow"    # 45 tests, ~25 s — mientras iteras
poetry run pytest                  # 48 tests, ~6 min — antes de commitear
```

Los `slow` son los **golden files** que congelan el camino legado bit a bit
(`rtol=1e-12`), en cadena G1 (funciones especiales) → G2 (`FermiPotentials`) →
G3 (matrices) → G4 (autovalores end-to-end). Conservan a propósito dos bugs del
C++ original —el factor de unidades en `EhtoGHz` y una asimetría de matriz—
para no perder fidelidad.

**Si un cambio mueve un golden, es un cambio de física: para y repórtalo, no
regeneres el golden.**

---

## Extensión

**Primero decide de qué sistema es** — ese es el punto de toda la estructura.

- **Polar (vigente)**: términos de `H_mol` en `rb_krb_polar/charge_dipole.py`;
  montaje de base y cero de energía en `bop_system.py`.
- **Perturbador neutro**: `rb_neutral_perturber/fermi_krb.py`. El camino legado
  (`trimer.py`, `fermi_potentials.py`) está congelado — no lo extiendas,
  replica en la capa moderna.
- **Matemáticas**: `mathlib/`. No metas contexto físico ahí.

Regla de capas: `mathlib` → `basis` → `systems`. Nada de `mathlib/` o `basis/`
debería conocer un sistema concreto; las excepciones actuales están listadas
como deuda en `.claude/ARCHITECTURE.md`.

Cada punto de R es independiente, así que el barrido se paraleliza con
`multiprocessing` o `joblib` sin más cuidado que mantener la salida ordenada.

---

## Documentación

Toda la investigación se documenta en `docs/` (índice en
[`docs/INDEX.md`](docs/INDEX.md)). `docs/archive/rb_neutral_perturber/`
conserva la saga del pseudopotencial: técnicamente correcta, pero de **otro
sistema físico**.

## Contacto

Migración a Python: Javier Aguilera — aguilerajavier58@gmail.com
