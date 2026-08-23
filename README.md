# trimero_mod

Herramientas numéricas para la estructura electrónica y la dinámica nuclear de
moléculas Rydberg de largo alcance basadas en Rb. El proyecto construye
Hamiltonianos en bases angulares acopladas, calcula curvas de potencial
adiabático (BOP), orientación molecular y observables no adiabáticos.

La física vigente está resumida en [`docs/STATUS.md`](docs/STATUS.md) y la
estructura interna en [`.claude/ARCHITECTURE.md`](.claude/ARCHITECTURE.md).

## Sistemas físicos

| Sistema | Paquete | Interacción principal | Script de producción |
|---|---|---|---|
| Rb\*+KRb | `polar_rydberg` + `KRB` | carga–dipolo | `compute_bop_curve.py` |
| Rb\*+RbCs | `polar_rydberg` + `RBCS` | carga–dipolo | `compute_bop_curve.py` |
| Rb\*+2Rb neutros | `rb_neutral_perturber` | Fermi s+p | `compute_trimer_curves.py` |
| Rb\*+Rb+RbCs | `hybrid_neutral_polar` | carga–dipolo + Fermi | `compute_hybrid_curves.py` |
| Dinámica nuclear | `nonadiabatic_dynamics` | canales acoplados | `analyze_nonadiabatic_*.py` |

El sistema polar puro nunca contiene pseudopotencial de Fermi. `BOPSystem` se
mantiene como compatibilidad histórica; la API de producción es
`PolarBOPSystem`.

## Instalación y validación

Requiere Python 3.13 o posterior y Poetry:

```bash
poetry install
poetry run pytest -m "not slow"
```

Antes de integrar cambios de física:

```bash
poetry run pytest
```

La suite actual contiene 156 tests. Los más lentos protegen mediante goldens
el camino migrado desde C++.

## Ejecución por sistema

Todos los comandos se ejecutan desde la raíz del repositorio.

### 1. Rb\*+KRb polar

```bash
# Curvas BOP, bloques M_J=0,1
poetry run python scripts/compute_bop_curve.py \
  --molecule krb --n-manifold 25 --n-max 6 --mj 0 1

# Orientación del dipolo
poetry run python scripts/compute_orientation_curve.py \
  --molecule krb --n-manifold 25 --n-max 6 --mj 0

# Comparaciones entre manifolds ya calculados
poetry run python scripts/compare_bop_curves_n.py \
  --molecule krb --n 24 25 26 27
poetry run python scripts/compare_orientation_n.py \
  --molecule krb --n 24 25 26 27
```

Salidas: `plots/rb_krb_polar/data/` y `plots/rb_krb_polar/figures/`.

### 2. Rb\*+RbCs polar

Usa `B=490.17 MHz`, `d=1.225 D`. Para `n=25`, `N_max=6` está convergido
frente a `N_max=8` dentro de 0.0023 GHz en los puntos de control.

```bash
# Convergencia rotacional
poetry run python scripts/check_polar_convergence.py \
  --molecule rbcs --n-manifold 25 --n-max 4 6 8 \
  --r 500 1000 1500

# Curvas completas
poetry run python scripts/compute_bop_curve.py \
  --molecule rbcs --n-manifold 25 --n-max 6 --mj 0 1 \
  --ymin -60

# Orientación
poetry run python scripts/compute_orientation_curve.py \
  --molecule rbcs --n-manifold 25 --n-max 6 --mj 0

# Figuras únicas de comparación para n=24,...,30 (requieren los NPZ previos)
poetry run python scripts/compare_bop_curves_n.py \
  --molecule rbcs --n 24 25 26 27 28 29 30 --mj 0
poetry run python scripts/compare_orientation_n.py \
  --molecule rbcs --n 24 25 26 27 28 29 30 --mj 0 --rmax-plot 800

# Figura de orientación y alineamiento tipo Fig. 3, comparando n=25 y n=29
# (los NPZ deben haberse calculado con la versión que incluye COS2)
poetry run python scripts/plot_orientation_alignment.py \
  --molecule rbcs --n 25 29 --mj 0 --n-max 6 --rmin 400 --rmax 1800

# Potenciales con campo DC paralelo a Z, un panel independiente por campo
poetry run python scripts/compute_field_curves.py \
  --molecule rbcs --n-manifold 25 --n-max 6 --mj 0 \
  --fields 0 100 300 500 --rmin 400 --rmax 1800 --step 5
```

Salidas: `plots/rb_rbcs_polar/data/` y `plots/rb_rbcs_polar/figures/`.
Resultados actuales: [`docs/analysis_rb_rbcs_curvas_n25.md`](docs/analysis_rb_rbcs_curvas_n25.md).

### 3. Perturbadores neutros: Rb–Rb\*–Rb

```bash
# Σ y campos 0, 100, 300 y 500 V/m
poetry run python scripts/compute_trimer_curves.py

# Π
poetry run python scripts/compute_trimer_curves.py --symmetry Pi

# Sólo onda s, dímero y campos elegidos
poetry run python scripts/compute_trimer_curves.py \
  --symmetry Sigma --fields 0 500 --s-wave-only --dimer
```

Opciones: `--n-manifold`, `--symmetry`, `--fields`, `--rmin`, `--rmax`,
`--s-wave-only`, `--dimer` y `--no-plot`.

Salidas: `plots/rb_neutral_perturber/data/` y
`plots/rb_neutral_perturber/figures/`.

El camino directo `trimer.py` es legado congelado. La física nueva usa
`linear_trimer.py` mediante `compute_trimer_curves.py`.

### 4. Híbrido Rb\*+Rb+RbCs

El Rb neutro está en `theta=pi` a distancia `R1`; RbCs está en `theta=0` a
distancia variable `R2`:

```bash
poetry run python scripts/compute_hybrid_curves.py \
  --n-manifold 35 --n-max 6 --mj 0 \
  --r1 600 900 1100 --rmin 500 --rmax 1500 --step 25
```

La curva se inicia desde la referencia polar sin Fermi y se sigue por máximo
solapamiento. Se guardan peso de manifold, pesos de 38s/37p/36d,
solapamientos y espectro auxiliar.

Salidas: `plots/hybrid_neutral_polar/data/` y
`plots/hybrid_neutral_polar/figures/`.

### 5. Dinámica no adiabática

El paquete implementa acoplamiento de derivada, normalización energética,
canales acoplados, estabilización, tasas de decaimiento y factores de
Franck–Condon. Los scripts actuales reproducen el caso de demostración `n=25`,
pero consumen un barrido electrónico previo con arrays `R`, `W` y `V`. Ese
barrido no está versionado; su ruta es el primer argumento posicional:

```bash
poetry run python scripts/analyze_nonadiabatic_n25_crossing.py \
  /ruta/fase6_sweep_n25_MJ0_states54_55.npz \
  plots/hybrid_neutral_polar/data/fase6_n25_results.npz

poetry run python scripts/analyze_nonadiabatic_n25_crossing_wide.py \
  /ruta/fase6b_sweep_n25_MJ0_states54_55.npz \
  plots/hybrid_neutral_polar/data/fase6b_n25_results.npz
```

Salidas por defecto:

```text
plots/hybrid_neutral_polar/data/fase6_n25_results.npz
plots/hybrid_neutral_polar/data/fase6b_n25_results.npz
```

No se debe aplicar este pipeline a una curva nueva antes de validar su
convergencia, referencia energética y seguimiento. Véase
[`docs/PLAN_nonadiabatic_dynamics.md`](docs/PLAN_nonadiabatic_dynamics.md).

### 6. Diagramas de geometría

```bash
poetry run python scripts/draw_geometry.py \
  --hybrid 600 900 \
  -o plots/geometry/hybrid_neutral_polar_esquema.png \
  --formats png pdf
```

Para cuerpos arbitrarios consulta `poetry run python scripts/draw_geometry.py --help`.

## Opciones comunes de los scripts polares

| Opción | Significado |
|---|---|
| `--molecule krb|rbcs` | especie molecular |
| `--n-manifold` | manifold cuasi-degenerado de Rb |
| `--n-max` | corte rotacional |
| `--mj` | bloque o bloques de proyección total |
| `--rmin`, `--rmax`, `--step` | malla radial en `a0` |
| `--weight` | peso mínimo de manifold |
| `--reuse` | reutiliza un `.npz` compatible (`compute_bop_curve.py`) |
| `--no-plot` | omite la figura (`compute_bop_curve.py`) |
| `--out` | cambia la ruta de la figura |

Los `.npz` polares contienen `R`, `E`, `K`, `W`, `spectrum` y metadatos de
especie, constantes moleculares, base y criterio de carácter.

## Organización del código

```text
src/trimero/
├── mathlib/                    primitivas angulares y especiales
├── basis/                      base cuántica y radiales
├── simulation/                 seguimiento de curvas
├── visualization/              diagramas geométricos
└── systems/
    ├── polar_molecule.py       catálogo KRb/RbCs
    ├── polar_rydberg/          motor polar genérico
    ├── rb_krb_polar/           operador carga–dipolo + legado compatible
    ├── rb_rbcs_polar/          configuración pública RbCs
    ├── rb_neutral_perturber/   Fermi moderno y legado
    ├── hybrid_neutral_polar/   sistema híbrido
    └── nonadiabatic_dynamics/  dinámica nuclear
```

`plots/` se organiza por sistema y luego por `data/` y `figures/`; consulta
[`plots/README.md`](plots/README.md). `archive/` contiene resultados históricos.

## Reglas científicas

- No añadir Fermi al sistema polar puro.
- Seleccionar curvas por carácter; no confiar en un índice fijo a través de
  cruces evitados.
- Verificar `N_max` para cada especie y manifold.
- No regenerar un golden para ocultar una regresión.
- Las unidades internas son atómicas; las conversiones se hacen en las
  fronteras de entrada y salida.

## Documentación

- [`docs/INDEX.md`](docs/INDEX.md): índice general.
- [`docs/STATUS.md`](docs/STATUS.md): estado científico vigente.
- [`docs/PLAN_rb_rbcs_polar.md`](docs/PLAN_rb_rbcs_polar.md): sistema RbCs.
- [`.claude/ARCHITECTURE.md`](.claude/ARCHITECTURE.md): arquitectura detallada.

## Contacto

Migración a Python: Javier Aguilera — aguilerajavier58@gmail.com
