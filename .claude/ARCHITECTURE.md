# Arquitectura del proyecto

**Última actualización**: 2026-08-20 · Python 3.13+, NumPy 2.3+, SciPy 1.16+

> Estado científico vigente: **[`docs/STATUS.md`](../docs/STATUS.md)**. Este
> documento describe el **código**; `STATUS.md` describe la **física**.

## El punto de partida: son DOS sistemas físicos

Lo más importante que hay que saber antes de tocar nada. El repositorio cubre
dos problemas distintos que durante varias rondas estuvieron mezclados en el
mismo espacio de nombres, hasta que la mezcla produjo resultados con premisa
equivocada. Desde la reorganización del 2026-08-20 el árbol de paquetes los
separa explícitamente.

| | `rb_krb_polar` | `rb_neutral_perturber` |
|---|---|---|
| perturbador | KRb, **polar** | átomo/molécula **neutra** |
| interacción | carga-dipolo, `−d·F_ryd` | dispersión de contacto (pseudopotencial de Fermi) |
| Hamiltoniano | `H_ad = H_A + H_mol` | `H = H_a + V_Fermi` |
| referencia | Aguilera-Fernández 2015 / González-Férez 2015 | Aguilera-Fernández 2016 |
| estado | **vigente** | verde y protegido, línea en pausa |

**El pseudopotencial de Fermi no interviene en el sistema polar.** Si te ves
añadiendo `V_Fermi` a una curva de Rb\*-KRb, para y lee `docs/STATUS.md`.

## Capas

```
src/trimero/
├── mathlib/          primitivas matemáticas, sin contexto físico
│   ├── angular.py        wigner_3j, gaunt            (lo usa el lado polar)
│   ├── special.py        Spherical, DRnl, DOlm, DPhilm, hydrogenicR  (lado legado)
│   └── laplacian.py      reexporta special.py
├── basis/            enumeración de la base, compartida
│   ├── quantum.py        CoupledBasis, QuantumBasisBlock — |l m_l> ⊗ |N M_N>,
│   │                     bloqueada por M_J
│   └── radial.py         RadialBasis — tablas radiales hidrogenoides
├── simulation/
│   └── bop_tracking.py   trace_curve — rastreo de curvas, agnóstico del sistema
└── systems/
    ├── rb_atom.py        Atom.E_Rb() — defectos cuánticos de Rb. COMPARTIDO:
    │                     es la fuente de verdad de los dos sistemas
    ├── rb_krb_polar/
    │   ├── charge_dipole.py   RydbergElectronField, ChargeDipoleHamiltonian,
    │   │                      rydberg_diagonal, constantes de KRb
    │   ├── bop_system.py      BOPSystem — monta base + H para un manifold dado
    │   └── rb_defects.py      DELTA0_NS_PAPER, neighbor_levels(), n_star_nl()
    └── rb_neutral_perturber/
        ├── fermi_krb.py         ScatteringLengths, FermiPseudopotential
        │                        (vectorizado, sobre CoupledBasis)
        ├── linear_trimer.py     SymmetricLinearTrimer — H = H₀ + F·r + V₁ + V₂
        │                        del trímero lineal simétrico. CAPA MODERNA
        ├── fermi_potentials.py  FermiPotentials  ── legado, golden files
        └── trimer.py            Trimer_energies_field  ── legado, golden files
```

Regla de dependencia: `mathlib` → `basis` → `systems`. Nada de `mathlib` o
`basis` debería conocer un sistema concreto (ver **Deuda** más abajo: hoy hay
dos excepciones).

### Grafo real de dependencias

```
basis/radial.py                       → systems.rb_krb_polar.rb_defects        ⚠️
mathlib/laplacian.py                  → mathlib.special
systems/rb_atom.py                    → mathlib.{laplacian,special}
rb_krb_polar/bop_system.py            → basis.{quantum,radial}
                                        rb_krb_polar.{charge_dipole,rb_defects}
                                        rb_neutral_perturber.fermi_krb          ⚠️
rb_krb_polar/charge_dipole.py         → basis.{quantum,radial} · mathlib.angular
                                        systems.rb_atom · rb_krb_polar.rb_defects
rb_krb_polar/rb_defects.py            → systems.rb_atom
rb_neutral_perturber/fermi_krb.py     → basis.{quantum,radial}
                                        rb_krb_polar.rb_defects                 ⚠️
rb_neutral_perturber/linear_trimer.py → basis.{quantum,radial} · systems.rb_atom
                                        rb_krb_polar.rb_defects                 ⚠️
                                        rb_neutral_perturber.fermi_krb
rb_neutral_perturber/fermi_potentials → mathlib.{laplacian,special}
rb_neutral_perturber/trimer.py        → systems.rb_atom
                                        rb_neutral_perturber.fermi_potentials
```

Las cuatro ⚠️ son la MISMA deuda (`rb_defects` mezcla física atómica de Rb
con la composición de la base del paper polar), conocida y documentada, no
descuidos. Ver abajo.

## El camino vigente: Rb\*-KRb

### `BOPSystem` es la pieza central

```python
from trimero.systems.rb_krb_polar.bop_system import BOPSystem
from trimero.systems.rb_krb_polar.rb_defects import DELTA0_NS_PAPER

sysm = BOPSystem(n_manifold=25, delta0_ns=DELTA0_NS_PAPER)
H = sysm.hamiltonian(R=800.0, M_J=0, fermi=False)   # 1113 × 1113
```

Construir cuesta unos segundos (tablas radiales); **reutiliza la instancia para
todo un barrido en R**. Todo lo que depende del manifold sale de un único
parámetro `n_manifold`:

| qué | de dónde sale |
|---|---|
| l máximo del manifold | `n_manifold − 1` → `CoupledBasis` |
| niveles vecinos (n+1)d, (n+2)p, (n+3)s | `rb_defects.neighbor_levels(n)` |
| funciones radiales | `RadialBasis` |
| energías diagonales de H_A | `rydberg_diagonal(n_manifold)` |
| cero de energía | `E_manifold = −0.5/n²`, exacto para l ≥ 3 |

**`fermi=False` no significa «apagar un término de este sistema»**: significa
que el pseudopotencial no forma parte del modelo polar. Que el flag exista es
deuda, no diseño.

### Identificación de la curva: por CARÁCTER, no por índice

`character_curve()` devuelve la curva adiabática **más baja con peso de manifold
> 50 %**, evaluado en cada R. Un índice fijo identificado en un extremo no vale:
los estados de (n+1)d y (n+2)p producen cruces evitados y el índice cambia de
objeto por el camino.

### Punto de entrada

```bash
poetry run python scripts/compute_bop_curve.py --n-manifold 25 --mj 0 1
```

**Único script de producción.** Los doce `run_*.py` / `analyze_*.py` de las
rondas de exploración están en `scripts/archive/`, con sus imports arreglados
pero sin garantía de que su física siga siendo la vigente.

## El camino legado: perturbador neutro

`Trimer_energies_field(n1, dc_field_au)` es la traducción directa del C++
original: lee `data/Wavefunction/*.dat`, construye la matriz `n1²×n1²` con la
lógica de casos A/B/C/D según el `l` de cada índice, diagonaliza por cada fila
de R y escribe `Trimer_R_sp_wave_*.dat`.

**Está congelado y protegido por golden files bit a bit** (`rtol=1e-12`), en
`tests/systems/rb_neutral_perturber/characterization/`. La cadena es g1
(funciones especiales) → g2 (`FermiPotentials`) → g3 (matrices) → g4
(autovalores end-to-end). Los tests documentan explícitamente dos bugs del
legado que **se conservan a propósito** para que el golden siga siendo fiel:
un factor de unidades en `EhtoGHz` y una asimetría de la matriz.

⚠️ **No refactorices este camino «de paso».** Cualquier cambio que altere un
golden es un cambio de física, y hay que tratarlo como tal.

**La capa moderna equivalente es `linear_trimer.py`**, que resuelve el mismo
sistema físico sobre `CoupledBasis(N_max=0)` y sin tres bugs de datos que el
legado arrastra (`rvsDR38s.dat` × 10⁻⁶, desalineación de mallas 780 vs 776
filas, `EhtoGHz` × 10³). Los tres están documentados en
`docs/analysis_trimero_lineal_campo_dc.md` §7. Física nueva va ahí, no aquí.

## Tests

```
tests/
├── conftest.py                    pone src/ en el path
├── basis/                         enumeración de la base (compartida)
└── systems/
    ├── rb_krb_polar/              charge_dipole, rydberg_field, bop_system,
    │                              regression_fig1
    └── rb_neutral_perturber/      fermi_krb, linear_trimer (T1-T9 analíticos),
                                   regression_trimer_2016 (anclas del paper),
                                   characterization/ (goldens del legado)
```

**77 tests.** `poetry run pytest` completo tarda ~6 min (los `slow` son los
goldens end-to-end del legado). Para iterar: `pytest -m "not slow"`, ~25 s.

`tests/systems/rb_krb_polar/test_regression_fig1.py` ancla los números
verificados de la Fig. 1 (n=25, M_J=0): profundidad −23.100 GHz, E(1800 a₀) =
−0.338 GHz, 8 mínimos locales, contra `plots/rb_krb_polar/fig1_ad_MJ0_n25.npz`.

## Datos y artefactos

| ruta | qué |
|---|---|
| `data/Wavefunction/` | `.dat` de entrada. `rvsAS`/`rvsAP` son longitudes de dispersión: **insumo exclusivo del perturbador neutro** |
| `plots/rb_krb_polar/` | lo vigente del sistema polar: `fig1_ad_*` y sus `.npz` |
| `plots/rb_neutral_perturber/` | lo vigente del neutro: `trimer_lineal_*` y sus `.npz` |
| `plots/archive/<sistema>/` | figuras y datos de rondas superadas, divididos por sistema |
| `docs/` | los 5 documentos de referencia activa + `STATUS.md` |
| `docs/archive/rb_neutral_perturber/` | la saga del pseudopotencial: correcta, pero de otro sistema |

⚠️ `plots/` **va siempre al commit**, nunca al `.gitignore`. El árbol de
`plots/` **espeja `src/trimero/systems/`**: un subdirectorio por sistema físico.
`graphify-out/` NO se commitea (ver CLAUDE.md §7).

## Invariantes

1. **Una sola definición de la base.** La composición manifold + (n+1)d +
   (n+2)p + (n+3)s vive en `rb_defects.neighbor_levels()`. No se escribe a mano
   en ningún otro sitio.
2. **`Atom.E_Rb()` es la fuente de verdad** de los defectos cuánticos de Rb.
   `DELTA0_NS_PAPER` es un override local para comparar con un paper concreto;
   `delta0_ns=None` delega exactamente en `Atom` sin cambiar nada.
3. **Unidades**: interno en unidades atómicas (a₀, E_h). La conversión a GHz es
   explícita y sólo en la frontera de salida (`GHZ_PER_HARTREE`).
4. **Determinismo**: mismo input → mismo output. Los goldens dependen de ello.
5. **Los golden files son ley.** Si un cambio los mueve, es un cambio de física.

## Deuda técnica

1. **`BOPSystem` sigue acoplado a `fermi_krb`.** `__post_init__` construye
   siempre un `FermiPseudopotential` (lee `rvsAS.dat`/`rvsAP.dat`) y
   `hamiltonian()` tiene `fermi=True` por defecto, aunque el modelo polar no lo
   use. Desacoplarlo **cambia números** y hay tres tests que dependen del
   comportamiento actual, así que se dejó para una ronda propia.
2. **`rb_defects.py` mezcla dos cosas**: física atómica de Rb (compartida) y
   composición de la base del paper polar. Por eso `basis/radial.py` y
   `rb_neutral_perturber/fermi_krb.py` importan de `rb_krb_polar/`. Separarlas
   cerraría las dos inversiones de capa restantes.

## Rendimiento

| operación | complejidad | nota |
|---|---|---|
| montar `BOPSystem` | — | ~segundos, **una sola vez por barrido** |
| construir H(R) | O(dim²) | dim = 1113 para n=25, M_J=0 |
| diagonalizar | O(dim³) | ~1.1 s/punto — el cuello de botella |
| barrido en R | O(n_R × dim³) | 281 puntos ≈ 5 min. Paralelizable: cada R es independiente |

`eigvalsh` es ~2× más rápido que `eigh`, pero el criterio de carácter necesita
los autovectores, así que la curva BOP no puede usarlo.

## Puntos de extensión

- **Otro manifold**: `--n-manifold`. No hay nada cableado a n=24 ni n=25.
- **Otro término en H_mol**: `ChargeDipoleHamiltonian` en `charge_dipole.py`.
- **Perturbador neutro, geometría lineal simétrica**: hecho y validado contra
  Aguilera-Fernández 2016. `scripts/compute_trimer_curves.py`. Ver
  `docs/analysis_trimero_lineal_campo_dc.md`.
- **Otras geometrías del trímero** (asimétrica §III.B, planar §IV del paper):
  toda la geometría está encapsulada en `SymmetricLinearTrimer.parity_factor`.
  Extenderla es sustituir ese factor por los armónicos esféricos evaluados en
  cada θᵢ, con `A_s`/`A_p` propios de cada `R_i` (ya no comparten valor).
  El resto —base, campo, cero de energía— no cambia.
- **`compute_bop_curve.py --system neutral`** sigue lanzando
  `NotImplementedError`: ese camino era la curva BOP tipo Rb*-KRb, no el
  trímero. El barrido del trímero es `compute_trimer_curves.py`.
