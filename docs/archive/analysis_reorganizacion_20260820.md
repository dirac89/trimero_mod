# Reorganización del repositorio por sistema físico

**Fecha**: 2026-08-20
**Autor**: Javier Aguilera
**Relevancia**: Registro de qué se movió y por qué. **No hay ningún cambio de
física en esta ronda.** Los 46 tests preexistentes pasan con exactamente los
mismos nombres y los mismos números antes y después.
**Tipo**: analysis

## Resumen

`docs/analysis_fig1_carga_dipolo_sin_fermi.md` estableció que el repositorio
contenía **dos sistemas físicos distintos** mezclados en el mismo espacio de
nombres:

- **Rb\*-KRb polar** — `H_ad = H_A + H_mol`, KRb como dipolo puntual en el campo
  del Rydberg (Aguilera-Fernández 2015 / González-Férez 2015). **Vigente.**
- **Perturbador neutro** — dispersión de contacto, pseudopotencial de Fermi s+p
  (Aguilera-Fernández 2016). Correcto, pero **otro problema**.

Mientras estuvieron mezclados, la maquinaria del segundo (remapeo `k(R)`,
ventana de exclusión de la resonancia de onda p, tope de dominio en `2n²a₀`) se
aplicaba por defecto al primero, donde no pinta nada. Esta reorganización separa
los dos en el árbol de paquetes, en los tests y en `docs/`, sin tocar una sola
fórmula.

## Palabras clave

- separación por sistema físico, `rb_krb_polar` / `rb_neutral_perturber`
- reorganización sin cambio de física
- `git mv`, historial preservado

## 1. Estructura resultante

```
src/trimero/
├── basis/          quantum.py (CoupledBasis), radial.py (RadialBasis)   [compartido]
├── mathlib/        angular.py (3j/Gaunt), special.py, laplacian.py      [compartido]
├── simulation/     bop_tracking.py (trace_curve, agnóstico del sistema) [compartido]
└── systems/
    ├── rb_atom.py                    Atom.E_Rb(): defectos cuánticos de Rb  [compartido]
    ├── rb_krb_polar/                 H_ad = H_A + H_mol                     [VIGENTE]
    │   ├── charge_dipole.py
    │   ├── bop_system.py
    │   └── rb_defects.py
    └── rb_neutral_perturber/         pseudopotencial de Fermi               [en pausa]
        ├── fermi_krb.py              vectorizado, sobre CoupledBasis
        ├── atom.py → (ver §2)
        ├── fermi_potentials.py       legado, protegido por golden files
        └── trimer.py                 legado, protegido por golden files
```

## 2. Tabla de movimientos

### 2.1 `src/`

| origen | destino | por qué |
|---|---|---|
| `hamiltonians/charge_dipole.py` | `systems/rb_krb_polar/charge_dipole.py` | `-d·F_ryd`: es el término característico del sistema polar |
| `simulation/bop_system.py` | `systems/rb_krb_polar/bop_system.py` | monta la base y H del sistema polar; no es genérico |
| `systems/rb_defects.py` | `systems/rb_krb_polar/rb_defects.py` | `DELTA0_NS_PAPER` y `neighbor_levels()` son de la comparación con el paper polar |
| `hamiltonians/fermi_krb.py` | `systems/rb_neutral_perturber/fermi_krb.py` | pseudopotencial de contacto = perturbador neutro |
| `hamiltonians/fermi.py` | `systems/rb_neutral_perturber/fermi_potentials.py` | ídem, versión legado. Renombrado a su nombre real (`FermiPotentials`) |
| `hamiltonians/trimer.py` | `systems/rb_neutral_perturber/trimer.py` | el trímero legado son dos perturbadores neutros alrededor del Rydberg |
| `systems/atom.py` | **`systems/rb_atom.py`** | ver más abajo |
| `hamiltonians/`, `io/` | *(eliminados)* | `hamiltonians/` quedó vacío; `io/` nunca tuvo contenido |

**`atom.py` no bajó a `rb_neutral_perturber/`**, aunque en el plan inicial iba
ahí como «código legado». `Atom.E_Rb()` es la fuente de verdad de los defectos
cuánticos de Rb y lo usan **los dos** sistemas (`rb_defects.py` y
`charge_dipole.py`, ambos del lado polar). Anidarlo bajo un subsistema habría
creado dos aristas cruzadas polar→neutro artificiales. Queda en `systems/` como
capa compartida, que es lo que físicamente es: el átomo de Rb.

### 2.2 `tests/` — mismo contenido, sólo ubicación

| origen | destino |
|---|---|
| `tests/hamiltonians/test_charge_dipole.py` | `tests/systems/rb_krb_polar/` |
| `tests/hamiltonians/test_rydberg_field.py` | `tests/systems/rb_krb_polar/` |
| `tests/simulation/test_bop_system.py` | `tests/systems/rb_krb_polar/` |
| `tests/hamiltonians/test_fermi_krb.py` | `tests/systems/rb_neutral_perturber/` |
| `tests/characterization/` (4 tests + `goldens/*.npz` + 2 generadores) | `tests/systems/rb_neutral_perturber/characterization/` |
| `tests/basis/test_basis_enumeration.py` | *(sin mover: la base es compartida)* |

`tests/mathlib/` **no se creó**: el único test que ejerce `mathlib/special.py`
es `test_g1_special.py`, y es un golden del legado que forma cadena con
g2→g3→g4. Separarlo habría roto una suite coherente para ganar una carpeta.

### 2.3 `docs/`

Se quedan en `docs/` los cinco de referencia activa, enlazados desde el nuevo
`docs/STATUS.md`: `analysis_fig1_carga_dipolo_sin_fermi.md`,
`analysis_base_correcta_3_vecinos.md`, `analysis_verificacion_tabla_I.md`,
`analysis_validacion_carga_dipolo.md`, `analysis_campo_electron_rydberg.md`.

A `docs/archive/rb_neutral_perturber/` van los siete de la saga del
pseudopotencial. **Su contenido técnico sigue siendo correcto**; lo que no
aplica es su uso en la comparación con el paper polar. Cada uno lleva ya su
aviso de premisa al principio:

| documento | qué sigue valiendo |
|---|---|
| `analysis_pseudopotencial_fermi_krb.md` | derivación de V_s/V_p y el remapeo k(R) |
| `analysis_resonancia_onda_p.md` | la resonancia ³Pᵒ a 24.8 meV, bien localizada |
| `analysis_interpolacion_polo_Ap.md` | la interpolación de 1/A_p a través del polo |
| `analysis_ventana_exclusion_resonancia.md` | el criterio de ventana (polo ± 2×FWHM) |
| `analysis_extension_dominio_fig1.md` | la medida de que V_Fermi no es despreciable en el borde |
| `analysis_curva_bop_MJ0.md` | la metodología adiabático/diabático |
| `analysis_procedencia_rvsAS_rvsAP.md` | procedencia de `rvsAS.dat`/`rvsAP.dat`, insumo **exclusivo** del perturbador neutro |

Se dejan donde estaban, sin clasificar, `AUDIT_MIGRACIÓN_PYTHON.md`,
`DISEÑO_HAMILTONIANO_KRBC.md` y `_v2.md`: no son ni referencia activa ni
material del perturbador neutro, así que ninguna de las dos reglas los alcanza.

### 2.4 `scripts/` y `plots/`

Doce `run_*.py` / `analyze_*.py` de las rondas de exploración → `scripts/archive/`,
con sus imports y sus rutas a `.npz` actualizados para que sigan ejecutándose.
Los sustituye un único script de producción, `scripts/compute_bop_curve.py`.

`plots/` conserva el PNG vigente `fig1_ad_MJ0_MJ1_n25.png` **y los dos `.npz`
que lo respaldan** (`fig1_ad_MJ0_n25.npz`, `fig1_ad_MJ1_n25.npz`): son los datos
verificados contra los que se comprueba la regresión, y archivarlos dejaría el
criterio de verificación fuera de sitio. Los otros cinco PNG, cinco `.npz` y
`plot_energies.py` → `plots/archive/`.

## 3. Las tres trampas del movimiento

Ninguna es de física, pero las tres habrían roto la suite en silencio:

1. **`fermi_krb.py` bajó un nivel y perdió `data/`.** Usaba
   `Path(__file__).resolve().parents[3] / "data" / "Wavefunction"`. En
   `systems/rb_neutral_perturber/` hace falta `parents[4]`. Sin el arreglo,
   `ScatteringLengths` revienta al construirse, y con ella
   `BOPSystem.__post_init__` — es decir, **toda la suite del lado polar**, que
   ni siquiera usa esas tablas.
2. **Los goldens también bajaron dos niveles.** `REPO = parents[2]` → `parents[4]`
   en `generate_goldens.py`, `generate_goldens_g3.py` y `test_g4_eigenvalues.py`.
3. **`hamiltonians.fermi` es prefijo de `hamiltonians.fermi_krb`.** Un `sed`
   ingenuo convierte `fermi_krb` en `fermi_potentials_krb`. Hubo que sustituir
   `fermi_krb` primero y anclar el patrón de `fermi` con `[^_]`.

## 4. Lo que NO se tocó, y por qué

**`BOPSystem` sigue acoplado a `fermi_krb`.** Su `__post_init__` construye
siempre un `FermiPseudopotential` y `hamiltonian()` tiene `fermi=True` por
defecto, aunque el modelo polar no lo use. Es la última arista
`rb_krb_polar → rb_neutral_perturber` y es exactamente la mezcla que esta ronda
quería deshacer — pero **desacoplarla cambia números**, y hay tres tests
(`test_s2_dominio_lo_fija_el_np`, `test_s4_resonancia_misma_energia_distinto_R`,
`test_s5_fermi_off_es_cero_exacto`) que dependen del comportamiento actual.
Queda anotada como deuda en `docs/STATUS.md`. El lado polar se protege pasando
`fermi=False` explícitamente.

Lo mismo con la partición de `rb_defects.py` en su mitad atómica (compartida) y
su mitad de composición de base (polar), que cerraría las dos inversiones de
capa restantes.

## 5. Verificación

- **`git mv` en todos los movimientos**, no copiar+borrar: el historial de cada
  fichero se conserva. Los ficheros que aparecen como añadidos y no como
  renombrados eran trabajo de la ronda anterior aún sin commitear, sin historial
  previo que preservar.
- **`pytest` completo (incluidos `slow` y golden files) después de cada
  movimiento significativo**, no sólo al final: pasos 2, 3, 4 y verificación
  final. En los cuatro cortes, los 46 tests preexistentes pasaron **con nombres
  idénticos al baseline**, comprobado con `diff` de la lista nominal.
- **Reproducción literal de los números verificados** con el script de
  producción nuevo: ver `docs/STATUS.md` y
  `tests/systems/rb_krb_polar/test_regression_fig1.py`.

## Conclusiones y aplicación al proyecto

El árbol de paquetes ahora dice de qué sistema físico es cada cosa, que es
justo la información cuya ausencia permitió aplicar el pseudopotencial de Fermi
al sistema equivocado durante varias rondas. `docs/STATUS.md` es el punto de
entrada; `docs/archive/rb_neutral_perturber/` conserva intacto el trabajo del
otro sistema por si se retoma esa línea.

## Referencias

- `docs/STATUS.md` — estado vigente tras la reorganización.
- `docs/analysis_fig1_carga_dipolo_sin_fermi.md` — la corrección de premisa que
  motivó esta separación.
- `docs/superpowers/specs/2026-08-18-refactor-estructura-design.md` — el refactor
  anterior, por capas; éste es el siguiente, por sistema físico.
