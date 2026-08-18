# Refactor de estructura: núcleo como paquete + scripts

**Fecha**: 2026-08-18
**Autor**: Javier Aguilera
**Estado**: diseño aprobado, pendiente de plan de implementación
**Relevancia**: la estructura plana actual de `src/` bloquea la integración del
pseudopotencial de Fermi con el Hamiltoniano carga-dipolo, que es el siguiente
objetivo científico del proyecto.

## Resumen

El repositorio tiene hoy 2650 líneas en un `src/` plano donde conviven
biblioteca, ejecutables y tests, sin ser un paquete importable. Este documento
define el refactor a un **híbrido**: un paquete instalable `trimero` con la
física reutilizable, más scripts sueltos que lo importan.

El refactor se ejecuta en once pasos verificables, cada uno respaldado por tests
de caracterización generados *antes* de mover nada. Ningún paso cambia física
salvo uno, que corrige un bug de unidades y lo hace de forma explícita.

## Motivación

### Problemas estructurales

1. **`src/` no es un paquete.** Los imports son `from atom import Atom`, así que
   todo sólo funciona ejecutando desde dentro de `src/`. De ahí que los tests se
   lancen con `cd src && python test_x.py`.
2. **Tres tipos de código mezclados**: biblioteca (`atom`, `charge_dipole`,
   `quantum_basis`…), ejecutables (`main`, `run_*`, `benchmark_*`) y tests
   (`test_*`), todos en el mismo directorio.
3. **`laplacian.py` es un shim vacío**: 7 líneas que re-exportan 4 nombres de
   `math_aux.py`, herencia de la migración desde C++. `atom.py` lo consume con
   `from laplacian import *`.
4. **`main.py` ejecuta al importarse**: la llamada a
   `test_trimer_energies_field()` está a nivel de módulo, fuera del
   `if __name__ == "__main__"`.
5. **Los `.dat` de salida se escriben en la raíz** del repositorio y estaban
   trackeados en git.
6. **Los tests usan un harness manual propio**; cada archivo reimplementa su
   runner y su resumen.

### Bug de unidades (medido)

```
src/charge_dipole.py:61   HZ_PER_HARTREE = 6.579683920502e15   # Hz, CODATA 2018
src/trimer.py:13          EhtoGHz        = 6.579683920729e9    # rotulado "GHz"
src/main.py:10            EhtoGHz        = 6.579683920729e9    # rotulado "GHz"
```

1 E_h = 6.5797e15 Hz = **6.5797e6 GHz**, pero `EhtoGHz` vale `6.5797e9`, mil
veces más: ese número es el factor Hartree→**MHz**. `trimer.py:234` lo aplica al
escribir `Trimer_R_sp_wave_N35_R_300_GHz.dat`, de modo que **el fichero rotulado
GHz contiene MHz**.

Aparte, `...729` frente a `...502` son revisiones CODATA distintas (diferencia
relativa ~3e-11, numéricamente irrelevante, pero delata que nadie las comparó).

### Bug de asimetría en la matriz de campo DC (hallado en el paso 0)

Descubierto al generar el golden G3, no durante el diseño. En
`Trimer_energies_field`, la matriz `field` se construye con:

```python
for i in range(n1 - 1):        # i llega sólo hasta n1-2
    for k1 in range(1, 2 * i + 2):
        for j in range(n1):    # j llega hasta n1-1
```

El par `(i=n1-1, j=n1-2)` nunca se rellena, mientras que `(i=n1-2, j=n1-1)`
sí. El acoplamiento Stark entre las dos últimas capas `l` queda por tanto sólo
en el triángulo superior de `spV`.

Como `np.linalg.eigh` lee por defecto el triángulo **inferior** (`UPLO='L'`),
ese acoplamiento **se descarta silenciosamente**: los autovalores se calculan
como si no existiera. Medido a `n1=5, dc=0.1`: 14 elementos asimétricos,
todos confinados a los bloques `(l=3, l=4)`, con magnitud 60.85.

Sin campo DC (`dc=0`) la matriz sí es exactamente simétrica, porque `field` es
idénticamente nula y la asimetría no llega a materializarse.

El golden G3 congela el bug tal cual (`test_g3_spv_asymmetry_is_the_known_legacy_bug`),
de modo que corregirlo obligue a regenerar el golden de forma consciente.
**La corrección está propuesta como paso 4b y requiere decisión del autor**,
por ser un cambio de física y no de estructura.

### Coste computacional (medido)

```
coste por elemento de matriz : 4.05 ms

  n1= 5   dim=  25    1 fila =     2.5 s     479 filas =   0.34 h
  n1=10   dim= 100    1 fila =    40.5 s     479 filas =   5.4  h
  n1=35   dim=1225    1 fila =  6077   s     479 filas = 808.6  h   (34 días)
```

**La simulación de producción nunca se ha ejecutado completa.** El
`range(297, 307)` que apareció en `trimer.py` no era un despiste: era la única
forma de obtener resultados en tiempo finito.

El origen está localizado:

```
hydrogenicR    187.4 us     <- genlaguerre(n-l-1, 2l+1) construye un objeto
DRnl           408.9 us     <- polinomio nuevo en CADA llamada
Spherical        1.8 us
DOlm/DPhilm     16.4 us
```

Con `lru_cache` sobre `genlaguerre`, `hydrogenicR` baja de 187 µs a 14.9 µs
(**13×**) devolviendo valores **idénticos bit a bit**, verificado con `==`
exacto sobre (n=35,l=4), (n=24,l=3) y (n=27,l=0).

## Decisiones de diseño

| Decisión | Alternativas descartadas | Razón |
|---|---|---|
| Híbrido: núcleo paquete + scripts | biblioteca instalable pura; colección de scripts | la física es reutilizable, pero los escaneos exploratorios no merecen ceremonia de entry point |
| Agrupación por capas de abstracción | por prefijo plano; por sistema físico | los dos sistemas comparten ~70% del código, así que agrupar por sistema vaciaría las carpetas; las capas hacen que un import ilegal sea visible |
| Caracterizar antes de descomponer | mover sin más; descomponer directamente | el legado no tiene un solo test; descomponerlo a ciegas es cambiar física sin red |

## Sección 1 — Arquitectura de capas y mapeo de módulos

**Regla estructural: cada capa sólo importa de capas inferiores.** Nunca lateral
entre paquetes hermanos, nunca hacia arriba.

```
CAPA 0  mathlib/       matemática pura          no importa nada de trimero
CAPA 1  basis/         espacios de Hilbert      importa: mathlib
CAPA 2  systems/       objetos físicos          importa: mathlib, constants
CAPA 3  hamiltonians/  operadores               importa: 0, 1, 2
CAPA 4  io/            frontera con disco       importa: nada de física
CAPA 5  simulation/    orquestación             importa: todo
```

### Árbol objetivo

```
trimero_mod/
├── src/trimero/
│   ├── __init__.py
│   ├── constants.py
│   ├── mathlib/
│   │   ├── special.py
│   │   └── angular.py
│   ├── basis/
│   │   ├── quantum.py
│   │   └── radial.py
│   ├── systems/
│   │   └── atom.py
│   ├── hamiltonians/
│   │   ├── base.py
│   │   ├── charge_dipole.py
│   │   ├── fermi.py
│   │   └── trimer.py
│   ├── io/
│   │   ├── wavefunctions.py
│   │   └── results.py
│   └── simulation/
│       └── scan.py
├── scripts/
│   ├── run_trimer.py
│   ├── run_charge_dipole.py
│   ├── run_full_hamiltonian.py
│   └── benchmark_diagonalization.py
├── tests/
│   ├── test_architecture.py
│   ├── mathlib/    basis/    hamiltonians/
│   └── characterization/
├── data/Wavefunction/
├── output/            (generado, ignorado por git)
└── docs/
```

### Mapeo archivo por archivo

| Destino | Origen | Contenido |
|---|---|---|
| `constants.py` | nuevo | `HZ_PER_HARTREE`, factores GHz/MHz con nombre honesto, `B_KRB_AU` y `D_KRB_AU` (hoy en `run_*`) |
| `mathlib/special.py` | `math_aux.py` | `Spherical`, `DRnl`, `DOlm`, `DPhilm`, `hydrogenicR` |
| `mathlib/angular.py` | `angular_algebra.py` | `wigner_3j`, `gaunt` |
| `basis/quantum.py` | `quantum_basis.py` | `QuantumBasisBlock`, `CoupledBasis` |
| `basis/radial.py` | `rydberg_radial.py` | `RadialBasis` |
| `systems/atom.py` | `atom.py` | `Atom` |
| `hamiltonians/base.py` | nuevo | ABC `Hamiltonian` |
| `hamiltonians/charge_dipole.py` | `charge_dipole.py` | `ChargeDipoleHamiltonian`, `RydbergElectronField`, `rydberg_diagonal` |
| `hamiltonians/fermi.py` | `fermi_potentials.py` | `FermiPotential`; `Vs` y `Vp` se conservan como métodos, `Vsp` pasa a `matrix_element` en el paso 7 |
| `hamiltonians/trimer.py` | `trimer.py` (parcial) | sólo el ensamblado de H |
| `io/wavefunctions.py` | `trimer.py` (parcial) | los `np.loadtxt` de `data/Wavefunction/` |
| `io/results.py` | `trimer.py` (parcial) | escritura de `.dat` |
| `simulation/scan.py` | `trimer.py` (parcial) | el bucle sobre R |
| **eliminado** | `laplacian.py` | shim de re-export |

### Cambios que no son mover archivos

1. **`laplacian.py` desaparece.** `atom.py` (hoy `from laplacian import *`) y
   `fermi_potentials.py` pasan a importar de `mathlib.special` por nombre.
2. **`trimer.py` se parte en cuatro** destinos, según las cuatro cosas que hoy
   hace en una sola función. Ocurre en el paso 8, después de los goldens.
3. **`quantum_basis.py:217 test_basis()`** sale del módulo de producción a
   `tests/`.
4. **`main.py` deja de ejecutar al importarse**: se elimina la llamada suelta y
   `scripts/run_trimer.py` queda sólo con `if __name__ == "__main__"`.

### Cumplimiento automático

`tests/test_architecture.py` recorre los imports de cada módulo del paquete y
falla si alguno apunta a una capa superior o lateral. Sin esta comprobación la
estructura se degrada en meses; la disciplina humana no basta.

## Sección 2 — Jerarquía POO de operadores

### El problema

`ChargeDipoleHamiltonian` ya es una clase-operador correcta: se construye con
sus constantes físicas (`B`, `d`) y expone `matrix_element(bra, ket, R)` y
`build(block, R)`.

`FermiPotentials` es lo contrario, un objeto-función:

```python
def __init__(self, s, n, n2, li, lj, mi, mj, r1, theta1,
             As1, wave1, wave2, Ap1, Dwave1, Dwave2):   # 15 posicionales
```

Se instancia para calcular **un solo elemento** y se descarta. Mezcla constantes
físicas (`As1`, `Ap1`, `s`), estado cuántico del par (`li, lj, mi, mj, n, n2`) y
datos tabulados de disco (`wave1, wave2, Dwave1, Dwave2`). Cualquier permutación
de argumentos se ejecuta sin error y devuelve un número plausible.

Además, `Vs` tiene cuatro ramas (`li>2 and lj>2`, `li<=2 and lj>2`,
`li>2 and lj<=2`, `else`) y `Vp` otras cuatro. **Las ocho calculan la misma
física**; lo único que cambia es el origen de la función radial: `hydrogenicR()`
si `l>2`, o el valor tabulado `wave1_`/`wave2_` si `l<=2`. Es una decisión de
fuente de datos incrustada ocho veces dentro de las fórmulas del potencial.

### Las dos interfaces

```python
class Hamiltonian(ABC):
    """Operador hermítico sobre la base acoplada, parametrizado por R."""

    @abstractmethod
    def matrix_element(self, bra: State, ket: State, R: float) -> float: ...

    def build_reference(self, block: QuantumBasisBlock, R: float) -> np.ndarray:
        """Barrido dim² completo. Genérico: sirve a toda subclase."""

    def build(self, block: QuantumBasisBlock, R: float) -> np.ndarray:
        """Por defecto delega en build_reference. Se sobrescribe para
        explotar las reglas de selección propias del operador."""
        return self.build_reference(block, R)

    def coupling_partners(self, state: State, block) -> Iterator[State]:
        """Estados con los que `state` puede acoplarse. Vacío = todos."""
```

```python
class WavefunctionSource(ABC):
    """De dónde sale u_nl(r): hidrogénica, tabulada de disco, o Coulomb."""

    @abstractmethod
    def radial(self, n: int, l: int, r: float) -> float: ...

    @abstractmethod
    def d_radial(self, n: int, l: int, r: float) -> float: ...
```

Con `WavefunctionSource`, las ocho ramas colapsan a una expresión por potencial.
La decisión `l>2 ? hidrogénica : tabulada` se muda a `HybridSource(l_threshold=2)`,
donde es **una** condición revisable en lugar de ocho copias que pueden divergir.

`FermiPotentials` pasa a:

```python
class FermiPotential(Hamiltonian):
    def __init__(self, a_s, a_p, source: WavefunctionSource, include_p=True): ...
    def matrix_element(self, bra, ket, R): ...     # antes: Vsp()
```

De 15 posicionales a 4 con nombre; el estado cuántico viaja en los argumentos de
`matrix_element`, que es donde le corresponde.

### Composición

```python
class CompositeHamiltonian(Hamiltonian):
    def __init__(self, terms: Sequence[Hamiltonian]): ...
    def matrix_element(self, bra, ket, R):
        return sum(t.matrix_element(bra, ket, R) for t in self.terms)
```

Con lo cual el objetivo pendiente —integrar Fermi con carga-dipolo— es:

```python
H = CompositeHamiltonian([
    ChargeDipoleHamiltonian(B, d, electron_field=field),
    FermiPotential(a_s, a_p, source),
])
```

**Unión, no intersección.** El `build()` optimizado de cada término visita sólo
los pares que *sus* reglas de selección permiten. Al componer, el patrón de
dispersión del compuesto es la **unión** de los patrones: si carga-dipolo acopla
ΔN=±1 y Fermi acopla Δl arbitrario, el compuesto debe visitar ambos conjuntos.
Para eso existe `coupling_partners()`: `CompositeHamiltonian.build()` une lo que
declara cada término, y un término que devuelve vacío fuerza el barrido
completo. Intersectar daría ceros silenciosos en elementos reales — el bug que
no se detecta hasta comparar con el paper.

### Preservación de la garantía de regresión

Existe un test que exige que `ChargeDipoleHamiltonian(electron_field=None)`
reproduzca su matriz **bit a bit** (`np.array_equal`). La ABC no lo pone en
riesgo:

- `ChargeDipoleHamiltonian` **conserva su `build()` actual sin tocar una línea**;
  sólo pasa a declararse como override del de la ABC.
- `cos_theta_element` mantiene `N_max = max(N1, N2)`, expresión única para ambas
  ramas, de modo que `||H - H.T||` sigue siendo exactamente 0 en el caso
  ion-solo. No se separa en dos ramas.
- La distinción `build()` / `build_reference()` sube a la ABC en vez de
  desaparecer; los tests que comparan ambos siguen valiendo sin cambios.

El único cambio de comportamiento real de esta sección ocurre dentro de
`FermiPotentials`, que es código sin tests: por eso los goldens van antes.

## Sección 3 — Caracterización del legado y plan de migración

### Golden files en cascada

Cuatro niveles de creciente integración. Si algo se rompe, el nivel más profundo
que falla localiza el paso culpable.

| Nivel | Captura | Tolerancia |
|---|---|---|
| **G1** funciones puras | `hydrogenicR`, `DRnl`, `DOlm`, `DPhilm`, `Spherical` sobre rejilla de (n,l,m,r,θ) | bit a bit |
| **G2** elemento Fermi | `Vs()`, `Vp()`, `Vsp()` en las 4 ramas de `li/lj` | bit a bit |

G2 se captura en el paso 0 contra la API actual (`Vsp()`). A partir del paso 7,
`Vsp()` ya no existe y el golden se compara contra `matrix_element(bra, ket, R)`
sobre los mismos estados; los valores esperados no cambian, sólo la llamada que
los produce.
| **G3** matrices | `field` y `spV` a `n1=5`, una fila | bit a bit |
| **G4** salida final | autovalores de 3 filas a `n1=5`: un golden del fichero en `au` y otro del fichero rotulado `_GHz` | `rtol=1e-12` |

G1–G3 exigen igualdad exacta por ser aritmética determinista pura. G4 afloja a
`rtol=1e-12` porque `np.linalg.eigh` delega en LAPACK y el último bit puede
variar entre versiones de BLAS; exigir exactitud daría fallos espurios en otra
máquina.

**G2 es el juez de la sección 2**: es lo que demuestra que colapsar las ocho
ramas de `Vs`/`Vp` a una sola expresión sigue dando el mismo número en las cuatro
combinaciones de `li/lj`.

Los goldens se generan a `n1=5`, no a 35. Con `n1=5` el `max_dim` es 25 y las
cuatro ramas `i<3/j<3`, `i<3/j>2`, `i>2/j<3`, `i>2/j>2` se visitan todas: la
cobertura estructural no depende de `n1` grande.

Coste real medido en el paso 0: un run completo de 479 filas a `n1=5` tarda
**213 s** (no los ~20 min estimados; con `n1=5` muchos pares caen en la rama
tabulada, más barata que la hidrogénica). Los cuatro goldens ocupan 56 KB en
total, así que se versionan en git sin reparos.

G2 incluye `(mi, mj) = (0, 0)` a propósito: con `m != 0` los armónicos se
anulan en θ=0 y θ=π, que son exactamente los dos ángulos que usa `trimer.py`
(`theta=0`, `theta1=pi`). Sin el caso `m=0`, `Vs` valdría cero en dos tercios
de los casos y el golden vigilaría ceros triviales. `test_g2_covers_all_four_branches`
comprueba que cada rama aporta valores no nulos.

### Tests existentes

Los 14 tests migran a pytest **sin reescribirse**: ninguna función `test_*` lleva
parámetros y todas usan `assert`, así que pytest las descubre tal cual. El
trabajo es añadir `pytest` a las dependencias, mover los archivos a `tests/`, y
borrar el bloque `if __name__ == "__main__"` con su runner manual. La lógica no
se toca: es la red de seguridad, no puede moverse a la vez que lo que vigila.

### Orden de migración

Cada paso es un commit; su verificación debe pasar antes de empezar el
siguiente.

| # | Paso | Verificación |
|---|---|---|
| 0 | Generar goldens G1–G4 **con el código actual, sin mover nada** | los goldens se reproducen dos veces seguidas |
| 1 | Añadir `pytest`; mover los 14 tests a `tests/` | 14/14 pasan vía `pytest` |
| 2 | Crear `src/trimero/`, mover módulos, arreglar imports | 14/14 + G1–G4 |
| 3 | Eliminar `laplacian.py`; imports explícitos | 14/14 + G1–G4 |
| 4 | `constants.py` y **arreglo del bug de unidades** | G1–G3 intactos; G4-au intacto; **G4-GHz cambia ×1000** |
| 4b | **Arreglo de la asimetría de `field`** (`range(n1-1)` -> `range(n1)`); pendiente de decisión | G1–G2 intactos; **G3 y G4 cambian en los casos con dc != 0** |
| 5 | `lru_cache` en `genlaguerre` | G1–G4 bit a bit + speedup medido |
| 6 | ABC `Hamiltonian`; `ChargeDipole` sólo declara el override | 14/14, incluido el `np.array_equal` de regresión |
| 7 | `WavefunctionSource`; colapsar las 8 ramas de Fermi | **G2** |
| 8 | Partir `trimer.py` en `io/`, `hamiltonians/`, `simulation/` | G3–G4 |
| 9 | `test_architecture.py` | el nuevo test pasa |
| 10 | `CompositeHamiltonian` + integrar Fermi con carga-dipolo | objetivo científico pendiente |

**El paso 4 es el único que cambia una salida a propósito.** Al corregir
`EhtoGHz`, el fichero rotulado GHz pasa a contener GHz reales, y su golden cambia
por un factor 1000. Se actualiza conscientemente y se documenta en el commit; el
golden en `au`, que no depende del factor, sirve de control de que no se coló
ningún otro cambio.

El paso 5 va antes del 8 a propósito: con el cacheo activo, regenerar goldens en
los pasos posteriores cuesta 13× menos.

## Fuera de alcance

- **Viabilidad de `n1=35`.** Incluso con el 13× del paso 5, las 479 filas siguen
  siendo ~62 h. Bajarlo de verdad exige vectorizar sobre la rejilla de `r` y
  explotar las reglas de selección en el bucle `i,j` (hoy se recorren los 1225²
  pares aunque la mayoría sean cero por simetría). Es un proyecto posterior,
  sobre la estructura y los tests ya establecidos.
- **La función radial del estado 27s**, que usa hidrogenoide con `n_eff=24` en
  lugar de Coulomb con `n*=23.869` (~1.1% en extensión radial). Documentado en
  `docs/analysis_campo_electron_rydberg.md`.
- **Procedencia de `rvsAS.dat` / `rvsAP.dat`**: el script generador sigue sin
  localizarse. Ver `docs/analysis_procedencia_rvsAS_rvsAP.md`.
- **Comparación con las figuras del paper**, que requiere el paso 10 completo.

## Referencias

- González-Férez, Sadeghpour & Schmelcher, *New J. Phys.* **17**, 013021 (2015).
  DOI: 10.1088/1367-2630/17/1/013021
- `docs/analysis_validacion_carga_dipolo.md` — validación de la ronda del ion
- `docs/analysis_campo_electron_rydberg.md` — campo del electrón Rydberg
- `docs/analysis_procedencia_rvsAS_rvsAP.md` — procedencia de los datos de entrada
