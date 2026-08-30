# Arquitectura del proyecto

**Última actualización:** 2026-08-23

**Entorno:** Python 3.13+, NumPy 2.3+, SciPy 1.16+

Este documento describe el código. El estado científico vigente está en
[`docs/STATUS.md`](../docs/STATUS.md); los comandos de usuario están en
[`README.md`](../README.md).

## Mapa de sistemas

El repositorio ya no contiene solamente “polar” y “neutro”. Hay un núcleo
polar compartido, dos configuraciones moleculares, un sistema neutro, su
composición híbrida y una capa posterior de dinámica nuclear.

```text
                         polar_molecule.{KRB,RBCS}
                                   │
                                   ▼
mathlib ──► basis ──► rb_krb_polar/charge_dipole
                    │              │
                    │              ▼
                    │       polar_rydberg/PolarBOPSystem
                    │          │                    │
                    │          ▼                    ▼
                    │     Rb*+KRb              Rb*+RbCs
                    │                               ├──► Rb*+2RbCs
                    │                               │
                    ▼                               │
          rb_neutral_perturber                      │
                    │                               │
                    └────────► hybrid_neutral_polar ◄┘
                                      │
                                      ▼
                         nonadiabatic_dynamics
```

### Hamiltonianos

| Sistema | Hamiltoniano |
|---|---|
| polar puro | `H_A + B N² - d·(F_ion + F_elec)` |
| polar doble | `H_A + Σ_i[B N_i² - d_i·F_ryd(R_i)] + V_dd + H_F` |
| neutro | `H_A + F_ext·r + V_Fermi` |
| híbrido | `H_A + H_mol(R2) + V_Fermi^pi(R1)` |
| dinámica | ecuación nuclear sobre BOP y acoplamientos derivados |

El pseudopotencial de Fermi no pertenece al sistema polar puro.

## Capas y responsabilidades

```text
src/trimero/
├── mathlib/
│   ├── angular.py               wigner_3j y Gaunt
│   ├── special.py               funciones especiales del legado
│   └── laplacian.py             compatibilidad con special.py
├── basis/
│   ├── quantum.py               CoupledBasis y bloques M_J
│   └── radial.py                radiales e integrales G^k/Z^k
├── simulation/
│   └── bop_tracking.py          seguimiento agnóstico por solapamiento
├── visualization/
│   └── geometry_diagram.py      diagramas; no construye Hamiltonianos
└── systems/
    ├── rb_atom.py               energías y defectos cuánticos de Rb
    ├── polar_molecule.py        PolarMolecule, KRB, RBCS
    ├── polar_rydberg/
    │   └── polar_system.py      PolarBOPSystem sin Fermi
    ├── rb_krb_polar/
    │   ├── charge_dipole.py     operador carga–dipolo compartido
    │   ├── rb_defects.py        niveles vecinos y override del paper
    │   └── bop_system.py        BOPSystem histórico compatible
    ├── rb_rbcs_polar/
    │   └── __init__.py          RbRbCsPolarSystem
    ├── double_polar_rydberg/
    │   ├── basis.py             base de dos rotores, bloques M_J
    │   ├── contracted.py        base pendular local contraída experimental
    │   └── system.py            Hamiltoniano disperso y solver selectivo
    ├── rb_neutral_perturber/
    │   ├── fermi_krb.py         FermiPseudopotential moderno
    │   ├── linear_trimer.py     sistema lineal moderno y validado
    │   ├── fermi_potentials.py  legado protegido por goldens
    │   └── trimer.py            entrada end-to-end del legado
    ├── hybrid_neutral_polar/
    │   ├── hybrid_system.py     HybridNeutralPolar
    │   └── parity.py            transformación theta=pi
    └── nonadiabatic_dynamics/
        ├── coupling.py          acoplamientos de derivada
        ├── energy_normalization.py
        ├── coupled_channels.py
        ├── stabilization.py
        ├── decay_rates.py
        └── franck_condon.py
```

Regla conceptual: `mathlib → basis → systems → scripts`. `visualization` es
paralela a la física y no debe importar sistemas.

## Núcleo polar

### Catálogo molecular

`PolarMolecule` contiene nombre, constante rotacional en Hz y dipolo en Debye,
con propiedades derivadas `B_au`, `B_ghz` y `d_au`.

```python
from trimero.systems.polar_molecule import KRB, RBCS
```

El catálogo es la fuente de configuración de los motores nuevos. El híbrido
consume `RBCS`; sus argumentos `B_mhz` y `d_debye` llegan realmente a
`ChargeDipoleHamiltonian`. `charge_dipole.py` conserva aliases numéricos de KRb
por compatibilidad con imports y tests históricos (véase Deuda).

### `PolarBOPSystem`

```python
from trimero.systems.polar_molecule import RBCS
from trimero.systems.polar_rydberg import PolarBOPSystem

system = PolarBOPSystem(
    molecule=RBCS,
    n_manifold=25,
    N_max=6,
)
H = system.hamiltonian(R=800.0, M_J=0)
E, V = system.solve(R=800.0, M_J=0)
```

Responsabilidades:

- construir `CoupledBasis` y `RadialBasis`;
- montar `H_A + H_mol` sin importar datos de dispersión;
- proporcionar máscara de manifold, umbrales y cero energético;
- seleccionar el estado más bajo cuyo peso de manifold supera el umbral.

La base electrónica es `(n,l>=3)+(n+1)d+(n+2)p+(n+3)s`. Su composición sale de
`neighbor_levels(n)`; no debe duplicarse en scripts.

### Compatibilidad `BOPSystem`

`rb_krb_polar.BOPSystem` conserva el comportamiento histórico, incluida la
construcción de `FermiPseudopotential` y el argumento `fermi`. Algunos tests de
caracterización dependen de él. Los scripts polares nuevos no lo usan: llaman a
`PolarBOPSystem`, donde Fermi no existe ni como interruptor.

### Identificación de curvas

El índice espectral no identifica un objeto físico a través de cruces evitados.
La selección puntual usa peso de manifold; cuando se necesita continuidad se
usa máximo solapamiento entre autovectores consecutivos. Los `.npz` conservan
`K`, `W` y, cuando aplica, `overlap` para auditar la asignación.

## Perturbador neutro

`SymmetricLinearTrimer` es la capa moderna para Rb–Rb\*–Rb. Usa las tablas de
`data/Wavefunction/`, los canales s/p y el campo DC. El buen número cuántico es
`m_l`; no hay rotor molecular.

`trimer.py` y `fermi_potentials.py` son traducciones del C++ y están congelados
por goldens G1→G4. Conservan deliberadamente peculiaridades numéricas del
legado. No deben refactorizarse como efecto lateral de otro cambio.

## Sistema híbrido

`HybridNeutralPolar` compone, sobre la misma `CoupledBasis`:

```text
diag(H_A) + ChargeDipoleHamiltonian(R2) + FermiPseudopotential_pi(R1)
```

- RbCs define el eje `+Z` y está en `theta=0`.
- El Rb neutro está en `theta=pi`.
- `phase_pi(l1,l2)=(-1)^(l1+l2)` transforma el elemento de Fermi.
- Ambos términos conservan `M_J=m_l+M_N` en la geometría axial.
- `fermi=False` reproduce bit a bit `PolarBOPSystem(molecule=RBCS)` cuando se
  usa la misma base radial.

La curva híbrida de producción parte de una referencia polar sin Fermi y se
sigue adiabáticamente. El autovalor más bajo del Hamiltoniano completo no es,
en general, la curva de ligadura buscada.

## Dinámica no adiabática

Esta capa consume curvas electrónicas ya validadas; no decide qué BOP es
física. Flujo:

```text
BOP + autovectores
      │
      ▼
coupling ─► energy_normalization ─► coupled_channels
                                         │
                       ┌─────────────────┼──────────────┐
                       ▼                 ▼              ▼
                 stabilization      decay_rates   franck_condon
```

Antes de usarla deben estar fijados el cero absoluto, la masa reducida, la
malla, la continuidad del estado electrónico y el dominio donde conserva el
carácter pertinente.

## Scripts de producción

| Script | Sistema | Salida |
|---|---|---|
| `compute_bop_curve.py` | KRb/RbCs polar | `plots/rb_<mol>_polar/{data,figures}` |
| `compute_orientation_curve.py` | KRb/RbCs polar | misma raíz |
| `compare_bop_curves_n.py` | comparación polar | `figures/` |
| `compare_orientation_n.py` | comparación polar | `figures/` |
| `plot_orientation_alignment.py` | orientación y alineamiento polar, más rotor de referencia | `figures/` |
| `compute_field_curves.py` | KRb/RbCs polar con campo DC paralelo a Z | misma raíz; un panel por campo |
| `compute_double_rbcs_curves.py` | dos RbCs, BOP y orientación | `plots/rb_rbcs_rbcs_polar/` |
| `check_double_rbcs_convergence.py` | convergencia del sistema polar doble | terminal |
| `check_double_rbcs_contracted.py` | convergencia pendular contraída | terminal |
| `check_polar_convergence.py` | convergencia `N_max` | terminal |
| `compute_trimer_curves.py` | perturbador neutro | `plots/rb_neutral_perturber/` |
| `compute_hybrid_curves.py` | híbrido | `plots/hybrid_neutral_polar/` |
| `analyze_nonadiabatic_*.py` | dinámica | datos de fases 6/6b |
| `draw_geometry.py` | diagramas | `plots/geometry/` |

`scripts/archive/` contiene exploraciones históricas y no es API de producción.

## Datos y artefactos

```text
data/Wavefunction/                 insumos del pseudopotencial neutro
plots/<sistema>/data/             resultados numéricos
plots/<sistema>/figures/          figuras
plots/geometry/                    esquemas compartidos
plots/archive/                     resultados históricos
graphify-out/                      grafo regenerable, no versionado
```

Los `.npz` polares nuevos incluyen metadatos de especie, `B`, `d`, `n`,
`N_max`, `M_J`, peso de carácter y versión de esquema. `--reuse` rechaza datos
sin metadatos compatibles.

## Tests e invariantes

El tamaño de la suite lo dice pytest, no este documento
(`poetry run pytest --collect-only -q | tail -1`; 177 tests el 2026-08-30, de
los cuales 16 marcados `slow`). Familias principales:

- base y reglas de selección;
- hermiticidad del operador carga–dipolo y campo electrónico;
- regresión de la Fig. 1 de KRb;
- equivalencia bit a bit entre núcleo genérico y compatibilidad KRb;
- equivalencia RbCs con el límite polar del híbrido;
- límites neutro/polar del híbrido y conservación de `M_J`;
- regresión del trímero neutro de 2016;
- dinámica no adiabática por módulo;
- goldens G1→G4 del legado;
- diagramas geométricos.

Invariantes:

1. `Atom.E_Rb()` es la fuente de verdad de energías atómicas; el override
   `DELTA0_NS_PAPER` es local a la comparación bibliográfica.
2. Unidades internas: `a0`, `E_h`, carga elemental. Conversiones sólo en
   entrada/salida.
3. Las matrices deben ser hermíticas dentro de la tolerancia documentada.
4. Los límites apagando un acoplamiento deben coincidir con objetos frescos e
   independientes, preferentemente bit a bit.
5. Un golden no se regenera para hacer desaparecer una regresión.
6. Cada especie/manifold debe demostrar convergencia en `N_max`.

## Deuda técnica

1. `ChargeDipoleHamiltonian`, `rydberg_diagonal` y `rb_defects` aún viven bajo
   el nombre histórico `rb_krb_polar`, aunque el motor genérico también los usa
   para RbCs. La física es genérica; el nombre de paquete no lo es todavía.
2. `rb_defects.py` mezcla física atómica compartida y composición de la base
   del paper polar. Esto crea imports desde `basis/radial.py` y desde el lado
   neutro hacia un paquete nombrado KRb.
3. `BOPSystem` histórico sigue acoplado a Fermi. Se conserva únicamente por
   compatibilidad y goldens; no debe reaparecer en scripts polares nuevos.
4. Las constantes `B_KRB_*`/`D_KRB_*` siguen exportadas desde
   `charge_dipole.py`, además del catálogo molecular. Deben converger hacia
   aliases del catálogo cuando se rompa la dependencia circular de imports.
5. Las radiales de los niveles vecinos usan hidrogenoides con `n_eff` entero;
   la energía sí usa el defecto cuántico exacto.

## Rendimiento

Para `n=25`, `N_max=6`, el bloque `M_J=0` tiene dimensión 1113. Construir las
radiales cuesta segundos; deben reutilizarse durante todo el barrido. La
diagonalización densa domina el coste. `eigvalsh` es más rápido, pero la
selección por carácter y el seguimiento requieren autovectores y, por tanto,
`eigh`.

Cada radio es matemáticamente independiente, aunque la lógica de seguimiento
por solapamiento debe consumir los resultados en orden radial.
