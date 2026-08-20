# Estado vigente del proyecto

**Fecha**: 2026-08-20 · **Autor**: Javier Aguilera

Documento corto y de entrada. Lo que no esté aquí, o no esté enlazado desde
aquí, no es referencia activa.

El repositorio contiene **dos sistemas físicos distintos** que durante varias
rondas estuvieron mezclados en el mismo espacio de nombres. Están separados
desde la reorganización del 2026-08-20 (`docs/archive/analysis_reorganizacion_20260820.md`).

| | `systems/rb_krb_polar/` | `systems/rb_neutral_perturber/` |
|---|---|---|
| perturbador | KRb, **polar** | átomo/molécula **neutra** |
| interacción | carga-dipolo, `-d·F_ryd` | dispersión de contacto, pseudopotencial de Fermi |
| referencia | Aguilera-Fernández 2015 / González-Férez 2015 | Aguilera-Fernández 2016 |
| estado | **vigente** | código verde, línea en pausa |

---

## Rb*-KRb: el modelo vigente

### Hamiltoniano

```
H_ad(R) = H_A + H_mol
        = diag(E_ryd) + [ B·N² − d·F_ion(R) − d·F_elec(R) ]
```

Ec. 1 de Aguilera-Fernández et al. 2015. KRb es un **dipolo puntual** en el
campo eléctrico del Rydberg.

**No lleva pseudopotencial de Fermi.** El de contacto modela un perturbador
neutro; aplicarlo aquí fue la premisa equivocada de varias rondas. Como
consecuencia **no existen** en este sistema: remapeo semiclásico `k(R)`,
ventana de exclusión de la resonancia de onda p, cola butterfly ni tope de
dominio en `2n²a₀`. `H_A` y `H_mol` se evalúan con la función de onda
hidrogenoide, definida en todo `R`, así que el rango R ∈ [400, 1800] a₀ sale
entero y continuo.

### Base electrónica

Manifold cuasi-degenerado **más tres niveles vecinos individuales**:

```
(n, l ≥ 3)  +  (n+1)d  +  (n+2)p  +  (n+3)s
```

Para n=25: manifold + 26d + 27p + 28s → dim(M_J=0) = 1113, dim(M_J=1) = 1106.
Son **tres** vecinos, no uno: usar sólo (n+3)s es una base incompleta. Única
definición en `rb_defects.neighbor_levels()`; no se escribe a mano en ningún
otro sitio.

### Identificación de la curva: por CARÁCTER, no por índice

La curva BOP es **la más baja con peso de manifold > 50 %**, evaluado en cada
`R`. Un índice fijo identificado en un extremo **no** vale: los estados de
(n+1)d y (n+2)p producen cruces evitados con las del manifold y el índice
cambia de objeto por el camino (medido: k=55 en R=1800 a₀ para M_J=0, distinto
a R menor). Esto sigue siendo necesario **aunque no haya V_Fermi**.

### δ₀(ns) = 3.13180

El defecto cuántico s por defecto del código (`rb_atom.Atom.E_Rb()`) es
δ₀ = 3.1311804, de Li, Mourachko, Noel & Gallagher, PRA 67, 052502 (2003). Es
correcto y es la fuente de verdad del código legacy, protegida por golden files.

Pero la Tabla I de González-Férez 2015 **no se reproduce con él**: p, d y f
concuerdan a ±0.004 GHz y los dos niveles s fallan por −0.299 y +0.265 GHz.
Ambos piden el mismo desplazamiento constante de δ₀ (+6.19e−4), lo que descarta
δ₂, la dependencia en n, la masa reducida y las unidades. El paper toma sus
energías de Marinescu, Sadeghpour & Dalgarno, PRA 49, 982 (1994), que es un
cálculo de **potencial modelo**, no un ajuste Rydberg-Ritz; que la discrepancia
esté confinada a l=0 —el canal más penetrante en el core— encaja con esa
diferencia de método.

`DELTA0_NS_PAPER = 3.13180` es por tanto un **override local para comparar con
este paper**, no una corrección general. `delta0_ns=None` delega exactamente en
`Atom.E_Rb()` y no cambia nada.

### Cómo se calcula

```bash
poetry run python scripts/compute_bop_curve.py --n-manifold 25 --mj 0 1
```

Único script de producción. Los doce `run_*.py` / `analyze_*.py` de las rondas
de exploración están en `scripts/archive/`.

---

## Deuda técnica conocida

1. **`BOPSystem` sigue acoplado a `fermi_krb`.** `__post_init__` construye
   siempre un `FermiPseudopotential` (lee `rvsAS.dat`/`rvsAP.dat`) y
   `hamiltonian()` tiene `fermi=True` por defecto, aunque el modelo polar no lo
   use. Es la última arista `rb_krb_polar → rb_neutral_perturber`. **No se tocó
   en la reorganización porque desacoplarlo cambia números** y hay tres tests
   (`test_s2`, `test_s4`, `test_s5`) que dependen del comportamiento actual.
   El lado polar se protege pasando `fermi=False` explícitamente.

2. **`rb_defects.py` mezcla dos cosas.** `mu_ns`/`n_star_ns`/`n_star_nl`/
   `energy_rb` son física atómica de Rb, compartida; `NEIGHBOR_DN`/
   `neighbor_levels`/`DELTA0_NS_PAPER` son la composición de la base del paper
   polar. Por eso `basis/radial.py` (capa compartida) y
   `rb_neutral_perturber/fermi_krb.py` importan de `rb_krb_polar/`. Separarlas
   cerraría las dos inversiones de capa que quedan.

---

## Referencia activa

Estos cinco documentos siguen vigentes. **Todo lo demás en `docs/archive/` es
material del otro sistema físico**, técnicamente correcto pero no aplicable a
Rb*-KRb.

| documento | qué aporta |
|---|---|
| [`analysis_fig1_carga_dipolo_sin_fermi.md`](analysis_fig1_carga_dipolo_sin_fermi.md) | **el cálculo bueno**: `H_ad = H_A + H_mol`, curvas M_J=0 y M_J=1 completas, y la corrección de premisa que ordena todo lo anterior |
| [`analysis_base_correcta_3_vecinos.md`](analysis_base_correcta_3_vecinos.md) | la base de tres vecinos y el criterio de carácter (§5.1) |
| [`analysis_verificacion_tabla_I.md`](analysis_verificacion_tabla_I.md) | δ₀(ns) = 3.13180 y su justificación completa |
| [`analysis_validacion_carga_dipolo.md`](analysis_validacion_carga_dipolo.md) | `B·N²` + campo del ion: derivación, 4 tests analíticos, escalado 1/R⁴ |
| [`analysis_campo_electron_rydberg.md`](analysis_campo_electron_rydberg.md) | campo del electrón Rydberg (Ec. A.6–A.10): expansión multipolar, validación contra cuadratura 2D |

Números de referencia verificados (n=25, M_J=0, R ∈ [400, 1800] a₀, paso 5 a₀),
que `tests/systems/rb_krb_polar/test_regression_fig1.py` protege:

| magnitud | valor |
|---|---|
| pozo más profundo | **−23.100 GHz** (R = 400 a₀) |
| E(R = 1800 a₀) | **−0.338 GHz** |
| mínimos locales | **8** |

## Referencias

- Aguilera-Fernández, J., Sadeghpour, H. R., Schmelcher, P. & González-Férez, R.,
  *J. Phys.: Conf. Ser.* **635**, 012023 (2015), arXiv:1507.07972 — Ec. 1, Fig. 1.
- González-Férez, R., Sadeghpour, H. R. & Schmelcher, P., *New J. Phys.* **17**,
  013021 (2015) — Tabla I, Fig. 2.
- Marinescu, M., Sadeghpour, H. R. & Dalgarno, A., *Phys. Rev. A* **49**, 982 (1994).
- Li, W., Mourachko, I., Noel, M. W. & Gallagher, T. F., *Phys. Rev. A* **67**,
  052502 (2003) — defectos cuánticos por defecto del código.
