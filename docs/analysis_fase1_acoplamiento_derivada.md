# Fase 1 — Acoplamiento no adiabático de derivada ⟨Ψᵢ(R)|d/dR|Ψⱼ(R)⟩

**Fecha**: 2026-08-22
**Autor**: Javier Aguilera (con Claude Code)
**Relevancia**: primer módulo de `docs/PLAN_nonadiabatic_dynamics.md`, la
réplica para Rb*-RbCs de la maquinaria de dinámica no adiabática de
Mellado-Alcedo et al., PRA 110, 013314 (2024). Sin este acoplamiento no hay
ecuación de canales acoplados (Fase 2) ni tasas de decaimiento (Fase 4).

## Resumen

Se implementó `trimero.systems.nonadiabatic_dynamics.coupling`, módulo
genérico (no conoce ningún sistema físico concreto) que calcula
⟨Ψᵢ(R)|d/dR|Ψⱼ(R)⟩ por diferencias finitas centradas a partir de cualquier
familia `solve(R) -> (w, V)` de Hamiltonianos ya diagonalizados. El paso
crítico —fijar la continuidad de fase de los autovectores entre R vecinos,
porque `numpy.linalg.eigh` no la garantiza— se resuelve maximizando el
solapamiento con el paso anterior de forma secuencial.

Se verificaron, en este orden y ANTES de aceptar el módulo:

1. Antisimetría exacta (hasta el orden de la discretización).
2. Diagonal nula (idem).
3. Acoplamiento pequeño y suave lejos de cruces, con un sistema real.
4. Pico de acoplamiento en un cruce evitado real, ya documentado en el
   proyecto (`plots/rb_krb_polar/data/fig1_ad_MJ0_n25.npz`).
5. Convergencia con el paso de diferencias finitas (orden 2, verificado).

Los 11 tests están en verde. Ningún test se relajó para acomodar un resultado
sospechoso: los dos ajustes de tolerancia que hicieron falta (ver §4) se
verificaron primero como error de discretización esperado, no como bug.

## Método

### Modelo de referencia analítico (2 niveles)

Para los tests (a), (b) y (e) se usa un sistema 2×2 real simétrico con
solución analítica exacta:

```
H(R) = [[a(R), c], [c, b(R)]],   a(R) = -k(R-R0),  b(R) = +k(R-R0)
```

Con ángulo de mezcla θ(R) = ½·atan2(2c, a(R)-b(R)), los autovectores son
Ψ_lower=(cosθ,sinθ), Ψ_upper=(-sinθ,cosθ), y el acoplamiento de derivada
exacto es

```
|A_12(R)| = |dθ/dR| = 2ck / [4k²(R-R0)² + 4c²]
```

(derivación por regla de la cadena de θ(R), en el docstring de
`tests/systems/nonadiabatic_dynamics/test_coupling.py`). Con k=0.01, c=0.05,
R0=10.0 el pico exacto en R0 es |A_12(R0)| = k/(2c) = 0.1.

### Cruce evitado real usado en el test (d)

`plots/rb_krb_polar/data/fig1_ad_MJ0_n25.npz` (barrido de producción,
`n_manifold=25`, `M_J=0`, sin Fermi, paso 5 a₀) registra que el índice k de
la curva de carácter de manifold cambia de 54 a 55 cerca de R=1300 a₀. Un
escaneo fino del hueco de energía w[55]-w[54] (paso 5 a₀ entre R=1270 y
1330 a₀, `BOPSystem(n_manifold=25, delta0_ns=DELTA0_NS_PAPER)`, `M_J=0`,
`fermi=False`) da:

```
R= 1270.0  gap(54,55)=3.804420e-07 Eh
R= 1290.0  gap(54,55)=1.550457e-07 Eh
R= 1300.0  gap(54,55)=4.724116e-08 Eh
R= 1305.0  gap(54,55)=4.982161e-09 Eh   <- mínimo real
R= 1310.0  gap(54,55)=5.595015e-08 Eh
R= 1320.0  gap(54,55)=1.538857e-07 Eh
```

El mínimo real está en R≈1305 a₀ (hueco ≈5×10⁻⁹ Eh), más estrecho de lo que
sugería la resolución de 5 a₀ del barrido de producción — coherente con la
nota de `docs/STATUS.md` sobre por qué el índice de carácter no es un
identificador fijo (§"Identificación de la curva: por CARÁCTER").

## Resultados numéricos reales

Suite completa (`poetry run pytest tests/systems/nonadiabatic_dynamics/ -v -s`,
incluye el test `slow`):

```
tests/systems/nonadiabatic_dynamics/test_coupling.py::test_fix_eigenvector_signs_flips_negative_overlap_columns PASSED
tests/systems/nonadiabatic_dynamics/test_coupling.py::test_fix_eigenvector_signs_leaves_positive_overlap_columns PASSED
tests/systems/nonadiabatic_dynamics/test_coupling.py::test_antisymmetry_toy_model PASSED
tests/systems/nonadiabatic_dynamics/test_coupling.py::test_diagonal_zero_toy_model PASSED
tests/systems/nonadiabatic_dynamics/test_coupling.py::test_antisymmetry_and_diagonal_zero_real_system_away_from_crossing PASSED
tests/systems/nonadiabatic_dynamics/test_coupling.py::test_smooth_away_from_crossings_real_system PASSED
tests/systems/nonadiabatic_dynamics/test_coupling.py::test_peak_at_real_documented_avoided_crossing
  |A(54,55)| en cruce (R=1305): 5.0672e-01
  |A(54,55)| lejos    (R=900) : 5.9795e-03
PASSED
tests/systems/nonadiabatic_dynamics/test_coupling.py::test_convergence_with_finite_difference_step
  errores vs h=(0.2,0.1,0.05): [2.2049825835077086e-05, 5.510357863433568e-06, 1.3774578919156788e-06]
PASSED
tests/systems/nonadiabatic_dynamics/test_coupling.py::test_eigenbasis_along_r_matches_analytic_coupling_magnitude PASSED
tests/systems/nonadiabatic_dynamics/test_coupling.py::test_derivative_coupling_rejects_nonuniform_grid PASSED
tests/systems/nonadiabatic_dynamics/test_coupling.py::test_derivative_coupling_requires_v_or_solve PASSED

============================= 11 passed in 12.25s ==============================
```

### (d) Pico en el cruce evitado real

|A(54,55)| en R=1305 a₀ (h=1 a₀): **0.5067**
|A(54,55)| en R=900 a₀ (h=1 a₀, lejos de cualquier cruce conocido de este par):
**5.980×10⁻³**

Razón: **84.75×**. El test exige >50× (margen sobre el valor medido, sin ser
tan laxo que un bug que reduzca el pico a la mitad pase desapercibido).

### (e) Convergencia con el paso h (orden 2)

Punto de prueba R=8.0 (modelo analítico 2 niveles, lejos del pico agudo en
R0=10, régimen donde domina el término O(h²)):

| h | error vs. exacto |
|---|---|
| 0.2 | 2.2050×10⁻⁵ |
| 0.1 | 5.5104×10⁻⁶ |
| 0.05 | 1.3775×10⁻⁶ |

Razones sucesivas: 4.002× y 4.001× — coincide con el orden 2 esperado de la
diferencia centrada (el error debería caer ~4× al reducir h a la mitad).

## Dos ajustes de tolerancia (documentados, no bugs)

Al ejecutar la suite por primera vez fallaron dos tests por miscalibración
del criterio, no por un error del módulo. Se investigó cada uno antes de
tocar nada:

1. **`test_diagonal_zero_toy_model`**: el primer `atol=1e-8` era
   injustificadamente estricto. La identidad A_ii=0 es exacta
   *analíticamente* (de la normalización ⟨Ψᵢ|Ψᵢ⟩=1 para todo R), pero la
   fórmula discreta de diferencias centradas sólo la respeta a O(h²). Se
   verificó explícitamente el escalado antes de tocar la tolerancia:

   ```
   h=0.200  max|diag|=2.0677e-05  max|diag|/h^2=5.1692e-04
   h=0.100  max|diag|=5.1735e-06  max|diag|/h^2=5.1735e-04
   h=0.050  max|diag|=1.2939e-06  max|diag|/h^2=5.1756e-04
   h=0.025  max|diag|=3.2349e-07  max|diag|/h^2=5.1759e-04
   ```

   El coeficiente h⁻² es estable (≈5.176×10⁻⁴) en cuatro valores de h — es
   error de discretización de orden 2, no un bug. `atol` se subió a `2e-6`
   (margen sobre el residuo medido con h=0.05, 1.294×10⁻⁶), documentando la
   razón en el propio test.

2. **`test_smooth_away_from_crossings_real_system`**: el criterio original
   comparaba el cociente `coupling[k+1]/coupling[k]` contra un piso absoluto
   `max(coupling[k], 1e-12)`. El acoplamiento real (BOPSystem n_manifold=6,
   N_max=2, estados 6 y 7, R∈[400,1800] a₀) decae **monótona y suavemente**
   de ~9×10⁻¹⁰ a ~6×10⁻¹⁵ Eh — nueve órdenes de magnitud en el rango, sin
   ningún salto real. El piso de 10⁻¹² se saturaba mucho antes del final del
   rango y producía cocientes artificialmente pequeños que no medían
   suavidad real, sino el propio piso. Se sustituyó por un chequeo de
   monotonía estricta (`np.diff(coupling) <= 0`, con margen numérico), que sí
   mide lo que el test pretende verificar (ausencia de picos o rebotes que
   delatarían un cruce oculto) sin depender de un piso arbitrario.

Ninguno de los dos ajustes cambia lo que el test verifica; ambos se
investigaron primero con el módulo ya implementado y datos reales antes de
decidir el nuevo criterio, según la disciplina de este documento maestro.

## Conclusiones y aplicación al proyecto

- El módulo `coupling.py` es correcto dentro del orden esperado de su
  discretización (verificado con referencia analítica independiente, no sólo
  contra sí mismo).
- El acoplamiento de derivada, aplicado a un cruce evitado REAL ya presente
  en los datos de producción del sistema Rb*-KRb, muestra el comportamiento
  físico esperado: pico agudo (∼0.5, casi divergente) en el mínimo real del
  hueco de energía, y valores 1-4 órdenes de magnitud menores lejos de él.
- Confirma que la disciplina de identificación de curvas "por CARÁCTER, no
  por índice fijo" de `docs/STATUS.md` tiene una contraparte cuantitativa
  directa: el cruce evitado que obliga a ese criterio es precisamente donde
  el acoplamiento no adiabático diverge — el punto donde la aproximación
  Born-Oppenheimer (seguir UNA curva adiabática) deja de ser válida.
- Con esto queda cerrada la Fase 1. La Fase 2 (ecuación radial de canales
  acoplados) puede alimentarse directamente de `derivative_coupling` sin
  cambios: acepta tanto `solve(R)` como una malla `V` precalculada, tal como
  pide el plan.

## Archivos tocados

- `src/trimero/systems/nonadiabatic_dynamics/__init__.py` (nuevo)
- `src/trimero/systems/nonadiabatic_dynamics/coupling.py` (nuevo)
- `tests/systems/nonadiabatic_dynamics/test_coupling.py` (nuevo, 11 tests)
- `docs/PLAN_nonadiabatic_dynamics.md` (nuevo, Fase 0)
- `docs/analysis_fase1_acoplamiento_derivada.md` (este documento)

## Referencias

- Mellado-Alcedo, Guttridge, Cornish, Sadeghpour & González-Férez,
  *Phys. Rev. A* **110**, 013314 (2024), arXiv:2401.09618 — Ec. 1-4 y 8, el
  paper que se está replicando para Rb*-RbCs.
- `plots/rb_krb_polar/data/fig1_ad_MJ0_n25.npz` — datos de producción donde se
  localizó el cruce evitado real usado en el test (d).
- `docs/STATUS.md` §"Identificación de la curva: por CARÁCTER, no por
  índice" — el precedente cualitativo que este documento cuantifica.
