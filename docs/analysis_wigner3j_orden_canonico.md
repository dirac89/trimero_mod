# Orden canónico en `wigner_3j`: hermiticidad exacta y no regresión

**Fecha**: 2026-08-21
**Autor**: ox-alpha (con dirección de Javier Aguilera)
**Relevancia**: Arregla la asimetría de hermiticidad (‖H−Hᵀ‖/‖H‖ ≈ 2.8e−13) del sistema híbrido neutro-polar a n=35 y de todo lo que usa `gaunt`/`RydbergElectronField`; deja la hermiticidad a nivel 3e−19.
**Tipo**: analysis

## Resumen

El test L-c (hermiticidad completa del híbrido) fallaba con ‖H−Hᵀ‖_F/‖H‖_F = 2.8e−13
frente al convenio del repo de 1e−14. La causa raíz NO estaba en el ensamblado sino
en `wigner_3j` (`src/trimero/mathlib/angular.py`): la fórmula de Racah evaluada en
punto flotante no es simétrica bajo permutación de columnas, así que
⟨a|O|b⟩ y ⟨b|O|a⟩, que pasan por símbolos 3j relacionados por permutación,
diferían a nivel 1e−9..1e−10 relativo por símbolo.

La solución es un **orden canónico de argumentos** dentro de `wigner_3j`: las 12
formas equivalentes del símbolo (6 permutaciones de columnas × 2 orientaciones
globales de los m) colapsan en UNA sola evaluación de Racah bit a bit, y las
relaciones entre formas se restauran con fases enteras exactas ±1.

## Palabras Clave

- símbolos 3j de Wigner
- fórmula de Racah
- cancelación catastrófica
- hermiticidad bit a bit
- orden canónico de argumentos

## Contenido Principal

### 1. Cadena causal verificada experimentalmente

| Eslabón | Asimetría medida (antes del arreglo) |
|---|---|
| `wigner_3j` swap columnas (l~30) | 2.7e−8 rel (caso k=37: …900973 vs …759456) |
| `gaunt` peor caso | 1.4e−9 |
| `field_element` swap | 7.0e−16 abs sobre elementos ~1.25e−6 (5.6e−10 rel) |
| H polar puro n=35, R=800 | ‖H−Hᵀ‖/‖H‖ = 1.3e−13 |
| H híbrido L-c | 2.8e−13 |

La cancelación de la suma-k quedó descartada como causa (sólo 8.3 % del peor par;
el ruido esperado era ~3e−23). El problema es el redondeo de gammaln en las
distintas ordenaciones equivalentes de la serie de Racah.

### 2. El arreglo: canonización dentro de `wigner_3j`

Regla implementada (todo en aritmética entera salvo la propia serie):

1. **Orientación global de los m**: se comparan lexicográficamente las secuencias
   ordenadas de pares `(j, m)` con m tal cual y con m volteado; la mayor gana.
   Voltear cuesta la fase exacta `(−1)^{j1+j2+j3}` (simetría m→−m).
2. **Orden ascendente de pares** `(j, m)`: el mayor j queda en la posición 3, lo
   que minimiza `j1+j2−j3`, acorta la serie alternante y reduce la cancelación.
3. **Fase de permutación**: inversiones ESTRICTAS de la secuencia ya orientada
   (pares idénticos no cuentan); impar ⇒ fase exacta `(−1)^{j1+j2+j3}`.

Con esto cualquier permutación de columnas o volteo simultáneo de los tres m
comparte la MISMA evaluación de Racah bit a bit. Las simetrías de Regge quedan
fuera a propósito (no se necesitan para hermiticidad).

Detalle importante descubierto al iterar: canonizar SÓLO columnas dejaba un
residuo (rel 1.9e−13 en L-c) porque ⟨i|V|j⟩ y ⟨j|V|i⟩ involucran símbolos
relacionados por permutación **y** volteo de m a la vez; sin fijar la orientación
global caen en M-vectores canónicos distintos. Con la regla completa, L-c baja a
**3.0e−19 relativo**.

### 3. Validación contra aritmética racional exacta

Referencia independiente: fórmula de Racah con la suma alternante evaluada en
`fractions.Fraction` (factoriales enteros exactos); sólo la raíz del prefactor es
flotante (~1 ulp). Muestra de 300 casos de los 3638 que cambian en la rejilla
m=0, l ≤ 26:

```text
error total legacy vs exacto: 5.405e-11
error total nuevo  vs exacto: 4.846e-11
ratio agregado: 1.12x a favor del nuevo
por caso: 166 mejora / 0 igual / 134 peora   (redondeo aleatorio caso a caso)
peor error relativo del nuevo frente al exacto: 2.099e-10
```

Es decir: la canonización garantiza CONSISTENCIA exacta y una mejora
estadística modesta; caso a caso el error flotante residual (~1e−10 rel en
l~26) es aleatorio entre formas equivalentes y muy inferior a cualquier
tolerancia física o de tests.

## Fórmulas Clave

```math
w_{3j}(j_2\,j_1\,j_3;\,m_2\,m_1\,m_3) = (-1)^{j_1+j_2+j_3}\, w_{3j}(j_1\,j_2\,j_3;\,m_1\,m_2\,m_3)
w_{3j}(j_1\,j_2\,j_3;\,-m_1\,-m_2\,-m_3) = (-1)^{j_1+j_2+j_3}\, w_{3j}(j_1\,j_2\,j_3;\,m_1\,m_2\,m_3)
```

Ambas se satisfacen ahora BIT A BIT (la segunda tras multiplicar por la fase
exacta; sólo puede diferir el signo de cero en símbolos nulos).

## Verificación (salidas reales)

### Simetría del símbolo

```text
=== (a) Simetría exacta bajo permutaciones ===
l<= 26, modo=m0:    1925 símbolos; peor desviación: None
l<= 26, modo=varios: 7700 símbolos; peor desviación: None
l<= 35, modo=m0:    4389 símbolos; peor desviación: None
=== (a2) Relación m->-m bit a bit ===
comprobados 15996; peor desviación: (0.0, ...)   # sólo signo de cero
=== Caso notorio (30,31,37): 6 permutaciones bit-idénticas ===
np.float64(-0.018544327145550165)  x6
```

### Hermiticidad L-c (n=35, híbrido)

```text
R1=900  R2=800  M_J=+0: dim=307, ||H-H^T||_F = 2.115e-21 (rel 3.0e-19)
R1=600  R2=1200 M_J=+3: dim=288, ||H-H^T||_F = 1.044e-21 (rel 1.5e-19)
R1=1100 R2=500  M_J=-2: dim=297, ||H-H^T||_F = 6.832e-22 (rel 9.7e-20)
```
(antes: 2.8e−13 relativo — mejora de ~6 órdenes de magnitud)

### No regresión

```text
tests/systems/hybrid_neutral_polar/            12 passed in 26.14s
suite rápida completa (-m "not slow")          75 passed, 14 deselected in 53.18s
test_regression_trimer_2016.py                 11 passed in 47.44s  (= baseline)
polar rápido (fig1 R1 + rydberg_field + cd)    14 passed, 1 failed (ver abajo)
```

## ⚠️ Desviación pendiente de decisión: fig1 R2

El test lento `test_r2_compute_bop_curve_reproduces_the_reference` exige
reproducir el golden `plots/rb_krb_polar/data/fig1_ad_MJ0_n25.npz` con
`rtol=1e-12` ("bit a bit"). El golden se generó con el `wigner_3j` VIEJO; el
nuevo es más preciso y mueve los autovalores por encima de esa tolerancia
ultra-estricta:

```text
R(a0)      E_ref(GHz)         E_nuevo(GHz)       diff_abs     diff_rel   k
400.0  -23.100158642148  -23.100158642151  2.85e-12  1.23e-13     53 = 53
600.0  -19.530216131768  -19.530216131778  1.07e-11  5.48e-13     54 = 54
800.0  -18.091318574598  -18.091318574616  1.85e-11  1.03e-12     54 = 54
1000.0 -19.057620440161  -19.057620440159  2.14e-12  1.12e-13     54 = 54
1200.0 -16.348166655673  -16.348166655672  7.14e-13  4.37e-14     54 = 54
1500.0  -2.126473245243   -2.126473245245  2.14e-12  1.01e-12     55 = 55
1800.0  -0.337614270502   -0.337614270500  2.14e-12  6.34e-12     55 = 55
```

Lectura física: desviaciones ≤ 2.1e−11 GHz ≈ 0.02 Hz; el índice de carácter k
es idéntico en los 7 puntos; la profundidad −23.100 GHz y E(1800)=−0.338 GHz
documentados no cambian (R1 sigue pasando). Es un cambio de PRECISIÓN, no de
física, pero excede la tolerancia declarada del test ⇒ según la instrucción de
trabajo, SE PARA Y SE REPORTA sin tocar el golden.

Opciones:
- **A. Regenerar el golden** con el código actual (la referencia física R1
  contra los números del documento seguiría pasando redondeada a 3 decimales).
- **B. Relajar rtol de R2** a 1e-9 documentando que el cambio proviene de una
  mejora de precisión de `wigner_3j` validada contra Racah racional exacto.
- **C. Revertir el arreglo** (no recomendado: reintroduce la asimetría de
  hermiticidad en TODO lo que usa gaunt).

## Conclusiones y Aplicación al Proyecto

- `wigner_3j` canónico elimina la única fuente sistemática de asimetría de la
  cadena angular; radiales (`dg_integrals`) y `ion_field_element` ya eran
  simétricos bit a bit; `field_element` cumple ⟨i|F_q|j⟩ = (−1)^q⟨j|F_{−q}|i⟩
  exactamente.
- El sistema neutro (trímero/Fermi) NO usa gaunt/wigner_3j: sus resultados son
  inalterables (11/11 idéntico a baseline).
- Pendiente: decisión A/B/C para fig1 R2; después, continuar con las Fases 1-2
  del plan del híbrido (barridos R2 y espectros), que quedaron en pausa.

## Referencias

- `src/trimero/mathlib/angular.py` — implementación canónica y docstring.
- Racah, G. *Phys. Rev.* **62**, 438 (1942) — fórmula original.
- [analysis_fig1_carga_dipolo_sin_fermi.md](analysis_fig1_carga_dipolo_sin_fermi.md)
  §5 — números documentados que R1 verifica.
- Scripts de verificación en
  `/var/folders/.../opencode/wigner_fix/` (verify_canonical.py,
  validate_exact.py, aggregate_stats.py, compare_grids.py, diag_chain.py).

## Notas Adicionales

- Observación colateral (pre-existente, no relacionada): `dg_integrals`
  devuelve NaN para algunos pares alto-l (p.ej. (25,26) con k=1, R=800) fuera
  de las combinaciones que produce el ensamblado; conviene auditarlo aparte.
- Los valores m=0 de la rejilla legacy vs nuevo difieren ≤ 1.44e−14 relativo;
  `gaunt` muestreada (10955 puntos) difiere ≤ 8.1e−16 absoluto.
