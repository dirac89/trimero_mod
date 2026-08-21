# Regeneración del golden fig1 (Opción A) y Fases 1-2 del híbrido Rb*-Rb-RbCs

**Fecha**: 2026-08-21
**Autor**: Javier Aguilera
**Relevancia**: Cierra la ronda del arreglo de `wigner_3j` (`analysis_wigner3j_orden_canonico.md`): regenera el golden afectado, verifica que los números publicados siguen siendo ciertos, descarta contaminación por NaN en producción y produce la primera física del sistema híbrido.
**Tipo**: analysis

## Resumen

Tras el arreglo de `wigner_3j` (orden canónico; hermiticidad 2.8e−13 → 3.0e−19
relativo), el test de regresión `test_r2` fallaba su rtol=1e−12 con desviaciones
de hasta 6.34e−12 relativo. Se ejecutó la **Opción A**: regenerar el golden
`plots/rb_krb_polar/fig1_ad_MJ0_n25.npz` con el wigner corregido. El cambio es
**mejora de precisión verificada contra aritmética racional exacta** (error total
5.405e−11 → 4.846e−11 frente al Racah en `fractions.Fraction`), **no una
corrección de física**: la tabla de los 7 valores de referencia se mueve como
mucho 1.8e−11 GHz en la cuarta-decimosexta cifra.

Además: (i) se reprodujeron **todos** los números de §5 de
`analysis_fig1_carga_dipolo_sin_fermi.md` con el golden nuevo; (ii) se diagnosticó
el NaN colateral de `dg_integrals` —es uso fuera de dominio de la API, no un bug
de producción—; (iii) se retomaron las Fases 1-2 del híbrido Rb*-Rb-RbCs con la
primera producción física (curvas E0(R2) a R1 ∈ {600, 900, 1100} a₀).

## Palabras Clave

- golden file
- regresión numérica
- orden canónico Wigner 3j
- hermiticidad
- sistema híbrido carga-dipolo + Fermi

---

## 1. Qué cambió exactamente

| elemento | antes | después |
|---|---|---|
| `wigner_3j` | evaluación dependiente del orden de entrada (asimetría ~1e−15..1e−13) | una única evaluación Racah para las 12 formas equivalentes |
| hermiticidad L-c n=35 | ‖H−Hᵀ‖_F = 2.8e−13 Eh | 2.115e−21 Eh (3.0e−19 rel) |
| golden `fig1_ad_MJ0_n25.npz` | precisión vieja | **regenerado** (243/281 puntos difieren bit a bit) |
| `test_regression_fig1.py` | R2 FALLABA (6.34e−12 > rtol=1e−12) | **2 passed in 11.33s** |

**Naturaleza del cambio**: mejora de redondeo. Validación contra el Racah racional
exacto (muestra de 300/3638 símbolos que cambian): error total legacy 5.405e−11,
nuevo 4.846e−11 (**1.12× mejor**); peor error relativo nuevo 2.099e−10. Caso a
caso, 166 mejoran y 134 empeoran (redondeo aleatorio, sin sesgo).

### Tabla de los 7 valores de referencia (old → new)

| R (a₀) | E_old (GHz) | E_new (GHz) | ΔE (GHz) |
|---|---|---|---|
| 400 | −23.100158642148 | −23.100158642151 | −2.4e−12 |
| 600 | …131768 | …131778 | +1.0e−11 |
| 800 | …574598 | …574616 | +1.8e−11 |
| 1000 | …440161 | …440159 | −2.0e−12 |
| 1200 | …655673 | …655672 | −1.0e−12 |
| 1500 | …245243 | …245245 | +2.0e−12 |
| 1800 | …270502 | …270500 | −2.0e−12 |

Máximo |ΔE| = 2.425e−11 GHz; máximo relativo 6.339e−12. **K idéntico bit a bit en
los 281 puntos** (el carácter no cambia); max|ΔW| = 3.93e−14; spectrum
max|Δ| = 2.639e−11 GHz.

Backup del golden viejo:
`/var/folders/.../opencode/wigner_fix/fig1_ad_MJ0_n25_OLD.npz` (temporal).

### Golden M_J=1 (`fig1_ad_MJ1_n25.npz`) — regenerado también

Mismo procedimiento (mismo script con M_J=1), backup viejo en temp
`fig1_ad_MJ1_n25_OLD.npz`. Resultados: **232/281 puntos difieren bit a bit**,
max|ΔE| = 2.640e−11 GHz (rel máx 7.006e−12), **K idéntico en todos los puntos**,
max|ΔW| = 1.64e−13, spectrum max|Δ| = 4.780e−11 GHz. Misma naturaleza: mejora de
redondeo sin cambio de carácter ni de física.

| R (a₀) | E_old (GHz) | E_new (GHz) | ΔE (GHz) |
|---|---|---|---|
| 400 | −19.014987950137 | −19.014987950147 | +1.0e−11 |
| 600 | −16.705046860188 | −16.705046860207 | +1.9e−11 |
| 800 | −16.191196929774 | −16.191196929796 | +2.2e−11 |
| 1000 | −15.967210321051 | −15.967210321051 | 0 |
| 1200 | −12.677384883122 | −12.677384883120 | −2.1e−12 |
| 1500 | −1.429645045877 | −1.429645045877 | 0 |
| 1800 | −0.276280619941 | −0.276280619939 | −1.4e−12 |

## 2. Verificación de los números documentados (§5 del fig1)

Script temporal `verify_doc_numbers.py` sobre el golden nuevo:

```
VEREDICTO GLOBAL: todos reproducidos
- pozo más profundo: -23.100159 GHz @ R=400 ✓ (doc: -23.100)
- mínimo interior:   -21.596890 GHz @ R=450 ✓
- E(1800):           -0.337614 GHz ✓ (doc: -0.338)
- umbral N=5:        -22.6437 GHz; curva encima salvo R<=405 (doc: ~420) ✓
- 8 mínimos locales: R/E/amplitudes (1.417,1.192,1.074,1.029,1.103,1.377,2.169)
  y separaciones (65,75,80,95,110,140,140) a₀ ✓
- cruces evitados:   -10 GHz @1295, -5 @1385, -2 @1510, -1 @1615 ✓
```

`docs/STATUS.md` cita −23.100 / −0.338 / 8 mínimos: sin cambios a la precisión
publicada, no requiere edición.

## 3. Diagnóstico del NaN en `dg_integrals`

Durante la verificación apareció un NaN al pedir integrales con l=25..34 sobre
una base por defecto. Causa raíz (`src/trimero/basis/radial.py`, `_hydrogenic_u`):

```python
gammaln(n-l) = gammaln(24-25) = gammaln(-1) = inf      # entero no positivo -> polo
eval_genlaguerre(-2, 51, rho) = 0                       # grado negativo
log_norm = inf                                          # exp(inf)*... -> NaN
u resultante: [nan nan nan nan nan]
```

Es decir: pedir un hidrogenoide con **l ≥ n** (estado ligado inexistente). No es
división por cero ni log de negativo en los integrales.

**Confirmación de no contaminación en producción** (`nan_diagnosis.py`):

| sistema | base real | l usados | u(l) finito | pares × R | no finitos |
|---|---|---|---|---|---|
| polar `BOPSystem(n=25)` | `RadialBasis(n_manifold=25, l_max=24)` | 0..24 (<25) | ✓ | 325 × 4 | **0** |
| híbrido `HybridNeutralPolar(n=35)` | `RadialBasis(n_manifold=35, l_max=34)` | 0..34 (<35) | ✓ | 630 × 6 | **0** |
| neutro (trímero) | `RadialBasis(n_manifold=35, l_max=34)` | sólo funciones radiales Fermi | ✓ | — | n/a |

Ambos sistemas construyen la base con `l_max = n_manifold − 1`
(`bop_system.py:113`, `hybrid_system.py:143`, `linear_trimer.py:221`), así que
ningún par del ensamblaje real puede caer fuera de dominio. El NaN exige llamar a
la API con argumentos que ningún sistema usa.

## 4. Fases 1-2 del híbrido: primera producción física

Sistema: `HybridNeutralPolar` — Rb*(n=35) + perturbador neutro Rb(5s) en θ=π
(V_Fermi^π con fase de paridad) + RbCs polar en θ=0 (rotor B=490.17 MHz,
d=1.225 D en campo del ion y del electrón).

```
H(R1,R2) = H_A(n=35) + H_mol^KRb(R2) + V_Fermi^π(R1)
```

Los 12 tests de límite pasan (construcción, θ→π, conservación de M_J, límites
polar puro / neutro puro). Script de producción nuevo:
`scripts/compute_hybrid_curves.py` — la matriz de Fermi sólo depende de R1, así
que se evalúa **una vez por R1** y se reutiliza en todo el barrido R2 (mismo
orden de sumas que `hamiltonian()` ⇒ números idénticos al método).

Parámetros: M_J=0, N_max=2 (dim 307), R2 ∈ [500,1500] a₀ paso 25 (41 puntos),
keep=40 autovalores.

Salida real:

```
sanidad ||H-H^T||_F en (R1=600, R2=500): 2.976e-21 Eh

--- R1 = 600 a0 ---  E0: [-110.9683, -102.7940] GHz; mínimo global -110.9683 @ R2=500
                     mínimos locales: 0     (41 puntos, 59.9 s primera pasada)
--- R1 = 900 a0 ---  E0: [-111.1080, -102.9517] GHz; mínimo global -111.1080 @ R2=500
                     mínimos locales: 0     (41 puntos, 2.2 s, caché caliente)
--- R1 =1100 a0 ---  E0: [-111.1022, -102.9174] GHz; mínimo global -111.1022 @ R2=500
                     mínimos locales: 0     (41 puntos, 2.2 s)
```

Lectura física:

- Curvas suaves y monótonas en [500,1500] a₀, sin mínimos locales: la ligadura
  viene del canal carga-dipolo (polar) y decrece al abrir R2.
- Dependencia en R1 débil (~0.14 GHz entre R1=600 y 900 en R2=500): el
  pseudopotencial π está al otro lado del core y su efecto decae con R1.
- El extremo R2→1500 converge hacia la cota polar pura (~−102.8 GHz), coherente
  con el test de límite L-b (E[0] = −102.54 GHz en su configuración).
- Datos: `plots/hybrid_neutral_polar/hybrid_curves_R1{600,900,1100}_n35.npz`;
  figura `hybrid_curves_MJ0_n35.png`.

## Conclusiones y Aplicación al Proyecto

1. La Opción A estaba justificada: el golden nuevo es objetivamente más preciso
   (validación racional exacta) y toda la suite queda verde con rtol=1e−12 intacto.
2. Los números publicados en docs siguen siendo ciertos tras la regeneración.
3. El NaN de `dg_integrals` es un asunto cerrado: uso fuera de dominio, cero
   impacto en producción.
4. El híbrido tiene infraestructura completa (sistema + 12 tests + script de
   producción) y sus primeras curvas. Siguientes pasos naturales: convergencia en
   N_max (>2), barridos con M_J≠0, y explorar R2 < 500 a₀ donde pueden aparecer
   estructuras.

## Referencias

- [`analysis_wigner3j_orden_canonico.md`](analysis_wigner3j_orden_canonico.md) — causa raíz y arreglo de `wigner_3j` (ronda anterior).
- [`analysis_fig1_carga_dipolo_sin_fermi.md`](analysis_fig1_carga_dipolo_sin_fermi.md) — §5, números verificados en §2 de este documento.
- Aguilera-Fernández et al., J. Phys. B **49**, 124002 (2016) — el trímero neutro; base del módulo híbrido.
