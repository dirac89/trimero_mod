# Pseudopotencial de Fermi para Rb*-KRb: remapeo k(R) e integración

**Fecha**: 2026-08-19
**Autor**: Javier Aguilera
**Relevancia**: Cierra `H = H_a + H_mol + V_Fermi` sobre `CoupledBasis`. Primer
punto de la migración cuyo espectro es potencialmente comparable con figuras de
González-Férez 2015.

## Resumen

Se implementó `trimero/hamiltonians/fermi_krb.py` con el pseudopotencial de Fermi
s+p sobre la base acoplada, incluyendo el remapeo semiclásico en `k(R)` de las
tablas `rvsAS.dat`/`rvsAP.dat` (opción (i), decidida por el usuario). 7 tests
nuevos; 37/37 en el repo, sin regresiones.

## 1. Por qué NO se reutilizó `fermi.FermiPotentials`

`mathlib.special.Spherical` es `lpmv` **crudo**: le falta la normalización
`√((2l+1)/4π·(l−m)!/(l+m)!)`. Sus términos `Vs` y `VpA` van sin normalizar
mientras `VpB`/`VpC` (vía `DOlm`/`DPhilm`, que sí la llevan) van normalizados.
Es internamente inconsistente, el factor depende de `l` y `trimer.py` suma
`Vsp()` sin corregirlo fuera. Además su lógica de ramas `l<=2` está atada a las
funciones de onda del manifold n=35 de 2016. Corregirla cambiaría `trimer.py`,
congelado por golden files y fuera de alcance. **Hallazgo reportado, no
arreglado.**

## 2. Física implementada

Con `Z ∥ R⃗`, el perturbador está en `θ=0`, donde `Y_lm(0,φ) = √((2l+1)/4π)δ_{m0}`
y `∇ψ_{lm}|_{θ=0} ≠ 0` sólo para `|m| ≤ 1`:

```
⟨l₁ 0|V_s|l₂ 0⟩   = (A_s/2)·R_{l₁}R_{l₂}·√((2l₁+1)(2l₂+1))
⟨l₁ 0|V_p|l₂ 0⟩   = (3/2)A_p·R'_{l₁}R'_{l₂}·√((2l₁+1)(2l₂+1))
⟨l₁±1|V_p|l₂±1⟩   = (3/4)A_p·R_{l₁}R_{l₂}/R²·√(l₁(l₁+1)l₂(l₂+1)(2l₁+1)(2l₂+1))
```

`V` es puramente electrónico: diagonal en `(N, M_N)` y en `m_l`, con `|m_l| ≤ 1`.
Conserva `M_J` **por ausencia de acoplamiento al rotor**, no por una regla de
selección de momento angular como `H_mol`. Mezcla `l` fuertemente: ahí está la
ligadura ULRM.

## 3. Aproximaciones explícitas

1. **KRb como perturbador puntual**: un único par `(A_s, A_p)`, sin estructura
   interna. El paper de 2015 no discute su validez para un dímero polar.
2. **Remapeo k(R)**: asume que la columna 2 de las tablas es `A(k)` calculada
   sólo con la energía semiclásica de n=35, sin dependencia adicional en `l`.
   Procedencia aún sin verificar. Sustituible por la opción (ii) reemplazando
   `remap_R` por la identidad y apuntando `data_dir` a tablas nuevas.
3. **Energía media para elementos fuera de la diagonal**: entre estados de `n*`
   distinto (27s vs manifold, 1 % de diferencia) se usa `E = (E_i+E_j)/2`, que
   es simétrica en `i↔j` y preserva la simetría de la matriz.
4. **Radial del 27s**: sigue con `n_eff = 24` (limitación previa). El remapeo sí
   usa el `n*` exacto 23.868513.

## 4. Dominio del remapeo — hay casos fuera de rango

`R' = 1/(E_n + 1/R − E₃₅)`, tabla en `R' ∈ [111, 2448] a₀`:

| R [a₀] | k(R), n*=24 | R' | ¿en tabla? |
|---|---|---|---|
| 300 | 0.070218 | 348.0 | sí |
| 600 | 0.039965 | 828.7 | sí |
| 900 | 0.022048 | 1535.6 | sí |
| 1100 | 0.009059 | 2226.2 | sí |
| 1152 | — | 2450.0 | **no** (clásicamente prohibido) |
| 1500 | — | 4836.2 | **no** (clásicamente prohibido) |

**`R = 1500 a₀` está fuera de dominio**: supera el punto de retorno clásico
externo `2n² = 1152 a₀` de n=24, `k² < 0`, no hay longitud de dispersión
definida. `scattering()` lanza `ValueError` en lugar de extrapolar. Ventana
válida para n=24: `R ∈ [105.6, 1151.5] a₀`.

## 5. Resultados de los tests (7/7)

| Test | Resultado |
|---|---|
| F0 dR/dr analítica vs diferencias finitas | peor rel. `8.1e-09` |
| F1 dominio del remapeo | tabla arriba; `R=1500` lanza `ValueError` ✓ |
| **F2 formas cerradas vs ψ,∇ψ numéricos en (0,0,R)** | **peor rel. `3.8e-07`**; converge con `h`: `2.4e-05 → 3.8e-07` |
| F3c conservación de `M_J` | 4800 pares con `M_J` distinto: `0.0` exacto; `Δl≠0` no nulo ✓; `|m_l|=2` nulo ✓ |
| F3b hermiticidad | 164836 pares del soporte: `max|V_ij−V_ji|/max|V| = 4e-16`; `‖H−Hᵀ‖/‖H‖ = 5.5e-16` |
| F3a límite `A_s=A_p=0` | `max|ΔE| = 0.0` exacto; control positivo 18.665 GHz |
| F3d magnitud | ver abajo |

Dos comprobaciones resultaron **vacuas** en la primera versión (el bloque
`M_J=−20` no contiene ningún estado con `|m_l| ≤ 1`, así que `max|V| = 0`). Se
corrigieron para barrer el soporte real de `V` y se añadió un assert que impide
que vuelvan a ser vacuas.

## 6. Magnitud: V_Fermi vs carga-dipolo

Desplazamiento máximo del espectro respecto del desacoplado, bloque `M_J=0`:

| R [a₀] | ΔE carga-dipolo [GHz] | ΔE Fermi [GHz] | cociente | A_p [a₀³] |
|---|---|---|---|---|
| 400 | 18.14 | 159.66 | 8.8 | — |
| 500 | 17.69 | 341.01 | 19.3 | — |
| 600 | 17.37 | **469.59** | **27.0** | −49327 |
| 800 | 17.46 | 41.79 | 2.4 | — |
| 900 | 18.25 | 23.62 | 1.3 | −12904 |
| 1100 | 16.34 | 14.72 | 0.90 | −19399 |

**V_Fermi domina por 9–27× en `R ≲ 800 a₀`**, como dice la literatura. Se
desvanece al acercarse al punto de retorno externo (1152 a₀) porque `V ∝ |ψ(R)|²`
y la densidad electrónica cae ahí — no es una señal de alarma sino la razón por
la que las ULRM se ligan en `R ≲ 2n²`. El pico en `R=600` sigue a `A_p`, que
atraviesa la resonancia de forma p (`−49327 a₀³`).

## 7. Archivos

- `src/trimero/hamiltonians/fermi_krb.py` — `ScatteringLengths`, `FermiPseudopotential`, radiales analíticas.
- `tests/hamiltonians/test_fermi_krb.py` — 7 tests con referencia por ψ/∇ψ numéricos.
- `scripts/run_full_with_fermi.py` — espectro completo.

**Diseño**: módulo nuevo en `hamiltonians/`, no extensión de `fermi.py` (ver §1).
`ScatteringLengths` se inyecta en `FermiPseudopotential`, igual que
`RydbergElectronField` en `ChargeDipoleHamiltonian`: permite desactivar el
acoplamiento (`enabled=False`) para el test del límite nulo y sustituir el
remapeo sin tocar la física. El ensamblado `H_a+H_mol+V` se hace en el script,
a la espera del `CompositeHamiltonian` que la sesión `refactor` tiene planeado.

## Referencias

- González-Férez, Sadeghpour & Schmelcher, *New J. Phys.* **17**, 013021 (2015).
- `docs/analysis_procedencia_rvsAS_rvsAP.md` — procedencia de las tablas (abierta).
