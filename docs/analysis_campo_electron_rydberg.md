# Campo del electrón Rydberg en `F_ryd`: expansión multipolar y validación

**Fecha**: 2026-08-18
**Autor**: Javier Aguilera
**Relevancia**: Completa `H_mol = B·N² - d·F_ryd(R,r)` (Ec. 3-4 de González-Férez
2015) con el segundo término de `F_ryd`, el más delicado. Cierra la parte
carga-dipolo del Hamiltoniano; queda fuera el pseudopotencial de Fermi.

## Resumen

Se implementó el término `-d·(e·(r-R)/|r-R|³)` mediante la expansión completa en
armónicos esféricos con integrales radiales dependientes de `r ≶ R`
(Ec. A.6-A.10). Se validó con 9 tests, el central de los cuales compara la
expansión contra una **cuadratura 2D directa del campo crudo** `(r-R)/|r-R|³`,
sin expansión, sin símbolos 3j y sin separación radial/angular. Acuerdo relativo
de `1e-11` a `R = 1500` y `2500 a₀`.

Ningún test falló por un error de la implementación: los cuatro fallos de la
primera ejecución fueron errores en los **valores de referencia que escribí en
los tests** (ver §6).

## 1. Derivación implementada

Con el ion Rb⁺ en el origen, KRb en `R⃗` y `Z ∥ R⃗` (mismo convenio que la ronda
del ion, no redefinido):

```
F⃗_elec/e = ∇_R (1/|r⃗-R⃗|) = ẑ ∂f/∂R - ρ̂ (1/R) ∂f/∂γ
f = Σ_k g_k(r,R) P_k(cos γ),    g_k = r_<^k / r_>^{k+1}
```

con `γ` el ángulo entre `r⃗` y `R⃗` y `ρ̂` el unitario radial cilíndrico. Usando
`dP_k(cosγ)/dγ = -sinγ P_k'` y `P_k^1 = -√(1-x²) P_k'` (Condon-Shortley):

```
F_z = Σ_k (∂g_k/∂R) P_k(cosθ)
F_ρ = -(1/R) Σ_k g_k P_k^1(cosθ)
```

En componentes esféricas `F_{±1} = ∓(F_x ± iF_y)/√2`:

```
⟨l₁m₁|F_0   |l₂m₂⟩ = Σ_k Z^k · √(4π/(2k+1))         · C(l₁m₁|k, 0|l₂m₂)
⟨l₁m₁|F_{±1}|l₂m₂⟩ = (1/(√2 R)) Σ_k G^k
                        · √(4π k(k+1)/(2k+1))        · C(l₁m₁|k,±1|l₂m₂)

G^k = ∫ u₁u₂ g_k dr ,   Z^k = ∫ u₁u₂ (∂g_k/∂R) dr ,   u = r·R_{nl}
```

`C` = coeficiente de Gaunt `∫Y*_{l₁m₁} Y_{kq} Y_{l₂m₂} dΩ`.

Acoplamiento al rotor mediante producto escalar en componentes esféricas:

```
⟨i|-d⃗·F⃗|j⟩ = -d·√(4π/3) Σ_{μ=-1,0,1} (-1)^μ ⟨N'M'|Y_{1μ}|N M⟩ · ⟨l'm'|F_{-μ}|l m⟩
```

**Reglas de selección resultantes**: `ΔM_N = +μ`, `Δm_l = -μ` ⇒ `ΔM_J = 0` con
`M_N` NO conservado; `ΔN = ±1` (d⃗ es rango 1 sobre el rotor); `Δl` libre dentro
del triángulo de Gaunt. Contraste con el término del ion, que es estrictamente
diagonal en `(l, m_l)` y en `M_N`.

**La suma en k no es una truncación**: el Gaunt anula todo `k` fuera de
`|l₁-l₂| ≤ k ≤ l₁+l₂` con `l₁+l₂+k` par, así que cada elemento de matriz es una
suma **finita y exacta**. `k_max` sólo existe para los tests (k≤1 = límite
dipolar, k=0 = monopolo).

Las integrales radiales se escriben con factores acotados (`(r/R)^k ≤ 1` dentro,
`(R/r)^k ≤ 1` fuera), y `R` se inserta como nodo exacto de la malla para que el
vértice de `g_k` no introduzca error de media celda.

## 2. Comprobación de consistencia con la ronda anterior

Poniendo `F_μ = δ_{μ0}/R²` en la fórmula general se recupera
`-d·⟨cosθ⟩/R²`, exactamente el término del ion ya verificado. El test 0c lo
comprueba numéricamente: `⟨N'M'|cosθ|NM⟩` obtenido vía Gaunt/3j coincide con
`cos_theta_element` de la ronda 1 en los 2401 pares (`max|diff| = 2.05e-15`).
Además, `ChargeDipoleHamiltonian` sin `electron_field` reproduce **bit a bit** la
matriz de la ronda anterior (`np.array_equal = True`).

## 3. Resultados de los tests

| Test | Criterio | Resultado |
|---|---|---|
| 0a | 3j/Gaunt vs analítico | 647 tripletes vs fórmula cerrada de Racah, `max|diff| = 8.3e-15`; ortogonalidad `1 - 4e-14` |
| 0b | radiales | `⟨u|u⟩ = 1` a 1e-12; `⟨r⟩` exacto a `2e-14` rel. |
| 0c | cierre con ronda 1 | `max|diff| = 2.05e-15` |
| **6** | **expansión vs cuadratura 2D bruta** | **`1.1e-11` rel. a R=2500 y R=1500** |
| 3 | `M_J` conservado, `M_N` no | elemento `ΔM_N=+1, Δm_l=-1` vale `2.06e-10` (≠0); sólo-ion da 0; `ΔN=±2` da 0; 4800 pares con `M_J` distinto dan `0.0` exacto |
| 2 | hermiticidad, ambos términos | `‖H-Hᵀ‖/‖H‖ = 1.1e-16`; `build == build_reference` bit a bit |
| 4 | límite dipolar 1er orden | diferencia relativa decrece monótonamente con R en los 3 elementos |
| 5 | cancelación del monopolo | k=0 vale `-1/R²` a `3.6e-15` rel.; residuo cae como `1/R⁴` (deriva 1.6e-3) |
| 1 | límite `R→∞` | `Δmax = 2.2e-15 E_h` a `R=20000` (tol. `1e-10`) |

A `R = 900 a₀` la cuadratura bruta pierde precisión porque el campo crudo `1/D³`
casi diverge en `r≈R`. Esto **se comprueba, no se afirma**: al refinar los nodos
(`n_x` = 400→3200) la fuerza bruta converge monótonamente hacia el valor de la
expansión (`2.1e-3 → 2.8e-5`). El residuo es de la cuadratura, no de la
expansión, que es analítica en el corte `r ≶ R`.

## 4. ¿Qué término domina?

**El cociente de normas `‖V_ele‖/‖V_ion‖` (1.1-3.4) NO responde la pregunta**,
porque los dos términos se cancelan parcialmente: el monopolo del electrón anula
el campo del ion (el átomo Rydberg es neutro). La medida correcta es el
desplazamiento del espectro:

| R [a₀] | sólo ion [MHz] | ion + electrón [MHz] | cociente |
|---|---|---|---|
| 600 | 1655.58 | 17372.96 | **10.5** |
| 900 | 430.31 | 18254.84 | **42.4** |
| 1200 | 147.77 | 11017.90 | **74.6** |
| 1500 | 62.18 | 1162.10 | **18.7** |
| 2000 | 19.94 | 81.36 | 4.1 |
| 3000 | 3.96 | 3.28 | 0.83 |
| 12000 | 0.0155 | 0.0003 | 0.02 |

**El término del electrón domina por un factor 10-75 en la región relevante**
(`R ≲ 1500 a₀`), lo que confirma lo que dice la literatura.

El vuelco del cociente por debajo de 1 para `R ≳ 3000 a₀` **no** significa que el
electrón deje de importar: significa que el modelo "sólo ion" es no físico a esas
distancias (es una carga desnuda). Al incluir el electrón, el campo neto pasa de
monopolar `1/R²` a dipolar `1/R³`, y el desplazamiento de segundo orden cae como
`~1/R⁶` (exponente medido entre R=6000 y 12000: **6.4**) en lugar de `1/R⁴`.

**Origen físico del efecto grande**: el acoplamiento del ion es estrictamente
diagonal en `l`; el del electrón conecta `l` distintos **dentro del manifold
cuasi-degenerado n=24**. A `R=900`, `max|V_ele|` con `l'≠l` vale `3.09e-7 E_h`,
un canal que el ion no tiene en absoluto.

⚠️ `R = 1500 a₀` está **más allá** del punto de retorno clásico externo
`2n² = 1152 a₀` para n=24: la molécula queda fuera de la órbita Rydberg. La
región ULRM representativa para este manifold es `R ≈ 600-1200 a₀`, donde la
dominancia del electrón es mayor (42-75×).

## 5. Limitaciones abiertas

1. **Falta el pseudopotencial de Fermi** `V(r)` (dispersión s y p del electrón
   Rydberg con KRb, `fermi_potentials.py`). Es el mecanismo de ligadura de las
   ULRM. Los números de arriba son **sólo** la parte carga-dipolo.
2. **Función radial del 27s aproximada**: se usa una hidrogenoide con
   `n_eff = 24` en lugar de la función de Coulomb con `n* = 27 - μ_s = 23.869`
   (error ~1.1 % en extensión radial; afecta a 7 de los 1016 estados del bloque).
   `scipy.special.hyperu` da NaN para `l=0` con `n*` no entero. La energía sí usa
   el defecto cuántico exacto.
3. **Sin comparación con el paper**: los tests validan estructura, límites
   asintóticos y consistencia interna, no magnitudes absolutas frente a ninguna
   figura de González-Férez 2015. Esa comparación requiere el punto 1.
4. `build` del bloque M_J=0 tarda ~6.3 s la primera vez (poblado de cachés) y
   ~0.34 s después para el mismo R.

## 6. Errores detectados por el proceso

Los 4 fallos de la primera ejecución fueron **errores en mis valores de
referencia**, no en la implementación — el propio proceso de escribir los tests
primero los expuso:

| Fallo | Causa |
|---|---|
| `(1 1 0;000) = 0` | `j₁+j₂+j₃ = 2` es **par**; el símbolo vale `-1/√3`. El caso de suma impar es `(1 1 1;000)` |
| `(2 1 1;000) = -√(2/15)` | es permutación **cíclica** de `(1 1 2;000)`, luego `+√(2/15)` |
| `(l 1 l+1;000)` | el signo es `(-1)^{l+1}`, no `-` constante |
| test `ΔN=±2` | el par que escribí tenía `ΔN=+1` |
| tolerancia del test 5 | la cancelación es exacta sólo en `k=0`; el residuo cuadrupolar es 3.9 % del ion a R=3000, no ≪1e-3 |

La regla de suma de ortogonalidad de los 3j que escribí inicialmente también era
incorrecta (sumaba sobre `m₃`, dando `2j₃+1` en vez de 1).

## 7. Archivos

- `src/angular_algebra.py` — `wigner_3j` (Racah con log-factoriales), `gaunt`.
- `src/rydberg_radial.py` — `RadialBasis`: malla uniforme en `√r`, `u_{nl}`, `G^k`, `Z^k`.
- `src/charge_dipole.py` — `RydbergElectronField` (nueva) + `ChargeDipoleHamiltonian(electron_field=…)`.
- `src/test_rydberg_field.py` — 9 tests, incluida la referencia por cuadratura 2D.
- `src/run_full_hamiltonian_scan.py` — escaneo en R y comparación de magnitudes.

**Decisión de diseño**: el campo del electrón va en una clase separada
(`RydbergElectronField`) que se inyecta opcionalmente en
`ChargeDipoleHamiltonian`, en vez de fundirse con ella. Motivos: (a) el término
del ion es analítico y sin estado, el del electrón necesita malla radial,
funciones de onda y cachés por R; (b) con `electron_field=None` el Hamiltoniano
verificado en la ronda anterior queda intacto bit a bit, lo que da una regresión
comprobable; (c) permite ejercitar el término nuevo en aislamiento en los tests.

## Referencias

- González-Férez, Sadeghpour & Schmelcher, *New J. Phys.* **17**, 013021 (2015). DOI: 10.1088/1367-2630/17/1/013021
- Ni et al., *PCCP* **11**, 9626 (2009) — B(KRb); Ni et al., *Science* **322**, 231 (2008) — d(KRb)
