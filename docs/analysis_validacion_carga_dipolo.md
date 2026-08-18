# Validación del acoplamiento carga-dipolo (término del ion Rb⁺)

**Fecha**: 2026-08-18
**Autor**: Javier Aguilera
**Relevancia**: Primer bloque de `H_mol = B·N² - d·F_ryd(R,r)` (Ec. 3-4 de
González-Férez 2015) implementado sobre la base `CoupledBasis`. Fija el convenio
de signos y las reglas de selección antes de añadir el término, mucho más
delicado, del campo del electrón Rydberg.

## Resumen

Se implementó `src/charge_dipole.py` con **dos** de los tres ingredientes de
`H_mol`: el rotor rígido `B·N²` y el campo del ion Rb⁺, `-d·(e·R/R³)`. El
término del campo del electrón Rydberg (`e·(r-R)/|r-R|³`, expansión Ec. A.6-A.10)
queda **explícitamente fuera** de esta ronda.

Se escribieron 4 tests analíticos **antes** de la implementación
(`src/test_charge_dipole.py`). El proceso detectó un bug real de indexación en el
primer intento: el término `B·N²` se aplicaba sin exigir `δ_{l,l'} δ_{m_l,m_l'}`,
de modo que conectaba estados electrónicos ortogonales. Sin el test del límite
`R→∞` el error habría pasado desapercibido: produce una mezcla espuria de
`1.5×10⁻⁴ E_h` (≈ 10⁶ MHz), del orden de la propia escala Rydberg.

## Física implementada

Con el ion Rb⁺ en el origen, KRb en `R⃗` y eje de cuantización `Z ∥ R⃗` (el mismo
que define `M_J = m_l + M_N`), el campo del ion en la molécula es
`F_ion = R⃗/R³ = (1/R²)Ẑ` (e=1 en u.a.), luego

```
-d⃗·F_ion = -(d/R²)·cos(θ_d)
```

Elemento de matriz en la base `{(l, m_l, N, M_N)}`:

```
⟨l' m_l' N' M_N'| -d·F_ion |l m_l N M_N⟩
    = -δ_{l l'} δ_{m_l m_l'} (d/R²) ⟨N' M_N'|cosθ|N M_N⟩

⟨N-1 M|cosθ|N M⟩ = ⟨N M|cosθ|N-1 M⟩ = sqrt( (N²-M²) / ((2N-1)(2N+1)) )
```

Se codifica con `N_max = max(N,N')` para que ambas ramas sean la **misma**
expresión: la simetría de la matriz es exacta bit a bit, no aproximada
(`||H - H.T||_F = 0.0` exactamente, `np.array_equal(H, H.T) == True`).

Reglas de selección resultantes: `Δl = 0`, `Δm_l = 0`, `ΔM_N = 0`, `ΔN = ±1`,
y por tanto `ΔM_J = 0`. El término del ion es puramente axial, así que **no**
mezcla `M_N`; el término del electrón Rydberg sí lo hará (tiene componentes
perpendiculares a `Ẑ`), conservando `M_J` pero no `M_N`.

## Parámetros

| Parámetro | Valor | En u.a. | Fuente |
|---|---|---|---|
| `B(KRb)` | 1.114 GHz | 1.693090449×10⁻⁷ E_h | Ni et al. 2009 |
| `d(KRb)` | 0.566 D | 0.2226815538 e·a₀ | Ni et al. 2008 |
| `E(n=24, l≥3)` | — | −8.680555556×10⁻⁴ E_h | `atom.Atom.E_Rb()` |
| `E(27s)` | — | −8.776457784×10⁻⁴ E_h | `atom.Atom.E_Rb()` |

El 27s queda 9.590×10⁻⁶ E_h ≈ 63.1 GHz **por debajo** del manifold n=24, lo que
es consistente con que el paper lo trate como estado vecino relevante.

## Resultados de los 4 tests analíticos

| Test | Criterio | Resultado medido |
|---|---|---|
| 3. `⟨N,M_N\|B·N²\|N,M_N⟩ = B·N(N+1)` | igualdad exacta | diferencia `0.0` para todos los 49 `(N,M_N)`; 2352 elementos fuera de diagonal, todos `0.0` |
| 4. `H_cd` no mezcla `M_J` | error o cero explícito | `matrix_element` lanza `ValueError`; 18000 pares con `M_J` distinto evaluados **sin** la comprobación dan `max\|elemento\| = 0.0` |
| 2. Hermiticidad a `R=1500 a₀` | `\|\|H-H.T\|\| ≈ 0` | `0.0` exacto; barrido completo de 40000 pares con `max\|H_ij - H_ji\| = 0.0` |
| 1. Límite `R→∞` (`R=10000 a₀`) | espectro = `E_Ryd + B·N(N+1)` | `Δmax = 4.881×10⁻¹² E_h` (tolerancia `10⁻¹⁰`) |

Cada test lleva un **control positivo** que impide que pase de forma vacua: se
verifica que a `R` grande el acoplamiento sigue siendo no nulo, que dentro de un
bloque `M_J` hay elementos no nulos, y que a `R=1500 a₀` el espectro **sí** se
aparta del desacoplado (62.2 MHz). Sin ellos, una implementación que devolviera
un acoplamiento idénticamente nulo pasaría los tests 1, 2 y 4.

## Escalado del desplazamiento: confirmación de segundo orden

Salida de `src/run_charge_dipole_block.py` (bloque `M_J=0`, dim 1016):

| R [a₀] | `max\|E − (E_Ryd + B·N(N+1))\|` [E_h] | shift·R⁴ [E_h·a₀⁴] |
|---|---|---|
| 500 | 4.307274×10⁻⁷ | 2.692×10⁴ |
| 1000 | 4.459249×10⁻⁸ | 4.459×10⁴ |
| 1500 | 9.449584×10⁻⁹ | 4.784×10⁴ |
| 3000 | 6.018476×10⁻¹⁰ | 4.875×10⁴ |
| 10000 | 4.881266×10⁻¹² | 4.881×10⁴ |

El producto `shift·R⁴` converge a ≈ 4.881×10⁴, es decir el desplazamiento va
como `R⁻⁴`. Es exactamente lo esperado para teoría de perturbaciones de **segundo
orden** en un acoplamiento que cae como `R⁻²` (no hay término de primer orden
porque `cosθ` no tiene elementos diagonales en `N`). La desviación a `R=500 a₀`
señala la entrada de órdenes superiores. Este escalado es una comprobación
independiente del convenio `d/R²` que ningún test unitario impone directamente.

Estimación analítica de control a `R=1500 a₀`: `V ≈ (d/R²)·⟨cosθ⟩ = 5.71×10⁻⁸ E_h`,
denominador `ΔE = 2B = 3.39×10⁻⁷ E_h`, luego `V²/ΔE ≈ 9.6×10⁻⁹ E_h` frente al
`9.45×10⁻⁹ E_h` medido. Coinciden.

## Conclusiones y aplicación al proyecto

1. El convenio de signos y las reglas de selección del término del ion quedan
   fijados y verificados. La siguiente ronda (campo del electrón Rydberg) debe
   **reutilizar** este convenio, no redefinirlo.
2. El bloque `M_J=0` (1016×1016) se construye en ~0.35 s y se diagonaliza en
   ~0.14 s: coherente con el benchmark previo, sin problema de coste.
3. `||H-H.T||` exactamente cero no es casualidad: depende de escribir el elemento
   de `cosθ` con `N_max = max(N,N')`. Conviene mantener ese patrón al añadir la
   expansión de las Ec. A.6-A.10.
4. **Pendiente**: los 4 tests validan estructura y límites asintóticos, no la
   magnitud absoluta frente al paper. No se ha comparado con ninguna figura de
   González-Férez 2015 todavía, y no se podrá hasta que el término del electrón
   Rydberg esté implementado (es el dominante a los `R` de interés).
5. **Efecto colateral**: se eliminó `sph_harm` de la lista de imports de
   `src/math_aux.py`. Era un import muerto (ningún uso en el repo) que SciPy 1.17
   ya no exporta y que impedía importar `atom.py`. Cambio de cero efecto sobre el
   comportamiento.

## Archivos

- `src/charge_dipole.py` — `ChargeDipoleHamiltonian`, `rydberg_diagonal`, constantes.
- `src/test_charge_dipole.py` — los 4 tests analíticos + controles positivos.
- `src/run_charge_dipole_block.py` — construcción/diagonalización del bloque y escalado en R.

## Referencias

- González-Férez, Sadeghpour & Schmelcher, *New J. Phys.* **17**, 013021 (2015). DOI: 10.1088/1367-2630/17/1/013021
- Ni et al., *Phys. Chem. Chem. Phys.* **11**, 9626 (2009) — B(KRb)
- Ni et al., *Science* **322**, 231 (2008) — d(KRb)
- CODATA 2018 (NIST): 1 E_h = 6.579683920502×10¹⁵ Hz; 1 D = 0.393430307 e·a₀
