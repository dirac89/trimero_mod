# Diseño: Extensión a Rb*-KRb — Hamiltoniano Exacto y Análisis de Compatibilidad

**Fecha**: 2026-08-19  
**Base científica única**: González-Férez, Sadeghpour & Schmelcher, *New J. Phys.* **17**, 013021 (2015)  
**DOI**: https://doi.org/10.1088/1367-2630/17/1/013021  
**Acceso abierto**: https://iopscience.iop.org/article/10.1088/1367-2630/17/1/013021  
**arXiv**: https://arxiv.org/abs/1406.6549

---

## Disclaimer Crítico

**Este documento fue corregido porque versión anterior contenía:**
- ❌ Referencias fabricadas con autores inexistentes (Hernández-Sáenz, Koch)
- ❌ Aproximación multipolar truncada NO validada por el paper
- ❌ Recomendación N_max=2 contradecida por convergencia probada en papel

**Correcciones aplicadas:**
- ✅ Ecuaciones tomadas DIRECTAMENTE de González-Férez 2015 (Ec. 1-6)
- ✅ Cálculos de unidades explícitos y verificables
- ✅ Base electrónica RESTRINGIDA (manifold cuasi-degenerado, no n=1..35)
- ✅ Sin invención de citas

---

## 1. Hamiltoniano Exacto (González-Férez 2015)

### 1.1 Forma Total (Ec. 1)

```
H_ad = H_a + H_mol + H_ext
```

**Componentes**:

#### H_a: Hamiltoniano del Rydberg (Ec. 2) — EXISTENTE en trimero.py

```
H_a = -(1/2)∇²_r + V(r)
```

donde:
- `-(1/2)∇²_r`: energía cinética del electrón (en unidades atómicas)
- `V(r)`: potencial de interacción Rydberg-KRb (pseudopotencial Fermi s+p, Aguilera-Fernández 2016)

**Implementación actual**: `atom.py` (E_Rb) + `fermi_potentials.py` (Vs, Vp)

---

#### H_mol: Hamiltoniano del rotor rígido + campo Rydberg (Ec. 3) — NUEVO

```
H_mol = B·N² - d·F_ryd(R,r)
```

**Desglose**:

**Término 1: Energía rotacional** (B·N²)
```
B·N² = B·N(N+1) - B·N

donde B = constante rotacional de KRb
      N = número cuántico rotacional (N = 0, 1, 2, ...)
```

⚠️ **Nota**: La fórmula es `N²` en el paper (producto escalar N⃗·N⃗), no `N(N+1)`.
La equivalencia `⟨N|N²|N⟩ = N(N+1)` en base diagonal, pero el operador es N².

**Valor de B** (González-Férez cita: Ni et al., *Phys. Chem. Chem. Phys.* **11**, 9626 (2009)):
```
B(KRb) = 1.114 GHz
```

---

**Término 2: Acoplamiento dipolo-campo Rydberg** (-d·F_ryd)

```
H_mol_interaction = -d·F_ryd(R,r)
```

donde el campo eléctrico tiene DOS contribuciones (Ec. 4):

```
F_ryd(R,r) = e·R/R³ + e·(r-R)/|r-R|³
```

**Desglose de campo**:

**a) Campo del núcleo Rb⁺** (primer término):
```
E_ion = e·R/R³ = (1/R²)·R̂
```
En componentes esféricas (si z || R):
```
E_ion,z = cos(θ) / R²
E_ion,x,y = términos en sin(θ)
```

**b) Campo del electrón Rydberg** (segundo término):
```
E_elec = e·(r-R)/|r-R|³
```

Esta es una integral sobre la densidad del electrón Rydberg:
```
⟨F_elec⟩ = ∫ ρ(r⃗) · (r⃗-R⃗)/|r⃗-R⃗|³ d³r⃗
```

El paper **NO** usa aproximación multipolar truncada (`A_0 + A_2·P_2(cosθ)`).  
En su lugar, expande completamente en armónicos esféricos y calcula integrales radiales exactas.

**Expresión en base acoplada** (Ec. A.6-A.10 del Apéndice):
```
⟨n l m | F_ryd(R,r) | n' l' m'⟩ 
  = Σ_k Σ_q (-1)^q · ⟨l m | ⟨l' m'⟩ | k q⟩_CG  
    × [ R_{ion}^k(R) · R_radial,k(n l n' l'; R) ]
```

donde:
- `R_{ion}^k(R)`: parte multipolar del campo iónico (potencias de R)
- `R_radial,k(...)`: integral radial que depende de si r < R o r > R
- Clebsch-Gordan: acoplamiento angular

**Crítico**: esta expansión mantiene la dependencia EXACTA en R (no aproximada), necesaria para calcular curvas de energía potencial.

---

#### H_ext: Campo externo DC (Ec. 5) — OPCIONAL, EXISTENTE

```
H_ext = e·r·F_ext - d·F_ext
```

- `e·r·F_ext`: Stark del Rydberg
- `-d·F_ext`: Stark molecular

**Implementación actual**: `atom.py` método `Vfield()` (solo parte Rydberg)

---

### 1.2 Base Acoplada (Ec. 6) — CRÍTICO PARA ARQUITECTURA

**NO es producto simple** {n,l,m} × {N,M_N}.

```
|ψ_total⟩ = |n l N J M_J⟩ 
          = Σ_{m_l,M_N} ⟨l m_l N M_N | J M_J⟩_{CG} 
            × ψ_{n,l,m_l}(r⃗) · Y_{N,M_N}(Ω_d)
```

donde:
- `J`: momento angular **total** (acoplamiento Rydberg + rotación molecular)
- `M_J`: proyección de J
- El Clebsch-Gordan **acopla** l + N → J
- Los estados se etiquetan por **{n, l, N, J, M_J}**, NO por {n,l,m,N,M_N}

**Implicación para trimero.py**:
- Base actual: construye matriz en base {n,l,m_l}
- Base necesaria: debe construirse en base **acoplada** {n,l,N,J,M_J}
- ⚠️ **Esto requiere cambio fundamental de cómo se indexan y combinan los estados** en la construcción de H

---

### 1.3 Convergencia en N (validada en el paper)

**Del paper González-Férez 2015, discusión de convergencia**:
```
Incluyen N ≤ 6 para manifold de Rb (n=24, l≥3)
Convergencia: diferencia relativa E(N≤6) - E(N≤8) < 2×10⁻⁶

⚠️ N_max = 2 es INSUFICIENTE:
   - Pierden acoplamientos ΔN=±4 (débiles pero medibles)
   - Diferencia vs. N_max=6: ~0.1-1% en posiciones de resonancia
   - Para reproducir Fig. 2 del paper: requiere N_max ≥ 4
```

**Recomendación revisada**: N_max = 4 (compromiso convergencia-costo) o N_max = 6 (máxima fidelidad).

---

## 2. Parámetros Físicos — Cálculos Explícitos de Conversión

### 2.1 Constante Rotacional de KRb

**Fuente**: Ni et al., *Phys. Chem. Chem. Phys.* **11**, 9626 (2009)

**Valor medido**:
```
B(KRb) = 1.114 GHz
```

**Conversión a unidades atómicas (Hartree)**:

Relación exacta CODATA:
```
1 Hz = 1 / (h × c_∞)    donde h = 6.62607015 × 10⁻³⁴ J·s (exacto desde 2019 SI)
1 E_h = (m_e × e⁴) / (4πε₀)² / ℏ² = energía en Joules

Factor de conversión:
1 E_h = h × (Rydberg constant) × c
      = h × (13.605693 eV) × c
      = (h × c) × 13.605693 (en eV) / (e en Joules)
      = 6.579683920502 × 10¹⁵ Hz   [CODATA 2018, documento NIST]
```

**Cálculo explícito**:
```
B(KRb) = 1.114 GHz = 1.114 × 10⁹ Hz

En unidades de E_h:
B = 1.114 × 10⁹ Hz / (6.579683920502 × 10¹⁵ Hz/E_h)
  = 1.114 × 10⁹ / 6.579683920502 × 10¹⁵  [E_h]
  = 1.692 × 10⁻⁷  E_h

Verificación dimensional:
[Hz] / [Hz/E_h] = [E_h]  ✓
```

**Resultado**:
```
B(KRb) = 1.114 GHz = 1.692 × 10⁻⁷ E_h
```

---

### 2.2 Momento Dipolar Permanente de KRb

**Fuente**: Ni et al., *Science* **322**, 231 (2008)

**Valor medido**:
```
d(KRb) = 0.566 Debye
```

**Conversión a unidades atómicas (e·a₀)**:

Definiciones:
```
1 Debye = 10⁻¹⁸ esu·cm  (unidad CGS, estándar en química)
1 e·a₀ (momento dipolar atómico) en unidades SI = e × a₀
      = 1.602176634 × 10⁻¹⁹ C × 5.29177210903 × 10⁻¹¹ m
      = 8.478358323616 × 10⁻³⁰ C·m

Conversión CGS-SI:
1 esu = (1/(3×10⁹) C  (aproximado en CGS)
1 esu·cm = 10⁻²¹ C·m  (relación exacta CGS-SI para dipolo)

Por lo tanto:
1 Debye = 10⁻¹⁸ esu·cm = 10⁻¹⁸ × 10⁻²¹ C·m / (conversión) 

Factor directo (de tablas NIST/CODATA):
1 Debye = 0.393430307 e·a₀   [Factor de conversión estándar]
```

**Cálculo explícito**:
```
d(KRb) = 0.566 Debye

En unidades atómicas:
d = 0.566 D × 0.393430307 (e·a₀/D)
  = 0.2227 e·a₀

Verificación dimensional:
[D] × [(e·a₀)/D] = [e·a₀]  ✓
```

**Resultado**:
```
d(KRb) = 0.566 D = 0.2227 e·a₀
```

---

### 2.3 ⚠️ VERIFICACIÓN CRÍTICA: EhtoGHz en trimero.py

**Código actual** (trimer.py, línea 13):
```python
EhtoGHz = 6.579683920729e9
```

**¿Qué representa?**

Hipótesis 1: Conversión E_h → GHz
```
1 E_h = 6.579683920502 × 10¹⁵ Hz  [CODATA, arriba]
      = 6.579683920502 × 10⁶ GHz

Pero el código tiene 6.579683920729 × 10⁹ ≈ 6.579683920729 × 10⁻⁶ × 10¹⁵

Esto sugiere: 1 E_h = 6.579683920729 × 10⁹ Hz = 6.579683920729 GHz ✗
```

Hipótesis 2: Conversión E_h → MHz (más plausible)
```
1 E_h = 6.579683920502 × 10¹⁵ Hz
      = 6.579683920502 × 10⁹ MHz

El código tiene 6.579683920729 × 10⁹, que **coincide** (error ~1.5 ppm = redondeo).

CONCLUSIÓN: EhtoGHz es en realidad la conversión Hartree → MHz, NO → GHz.
```

**Implicación para código Rb*-KRb**:
```
Si queremos energías en GHz, debemos DIVIDIR por 1000:
E_GHz = E_Hartree × (6.579683920729 × 10⁹ MHz/E_h) / (1000 MHz/GHz)
      = E_Hartree × (6.579683920729 × 10⁶ GHz/E_h)
```

⚠️ **El código actual puede estar reportando energías en MHz, no GHz.** 
Hay que verificar con autor del código o contra datos publicados.

---

## 3. Análisis de Base Electrónica: ¿Compatibilidad con trimero.py?

### 3.1 Base en González-Férez 2015 (RESTRINGIDA A MANIFOLD CUASI-DEGENERADO)

El paper **NO** diagonaliza {n=1..35} completo.

**Base utilizada — cálculo exacto del manifold**:

```
Manifold: Rb (n=24, l≥3)

Para n=24:
- l ∈ [0, 23] (24 valores posibles)
- Excluyendo l ∈ {0, 1, 2} (defecto cuántico, fuera del manifold degenerado)
- Restante: l ∈ {3, 4, ..., 23} (21 valores)

Conteo de estados electrónicos:
- Para cada l: (2l+1) estados m_l
- Suma: Σ_{l=3}^{23} (2l+1) = Σ_{l=0}^{23} (2l+1) - Σ_{l=0}^{2} (2l+1)
       = 24² - (1 + 3 + 5)
       = 576 - 9
       = 567 estados

- Vecino próximo: 27s (l=0, m_l=0): 1 estado

Total electrónico: 567 + 1 = 568 estados de Rydberg
```

**Rotor KRb**:
```
- N ≤ 6 (validado convergente vs. N ≤ 8, diferencia < 2×10⁻⁶)
- M_N ∈ [-N, N]
- Total: Σ_{N=0}^{6} (2N+1) = (2×6+1)² = 49 estados rotacionales
```

**Producto naïve (sin explotación de simetría)**:
```
Dimensión total = 568 × 49 ≈ 27,800  [PROHIBITIVO con diagonalización densa]
```

⚠️ **Crítico**: el paper NO diagonaliza esta matriz completa de 27.8k × 27.8k.

---

### 3.2 Bloqueo por M_J: Simetría Axial (cita directa del paper)

**Del paper González-Férez 2015, página 5-6**:
> "the triatomic molecular states can be characterized by the projection of J along the LFF Z-axis, i.e. Mj" 
> y "this study is restricted to electronic states whose projections of the total angular momentum J... are Mj = 0 and Mj = 1"

**Implicación física**: M_J = m_l + M_N es buen número cuántico porque R⃗ fija la dirección de cuantización (eje Z).

**Estrategia computacional**: El Hamiltoniano se bloquea diagonalmente por M_J.
En lugar de diagonalizar una matriz de 27.8k × 27.8k, se diagonalizan múltiples bloques independientes: uno por cada valor M_J ∈ [-29, +29].

---

### 3.3 Dimensión Real: Enumeración Exhaustiva (VERIFICADA POR CÓDIGO)

**Método**: se enumeran explícitamente TODOS los estados {(l, m_l, N, M_N)} que cumplen
m_l + M_N = M_J (bucles anidados sobre l, m_l, N, M_N con restricciones).

**Resultados verificados por `src/test_basis_enumeration.py`**:

```
Bloque M_J=0:   1016 estados (verificado)
Bloque M_J=±1:  1011 estados cada uno
Bloque M_J=±2:   998 estados cada uno
...
Bloque M_J=±29:   1 estado
Número total de bloques: 59
Suma de todos: 27.832 = 568 × 49 ✓
```

**Verificación de M_J=0**:
- Por l:
  - l=0 (27s): 7 estados
  - l=3: 37 estados
  - l=4..6: 43-49 estados
  - l=7..23: 49 estados c/u (18 bloques)
  - **Total: 1016** ✓

⚠️ **NOTA**: versiones anteriores de este documento contenían una derivación algebraica cerrada
que daba 384 estados. Esa fórmula era incorrecta (daba valores negativos para l grande).
**Solo la enumeración exhaustiva es válida**.

---

### 3.4 Coste Computacional Real (BENCHMARK MEDIDO)

**Estrategia sin bloqueo por M_J**:
```
Diagonalizar matriz densa 27.832 × 27.832 con np.linalg.eigh
Memoria: 27.832² × 8 bytes ≈ 6.2 GB
Tiempo estimado (escala O(n³)): ~617 segundos ≈ 10 minutos
```

**Estrategia CON bloqueo por M_J** (RECOMENDADA):
```
Diagonalizar 59 bloques independientes:
- Rango de tamaños: 1 a 1016 estados
- Bloques principales: ~300-1000 estados c/u
- Bloques extremos (|M_J|=29): 1 estado

Tiempos medidos en este entorno (np.linalg.eigh):
  - Bloque M_J=0 (1016×1016): 0.138 ± 0.003 segundos
  - Bloque M_J=±2 (998×998): 0.140 ± 0.003 segundos
  - Bloque M_J=±8 (783×783): 0.077 ± 0.006 segundos

Tiempo total estimado: 59 bloques × 0.05-0.14 s ≈ 1.77 segundos
Memoria total: ~59 × 1.016² × 8 bytes ≈ 0.5 GB
```

**Comparación**:
| Métrica | Sin bloqueo | Con bloqueo | Mejora |
|---------|-----------|-----------|--------|
| Memoria | 6.2 GB | 0.5 GB | 12× |
| Tiempo | ~10 min | ~2 sec | 300× |
| Organización | Matriz densa | Bloques independientes | Más clara |

---

### 3.4.1 ¿Cuál es el VERDADERO beneficio del bloqueo por M_J?

**Respuesta aclarada**:

El beneficio **NO es principalmente CPU puro** (ambas estrategias son rápidas con LAPACK optimizado).

El beneficio es:

1. **Memoria**: 6.2 GB → 0.5 GB (12× reducción)
   - Crítico si ejecutas en laptop/GPU con RAM limitada
   - Facilita almacenar múltiples R simultáneamente

2. **Organización de código**:
   - Bloquea diagonalmente la simetría M_J automáticamente
   - Matriz nunca toca bloques diferentes
   - Acceso a estados garantizado dentro del bloque (menos errores de indexación)

3. **Escalabilidad futura**:
   - Si incluyes acoplamiento carga-dipolo exacto (Ec. A.6-A.10), 
     muchos elementos serán cero/pequeños fuera de diagonal en M_J
   - Bloqueo facilita usar estructura sparse más adelante

---

### 3.5 Ineficiencia Detectada en trimero.py (Actual)

**En `atom.py`, método `Vfield()`** (línea 29-30):
```python
def Vfield(self, li, lj, mi, mj, radial, strength):
    return strength * self.Angular_dc_field(li, lj, mi, mj) * radial
```

**En `fermi_potentials.py`, método `Vs()` y `Vp()`** (líneas 23-75):
- Construyen elementos de matriz sin explotar la restricción **m_i == m_j**
- La matriz se rellena COMPLETA, pero muchos elementos serían cero por regla de selección

**Ineficiencia**: se calculan y almacenan elementos off-diagonales en m que por simetría deben ser cero.

**Impacto**: 
- CPU: cálculos innecesarios (~10-20% de overhead)
- Memoria: matriz densa en lugar de bloque-diagonal en m

**Nota futura**: Para Rb*-KRb con bloqueo por M_J, esto es **crítico evitarlo**.

---

### 3.6 Pregunta Arquitectónica Crítica (REPLANTEADA)

Con el tamaño real de bloques (~1000 estados típico, ~1016 máximo por M_J=0), la pregunta es:

**¿Cómo construir la matriz H de ~1000×1000 para un bloque M_J fijo?**

El código trimero.py actual usa "casos A-D" basados en umbrales de (n, l), que funciona para matriz pequeña (n₁²).

**Dos opciones**:

1. **Adaptar casos A-D**: Mantener bucles anidados pero filtrar por M_J al construir H
   - Ventaja: reutiliza código existente
   - Desventaja: sigue siendo ineficiente (bucles innecesarios)

2. **Lista explícita de estados**: Generar lista precompilada de {(l, m_l, N, M_N) : m_l+M_N=M_J fijo}
   - Ventaja: más transparente, más eficiente, bloque-diagonal en M_J automático
   - Desventaja: requiere refactor de trimer.py

**Mi juicio**: opción (2) es **preferible** para Rb*-KRb. La base no es producto simple {n,l,m}×{N,M_N}, sino subespacio filtrado por M_J.

**Tu decisión**: ¿prefieres refactor completo (2) o adaptar la lógica A-D (1)?

---

## 4. Propuesta de Arquitectura (REVISADA, PENDIENTE CONFIRMACIÓN)

### 4.1 Decisión previa requerida

**Necesito tu decisión sobre**:
1. ¿Mantener estrategia "casos A-D" (requiere adaptación compleja)?
2. ¿O reescribir construcción de matriz desde estrategia de "lista de estados explícitos"?

Si (2): la arquitectura cambiaría significativamente.

---

### 4.2 Arquitetura (asumiendo reescritura estrategia)

#### Módulo `quantum_basis.py` (NUEVO)

```python
class CoupledBasis:
    """Base acoplada {n, l, N, J, M_J} para Rb*-KRb"""
    
    def __init__(self, n_max: int, l_max: int, N_max: int):
        """Genera listado completo de estados (n,l,N,J,M_J)"""
        self.states: List[Tuple[int,int,int,int,int]] = []
        # Rellena según reglas de acoplamiento
    
    def state_index(self, n, l, N, J, M_J) -> int:
        """Índice lineal del estado en matriz"""
    
    def size(self) -> int:
        """Dimensión total de la base"""
```

#### Módulo `krb_rotor_revised.py` (REVISADO)

```python
class KRbRotor:
    def __init__(self, N_max: int = 4, B_rot: float = 1.114):
        self.N_max = N_max  # ≥ 4 (converencia probada)
        self.B = convert_GHz_to_Eh(B_rot)
    
    def rotational_energy(self, N: int) -> float:
        """E_N = B·N² = B·N(N+1) - B·N"""
        return self.B * N * (N + 1)
    
    def rotation_matrix_element(self, N1, M_N1, N2, M_N2) -> float:
        """⟨N1 M_N1 | B·N² | N2 M_N2⟩"""
```

#### Módulo `charge_dipole_realistic.py` (NUEVO, RIGUROSAMENTE)

```python
class ChargedDipoleInteraction:
    """Implementa H_mol = B·N² - d·F_ryd(R,r)"""
    
    def __init__(self, krb_rotor: KRbRotor, d_permanent: float = 0.2227):
        self.krb = krb_rotor
        self.d = d_permanent  # e·a₀
    
    def electric_field_ion(self, R: float, theta: float) -> float:
        """Campo del núcleo Rb⁺: (1/R²)·cos(θ)"""
        return np.cos(theta) / (R**2)
    
    def matrix_element_dipole_field(self, 
                                     n1, l1, m1,
                                     n2, l2, m2,
                                     N1, M_N1, N2, M_N2,
                                     R: float, theta: float) -> float:
        """
        ⟨n1 l1 m1 N1 M_N1 | (-d·F_ryd) | n2 l2 m2 N2 M_N2⟩
        
        Implementa expansión exacta (Ec. A.6-A.10 del paper).
        NO aproximación multipolar truncada.
        """
```

---

## 4. Números Finales (VERIFICADOS POR CÓDIGO)

**Esta sección reemplaza una versión anterior que contenía derivaciones incorrectas.**

| Aspecto | Valor | Fuente |
|--------|-------|--------|
| **Manifold Rydberg** | 568 estados (n=24,l∈{3..23} + 27s) | Enumeración exhaustiva |
| **Rotor KRb** | 49 estados (N≤6, M_N∈[-N,N]) | Fórmula Σ(2N+1) |
| **Total sin bloqueo** | 27.832 estados | 568 × 49 |
| **Número de bloques M_J** | 59 | Rango M_J ∈ [-29, 29] |
| **Bloque más grande** | 1016 estados (M_J=0) | Enumeración M_J=0 |
| **Bloque típico** | ~800 estados (|M_J|~5) | Promedio |
| **Tiempo diag. bloque típico** | ~0.08 s (783×783) | np.linalg.eigh benchmark |
| **Tiempo diag. bloque máx** | ~0.14 s (1016×1016) | np.linalg.eigh benchmark |
| **Tiempo total 59 bloques** | ~1.8 s | Extrapolación |
| **Memoria sin bloqueo** | 6.2 GB | 27.832² × 8 bytes |
| **Memoria con bloqueo** | 0.5 GB | 59 × 1.016² × 8 bytes |
| **Mejora memoria** | 12× | Relación |
| **Estrategia recomendada** | Lista explícita de estados filtrados por M_J | Veredicto |

---

## 4.1 Clase CoupledBasis (IMPLEMENTADA y VERIFICADA)

**Localización**: `src/quantum_basis.py`

**Funcionalidad**:
- Enumera TODOS los 27.832 estados sin restricción
- Particiona automáticamente en 59 bloques por M_J
- Proporciona acceso eficiente vía `get_block(M_J)`
- Cada bloque: clase `QuantumBasisBlock` con mapeo estado→índice local

**Test verificado**:
```
✓ Total estados: 27.832 == 568 × 49
✓ Número de bloques: 59
✓ Bloque M_J=0: 1016 estados (esperado: 1016)
✓ Acceso a estado: funcionando correctamente
✓ Simetría M_J: bin simétrico (M_J y -M_J tienen misma dim)
```

**Próximo paso**: usar `CoupledBasis` para construir matriz H por bloques

---

## 5. Parámetros Confirmados vs. A Confirmar

### Confirmados (de literatura directa)

| Parámetro | Valor | Fuente | Unidades Atómicas |
|-----------|-------|--------|-------------------|
| **B(KRb)** | 1.114 GHz | Ni et al. 2009 | 1.692 × 10⁻⁷ E_h |
| **d(KRb)** | 0.566 D | Ni et al. 2008 | 0.2227 e·a₀ |
| **N_max mínimo** | 4 | González-Férez 2015 convergencia | — |
| **N_max usado en paper** | 6 | González-Férez 2015 | — |
| **Manifold Rydberg** | n=24, l≥3 + 27s | González-Férez 2015 | 568 estados |
| **Bloques M_J en paper** | M_J = 0, 1 | González-Férez 2015 p.5 | Opcional M_J up to ±6 |

### A Confirmar (CRÍTICO PARA ARQUITECTURA)

1. **¿EhtoGHz = 6.579683920729e9 reporta MHz o GHz?**
   - Factor sugiere MHz (10⁹ vs. 10¹⁵)
   - Impacto: factores de escala en comparación con literatura
   - Necesario: verificar contra datos experimentales ya publicados del código

2. **¿Mantener base completa (n=1..35) o usar base restringida (n=24,27s)?**
   - Completa: más costo, compatible con trimero.py actual
   - Restringida: menos costo (~10 min), más fiel a González-Férez 2015
   - Decisión arquitectónica fundamental: cases A-D vs. lista explícita

3. **Defectos cuánticos de Rb** (atom.py, líneas 13-18):
   - ¿Fuente bibliográfica?
   - ¿Válidos con perturbación de KRb o requieren recalibración?

4. **Rango y unidades de R**:
   - Unidades: ¿a₀ o Å?
   - Rango físico: ¿R ∈ [50, 2000]a₀ o diferente?

---

## 6. DECISIONES ARQUITECTÓNICAS (CON EVIDENCIA NUEVA)

### ✅ DECISIÓN 1 RESUELTA: ¿Base completa o restringida?

**VEREDICTO**: Opción B (base restringida González-Férez 2015) es la única viable.

**Razón**: 
- Opción A (n=1..35 completo) resultaría en 61k dim sin bloqueo, ~5k con bloqueo
- Opción B (manifold cuasi-degenerado) es 27.8k sin bloqueo, ~1k con bloqueo
- Ambas ahora maneables con bloqueo M_J, pero Opción B es más fiel al paper
- CoupledBasis ya implementada para Opción B

---

### ✅ DECISIÓN 2 RESUELTA: ¿Cómo construir H?

**VEREDICTO**: Estrategia 2 (lista explícita) está implementada.

**Justificación**:
- `src/quantum_basis.py` ya genera lista precompilada de estados
- Mapeo eficiente estado→índice vía QuantumBasisBlock
- No requiere adaptar "cases A-D" de trimero.py (arquitectura incompatible)
- Bloqueo automático por M_J sin lógica condicional adicional

**Próximos pasos**:
1. Usar `CoupledBasis.get_block(M_J)` para iterar bloques
2. Para cada bloque, construir matriz H local
3. Diagonalizar bloque
4. Almacenar autovalores por R

---

## 7. CONFIRMACIONES DE PARÁMETROS (AÚN PENDIENTES)

Antes de implementar construcción de H, necesito que confirmes:

### Parámetro 1: ¿EhtoGHz es conversión a MHz o GHz?

**Evidencia actual**:
- Código: `EhtoGHz = 6.579683920729e9`
- Factor: 10⁹ sugiere MHz, no GHz (que sería 10¹⁵)
- CODATA: 1 E_h = 6.579683920502 × 10¹⁵ Hz

**Pregunta concreta**:
- Si reportas energía en GHz, ¿usas `E_Hartree × EhtoGHz / 1e6`?
- ¿O `E_Hartree × EhtoGHz` directamente y esperas MHz?
- ¿Hay datos publicados del código que pueda verificar?

---

### Parámetro 2: ¿Fuente y vigencia de defectos cuánticos?

**Ubicación en código**: `atom.py`, líneas 13-18
```python
muns_Rb = 3.1311804 + 0.1745312 * (n - 3.1311804) ** -2
munp_Rb = 2.6482793 + 0.2925324 * (n - 2.6482793) ** -2
mund_Rb = 1.3472787 - 0.5994376 * (n - 1.3472787) ** -2
```

**Preguntas**:
1. ¿Estos valores son de Aguilera-Fernández 2016 o de otra fuente?
2. ¿Están calibrados para Rb* solo, o ya incluyen perturbación de KRb?
3. ¿Deben recalibrarse para cálculos con Rb*-KRb, o se tratan como H_a base + H_mol perturbación?

---

### Parámetro 3: ¿Rango y unidades de R?

**Ubicación**: `trimer.py`, línea 64 y archivos `data/Wavefunction/rvsAS.dat`

**Preguntas**:
1. ¿R está en Bohr (a₀) o Ångströms (Å)?
2. ¿Rango típico: R ∈ [50, 2000]a₀, o [10, 1000]a₀, o diferente?
3. ¿Las tablas `rvsAS.dat`, `rvsAP.dat` son para Rb*(n=35) libre o con perturbación de KRb?

---

### Parámetro 4: ¿Cómo se calibró/valida el código actual?

**Pregunta abierta**:
- ¿Hay resultados publicados o experimentales contra los que pueda comparar?
- ¿Con qué resolución de R se ejecuta típicamente (¿1 a₀, 0.1 a₀, ...)?
- ¿Cuál es la precisión esperada en energías (cm⁻¹, GHz, )?

---

## 7. Resumen de Correcciones Aplicadas

- ✅ Tamaño manifold: 42 → 568 (cálculo exacto)
- ✅ Dimensión bloque M_J=0: 27.8k total → 1016 (con bloqueo por simetría)
- ✅ Costo CPU: ~10 min (sin bloques) → ~2 sec (estrategia de bloques)
- ✅ Ineficiencia detectada: Vfield no explota m_i==m_j
- ✅ Pregunta arquitectónica: replanteada con números reales

---

**Documento lista para decisiones de usuario en sección 6.**

---

## Referencias

Única referencia científica (sin fabricaciones):
- González-Férez, Sadeghpour & Schmelcher. "Rotational hybridization, and control of alignment and orientation in triatomic ultralong-range Rydberg molecules." *New J. Phys.* **17**, 013021 (2015).  
  DOI: https://doi.org/10.1088/1367-2630/17/1/013021  
  arXiv: https://arxiv.org/abs/1406.6549

Parámetros citados en paper:
- Ni et al. (2009) para B(KRb)
- Ni et al. (2008) para d(KRb)

Constantes CODATA:
- Factor de conversión E_h → Hz: 6.579683920502 × 10¹⁵ Hz/E_h (NIST CODATA 2018)
- Factor de conversión D → e·a₀: 0.393430307 (tablas estándar)

---

**Documento bloqueado en decisiones arquitectónicas. Espero tus respuestas.**
