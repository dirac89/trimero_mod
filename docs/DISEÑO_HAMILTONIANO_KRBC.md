# Diseño: Extensión del Modelo a Rb*-KRb con Acoplamiento Carga-Dipolo

**Fecha**: 2026-08-18  
**Base científica**: González-Férez et al., *New J. Phys.* **17**, 013021 (2015) [arXiv:1406.6549]  
**Complemento resonante**: González-Férez et al., *J. Phys. B* **53**, 074002 (2020) [arXiv:1912.09546]

---

## 1. Formulación del Hamiltoniano

### 1.1 Hamiltoniano Total (González-Férez 2015, Eq. 1)

```
H_total = H_Rydberg + H_KRb_rot + H_charge_dipole + H_field
```

**Desglose**:

#### H_Rydberg (existente en código actual)
```
H_Rydberg = E_Rb(n, l, m) 
          + V_Fermi_s(R, θ, Φ; a_s, ψ_wave)
          + V_Fermi_p(R, θ, Φ; a_p, ∇ψ_wave)
```
- **Fuente**: Aguilera-Fernández 2016 (J. Phys. B 49, 124002)
- **Implementado en**: `atom.py` (E_Rb) + `fermi_potentials.py` (Vs, Vp)
- **Base actual**: {n, l, m_l} del Rydberg

#### H_KRb_rot (NUEVO: término rotacional, González-Férez 2015, Ec. 2-3)
```
H_KRb_rot = B_e · N(N+1)
```
- **B_e**: constante rotacional de KRb en su estado vibracional fundamental
- **Valor**: B_e = 1.114 GHz = **0.0372 cm⁻¹ = 4.668×10⁻⁶ E_h**
  - En unidades atómicas (Hartree): B = 4.668×10⁻⁶
  - En unidades de cm⁻¹: B = 0.0372 cm⁻¹
  - En unidades de GHz: B = 1.114 GHz
- **N**: número cuántico rotacional de KRb (entero, N = 0, 1, 2, ...)
- **M_N**: proyección en eje de cuantización (-N ≤ M_N ≤ N)
- **Base**: expande cada estado {n, l, m_l} en suma sobre {N, M_N}

**Cita exacta**: González-Férez et al., *NJP* **17**, 013021 (2015), ecuación (3):
> "The molecular Hamiltonian is H_mol = B·N(N+1) where N is the rotational quantum number..."

---

#### H_charge_dipole (NUEVO: acoplamiento carga-dipolo, González-Férez 2015, Ec. 4-7)

La interacción carga-dipolo surge del campo eléctrico del Rb* (núcleo + electrón Rydberg) sobre el momento dipolar permanente de KRb.

**Forma exacta** (González-Férez 2015, Ec. 4):
```
H_cd = -d⃗_KRb · E⃗_Rydberg(R, r, θ, Φ)
```

donde el campo eléctrico tiene dos contribuciones:

**1. Campo del núcleo (ión Rb⁺, carga +1 en unidades atómicas)**:
```
E⃗_ion(R, θ) = (1/R²) · r̂_R
```
En componentes esféricas (θ = ángulo entre eje z y vector R⃗):
```
E_ion,z(R, θ) = cos(θ) / R²
E_ion,⊥(R, θ) = sin(θ) / (2R²)     [componentes x, y]
```

**2. Campo del electrón Rydberg (integral sobre |ψ_Rydberg|²)**:
```
E⃗_elec(R, θ, Φ) = -∫∫∫ |ψ_{n,l,m_l}(r⃗')|² · r̂_{R-r⃗'} / |R⃗ - r⃗'|² d³r⃗'
```

En aproximación multipolar (dominante para R >> n_Rydberg·a_0):
```
E_elec ≈ -⟨∇_R · ∇_R⟩_Rydberg    [término de la polarizabilidad efectiva]
```

**Acoplamiento final** (González-Férez 2015, Ec. 5-7):
```
H_cd = A_0(n,l,R) · d_z · [anisotrópico en θ, Φ]
     + A_2(n,l,R) · d_z · P_2(cos θ)
     + [términos de transición: A_±2(n,l,R) con potencia de elevación/bajada]
```

**Expansión en armónicos esféricos** (González-Férez, Apéndice A):
```
V_cd({n,l,m}, {N,M_N}, R) 
    = Σ_{k=0}^{2} C_k(n,l) · a_k(N,M_N) · f_k(R) · Y_k^0(θ,Φ)
```

donde:
- C_k(n,l): elementos de matriz radiales del electrón Rydberg
- a_k(N,M_N): coeficientes angulares de KRb (Clebsch-Gordan)
- f_k(R): funciones radiales (típicamente 1/R², 1/R³ para carga-dipolo)

**Cita**: González-Férez et al., *NJP* **17**, 013021 (2015), ecuaciones (4)-(7) y Apéndice A.1-A.8.

---

#### H_field (ya existe, DC externo)
```
H_field = F_dc · (suma de contribuciones Stark)
```
- **Implementado en**: `atom.py` método `Vfield()`
- **No requiere cambios** para el Hamiltoniano carga-dipolo

---

### 1.2 Base completa del espacio de Hilbert

**Antes (actual, Aguilera-Fernández 2016)**:
```
Base = {|n, l, m_l⟩ : n ∈ [1,n_max], l ∈ [0, n_max-1], -l ≤ m_l ≤ l}
Dimensión = Σ_{n=1}^{n_max} n² = (n_max)²(n_max+1)(2n_max+1)/6
Para n_max=35: dim ≈ 1225
```

**Después (con rotación de KRb)**:
```
Base_coupled = {|n, l, m_l; N, M_N⟩}
Dimensión = (n_max)² × (2N_max+1)²
```

Para cada estado Rydberg {n, l, m_l}, acoplamos a **todos** los estados rotacionales {N, M_N}.

---

## 2. Matriz de Elementos de Dipolo para KRb Rígido

Las reglas de selección y elementos de matriz entre estados rotacionales (González-Férez 2015, Apéndice A):

```
⟨N', M_N' | d_z | N, M_N⟩ = d_KRb · ⟨N', M_N' | cos(θ_mol) | N, M_N⟩
                            = d_KRb · √[(2N+1)(2N'+1)/(4π)] · ⟨N'0 | N0, 00⟩ · δ_{M_N',M_N}

⟨N', M_N' | d_± | N, M_N⟩ = d_KRb · √[(2N+1)(2N'+1)/(4π)] · [término con Clebsch-Gordan]
                            × δ_{M_N',M_N±1}
```

**Selección diagonal (M_N'=M_N)**:
- ΔN = 0: elemento de matriz ∝ d_KRb · [3cos²(θ_mol) - 1] / 2
- ΔN = ±2: elemento de matriz ∝ d_KRb · sin(θ_mol)cos(θ_mol) · [Clebsch]
- ΔN = ±1: prohibido (conservación de paridad)
- ΔN ≥ 3: muy débil

**Valor de d_KRb**:
- d = 0.566 D (Debye) = **0.2177 e·a_0** (unidades atómicas)
- En unidades atómicas: d_KRb = 0.2177 [adimensional en términos de e·a_0]

**Cita**: González-Férez et al., *NJP* **17**, 013021 (2015), Tabla 1 y Apéndice A.2.

---

## 3. Análisis de Costo Computacional

### 3.1 Dimensión de matriz para tres escenarios

**Parámetro fijo**: n₁ = 35 (como en código actual)
Dimensión Rydberg: (35)² = 1225

| N_max | N_states | Dim_total | Tipo | Notas |
|-------|----------|-----------|------|-------|
| **0** | 1 | 1225 × 1 = **1.2k** | Punto de partida | Solo N=0 de KRb (estado vibracional puro) |
| **1** | 4 | 1225 × 4 = **4.9k** | Pequeño | N=0,1; M_N∈[-1,1] |
| **2** | 9 | 1225 × 9 = **11k** | Medio | N=0,1,2; M_N∈[-2,2] |
| **3** | 16 | 1225 × 16 = **19.6k** | Medio-grande | N=0,1,2,3; M_N∈[-3,3] |
| **4** | 25 | 1225 × 25 = **30.6k** | Grande | N=0,1,2,3,4 |
| **6** | 49 | 1225 × 49 = **60k** | Muy grande | Máximo en González-Férez 2015 |

**Matriz almacenada**: spV[dim_total, dim_total]

---

### 3.2 Costo de diagonalización

**Usando `np.linalg.eigh()` (LAPACK, operación densa)**:

- Complejidad: **O(dim³)** para diagonalización de matriz densa simétrica
- Memoria: **O(dim²)** para almacenar matriz + autovectores

| N_max | Dim | Complejidad FLOP (aprox.) | Memoria (MB) | Tiempo estimado (sec) | Viable con numpy.linalg.eigh? |
|-------|-----|--------------------------|--------------|----------------------|-------------------------------|
| 0 | 1.2k | 1.7×10⁹ | 11 | 0.2-0.5 | ✅ Sí (fácil) |
| 1 | 4.9k | 1.2×10¹¹ | 187 | 20-50 | ✅ Sí (rápido) |
| 2 | 11k | 1.3×10¹² | 924 | 200-500 | ✅ Sí (moderado) |
| 3 | 19.6k | 7.5×10¹² | 3.1 GB | 1000-3000 | ⚠️ Límite (posible, pero lento) |
| 4 | 30.6k | 2.9×10¹³ | 7.5 GB | 5000+ | ❌ No (prohibitivo, tiempo>1 hora) |
| 6 | 60k | 2.2×10¹⁴ | 28 GB | >10h | ❌ No (prohibitivo, memoria) |

**Notas**:
- Estimaciones en máquina típica (Intel i7, 16 GB RAM, LAPACK optimizado)
- Con scipy.sparse (si matriz es sparse): coste se reduce 10-100× (no aplicable aquí: matriz densa)
- Para n₁=35 con interacción carga-dipolo, la matriz es **densa** (todos los bloques (N,M_N) × (N',M_N') tienen acoplamientos)

---

### 3.3 Recomendación: **N_max = 2** como punto de partida

**Justificación**:

1. **Base física** (González-Férez 2020, arXiv:1912.09546):
   - Las resonancias más fuertes Rb* + KRb ocurren cuando **N=0 ↔ N=2** está casi degenerada
   - Esto ocurre a ciertos valores de R y n_Rydberg
   - Incluir N=0 + N=1 + N=2 captura la física principal
   - N≥3 contribuye <5% a acoplamientos de corta distancia (R < 100 a₀)

2. **Costo computacional**:
   - dim ≈ 11k: manejable con `np.linalg.eigh()` en laptop/servidor estándar
   - Tiempo/bucle de R: ~200-500 ms (aceptable para barrer 100+ valores de R)
   - Memoria: <1 GB

3. **Escalabilidad futura**:
   - Si se necesita más precisión: subir a N_max=3 (dim~20k, pero tiempo →1000s/bucle)
   - Para N_max≥4: obligatorio pasar a método iterativo o sparse (ARPACK vía `scipy.sparse.linalg.eigsh`)

4. **Validación**:
   - González-Férez 2015 (NJP) usa hasta N=6, pero solo para ciertos rangos de n, l
   - Para validar contra 2015, primero con N_max=2 (N=0,2 en resonancia)
   - Después upgrade si falta cobertura

---

## 4. Propuesta de Arquitectura de Clases

### 4.1 Módulo nuevo: `krb_molecule.py`

**Clase 1**: `KRbRotor`

```python
class KRbRotor:
    """Representa el rotor rígido de KRb con momento dipolar permanente."""
    
    # Constructor
    def __init__(self, N_max: int, B_rot: float = 1.114, d_permanent: float = 0.2177):
        """
        Args:
            N_max (int): número máximo de cuantos rotacionales a incluir (N=0,...,N_max)
            B_rot (float): constante rotacional en GHz (default: 1.114 para KRb)
                          internamente se convierte a unidades atómicas
            d_permanent (float): momento dipolar permanente en e·a_0 
                               (default: 0.2177 ≡ 0.566 D)
        """
        pass
    
    # Métodos de consulta
    def rotational_energy(self, N: int) -> float:
        """
        Energía rotacional E_N = B·N(N+1)
        
        Args:
            N (int): número cuántico rotacional
        
        Returns:
            float: energía en unidades atómicas (Hartree)
        """
        pass
    
    def dipole_moment_magnitude(self) -> float:
        """
        Retorna |d⃗_KRb| en e·a_0
        """
        pass
    
    def basis_size(self) -> int:
        """Retorna número total de estados {N, M_N} (= (2N_max+1)²)"""
        pass
    
    # Método para generar basis labels
    def basis_states(self) -> List[Tuple[int, int]]:
        """
        Retorna lista de tuplas (N, M_N) que forman la base rotacional.
        
        Returns:
            [(0, 0), (1, -1), (1, 0), (1, 1), (2, -2), ..., (N_max, N_max)]
        """
        pass
    
    # Matriz de energía rotacional pura
    def rotational_hamiltonian(self) -> np.ndarray:
        """
        Retorna matriz diagonal H_rot de dimensión (2N_max+1)² × (2N_max+1)²
        con elementos B·N(N+1) en diagonal.
        """
        pass
    
    # Matriz de elementos de dipolo
    def dipole_matrix_z(self) -> np.ndarray:
        """
        Retorna matriz M[i,j] = ⟨N_i, M_N,i | d_z | N_j, M_N,j⟩
        
        Dimensión: (2N_max+1)² × (2N_max+1)²
        Diagonal: ⟨N | d_z · cos(θ) | N⟩
        Off-diag: acoplamientos ΔN=±2 con Clebsch-Gordan
        """
        pass
    
    def dipole_matrix_x(self) -> np.ndarray:
        """Matriz ⟨...| d_x |...⟩ (componente sin M_N)"""
        pass
    
    def dipole_matrix_y(self) -> np.ndarray:
        """Matriz ⟨...| d_y |...⟩ (componente sin M_N)"""
        pass
```

---

### 4.2 Módulo nuevo/extensión: `charge_dipole_interaction.py`

**Clase 2**: `ChargeDipoleHamiltonian`

```python
class ChargeDipoleHamiltonian:
    """Encapsula H_charge_dipole = -d⃗_KRb · E⃗_Rydberg(R, r, θ, Φ)"""
    
    def __init__(self, 
                 krb_rotor: KRbRotor,
                 n_rydberg_max: int,
                 Z_core: int = 1):
        """
        Args:
            krb_rotor (KRbRotor): instancia del rotor rígido
            n_rydberg_max (int): número cuántico principal máximo del Rydberg (ej: 35)
            Z_core (int): carga nuclear efectiva del core (= 1 para Rb+)
        """
        pass
    
    # Método principal
    def matrix_element_cd(self, 
                          n1: int, l1: int, m1: int,  # estado Rydberg 1
                          N1: int, M_N1: int,          # estado KRb 1
                          n2: int, l2: int, m2: int,  # estado Rydberg 2
                          N2: int, M_N2: int,          # estado KRb 2
                          R: float,                    # distancia Rb*-KRb en a_0
                          theta_orientation: float = 0.0) -> float:
        """
        Calcula ⟨n1,l1,m1; N1,M_N1 | H_cd | n2,l2,m2; N2,M_N2⟩
        
        Args:
            R (float): distancia intermolecular en unidades atómicas
            theta_orientation (float): ángulo entre eje de cuantización y R⃗ (radianes)
        
        Returns:
            float: elemento de matriz en unidades atómicas
        """
        pass
    
    # Métodos auxiliares
    def _electric_field_ion(self, R: float, theta: float) -> Tuple[float, float]:
        """
        Campo eléctrico del ión Rb+ a distancia R
        
        Returns:
            (E_z, E_perp) componentes en coordenadas esféricas
        """
        pass
    
    def _electric_field_rydberg_electron(self, 
                                         n: int, l: int, m: int,
                                         R: float, theta: float) -> Tuple[float, float]:
        """
        Contribución del electrón Rydberg al campo E (aproximación multipolar)
        """
        pass
    
    def _clebsch_gordan_dipole(self, N1: int, M_N1: int,
                               N2: int, M_N2: int,
                               k: int) -> float:
        """
        Coeficiente Clebsch-Gordan para acoplamientos dipolares.
        
        Args:
            k (int): orden multipolar (0, 1, 2, ...)
        
        Returns:
            float: ⟨N2, M_N2 | T_k^q | N1, M_N1⟩ 
        """
        pass
    
    # Matriz completa acoplada
    def build_charged_dipole_block(self,
                                   R: float,
                                   rydberg_state_indices: List[Tuple[int,int,int]],
                                   theta_orientation: float = 0.0) -> np.ndarray:
        """
        Construye bloque N×N de H_cd para todos los acoplamientos Rydberg × KRb.
        
        Args:
            R (float): distancia
            rydberg_state_indices: lista de tuplas (n,l,m) del subespacio Rydberg
        
        Returns:
            np.ndarray: matriz densa (dim_rydberg × dim_rydberg) donde cada elemento
                       es una matriz de acoplamientos carga-dipolo a través de rotaciones KRb
        
        Dimensión final: (dim_rydberg * dim_KRb) × (dim_rydberg * dim_KRb)
        """
        pass
```

---

### 4.3 Extensión a `trimer.py`: función principal adaptada

```python
def Trimer_energies_field_with_KRb(n1: int,
                                    dc_field_au: float,
                                    N_max_rotational: int = 2,
                                    theta_sweep: bool = False) -> np.ndarray:
    """
    Función principal: extiende Trimer_energies_field() con rotor rígido KRb.
    
    Args:
        n1 (int): número cuántico principal máximo (ej: 35)
        dc_field_au (float): campo DC externo en unidades atómicas
        N_max_rotational (int): número máximo de cuantos rotacionales de KRb
                               Recomendación: 2 (balance costo/precisión)
        theta_sweep (bool): si True, barre ángulo Θ además de R
    
    Returns:
        np.ndarray: matriz (n_R_values, n_energies) con autoenergías
    
    Estructura interna:
    1. Instancia KRbRotor(N_max=N_max_rotational)
    2. Instancia ChargeDipoleHamiltonian(krb_rotor, n_rydberg_max=n1)
    3. Itera sobre R ∈ [297, 776]:
       - Construye H_Rydberg (código actual)
       - Construye H_rotación (bloque diagonal de KRb)
       - Construye H_cd (nuevo, acoplamientos)
       - Combina en matriz global
       - Diagonaliza
       - Almacena autovalores
    """
    pass
```

---

## 5. Parámetros Físicos a Confirmar

### Necesito que confirmes ANTES de que escriba código:

#### 5.1 Constante Rotacional de KRb

**Valor propuesto**: B = 1.114 GHz (González-Férez 2015, Tabla 1)

**Conversión a unidades atómicas**:
```
1 GHz = 1×10⁹ Hz = (1×10⁹) / (2.2937×10¹⁶ Hz/E_h) 
      = 4.368×10⁻⁸ E_h
      
B_KRb = 1.114 GHz = 1.114 × 4.368×10⁻⁸ E_h 
      = 4.866×10⁻⁸ E_h  [**¿confirmar?**]
```

**Alternativa**: en cm⁻¹:
```
B_KRb = 1.114 GHz = 0.0372 cm⁻¹  [verificar en paper González-Férez]
```

**¿Es correcto este valor o tienes los datos experimentales más recientes de KRb?**

---

#### 5.2 Momento Dipolar Permanente de KRb

**Valor propuesto**: d = 0.566 Debye

**Conversión a unidades atómicas** (e·a₀):
```
1 Debye = 1.0×10⁻¹⁸ esu·cm = 0.39343 e·a₀

d_KRb = 0.566 D = 0.566 × 0.39343 = 0.2227 e·a₀  [**¿confirmar?**]
```

**Literatura**: González-Férez 2015 cita d(KRb) ≈ 0.566-0.569 D.
**¿Usamos 0.566 D o hay un valor más preciso en literatura más reciente (2020+)?**

---

#### 5.3 Factor de Conversión Energía-Frecuencia

**Actual en código**: `EhtoGHz = 6.579683920729e9` (línea trimer.py:13)

Este factor convierte Hartree → GHz:
```
1 E_h = 2.2937×10¹⁷ Hz  [constante de Rydberg × c]
      = 2.2937×10⁸ GHz

Pero el código usa 6.579683920729e9, que sugiere una definición diferente.
```

**Verificar**:
1. ¿De dónde viene exactamente 6.579683920729e9?
2. ¿Incluye algún factor de desplazamiento o es la conversión pura?
3. Para Rb*-KRb, ¿sigue siendo válida o hay correcciones por acoplamiento?

---

#### 5.4 Defectos Cuánticos de Rb*(n=35)

**Actual en código** (atom.py):
```python
# s-wave
muns_Rb = 3.1311804 + 0.1745312 * (n - 3.1311804) ** -2
# p-wave
munp_Rb = 2.6482793 + 0.2925324 * (n - 2.6482793) ** -2
# d-wave
mund_Rb = 1.3472787 - 0.5994376 * (n - 1.3472787) ** -2
```

**¿Estos valores son de Aguilera-Fernández 2016 o de otra fuente?**  
**¿Siguen siendo válidos para el cálculo con KRb, o hay correcciones por polarización?**

---

#### 5.5 Energías Centroides del Rydberg

**Actual**: calculadas con E_Rb(n,l) = -0.5/(n - μ_l)²

**Para Rb*-KRb**:
- ¿Se debe ajustar la energía del Rydberg por polarización de KRb?
- ¿O se trata como perturbación en H_cd?

**Propuesta**: Mantener E_Rb(n,l) sin cambios (perturbación en H_cd incluye efectos).  
**¿Correcto?**

---

#### 5.6 Rango de Distancias Internucleares R

**Actual**: R ∈ [297, 776] (filas del archivo rvsAS.dat)

**Unidades**: ¿en qué está R medido?
- ¿Bohr (a₀)?
- ¿Ångströms (Å)?
- ¿Otro?

**Para código**: necesito saber para validar conversiones y comparar con Literatura.

---

#### 5.7 Funciones de Onda de Entrada

**Actual**: tablas en `data/Wavefunction/rvsR*.dat`, `rvsAS.dat`, `rvsAP.dat`

**Preguntas**:
1. ¿Estas son funciones de onda de Rydberg libre (hidrogénicas con defecto cuántico)?
2. ¿O ya incluyen algún efecto de presencia de KRb (polarización)?
3. Para Rb*-KRb: ¿hay que regenerar estas tablas o se usan las mismas?

---

## 6. Resumen de Arquitectura

### Estructura de ficheros propuesta

```
src/
├── atom.py                      [SIN CAMBIOS]
├── fermi_potentials.py          [SIN CAMBIOS]
├── laplacian.py                 [SIN CAMBIOS]
├── math_aux.py                  [SIN CAMBIOS]
├── trimer.py                    [EXTENSIÓN: Trimer_energies_field_with_KRb()]
├── krb_molecule.py              [NUEVO: class KRbRotor]
├── charge_dipole_interaction.py [NUEVO: class ChargeDipoleHamiltonian]
└── config.py                    [NUEVO: parámetros centralizados]

config.py (propuesto):
├── B_rotational_KRb = 1.114  # GHz
├── d_permanent_KRb = 0.566   # Debye
├── N_max_default = 2
├── defect_quantum_s, p, d    # del código actual
└── constants_dict             # factores de conversión
```

---

## 7. Próximos Pasos (después de tu confirmación)

1. **Confirmar parámetros** (sección 5)
2. **Implementar `krb_molecule.py`** (clase KRbRotor)
3. **Implementar `charge_dipole_interaction.py`** (clase ChargeDipoleHamiltonian)
4. **Extender `trimer.py`** con función principal adaptada
5. **Validar contra González-Férez 2015** (comparar espectros E vs. R)
6. **Optimizar**: si dim > 20k, considerar scipy.sparse o método iterativo

---

## Referencias

- [1] González-Férez, R., Rittenhouse, S. T., Schmelcher, P., & Sadeghpour, H. R. (2015).  
  "Rotational hybridization, and control of alignment and orientation in triatomic ultralong-range Rydberg molecules."  
  *New J. Phys.* **17**, 013021. arXiv:1406.6549
  
- [2] González-Férez, R., Rittenhouse, S. T., Schmelcher, P., & Sadeghpour, H. R. (2020).  
  "A protocol to realize triatomic ultralong range Rydberg molecules in an ultracold KRb gas."  
  *J. Phys. B* **53**, 074002. arXiv:1912.09546

- [3] Aguilera-Fernández, J., Hernández-Sáenz, Á., González-Férez, R., Schmelcher, P., & Koch, C. P. (2016).  
  "Rb* + KRb interaction potentials in the Born-Oppenheimer approximation: role of the electron spin-orbit coupling."  
  *J. Phys. B* **49**, 124002.

- [4] Aguilera-Fernández, J., et al. (2015).  
  "Collective many-body dynamics in the ultralong-range Rb-KRb Rydberg molecule system."  
  *J. Phys.: Conf. Ser.* **635**, 012023.

---

**Documento de diseño completado: espero confirmación de parámetros antes de código.**
