# Auditoría de Migración Python y Análisis del Modelo Rb*-KRb

**Fecha**: 2026-08-18  
**Rama auditoría**: `migrate-python` (actual) vs. `master` (C++ original)  
**Conclusión general**: Migración C++→Python es **correcta**, pero modelo actual es **Aguilera-Fernández 2016** (Rb* + neutros), **no** Rb*-KRb con carga-dipolo (2015/2017).

---

## 1. Veredicto sobre Hamiltoniano Carga-Dipolo

### Búsqueda en todas las ramas

```bash
# Ramas disponibles:
# - master (C++ original)
# - migrate-python (Python actual)
# - remotes/origin/master
# - remotes/origin/migrate-python
```

**Resultado de búsqueda**: ❌ **NO EXISTE** implementación del Hamiltoniano carga-dipolo

- **Búsquedas realizadas**:
  - `grep -r "dipole\|carga-dipolo\|KRb\|rotor\|M_N\|rotational constant"` → sin resultados
  - `grep -r "permanent.*dipole\|dipolo.*permanente"` → sin resultados
  - Revisión manual de todos los archivos `.cpp` y `.h` en `master`
  - Revisión de todos los módulos Python en `migrate-python`

- **Archivos relevantes revisados**:
  - `trimero.cpp` (líneas 49-599): θ hardcodeado (56: `theta=0.0`, 58: `theta1=pi`)
  - `FermiPotentials.cpp/h`: solo pseudopotenciales Fermi s+p puro
  - `Atom.cpp/h`: solo energías del Rydberg + campo DC externo (`Vfield`, `Angular_dc_field`)
  - `trimer.py`, `atom.py`, `fermi_potentials.py`: réplicas Python fiel de C++

### Modelo implementado: Aguilera-Fernández et al. 2016

**Configuración actual**:
- Sistema: Rb*(n=35) + 2 átomos **neutros** (sin KRb)
- Interacción: Pseudopotenciales Fermi s+p puro
- Base de Rydberg: defectos cuánticos para Rb (s, p, d)
- Campo externo: solo DC field (opcional, línea 213-220 en C++; línea 52-55 en Python)
- **SIN**: 
  - Constante rotacional B de rotor rígido
  - Momento dipolar permanente
  - Número cuántico rotacional N, M_N
  - Acoplamiento carga-dipolo ion-dipolo

**Referencias publicadas (modelo actual)**:
- Aguilera-Fernández et al., *J. Phys. B*, **49**, 124002 (2016)
- Rubidio + 2 átomos neutros (sin KRb)

**Modelos que FALTA implementar (2015/2017)**:
- González-Férez et al., *New J. Phys.*, **17**, 013021 (2015)
- Aguilera-Fernández et al., *J. Phys.: Conf. Ser.*, **635**, 012023 (2015)
- **Requieren**: Rb* + KRb (rotor rígido) con acoplamiento carga-dipolo

---

## 2. Comparación Estructural: C++ vs. Python

### Migración de clases y funciones

| Componente | C++ | Python | Estado |
|-----------|-----|--------|--------|
| `Atom` class | `Atom.cpp/h` | `atom.py` | ✅ Idéntico |
| `FermiPotentials` class | `FermiPotentials.cpp/h` | `fermi_potentials.py` | ✅ Idéntico |
| `Laplacian` utils | `Laplacian.cpp/h` | `laplacian.py` | ✅ Idéntico |
| `Trimer_energies_field()` | `trimero.cpp:49-599` | `trimer.py:5-242` | ✅ Lógica equivalente |
| Diagonalización | `gsl_eigen_symmv()` | `np.linalg.eigh()` | ✅ Compatible |

### Análisis de índices y bucles

**C++ (línea 186-191 en trimero.cpp)**:
```cpp
for(i = 0; i <= n1-1-1 ; i++)
{
    infile >> value;
    gsl_matrix_set (Radial, i, i+1, value);
}
```

**Python (línea 35-37 en trimer.py)**:
```python
for i in range(n1 - 1):
    value = float(infile.readline().strip())
    Radial[i, i + 1] = value
```

✅ **Correcto**: `range(n1-1)` ≡ `[0, n1-2]` ≡ C++ `i <= n1-2`

**Bucles anidados (construcción de matriz H)**:
- C++ (línea 243-550): 4 niveles anidados (i, k1, j, k2)
- Python (línea 75-225): 4 niveles anidados, **estructura idéntica**
- Deltas Kronecker, condiciones (i==j && m1==m2): **replicadas correctamente**

✅ **Conclusión**: La migración **respeta la lógica original**. No hay off-by-one errors ni reorganización de bucles.

---

## 3. Auditoría de Calidad del Código Python (migrate-python/src)

### 3.1 Errores y bugs

#### 🔴 BUG CRÍTICO: Bucle principal limitado

**Ubicación**: `trimer.py:64`

```python
for row in range(297, 307):  # Prueba corta: solo 10 iteraciones
```

**Problema**: Solo itera 10 filas (~2.6 GHz), cuando debería iterar **776 filas** como en C++.

**Impacto**: Resultados incompletos. Datos de salida truncados.

**Fix**:
```python
for row in range(297, rows):  # rows = 776
```

---

#### 🟡 ADVERTENCIA: Precisión numérica y escalas de energía

- C++ usa `double` (precisión IEEE 754, ~15-17 dígitos decimales)
- Python usa `numpy.float64` (equivalente)
- **Conversión de unidades**: `EhtoGHz = 6.579683920729e9` (línea 13)
  - Verificar que coincida exactamente con literatura (Rb*-KRb papers)

---

### 3.2 Defectos estructurales

#### **1. Ángulos hardcodeados sin parametrización**

**Líneas 10-12**:
```python
theta = 0.0
theta1 = pi
```

**Problema**: 
- θ₁ = π implica KRb **perpendicular** al eje de cuantización
- θ₂ = 0 implica Rb* **paralelo**
- No hay barrido angular ni exploración de configuraciones
- Para Rb*-KRb será crítico barrer θ(Θ) en toda la esfera

**Impacto**: Imposible extender a cálculos de potencial de superficie E(R, Θ)

---

#### **2. Estructura de casos anidados (Cases A, B, C, D)**

**Líneas 88-216**: Lógica de 4 casos con `if/elif` anidados

```python
if i < 3 and j < 3:          # CASO A
    # ...
elif i < 3 and j > 2:        # CASO B
    # ...
elif i > 2 and j < 3:        # CASO C
    # ...
elif i > lc and j > lc:      # CASO D
    # ...
else:
    # ...
```

**Problemas**:
- Muy frágil: cambiar un threshold (ej. `lc=2`) propaga regresiones
- Duplicación de código (~90 líneas) sin abstracción
- Difícil de depurar y validar
- No escalable si se añaden nuevas bases cuánticas

**Refactor sugerido**: Matriz de mapeos o generador de parámetros (n11, n21, n12, n22, wave_*, Dwave_*)

---

#### **3. Falta de parametrización de datos de entrada**

**Líneas 16-27**: Rutas hardcodeadas

```python
data_dir = "data/Wavefunction/"
As = np.loadtxt(data_dir + "rvsAS.dat")
# ...
```

**Problemas**:
- Sin validación de que archivos existan
- Sin checksum o versión de datos
- Imposible cambiar a funciones de onda KRb sin hard-editing

**Sugerencia**: Config file (YAML/JSON) con rutas y metadata

---

### 3.3 Ausencia de validaciones

- ❌ No se valida rango de R (¿R ∈ [1, ∞)? ¿unidades?¿ nm, a₀?¿)
- ❌ No se verifica convergencia de matriz H (autovalores → ∞?)
- ❌ No hay warnings si dim(H) es muy grande (matriz densa, ∝ n⁴)
- ❌ No se compara contra datos experimentales o publicados

---

## 4. Qué Falta Portar para Rb*-KRb Funcional

### 4.1 Nuevas estructuras de datos

```python
# Actual: base {n, l, m}
# Necesario: base {N, M_N, n, l, m}

class RydbergRotationalState:
    """Estado cuántico de Rydberg + rotor rígido"""
    N: int          # número cuántico rotacional de KRb
    M_N: int        # proyección en eje de cuantización
    n: int          # número cuántico principal de Rb*
    l: int          # momento angular orbital de Rb*
    m: int          # proyección de l
```

**Dimensión de espacio base**: 
- Actual: n₁² (si n₁=35 → 1225)
- Rb*-KRb: n₁² × (2N_max+1)² (si N_max=3 → ~1225×49 ≈ 60k)

---

### 4.2 Nuevos términos en Hamiltoniano

**Matriz H (forma bloque)**:

```
H = H_Rydberg + H_KRb_rot + H_charge_dipole + H_ion_dipole

H_Rydberg      → E_Rb(n,l) + pseudopotenciales Fermi s+p  [EXISTE]
H_KRb_rot      → B·N(N+1) + rotación rígida               [FALTA]
H_charge_dipole → -d·E (carga-dipolo Rb*-KRb)            [FALTA]
H_ion_dipole    → -d_KRb·E_ion (dipolo-ion)              [FALTA]
```

**Términos específicos a implementar** (González-Férez 2015, Ec. 1-15):

1. **Rotación rígida de KRb**: 
   ```
   B_e(R) · N(N+1)  donde B_e(R) = constante rotacional dependiente de R
   ```

2. **Acoplamiento Coulombiano Rb⁺-d(KRb)**:
   ```
   V_ion-dipole = (k_e·q / r³) · [3(d·r̂)r̂ - d]
   Términos ∝ cos(Θ), sin(Θ)
   ```

3. **Acoplamiento electrón Rydberg-dipolo KRb**:
   ```
   V_elec-dipole ∝ (función de onda Rydberg) × (momento dipolar KRb)
   ```

---

### 4.3 Funciones de onda nuevas

**Necesarios datos de entrada**:
- Función de onda rotacional de KRb: `R_N(R)` (en lugar de `R37p`, `R38s`)
- Matriz de elementos de dipolo: `<N,M_N|d_z|N',M_N'>`
- Constante rotacional: `B_e(R)` (tabla vs. R o fórmula)

**Archivos equivalentes a generar**:
```
data/Wavefunction/
├── rvsR_N0_N1.dat      # R(R) para KRb vibracional
├── rvsDR_N0_N1.dat     # dR/dR
├── dipole_matrix.dat   # elementos de matriz dipolar
└── B_rotational.dat    # B_e(R)
```

---

## 5. Riesgos para Reproducibilidad Numérica

### 5.1 Riesgos críticos (🔴)

| Riesgo | Descripción | Impacto | Mitigación |
|--------|-------------|--------|-----------|
| **Cambio de base** | Pasar de dim=1225 a dim≈60k | Matriz densa → requiere precondicionamiento | Usar sparse matrix (scipy.sparse) |
| **Escalas de energía** | Actual ~10⁰ mK; carga-dipolo ~10² mK | Pérdida de precisión si no se re-escala | Usar energía relativa (resta del mínimo) |
| **Nuevos parámetros** | B_e(Rb), d(KRb), α (polarizabilidad) | Si no coinciden literatura → divergencia | Citar Aguilera-Fernández 2015 tab. 1 |

### 5.2 Riesgos medios (🟡)

| Riesgo | Descripción | Impacto |
|--------|-------------|--------|
| **Ángulo Θ sin barrido** | C++ actual NO barre Θ | No hay validación vs. E(R,Θ) publicados |
| **Entrada de datos** | rvsR*.dat son para Rb*(n=35) + neutros | ¿Compatible con KRb rotacional? Verificar |
| **Diagonalización** | numpy.linalg.eigh usa LAPACK | Si matriz muy grande (~60k×60k) → timeout |

### 5.3 Verificación de convergencia

**Antes de usar en publicación**:
1. Calcular espectro E(R, Θ=0) a varios n₁ (n=30, 35, 40) → debe converger
2. Comparar potenciales enlazantes con Fig. 2 de González-Férez (2015)
3. Comparar números cuánticos {N, M_N, n} vs. asignaciones publicadas
4. Check: dispersión de estados ~1 GHz (vs. ~10 GHz para Rb* solo)

---

## 6. Lista Priorizada de Modificaciones

### Prioridad 1: Bugs críticos (semana 1)

- [ ] **Fix bucle principal** (trimer.py:64): `range(297, 307)` → `range(297, rows)`
- [ ] **Parametrizar ángulos**: pasar `theta`, `theta1` como argumentos (no constantes)
- [ ] **Verificar unidades**: R en a₀?, energía en E_h?, EhtoGHz correcto?

### Prioridad 2: Modelo carga-dipolo (semana 2-3)

- [ ] **Leer papers**: González-Férez 2015 (NJP 17) Ec. 1-20, Aguilera-Fernández 2015 (JPCS 635)
- [ ] **Clase RydbergRotationalState**: extender base cuántica a {N, M_N, n, l, m}
- [ ] **Hamiltoniano carga-dipolo**: codificar H_rot, H_ion-dipole, H_elec-dipole
- [ ] **Datos de entrada KRb**: preparar tablas B_e(R), funciones de onda rotacionales, matriz dipolar

### Prioridad 3: Refactor estructura (semana 2, paralelo)

- [ ] **Simplificar Cases A-D**: usar tabla o generador de parámetros
- [ ] **Config file (YAML)**: rutas datos, parámetros físicos, tolerancias
- [ ] **Logging**: reemplazar `print()` por `logging` estándar

### Prioridad 4: Validación (semana 3)

- [ ] **Benchmark vs. literatura**: E(R) comparar con tablas 2015/2017
- [ ] **Análisis de convergencia**: E(R) vs. n₁
- [ ] **Matriz de dispersión**: verificar estados del continuo

---

## 7. Resumen Ejecutivo

### Preguntas respondidas

**¿Existe Hamiltoniano carga-dipolo en C++?**  
❌ **NO**. Ni en master ni en ninguna rama.

**¿Qué falta para tener Rb*-KRb funcional?**  
Escribir **desde cero**:
1. Base cuántica {N, M_N, n, l, m}
2. Hamiltonianos de rotor rígido + carga-dipolo
3. Datos de entrada KRb (B_e, d, funciones onda)

**¿Es la migración Python correcta?**  
✅ **SÍ**. Replica fielmente C++. Un bug (bucle truncado).

**¿Cuánto tiempo?**  
2–3 semanas (desarrollo + validación contra literatura).

**¿Riesgos numéricos?**  
- Escalas de energía nuevas (~10² mK vs. 10⁰)
- Matriz 50× más grande (sparse matrix recomendado)
- Parámetros B_e, d deben coincidir literatura

---

## Apéndice A: Archivos auditoría

### Archivos C++ (master)
```
src/trimero.cpp          → líneas 49-599 (función principal)
src/Atom.cpp/h           → energías Rb* + campo DC
src/FermiPotentials.cpp/h → pseudopotenciales s+p
src/Laplacian.cpp/h      → funciones angulares (Y_lm, etc.)
```

### Archivos Python (migrate-python)
```
src/trimer.py            → función principal E(R)
src/atom.py              → clase Atom (replica C++)
src/fermi_potentials.py  → clase FermiPotentials
src/laplacian.py         → funciones auxiliares
src/math_aux.py          → hidrogenoide, especiales GSL
```

---

## Apéndice B: Referencias clave

- [1] González-Férez, R., et al. *New J. Phys.* **17**, 013021 (2015)
  - **Modelo**: Rb* + KRb (rotor rígido), Hamiltonian carga-dipolo
  - **Figuras críticas**: Fig. 2 (curvas E(R) con Θ)

- [2] Aguilera-Fernández, J., et al. *J. Phys.: Conf. Ser.* **635**, 012023 (2015)
  - **Parámetros**: B_e(Rb), d(KRb), defectos cuánticos
  - **Tablas**: valores numéricos para validación

- [3] Aguilera-Fernández, J., et al. *J. Phys. B* **49**, 124002 (2016)
  - **Modelo actual**: Rb* + 2 átomos neutros (sin carga-dipolo)
  - **Pseudopotenciales**: Fermi s+p, defecto cuántico Rb

---

**Fin del análisis.**

Última actualización: 2026-08-18  
Auditor: Claude Code AI
