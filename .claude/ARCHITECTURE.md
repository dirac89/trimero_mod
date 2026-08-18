# Arquitectura del Proyecto Trimero

## Visión General

Este proyecto implementa una **simulación cuántica de trimero atómico** que:
1. Carga datos de funciones de onda atómicas
2. Calcula potenciales de Fermi entre átomos
3. Construye una matriz Hamiltoniana para el sistema
4. Diagonaliza la matriz para obtener niveles de energía
5. Exporta autovalores en función del radio y campo eléctrico

## Capas de Arquitectura

### 1. **Capa de Entrada de Datos** (`data/Wavefunction/`)
- Archivos `.dat` con datos de funciones de onda: `rvsAS.dat`, `rvsAP.dat`, `rvsR38s.dat`, etc.
- Archivos `.txt` con valores esperados: `exp_val_r.txt`
- Todos los archivos se cargan en memoria al inicio en `Trimer_energies_field()`

### 2. **Capa Física** (módulos en `src/`)

#### `atom.py` - Clase `Atom`
Encapsula la física de un átomo individual:
- Propiedades: momento angular, paridad, índice de estado
- Métodos: cálculo de overlaps, funciones de estado

```python
class Atom:
    def __init__(self, name, l_value, parity, ...):
        # Define el estado atómico
```

**Responsabilidad**: Definir estados atómicos base sin interacción.

#### `fermi_potentials.py` - Clase `FermiPotentials`
Calcula potenciales de Fermi como función de la distancia:
- Métodos: `get_potential(r)`, `set_parameters()`
- Parametrización: exponenciales decrecientes típicas de Fermi

```python
class FermiPotentials:
    def __init__(self, a0_length):
        self.a0 = a0_length
    
    def get_potential(self, r_au):
        # Retorna V(r) en unidades atómicas
```

**Responsabilidad**: Proveer interacciones de dos cuerpos entre átomos.

#### `math_aux.py` - Funciones Matemáticas Especiales
Contiene utilidades numéricas:
- **Armónicos esféricos**: Y_lm(θ, φ)
- **Derivadas radiales**: dψ/dr
- **Integrales**: overlaps entre funciones de onda

```python
def spherical_harmonics(l, m, theta, phi):
    # Y_lm(θ, φ)

def radial_derivative(psi, r_grid):
    # Calcula dψ/dr numéricamente
```

**Responsabilidad**: Proporcionar primitivas matemáticas reutilizables.

#### `laplacian.py` - Interfaz Pública de Matemáticas
Expone funciones de `math_aux.py` de forma limpia:
```python
from laplacian import spherical_harmonics, radial_derivative
```

**Responsabilidad**: Mantener una API coherente y versátil.

### 3. **Capa de Simulación** (`trimer.py`)

#### Función Principal: `Trimer_energies_field(n1, dc_field_au, ...)`

**Responsabilidad Central**: Orquestar toda la simulación.

```python
def Trimer_energies_field(n1, dc_field_au, ...):
    # 1. Cargar datos
    data = load_all_data()
    
    # 2. Construir matriz de campo eléctrico
    field_matrix = build_dc_field_matrix(n1, dc_field_au)
    
    # 3. Bucle principal: para cada radio R
    for R in radius_points:
        # Construir H(R) con lógica de casos A, B, C, D
        H = build_hamiltonian(R, data, field_matrix)
        
        # Diagonalizar
        eigenvalues = scipy.linalg.eigvalsh(H)
        
        # Guardar
        save_eigenvalues(R, eigenvalues)
```

**Lógica de Casos Físicos**:
- **Caso A/B/C/D**: Diferentes configuraciones de solapamiento angular
  - Determinan qué elementos de matriz Hamiltoniana están activos
  - Afectan la dimensión de la matriz

**Matriz Hamiltoniana**:
```
H = H0 + H_interaction + H_field

H0: Energía de átomos individuales (diagonal)
H_interaction: Potenciales de Fermi (off-diagonal)
H_field: Interacción con campo eléctrico externo
```

### 4. **Capa de Entrada-Salida** (en `trimer.py` y `main.py`)

#### Lectura
- `numpy.loadtxt()` para `.dat`, `.txt`
- Manejo de rutas relativas a `data/Wavefunction/`

#### Escritura
- Formato: ASCII, una línea por punto de R
- Estructura: `R eigenvalue_1 eigenvalue_2 ... eigenvalue_N`
- Archivos: `Trimer_R_sp_wave_N{n1}_R_{int(100*radius)}_GHz.dat`, etc.

### 5. **Capa de Control** (`main.py`)

Punto de entrada que:
- Define parámetros físicos (n1, dc_field_au, radio_min, radio_max)
- Llama `Trimer_energies_field()`
- Incluye `test_trimer_energies_field()` para validación rápida

## Dependencias Entre Módulos

```
main.py
  └─> trimer.py (Trimer_energies_field)
       ├─> atom.py (Atom)
       ├─> fermi_potentials.py (FermiPotentials)
       ├─> math_aux.py (spherical_harmonics, derivadas, etc.)
       ├─> laplacian.py (interfaz de math_aux)
       └─> data/ (lectura de archivos .dat)
```

## Invariantes de Diseño

1. **Separación de Responsabilidades**:
   - `atom.py`: Define estados, no resuelve ecuaciones
   - `fermi_potentials.py`: Calcula V(r), no construye Hamiltonianos
   - `math_aux.py`: Primitivas matemáticas, no contexto físico
   - `trimer.py`: Orquesta y ejecuta la simulación

2. **Inmutabilidad de Datos**:
   - Los datos cargados (`rvsAS`, etc.) no se modifican después de la carga
   - Los resultados se escriben una única vez

3. **Unidades Consistentes**:
   - Entrada: Unidades atómicas (Bohr, Hartree)
   - Salida: Hartree (energía) o GHz (si se convierte)
   - Conversión explícita en funciones específicas

4. **Determinismo**:
   - Mismo input → Mismo output (sin aleatoriedad)
   - Útil para validación y reproducibilidad

## Flujo de Ejecución

```
1. main.py → Leer parámetros
2. Trimer_energies_field() → Cargar data/
3. loop R:
   a. Construir H(R) con casos A/B/C/D
   b. Diagonalizar H(R)
   c. Guardar autovalores
4. Salida: Archivos *.dat con resultados
```

## Puntos de Extensión

### Agregar Nueva Física
- **Nuevos términos en H**: Edita `build_hamiltonian()` en `trimer.py`
- **Nuevos potenciales**: Extiende `FermiPotentials`
- **Nuevos estados atómicos**: Añade instancias de `Atom`

### Mejora de Rendimiento
- **Paralelizar bucle de R**: Usa `multiprocessing.Pool` o `joblib`
- **Cachear diagonalizaciones**: Almacena H para R idénticos
- **Vectorizar operaciones**: Reemplaza loops con numpy donde sea posible

### Validación
- **Comparar con C++**: Los archivos `.dat` deben coincidir en precisión
- **Tests unitarios**: Crea en `tests/` para funciones de `math_aux.py`
- **Visualización**: Usa matplotlib para graficar espectros

## Consideraciones de Rendimiento

| Operación | Complejidad | Nota |
|-----------|------------|------|
| Cargar datos | O(n_datos) | Una única vez |
| Construir H | O(n_basis²) | Depende de n1 |
| Diagonalizar H | O(n_basis³) | Bottleneck principal |
| Bucle de R | O(n_radius × n_basis³) | Paralelizable |

Para n1 pequeño (5–10), es rápido. Para n1 > 30, considera paralelización.

---

**Última actualización**: 2026-08-18  
**Versión de migración**: Python 3.13+, NumPy 2.3.1+, SciPy 1.16.0+
