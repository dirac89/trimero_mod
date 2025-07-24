# Proyecto Trimero (Migración a Python)

Este proyecto simula la física de un trimero atómico usando potenciales de Fermi y diagonalización de matrices Hamiltonianas. Originalmente desarrollado en C++, ha sido completamente migrado a Python para facilitar su uso, mantenimiento y extensión.

## Estructura del Proyecto

```
trimero_mod/
├── pyproject.toml         # Configuración de Poetry y dependencias
├── README.md              # Este archivo
├── data/
│   └── Wavefunction/      # Archivos de datos necesarios para la simulación
└── src/
    ├── atom.py
    ├── fermi_potentials.py
    ├── laplacian.py
    ├── math_aux.py
    ├── main.py
    └── trimer.py
```

## Dependencias

- Python >= 3.11
- numpy
- scipy
- poetry (para gestión de entorno y dependencias)

Instala las dependencias y crea el entorno virtual con:
```sh
poetry install
```

## Descripción de los módulos principales

- **src/main.py**: Punto de entrada del proyecto. Llama a la función principal de simulación y permite ejecutar pruebas rápidas.
- **src/trimer.py**: Implementa la función principal `Trimer_energies_field`, que realiza la simulación física, lee los datos, construye y diagonaliza la matriz Hamiltoniana, y guarda los resultados.
- **src/atom.py**: Implementa la clase `Atom`, que encapsula la física atómica relevante.
- **src/fermi_potentials.py**: Implementa la clase `FermiPotentials`, que calcula los potenciales de Fermi para la simulación.
- **src/math_aux.py**: Contiene funciones matemáticas especiales (armónicos esféricos, derivadas radiales, etc.) necesarias para los cálculos físicos.
- **src/laplacian.py**: Expone las funciones matemáticas principales de `math_aux.py` para mantener una interfaz modular.

## Archivos de datos

Todos los archivos de datos necesarios para la simulación deben estar en `data/Wavefunction/`. Ejemplo de archivos requeridos:
- rvsAS.dat
- rvsAP.dat
- rvsR38s.dat
- rvsR36d.dat
- rvsR37p.dat
- rvsDR38s.dat
- rvsDR36d.dat
- rvsDR37p.dat
- exp_val_r.txt

## Ejecución del proyecto

Desde el directorio raíz del proyecto (`trimero_mod/`), ejecuta:
```sh
poetry run python src/main.py
```
Esto ejecutará la simulación principal con los parámetros definidos en `src/main.py`.

### Ejemplo de prueba rápida
En `src/main.py` hay una función de prueba que ejecuta la simulación con parámetros pequeños para validar el flujo:
```python
def test_trimer_energies_field():
    ...
```
Puedes descomentar la llamada a esta función para probar el flujo con `n1=5` y `dc_field_au=0.1`.

## Flujo principal del código
1. **Carga de datos**: Se leen todos los archivos de datos necesarios usando `numpy`.
2. **Construcción de la matriz de campo**: Se calcula la matriz de interacción de campo eléctrico.
3. **Bucle principal**: Para cada posición radial relevante, se construye la matriz Hamiltoniana usando la lógica de casos físicos (A, B, C, D), se diagonaliza y se guardan los autovalores.
4. **Salida**: Los resultados se guardan en archivos `.dat` en el directorio de ejecución.

## Ejemplo de resultados

Tras la ejecución, se generan archivos como:
- `Trimer_R_sp_wave_N35_R_300_GHz.dat`
- `Trimer_R_sp_wave_N35_R_300_au.dat`

Cada línea de estos archivos contiene el valor de R y los autovalores de la matriz Hamiltoniana para ese punto, por ejemplo:
```
1.234567890123456	-0.000123	-0.000456	...
1.345678901234567	-0.000234	-0.000567	...
...
```
Puedes analizar estos archivos con Python, Excel, gnuplot, etc.

## Trazabilidad y depuración
El código incluye mensajes `print` en los puntos clave del flujo para facilitar la trazabilidad y depuración.

## Extensión y mantenimiento
- Puedes modificar los parámetros físicos en `src/main.py` o `src/trimer.py`.
- Para agregar nuevos potenciales o física, extiende las clases en `src/atom.py` y `src/fermi_potentials.py`.
- Las funciones matemáticas pueden ampliarse en `src/math_aux.py` y exponerse en `src/laplacian.py`.

## Cómo contribuir
1. Haz un fork del repositorio y crea una rama para tu mejora o corrección.
2. Asegúrate de que tu código siga la estructura y estilo del proyecto.
3. Añade pruebas o ejemplos si es relevante.
4. Haz un pull request describiendo claramente tu aporte.

## Preguntas frecuentes (FAQ)

**¿Qué hago si falta un archivo de datos?**
- Verifica que el archivo esté en `data/Wavefunction/` y que el nombre sea exactamente igual (mayúsculas/minúsculas).

**¿Puedo usar mis propios datos?**
- Sí, solo asegúrate de que el formato y dimensiones sean compatibles con los archivos de ejemplo.

**¿Cómo cambio los parámetros físicos?**
- Modifica los valores en `src/main.py` o directamente en la llamada a `Trimer_energies_field`.

**¿Qué hago si obtengo un error de memoria?**
- Prueba con valores más pequeños de `n1` o ejecuta el código en una máquina con más RAM.

**¿Puedo paralelizar el cálculo?**
- El código actual es secuencial, pero puedes paralelizar el bucle principal usando multiprocessing o herramientas similares de Python.

## Contacto
Para dudas o mejoras, contacta al autor original o al responsable de la migración.
