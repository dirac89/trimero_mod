---
name: quick-test
description: Ejecución rápida de prueba con parámetros pequeños para validación
---

# Skill: quick-test

## Propósito
Validar que la simulación funciona correctamente con parámetros pequeños (n1=5, pocos puntos).

## Trigger
- Usuario dice: "prueba rápida", "test", "valida el flujo"
- O: `/quick-test`

## Flujo Automático

1. **Setup**:
   - n1 = 5 (muy pequeño, rápido)
   - dc_field_au = 0.1 (campo débil)
   - radio_min = 1.5 Bohr
   - radio_max = 3.0 Bohr
   - num_points = 5 (solo 5 puntos)

2. **Ejecución**:
   ```bash
   poetry run python src/main.py
   ```
   Con los parámetros pequeños comentados en `main.py`, o:
   ```python
   from src.trimer import Trimer_energies_field
   Trimer_energies_field(n1=5, dc_field_au=0.1, radius_min=1.5, radius_max=3.0, num_points=5)
   ```

3. **Validación**:
   - ✓ Datos cargados correctamente
   - ✓ Matriz Hamiltoniana construida (dimensión esperada)
   - ✓ Diagonalización exitosa
   - ✓ Autovalores reales y ordenados
   - ✓ Archivo `.dat` generado con formato correcto

4. **Reporte**:
   - Tiempo total de ejecución
   - Tamaño de matriz Hamiltoniana
   - Número de autovalores obtenidos
   - Primeros 3 niveles de energía

## Ejemplo de Salida

```
=== Quick Test ===
Parámetros: n1=5, dc_field=0.1 au, R=[1.5, 3.0] Bohr, 5 puntos
Tiempo: 2.34s

✓ Datos cargados
✓ Matriz H: 45x45
✓ Diagonalizaciones: 5/5 exitosas

Primeros autovalores (en Hartree):
  E1 = -0.523
  E2 = -0.521
  E3 = -0.515

✓ Test passou: archivo generado 'Trimer_R_sp_wave_N5_R_150_au.dat'
```

## Cuándo Usar

- Después de cambios en **física** o **matemáticas**
- Para validar que el proyecto está configurado correctamente
- Antes de ejecutar simulaciones grandes

## Tiempo Típico

< 5 segundos en máquina estándar.
