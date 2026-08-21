---
name: run-simulation
description: Ejecuta la simulación del trimero con parámetros configurables
---

# Skill: run-simulation

## Propósito
Ejecutar un barrido con parámetros específicos de física, en el sistema correcto.

## Trigger
- Usuario dice: "ejecuta la simulación", "run simulation con...", "simula con n1=X"
- O explícitamente: `/run-simulation`

## Parametrización

Solicita al usuario (si no están claros en la solicitud):
1. **n1**: Número máximo de quanta en el estado (ej: 5, 10, 35)
2. **dc_field_au**: Campo eléctrico en unidades atómicas (ej: 0.1, 0.5)
3. **radio_min, radio_max**: Rango de distancias a simular
4. **num_points**: Número de puntos de radio en los que evaluar

## Flujo

1. Valida que existan los archivos en `data/Wavefunction/`:
   - `rvsAS.dat`, `rvsAP.dat`, `rvsR38s.dat`, `rvsR36d.dat`, `rvsR37p.dat`
   - Archivos de derivadas: `rvsDR*.dat`
   - `exp_val_r.txt`

2. **Decide primero de qué sistema físico se trata** (ver `.claude/ARCHITECTURE.md`):

   **Rb*-KRb polar — el camino VIGENTE**:
   ```bash
   poetry run python scripts/compute_bop_curve.py \
       --n-manifold 25 --mj 0 1 --rmin 400 --rmax 1800 --step 5
   ```
   Los `.dat` de `data/Wavefunction/` NO hacen falta aquí: este sistema no usa
   longitudes de dispersión.

   **Perturbador neutro — camino legado, congelado**:
   ```python
   from trimero.systems.rb_neutral_perturber.trimer import Trimer_energies_field
   Trimer_energies_field(n1=..., dc_field_au=...)
   ```

3. Monitorea la ejecución y reporta:
    - Archivos generados (`plots/rb_krb_polar/fig1_ad_MJ*_n*.npz` + PNG, o
     `Trimer_R_sp_wave_*.dat` en el legado)
   - Número de puntos computados
   - Rango de energías resultantes

4. Valida salida:
   - Verifica que los autovalores sean reales
   - Comprueba que estén ordenados ascendentemente
   - Reporta cualquier NaN o Inf

## Salida Esperada

Archivos `.dat` en el directorio actual con formato:
```
R_value eigenvalue_1 eigenvalue_2 ... eigenvalue_N
1.234   -0.0001     -0.0002     ... -0.5000
1.345   -0.0001     -0.0003     ... -0.5001
...
```

## Ejemplo de Uso

```
Usuario: Ejecuta la simulación con n1=10 y dc_field_au=0.2
> /run-simulation

Claude: Voy a ejecutar la simulación con:
- n1: 10
- dc_field: 0.2 (unidades atómicas)
- Radio: [1.5, 5.0] Bohr, 50 puntos

[Validando datos...]
[Ejecutando...]
[Guardando resultados...]

✓ Completado: Trimer_R_sp_wave_N10_R_200_au.dat (50 puntos, 20 niveles de energía)
```

## Notas

- El tiempo de ejecución depende de n1 y num_points
- Para n1 > 30, considera paralelización
- Si hay error de memoria, reduce n1 o num_points
