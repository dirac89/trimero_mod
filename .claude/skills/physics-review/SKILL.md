---
name: physics-review
description: Revisa cambios de física antes de validarlos mediante simulación
---

# Skill: physics-review

## Propósito
Analizar cambios en la física del código (matriz Hamiltoniana, potenciales, etc.) para detectar inconsistencias antes de ejecutar simulaciones costosas.

## Trigger
- Usuario dice: "revisa la física", "physics check", "valida el cambio de potencial"
- O: `/physics-review`

## Checklist de Revisión

### 1. **Matriz Hamiltoniana**
- [ ] H es Hermitiana (H† = H)
- [ ] Dimensión correcta para los casos A/B/C/D
- [ ] Elementos diagonales son reales (energías individuales)
- [ ] Elementos off-diagonal corresponden a interacciones

### 2. **Potenciales de Fermi**
- [ ] V(r) es real y decreciente con r
- [ ] V(r → ∞) → 0 (decaimiento correcto)
- [ ] Parámetros a0 coherentes con unidades atómicas
- [ ] Simetría correcta para casos de solapamiento

### 3. **Funciones de Onda y Armónicos Esféricos**
- [ ] Y_lm(θ, φ) normalizadas correctamente
- [ ] Derivadas radiales dψ/dr bien calculadas (diferencias finitas)
- [ ] Overlaps entre funciones de onda < 1
- [ ] Integración numérica con precisión suficiente

### 4. **Campo Eléctrico**
- [ ] Matriz de campo diagonal (interacción simple)
- [ ] Elementos escalados correctamente por n1 y amplitud del campo
- [ ] Coherencia con convención de signos

### 5. **Unidades y Conversiones**
- [ ] Entrada en Bohr y Hartree (unidades atómicas)
- [ ] Conversión a GHz solo en salida (si se aplica)
- [ ] Factores de conversión verificados (1 Hartree ≈ 27.2 eV)

### 6. **Estabilidad Numérica**
- [ ] Matriz Hamiltoniana bien condicionada (número de condición < 10^10)
- [ ] Sin términos divergentes para r → 0
- [ ] Precisión de máquina suficiente para eigenvalue solver

## Flujo

1. **Identificar cambios**:
   - Inspecciona `git diff` o archivos modificados
   - Enfoca en: `trimer.py`, `fermi_potentials.py`, `math_aux.py`, `atom.py`

2. **Verificar matemática**:
   - Revisa derivaciones en comentarios
   - Compara con referencias teóricas

3. **Inspeccionar código**:
   - Valida operaciones vectoriales (shapes de numpy arrays)
   - Verifica ciclos de sumación correctos

4. **Proponer mejoras**:
   - Simplificar expresiones redundantes
   - Optimizar operaciones costosas

## Salida

Un reporte estructurado:
```
✓ Matriz Hamiltoniana: Verificado Hermitiano, dimensión correcta
✓ Potenciales: Decrecimiento correcto, parámetros coherentes
⚠ Armónicos esféricos: Revisar normalización en Y_22
✗ Campo eléctrico: Elemento (1,5) parece incorrecto, verificar cálculo

Recomendación: Ejecutar /quick-test antes de /run-simulation
```

## Ejemplo de Revisión

Si el usuario hace cambio en `fermi_potentials.py`:

```python
# Cambio propuesto:
def get_potential(self, r_au):
    # Antes: return self.a * np.exp(-r_au / self.a0)
    # Ahora:
    return self.a * np.exp(-r_au / self.a0) * (1 + self.b * r_au)
```

Claude verifica:
1. ¿Sigue siendo decreciente? ✓ (si self.b es pequeño)
2. ¿Límite r → 0? Sigue siendo finito (buen signo)
3. ¿Límite r → ∞? Sigue decayendo (correcto)
4. ¿Impacto en matriz H? Se recalculan elementos off-diagonal

Recomendación: **Válido**, ejecutar `/quick-test` para verificar impacto numérico.

## Notas

- No ejecuta código, solo revisa
- Útil para cambios teóricos grandes
- Acelera debugging evitando simulaciones inútiles
