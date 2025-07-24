import numpy as np
from atom import Atom
from fermi_potentials import FermiPotentials
from trimer import Trimer_energies_field

# Parámetros globales (deben definirse según el caso)
pqn = 35  # Ejemplo, ajustar según necesidad
dc_field_au = 0.0  # Ejemplo, ajustar según necesidad
pi = np.pi
EhtoGHz = 6.579683920729e9  # Ajustar si es necesario

if __name__ == "__main__":
    Trimer_energies_field(pqn, dc_field_au)

# Ejemplo de prueba (no sobrescribe el main real)
def test_trimer_energies_field():
    try:
        n1 = 5
        dc_field_au = 0.1
        print("Ejecutando prueba con n1=5 y dc_field_au=0.1...")
        Trimer_energies_field(n1, dc_field_au)
        print("Prueba completada exitosamente.")
    except Exception as e:
        print(f"Error durante la prueba: {e}")

# Descomenta la siguiente línea para ejecutar la prueba:
test_trimer_energies_field() 