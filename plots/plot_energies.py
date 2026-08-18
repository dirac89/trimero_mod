import numpy as np
import matplotlib.pyplot as plt
import sys
import os

# Parámetros configurables
output_file = 'Trimer_R_sp_wave_N35_R_300_GHz.dat'  # Ruta relativa desde la raíz del proyecto
num_levels = 5  # Número de niveles de energía a graficar

# Permitir pasar el archivo y número de niveles por argumentos
if len(sys.argv) > 1:
    output_file = sys.argv[1]
if len(sys.argv) > 2:
    num_levels = int(sys.argv[2])

# Cargar datos
try:
    data = np.loadtxt(output_file)
    # Si solo hay una fila, convertir a 2D
    if data.ndim == 1:
        data = data[np.newaxis, :]
except Exception as e:
    print(f"Error al leer el archivo {output_file}: {e}")
    sys.exit(1)

R = data[:, 0]
energies = data[:, 1:]

plt.figure(figsize=(8, 6))
for i in range(min(num_levels, energies.shape[1])):
    plt.plot(R, energies[:, i], label=f'Nivel {i+1}')

plt.xlabel('R (a.u.)')
plt.ylabel('Energía (GHz o a.u.)')
plt.title(f'Niveles de energía del trimero vs. distancia radial\nArchivo: {os.path.basename(output_file)}')
plt.legend()
plt.grid(True)
plt.tight_layout()

# Guardar la figura en el mismo directorio que este script
plot_path = os.path.join(os.path.dirname(__file__), 'energies_vs_R.png')
plt.savefig(plot_path)
print(f"Gráfico guardado en {plot_path}")
plt.show() 