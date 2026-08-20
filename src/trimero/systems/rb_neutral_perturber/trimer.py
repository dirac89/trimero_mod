import numpy as np
from trimero.systems.rb_atom import Atom
from trimero.systems.rb_neutral_perturber.fermi_potentials import FermiPotentials

def Trimer_energies_field(n1, dc_field_au):
    rubidium = Atom(n1, 0)  # oqn no está definido, ajustar si es necesario
    lc = 2
    cols = 2
    rows = 776
    theta = 0.0
    pi = np.pi
    theta1 = pi
    EhtoGHz = 6.579683920729e9  # Ajustar si es necesario

    # Leer matrices desde archivos (rutas y nombres corregidos y consistentes)
    data_dir = "data/Wavefunction/"
    print("Cargando archivos de datos...")
    As = np.loadtxt(data_dir + "rvsAS.dat"); print("Cargado rvsAS.dat")
    Ap = np.loadtxt(data_dir + "rvsAP.dat"); print("Cargado rvsAP.dat")
    As1 = np.loadtxt(data_dir + "rvsAS.dat"); print("Cargado rvsAS.dat (As1)")
    Ap1 = np.loadtxt(data_dir + "rvsAP.dat"); print("Cargado rvsAP.dat (Ap1)")
    R38s = np.loadtxt(data_dir + "rvsR38s.dat"); print("Cargado rvsR38s.dat")
    R36d = np.loadtxt(data_dir + "rvsR36d.dat"); print("Cargado rvsR36d.dat")
    R37p = np.loadtxt(data_dir + "rvsR37p.dat"); print("Cargado rvsR37p.dat")
    DR38s = np.loadtxt(data_dir + "rvsDR38s.dat"); print("Cargado rvsDR38s.dat")
    DR36d = np.loadtxt(data_dir + "rvsDR36d.dat"); print("Cargado rvsDR36d.dat")
    DR37p = np.loadtxt(data_dir + "rvsDR37p.dat"); print("Cargado rvsDR37p.dat")

    print(f"Manifold n1 = {n1} plots of E(R) type are included!")

    # Leer matriz radial
    print("Cargando exp_val_r.txt...")
    Radial = np.zeros((n1, n1))
    with open(data_dir + "exp_val_r.txt", "r") as infile:
        for i in range(n1 - 1):
            value = float(infile.readline().strip())
            Radial[i, i + 1] = value
    print("done reading radial field")

    max_dim = n1 * n1
    field = np.zeros((max_dim, max_dim))

    print("Construyendo matriz de campo...")
    for i in range(n1 - 1):
        for k1 in range(1, 2 * i + 2):
            m1 = k1 - i - 1
            for j in range(n1):
                for k2 in range(1, 2 * j + 2):
                    m2 = k2 - j - 1
                    if (i == j + 1) and (m1 == m2):
                        r_exp = Radial[j, i]
                        field[i * i + k1 - 1, j * j + k2 - 1] = rubidium.Vfield(i, j, m1, m2, r_exp, dc_field_au)
                    if (i == j - 1) and (m1 == m2):
                        r_exp = Radial[i, j]
                        field[i * i + k1 - 1, j * j + k2 - 1] = rubidium.Vfield(i, j, m1, m2, r_exp, dc_field_au)
    print("done field matrix trimer full nlm")

    # Abrir archivos de salida
    print("Abriendo archivos de salida...")
    store2DShift = open("Trimer_R_sp_wave_N35_R_300_GHz.dat", "w")
    store2DShift_au = open("Trimer_R_sp_wave_N35_R_300_au.dat", "w")

    print("Iniciando bucle principal sobre filas...")
    for row in range(297, rows):
        row1 = row
        R = As[row, 0]
        AS = As[row, 1]
        AP = Ap[row, 1]
        R1 = As1[row1, 0]
        AS1 = As1[row1, 1]
        AP1 = Ap1[row1, 1]
        print(f"\n[Row {row}] Construyendo la matriz Hamiltoniana...")
        spV = np.zeros((max_dim, max_dim))
        s = 1  # s=0 solo s-wave, s=1 s-wave + p-wave
        for i in range(n1):
            for k1 in range(1, 2 * i + 2):
                m1 = k1 - i - 1
                for j in range(n1):
                    for k2 in range(1, 2 * j + 2):
                        m2 = k2 - j - 1
                        delta_ij = 1.0 if (i == j and m1 == m2) else 0.0
                        field_contri = 0.0
                        if (i == j + 1) and (m1 == m2):
                            field_contri = field[i * i + k1 - 1, j * j + k2 - 1]
                        if (i == j - 1) and (m1 == m2):
                            field_contri = field[i * i + k1 - 1, j * j + k2 - 1]
                        # CASO A: i < 3 y j < 3
                        if i < 3 and j < 3:
                            if i == 1 and j == 1:
                                wave_R = R37p[row, 1]
                                Dwave_R = DR37p[row, 1]
                                wave_R1 = R37p[row1, 1]
                                Dwave_R1 = DR37p[row1, 1]
                            elif i == 2 and j == 2:
                                wave_R = R36d[row, 1]
                                Dwave_R = DR36d[row, 1]
                                wave_R1 = R36d[row1, 1]
                                Dwave_R1 = DR36d[row1, 1]
                            elif i == 0 and j == 0:
                                wave_R = R38s[row, 1]
                                Dwave_R = DR38s[row, 1]
                                wave_R1 = R38s[row1, 1]
                                Dwave_R1 = DR38s[row1, 1]
                            else:
                                wave_R = 0.0
                                Dwave_R = 0.0
                                wave_R1 = 0.0
                                Dwave_R1 = 0.0
                            n11 = n1 + 3 - i
                            n21 = n1 + 3 - i
                            n12 = n1 + 3 - j
                            n22 = n1 + 3 - j
                            wave_1 = wave_R
                            wave1_2 = 0.0
                            Dwave_1 = Dwave_R
                            Dwave1_2 = 0.0
                            wave_2 = wave_R1
                            wave2_1 = 0.0
                            Dwave_2 = Dwave_R1
                            Dwave2_1 = 0.0
                        # CASO B: i < 3 y j > 2
                        elif i < 3 and j > 2:
                            if i == 1:
                                wave_R = R37p[row, 1]
                                Dwave_R = DR37p[row, 1]
                                wave_R1 = R37p[row1, 1]
                                Dwave_R1 = DR37p[row1, 1]
                            elif i == 2:
                                wave_R = R36d[row, 1]
                                Dwave_R = DR36d[row, 1]
                                wave_R1 = R36d[row1, 1]
                                Dwave_R1 = DR36d[row1, 1]
                            elif i == 0:
                                wave_R = R38s[row, 1]
                                Dwave_R = DR38s[row, 1]
                                wave_R1 = R38s[row1, 1]
                                Dwave_R1 = DR38s[row1, 1]
                            else:
                                wave_R = 0.0
                                Dwave_R = 0.0
                                wave_R1 = 0.0
                                Dwave_R1 = 0.0
                            n11 = n1 + 3 - i
                            n21 = n1 + 3 - i
                            n12 = n1
                            n22 = n1
                            wave_1 = wave_R
                            wave1_2 = 0.0
                            Dwave_1 = Dwave_R
                            Dwave1_2 = 0.0
                            wave_2 = wave_R1
                            wave2_1 = 0.0
                            Dwave_2 = Dwave_R1
                            Dwave2_1 = 0.0
                        # CASO C: i > 2 y j < 3
                        elif i > 2 and j < 3:
                            if j == 0:
                                wave_R = R38s[row, 1]
                                Dwave_R = DR38s[row, 1]
                                wave_R1 = R38s[row1, 1]
                                Dwave_R1 = DR38s[row1, 1]
                            elif j == 1:
                                wave_R = R37p[row, 1]
                                Dwave_R = DR37p[row, 1]
                                wave_R1 = R37p[row1, 1]
                                Dwave_R1 = DR37p[row1, 1]
                            elif j == 2:
                                wave_R = R36d[row, 1]
                                Dwave_R = DR36d[row, 1]
                                wave_R1 = R36d[row1, 1]
                                Dwave_R1 = DR36d[row1, 1]
                            else:
                                wave_R = 0.0
                                Dwave_R = 0.0
                                wave_R1 = 0.0
                                Dwave_R1 = 0.0
                            n11 = n1
                            n21 = n1
                            n12 = n1 + 3 - j
                            n22 = n1 + 3 - j
                            wave_1 = 0.0
                            wave1_2 = wave_R
                            Dwave_1 = 0.0
                            Dwave1_2 = Dwave_R
                            wave_2 = 0.0
                            wave2_1 = wave_R1
                            Dwave_2 = 0.0
                            Dwave2_1 = Dwave_R1
                        # CASO D: i > lc y j > lc
                        elif i > lc and j > lc:
                            n11 = n1
                            n21 = n1
                            n12 = n1
                            n22 = n1
                            wave_1 = 0.0
                            wave1_2 = 0.0
                            Dwave_1 = 0.0
                            Dwave1_2 = 0.0
                            wave_2 = 0.0
                            wave2_1 = 0.0
                            Dwave_2 = 0.0
                            Dwave2_1 = 0.0
                        # Otros casos (por defecto)
                        else:
                            n11 = n1 + 3 - i
                            n21 = n1 + 3 - i
                            n12 = n1 + 3 - j
                            n22 = n1 + 3 - j
                            wave_1 = 0.0
                            wave1_2 = 0.0
                            Dwave_1 = 0.0
                            Dwave1_2 = 0.0
                            wave_2 = 0.0
                            wave2_1 = 0.0
                            Dwave_2 = 0.0
                            Dwave2_1 = 0.0
                        rubidium_mod = Atom(n11, i)
                        fermi_1 = FermiPotentials(s, n11, n12, i, j, m1, m2, R, theta, AS, wave_1, wave1_2, AP, Dwave_1, Dwave1_2)
                        fermi_2 = FermiPotentials(s, n21, n22, i, j, m1, m2, R1, theta1, AS1, wave_2, wave2_1, AP1, Dwave_2, Dwave2_1)
                        spV[i * i + k1 - 1, j * j + k2 - 1] = (
                            rubidium_mod.E_Rb() * delta_ij +
                            fermi_1.Vsp() +
                            fermi_2.Vsp() +
                            field_contri
                        )
        print(f"[Row {row}] Antes de la diagonalización...")
        evalues, evectors = np.linalg.eigh(spV)
        print(f"[Row {row}] Después de la diagonalización.")
        VShift = evalues[-1]
        print(f"[Row {row}] Guardando resultados...")
        store2DShift.write(f"{R:.15f}\t")
        store2DShift_au.write(f"{R:.15f}\t")
        for i in range(max_dim):
            store2DShift.write(f"{(evalues[max_dim - 1 - i] - rubidium.E_Rb()) * EhtoGHz}\t")
        store2DShift.write("\n")
        for i in range(max_dim):
            store2DShift_au.write(f"{(evalues[max_dim - 1 - i] - rubidium.E_Rb())}\t")
        store2DShift_au.write("\n")
    print("Cerrando archivos de salida...")
    store2DShift.close()
    store2DShift_au.close()
    print("Ejecución de Trimer_energies_field finalizada.") 