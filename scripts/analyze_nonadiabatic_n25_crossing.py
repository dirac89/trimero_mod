#!/usr/bin/env python3
"""
Fase 6 de docs/PLAN_nonadiabatic_dynamics.md: primera aplicación del
pipeline completo (Fases 1-5) a un dato real — el cruce evitado de
BOPSystem(n_manifold=25), M_J=0, estados 54/55, ya caracterizado en
docs/analysis_fase1_acoplamiento_derivada.md (R≈1305 a0).

Decisión explícita (ver docs/analysis_fase6_aplicacion_n25.md §1): se usa
n=25 (ya validado en las Fases 1-2) como DEMOSTRACIÓN del pipeline
completo, NO como predicción cuantitativa para el experimento real de
Ruttley et al. 2023 (que usa Rb(52s), no n=25).

Lee el barrido electrónico precalculado (R, W, V para los estados 54/55) y
aplica, en orden: Fase 1 (acoplamiento de derivada) -> Fase 2 (canales
acoplados) -> Fase 3 (estabilización, canal desacoplado) -> Fase 4 (tasas
de decaimiento) -> Fase 5 (Franck-Condon).

`run_pipeline` es la lógica reutilizable (la usa también
`analyze_nonadiabatic_n25_crossing_wide.py`, Fase 6b, con una ventana más
ancha y el paso extra de comparar contra estos resultados originales).
"""
import numpy as np

from trimero.systems.nonadiabatic_dynamics.coupling import (
    derivative_coupling, fix_eigenvector_signs,
)
from trimero.systems.nonadiabatic_dynamics.coupled_channels import (
    second_derivative_coupling, build_coupled_hamiltonian,
)
from trimero.systems.nonadiabatic_dynamics.stabilization import (
    stabilization_scan, track_trajectories, classify_stability, stable_state_summary,
)
from trimero.systems.nonadiabatic_dynamics.decay_rates import nonadiabatic_decay_rate
from trimero.systems.nonadiabatic_dynamics.franck_condon import (
    gaussian_wavepacket, franck_condon_factor,
)

GHZ_PER_HARTREE = 6.579683920502e6

# masa reducida Rb-87 + RbCs, ver docs/analysis_fase6_aplicacion_n25.md §3
MU_RB_RBCS = 113536.332476  # m_e

# cero de energia de BOPSystem(n_manifold=25): E(manifold)+KRb(N=0) = -0.5/n^2,
# exacto para l>=3 (rb_krb_polar/bop_system.py). Restar esto antes de
# convertir a GHz es lo que hace compute_bop_curve.py para reportar "V(R)"
# -- sin restarlo, las energias reportadas aqui abajo son solo el offset
# absoluto del Rydberg (~ -5264 GHz), no la profundidad de pozo relevante.
E_MANIFOLD_N25 = -0.5 / 25.0**2  # Eh

R0_CROSSING = 1305.0   # a0, el propio cruce evitado (Fase 1)
SIGMA_REAL_RUTTLEY = 945.0  # a0 (~50 nm), Ruttley et al. 2023


def load_sweep(path):
    d = np.load(path)
    return d["R"], d["W"], d["V"]


def run_pipeline(R, W, V, R0=R0_CROSSING, sigma_real=SIGMA_REAL_RUTTLEY,
                  L_min=60.0, n_L=20, n_keep=15, label=""):
    """
    Pipeline Fases 1-5 completo sobre un barrido (R, W, V) ya calculado.

    Returns:
        dict con R_mid, A, B, W_mid, bound (resumen de Fase 3), gammas
        (Fase 4), fc (Fase 5) -- todo lo necesario para inspección o para
        el test de consistencia de la Fase 6b.
    """
    h = R[1] - R[0]
    tag = f"[{label}] " if label else ""
    print(f"{tag}malla: {len(R)} puntos, R=[{R[0]:.1f},{R[-1]:.1f}] a0, h={h:.2f} a0")

    # --- Fase 1 ---
    V_fixed = V.copy()
    for k in range(1, len(R)):
        V_fixed[k] = fix_eigenvector_signs(V_fixed[k], V_fixed[k - 1])
    R_mid, A = derivative_coupling(R, V=V_fixed)
    _, B = second_derivative_coupling(R, V=V_fixed)
    W_mid = W[1:-1]
    idx_peak = np.argmax(np.abs(A[:, 0, 1]))
    print(f"{tag}Fase 1: |A(54,55)| max = {np.abs(A[idx_peak,0,1]):.4e} "
          f"en R={R_mid[idx_peak]:.1f} a0")

    # --- Fase 2 ---
    H_coupled = build_coupled_hamiltonian(R_mid, W_mid, A=A, B=B, mu=MU_RB_RBCS)
    print(f"{tag}Fase 2: H acoplada ensamblada, forma {H_coupled.shape}")

    # --- Fase 3: canal desacoplado d (=54) ---
    def hamiltonian_builder_d(L):
        n_box = max(3, int(round(L / h)))
        n_box = min(n_box, len(R_mid))
        Rb_ = R_mid[:n_box]
        Vb_ = W_mid[:n_box, 0][:, None]
        H = build_coupled_hamiltonian(Rb_, Vb_, mu=MU_RB_RBCS)
        return Rb_, H

    L_values = np.linspace(L_min, R_mid[-1] - R_mid[0], n_L)
    L_values, E_box = stabilization_scan(L_values, hamiltonian_builder_d, n_keep=n_keep)
    traj = track_trajectories(E_box)

    # CORRECCION IMPORTANTE (ver docs/analysis_fase6b_ventana_ampliada.md
    # §2): classify_stability compara |dE/dL| contra 2|E|/L, una formula
    # derivada asumiendo E referenciada al umbral de disociacion (E->0
    # lejos del pozo). Las energias de BOPSystem son ABSOLUTAS (~-8e-4 Eh,
    # dominadas por el offset del manifold Rydberg, ~-5264 GHz), no por la
    # escala de GHz relevante -- sin corregir esto, 2|E|/L es una
    # referencia absurdamente grande y (casi) cualquier estado sale
    # "estable". Se resta el umbral asintotico REAL del canal d, medido
    # directamente de la cola plana de W_mid[:,0] (ultimo 10% de la
    # ventana), antes de clasificar -- el desplazamiento es una constante,
    # no afecta a dE/dL, solo corrige la escala de comparacion.
    tail_n = max(50, len(W_mid) // 10)
    E_threshold_d = float(np.mean(W_mid[-tail_n:, 0]))
    L_mid_stab, stable, slope = classify_stability(
        L_values, traj - E_threshold_d, threshold_ratio=0.1)
    summary = stable_state_summary(L_mid_stab, traj[1:-1], stable, min_fraction=0.5)
    bound = sorted([s for s in summary if s["is_bound"]], key=lambda s: s["E_mean"])

    print(f"{tag}Fase 3: estados clasificados ligados: {len(bound)} (de n_keep={n_keep})")
    for s in bound:
        V_rel_ghz = (s["E_mean"] - E_MANIFOLD_N25) * GHZ_PER_HARTREE
        print(f"{tag}    E={s['E_mean']:.8e} Eh  V-E_manifold={V_rel_ghz:9.4f} GHz"
              f"  frac_estable={s['fraction_stable']:.2f}")
    if bound:
        depths = [(s["E_mean"] - E_MANIFOLD_N25) * GHZ_PER_HARTREE for s in bound]
        print(f"{tag}    rango V-E_manifold: [{min(depths):.4f}, {max(depths):.4f}] GHz")

    # --- Fase 4 ---
    gammas = []
    if bound:
        V_d, V_u = W_mid[:, 0], W_mid[:, 1]
        A_du, B_du = A[:, 0, 1], B[:, 0, 1]
        H_d_full = build_coupled_hamiltonian(R_mid, V_d[:, None], mu=MU_RB_RBCS)
        E_d_full = np.linalg.eigvalsh(H_d_full)
        for s in bound:
            idx = int(np.argmin(np.abs(E_d_full - s["E_mean"])))
            V_rel_ghz = (s["E_mean"] - E_MANIFOLD_N25) * GHZ_PER_HARTREE
            try:
                gamma, diag = nonadiabatic_decay_rate(
                    R_mid, V_d, V_u, A_du, B_du, idx, mu=MU_RB_RBCS)
                gammas.append({"E_mean": s["E_mean"], "gamma": gamma, "diag": diag})
                print(f"{tag}Fase 4: V-E_manifold={V_rel_ghz:.4f} GHz  "
                      f"Gamma={gamma*GHZ_PER_HARTREE*1e3:.4f} MHz")
            except ValueError:
                print(f"{tag}Fase 4: estado V-E_manifold={V_rel_ghz:.4f} GHz -- "
                      f"por debajo del mínimo de V_u en esta ventana, no evaluable")

    # --- Fase 5 ---
    window_span = R_mid[-1] - R_mid[0]
    fits = R0 - 5 * sigma_real >= R_mid[0] and R0 + 5 * sigma_real <= R_mid[-1]
    if not fits:
        print(f"{tag}Fase 5: AVISO -- paquete gaussiano REAL (R0={R0}, "
              f"sigma={sigma_real} a0) no cabe con margen 5sigma en la ventana "
              f"de {window_span:.0f} a0.")
    fc = []
    if bound:
        u_scat_real = gaussian_wavepacket(R_mid, R0, sigma_real)
        H_d_full = build_coupled_hamiltonian(R_mid, W_mid[:, 0][:, None], mu=MU_RB_RBCS)
        E_d_full, chi_d_full = np.linalg.eigh(H_d_full)
        for s in bound:
            idx = int(np.argmin(np.abs(E_d_full - s["E_mean"])))
            V_rel_ghz = (s["E_mean"] - E_MANIFOLD_N25) * GHZ_PER_HARTREE
            F_real = franck_condon_factor(R_mid, chi_d_full[:, idx], u_scat_real)
            fc.append({"E_mean": s["E_mean"], "F_real": F_real})
            print(f"{tag}Fase 5: V-E_manifold={V_rel_ghz:9.4f} GHz  F(sigma_real)={F_real: .4e}")

    return {"R_mid": R_mid, "A": A, "B": B, "W_mid": W_mid, "bound": bound,
            "gammas": gammas, "fc": fc, "sigma_fits": fits, "window_span": window_span}


def main(npz_path, out_npz):
    R, W, V = load_sweep(npz_path)
    results = run_pipeline(R, W, V)
    np.savez(out_npz, R_mid=results["R_mid"], A=results["A"], B=results["B"],
             W_mid=results["W_mid"], bound_summary=results["bound"])
    print(f"\nresultados guardados en {out_npz}")
    return results


if __name__ == "__main__":
    import sys
    npz_path = sys.argv[1] if len(sys.argv) > 1 else "fase6_sweep_n25_MJ0_states54_55.npz"
    out = sys.argv[2] if len(sys.argv) > 2 else "plots/hybrid_neutral_polar/fase6_n25_results.npz"
    main(npz_path, out)
