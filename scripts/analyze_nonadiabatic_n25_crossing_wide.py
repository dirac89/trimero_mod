#!/usr/bin/env python3
"""
Fase 6b de docs/PLAN_nonadiabatic_dynamics.md: cierra el pendiente §7.1 de
docs/analysis_fase6_aplicacion_n25.md — ampliar la ventana de R del cruce
evitado real (n=25, M_J=0, estados 54/55) más allá de los 250 a0
(R∈[1180,1430]) de la Fase 6 original.

Ventana nueva: R∈[500,1800] a0 (1300 a0, 5.2x más ancha), el mismo dominio
ya validado de la Fig. 1 completa (docs/STATUS.md). Coste estimado y
reportado ANTES de ejecutar: ~24.7 min a 1.138 s/punto (medido en la
Fase 6); ejecutado tras confirmar que estaba bajo el umbral de 30-40 min
acordado, sin pedir confirmación adicional.

Corre el pipeline completo (`run_pipeline` de
analyze_nonadiabatic_n25_crossing.py) sobre la ventana nueva, Y TAMBIÉN
sobre el SUBCONJUNTO de la ventana nueva que coincide con la ventana
original [1180,1430] -- comparación de consistencia barata antes de
aceptar los resultados nuevos como extensión válida (no una recomputación
independiente: son literalmente los mismos R, W, V, sólo recortados de la
malla nueva en vez de recalculados).
"""
import numpy as np

from analyze_nonadiabatic_n25_crossing import (
    run_pipeline, GHZ_PER_HARTREE, E_MANIFOLD_N25, SIGMA_REAL_RUTTLEY, R0_CROSSING,
)


def main(npz_path, out_npz, old_window=(1180.0, 1430.0)):
    d = np.load(npz_path)
    R, W, V = d["R"], d["W"], d["V"]

    print("=" * 78)
    print("PIPELINE SOBRE LA VENTANA NUEVA (500-1800 a0)")
    print("=" * 78)
    wide = run_pipeline(R, W, V, L_min=350.0, n_L=30, n_keep=300, label="ancho")

    print("\n" + "=" * 78)
    print("TEST DE CONSISTENCIA: mismo recorte que la Fase 6 original (1180-1430 a0)")
    print("=" * 78)
    mask = (R >= old_window[0]) & (R <= old_window[1])
    R_old, W_old, V_old = R[mask], W[mask], V[mask]
    redo = run_pipeline(R_old, W_old, V_old, n_L=20, n_keep=15, label="recorte-viejo")

    # --- comparacion explicita contra los numeros ya publicados de la Fase 6 ---
    print("\n" + "=" * 78)
    print("COMPARACION EXPLICITA vs. Fase 6 original (docs/analysis_fase6_aplicacion_n25.md)")
    print("=" * 78)
    orig_peak, orig_R = 5.0672e-01, 1305.0
    idx_peak = np.argmax(np.abs(redo["A"][:, 0, 1]))
    new_peak = np.abs(redo["A"][idx_peak, 0, 1])
    new_R = redo["R_mid"][idx_peak]
    print(f"pico |A(54,55)|: original={orig_peak:.4e} en R={orig_R:.1f}  "
          f"recorte-nuevo={new_peak:.4e} en R={new_R:.1f}  "
          f"dif_rel={abs(new_peak-orig_peak)/orig_peak:.2%}")

    orig_depths_ghz = [-16.2750, -15.6131, -15.0026, -14.4235, -13.8680,
                       -13.3298, -12.8021, -12.2792, -11.7567, -10.8593,
                       -10.7624, -9.8184, -9.8081, -9.0722, -8.9127]
    new_depths_ghz = sorted((s["E_mean"] - E_MANIFOLD_N25) * GHZ_PER_HARTREE
                            for s in redo["bound"])
    print(f"\nestados ligados: original={len(orig_depths_ghz)}  "
          f"recorte-nuevo={len(new_depths_ghz)}")
    print(f"  original      : {['%.4f' % v for v in orig_depths_ghz]}")
    print(f"  recorte-nuevo : {['%.4f' % v for v in new_depths_ghz]}")
    if len(orig_depths_ghz) == len(new_depths_ghz):
        diffs = [abs(a - b) for a, b in zip(orig_depths_ghz, new_depths_ghz)]
        print(f"  diferencia máxima por estado: {max(diffs):.4f} GHz")

    # --- resumen de la ventana ancha: cuantos estados con Gamma evaluable ---
    print("\n" + "=" * 78)
    print("RESUMEN VENTANA ANCHA")
    print("=" * 78)
    print(f"estados ligados totales: {len(wide['bound'])}")
    print(f"estados con Gamma evaluable: {len(wide['gammas'])}")
    for g in wide["gammas"]:
        v = (g["E_mean"] - E_MANIFOLD_N25) * GHZ_PER_HARTREE
        print(f"    V-E_manifold={v:9.4f} GHz   Gamma={g['gamma']*GHZ_PER_HARTREE*1e3:.4f} MHz")
    print(f"sigma_real={SIGMA_REAL_RUTTLEY} a0 cabe con margen 5sigma "
          f"en esta ventana ({wide['window_span']:.0f} a0): {wide['sigma_fits']}")

    np.savez(out_npz, R_mid=wide["R_mid"], A=wide["A"], B=wide["B"], W_mid=wide["W_mid"],
             bound_summary=wide["bound"], gammas=wide["gammas"])
    print(f"\nresultados guardados en {out_npz}")
    return wide, redo


if __name__ == "__main__":
    import sys
    npz_path = sys.argv[1] if len(sys.argv) > 1 else "fase6b_sweep_n25_MJ0_states54_55.npz"
    out = sys.argv[2] if len(sys.argv) > 2 else "plots/hybrid_neutral_polar/data/fase6b_n25_results.npz"
    main(npz_path, out)
