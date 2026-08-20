#!/usr/bin/env python3
"""
Curva BOP adiabática de M_J=0 para un manifold cualquiera de Rb*-KRb.

Generaliza `run_bop_curve.py` (que era n=24 + 27s) vía `BOPSystem`:

    poetry run python scripts/run_bop_manifold.py --n-manifold 25

Base electrónica: manifold (n, l≥3) MÁS los tres vecinos individuales
(n+1)d, (n+2)p, (n+3)s, tal como los define Aguilera-Fernández, Sadeghpour,
Schmelcher & González-Férez, J. Phys.: Conf. Ser. 635, 012023 (2015),
arXiv:1507.07972. `--neighbors 0` reproduce la base incompleta anterior.

H = H_a + H_mol + V_Fermi, delta0_ns = 3.13180, p_interpolation = "inverse".

MÉTODO: curvas ADIABÁTICAS directas, sin seguimiento por solapamiento. Dentro
de un bloque M_J las curvas no se cruzan (regla de no cruce), así que el
k-ésimo autovalor ORDENADO ya ES la k-ésima curva adiabática.

⚠️ Con la base completa NO basta con fijar el índice k que tiene carácter de
manifold en el borde del dominio: los estados de (n+1)d y (n+2)p caen justo en
el rango de energía relevante, y la curva de índice fijo cambia de carácter al
atravesar sus cruces evitados. Se calcula por tanto la curva de CARÁCTER —la más
baja con peso de manifold > 50 % en cada R— y se dibujan las dos para que se vea
dónde se separan. Con la base incompleta las dos coinciden, que es por lo que el
método de índice fijo funcionaba en docs/analysis_curva_bop_MJ0.md §2.

Hace, en este orden:
  1. Monta el sistema y calcula el dominio válido (todos los pares (l1,l2)).
  2. Barre R con `eigh` y guarda las N_CONTEXT curvas más bajas y el peso
     de manifold de cada una (hace falta para identificar la curva por carácter).
  3. Recalcula la ventana de exclusión de la resonancia de forma p PARA ESTE
     manifold (el polo está a la misma ENERGÍA pero a distinto R).
  4. Umbrales asintóticos ns + KRb(N) y check de residuo en el borde.
  5. Características de la curva de manifold, excluyendo la ventana.
  6. Figura con el mismo estilo que plots/archive/bop_curve_MJ0_excluded.png.
"""
import argparse
import os
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from trimero.systems.rb_krb_polar.charge_dipole import B_KRB_GHZ
from trimero.systems.rb_krb_polar.bop_system import GHZ_PER_HARTREE, BOPSystem
from trimero.systems.rb_krb_polar.rb_defects import DELTA0_NS_PAPER as D0

N_CONTEXT = 60
RULE = "=" * 84
DOC = "docs/analysis_interpolacion_polo_Ap.md"


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-manifold", type=int, default=25)
    ap.add_argument("--n-s", type=int, default=None,
                    help="n del estado ns vecino (por defecto n_manifold+3)")
    ap.add_argument("--neighbors", type=int, nargs="+", default=[0, 1, 2],
                    help="l de los vecinos individuales: 0 1 2 = (n+3)s, "
                         "(n+2)p, (n+1)d (base del paper). '0' sola = base "
                         "incompleta de las rondas anteriores")
    ap.add_argument("--tag", default=None, help="sufijo de los ficheros de salida")
    ap.add_argument("--step", type=float, default=5.0, help="paso en R [a0]")
    ap.add_argument("--ymin", type=float, default=-120.0)
    ap.add_argument("--ymax", type=float, default=25.0)
    ap.add_argument("--margin-factor", type=float, default=2.0,
                    help="semianchura de la ventana de exclusión en FWHM")
    ap.add_argument("--reuse", action="store_true",
                    help="reutiliza el .npz en vez de rebarrer R")
    ap.add_argument("--npz", default=None)
    ap.add_argument("--out", default=None)
    return ap.parse_args()


def segment_features(R, y, label):
    """Mínimos/máximos locales estrictos dentro de UN segmento contiguo."""
    mins = [i for i in range(1, len(y) - 1) if y[i] < y[i - 1] and y[i] < y[i + 1]]
    maxs = [i for i in range(1, len(y) - 1) if y[i] > y[i - 1] and y[i] > y[i + 1]]
    rows = []
    for i in mins:
        nxt = [m for m in maxs if m > i]
        rows.append({"segmento": label, "R": R[i], "E": y[i],
                     "depth": (y[nxt[0]] - y[i]) if nxt else float("nan")})
    return rows


def main():
    args = parse_args()
    sysm = BOPSystem(n_manifold=args.n_manifold, n_s=args.n_s, delta0_ns=D0,
                     neighbors=tuple(args.neighbors))
    tag = args.tag or (f"n{sysm.n_manifold}"
                       + ("" if len(args.neighbors) == 3 else "_incompleta"))
    npz_path = args.npz or f"plots/archive/bop_curve_MJ0_{tag}.npz"
    out_path = args.out or f"plots/archive/bop_curve_MJ0_{tag}_excluded.png"
    L_NAME = {0: "s", 1: "p", 2: "d"}

    blk = sysm.block(0)
    vecinos = " + ".join(f"{sysm.levels[l]}{L_NAME[l]}"
                         for l in sorted(sysm.levels, reverse=True))
    print(RULE)
    print(f"CURVA BOP  M_J = 0   —   Rb*-KRb, manifold n={sysm.n_manifold} "
          f"(l>={sysm.l_min}) + {vecinos}")
    print(RULE)

    # ------------------------------------------------------ 1. montaje
    R_lo, R_max = sysm.domain_bounds()
    R_min = max(110.0, float(np.ceil(R_lo)))
    print(f"\n[1] MONTAJE Y DOMINIO")
    print(f"    dim(bloque M_J=0) = {len(blk)}   "
          f"(l = {sorted(sysm.levels)} + manifold {sysm.l_min}..{sysm.l_max})")
    print(f"    delta0_ns = {D0}")
    print(f"    n* del manifold = {sysm.fermi.n_star_manifold():.6f}")
    for l in sorted(sysm.levels, reverse=True):
        print(f"      l={l} -> {sysm.levels[l]}{L_NAME[l]}   "
              f"n* = {sysm.n_star_of_l(l):.6f}   "
              f"radial hidrogenoide n = {sysm.radial.n_of_l(l)}   "
              f"2n*² = {2*sysm.n_star_of_l(l)**2:.2f} a0")
    print(f"    dominio calculado sobre TODOS los pares (l1,l2): "
          f"R in [{R_lo:.4f}, {R_max:.4f}] a0")
    turns = {l: 2.0 * sysm.n_star_of_l(l) ** 2 for l in sysm.levels}
    l_bind = min(turns, key=turns.get)
    print(f"    par más restrictivo: {sysm.levels[l_bind]}{L_NAME[l_bind]}-"
          f"{sysm.levels[l_bind]}{L_NAME[l_bind]}, "
          f"retorno clásico 2n*² = {turns[l_bind]:.2f} a0")
    print(f"    R_max queda {turns[l_bind]-R_max:.2f} a0 por debajo: lo fija el borde")
    print(f"    de la tabla (R'_max = {sysm.scattering.R_table.max():.0f} a0 frente a "
          f"2·35² = 2450 a0).")
    print("    ⚠️ Ese tope NO es 'donde se acaba la física': es donde deja de estar")
    print("    definido el REMAPEO SEMICLÁSICO k(R) con el que leemos las tablas de")
    print("    dispersión. H_mol no tiene esa restricción.")
    print(f"    se traza desde R = {R_min:.1f} a0")
    print(f"    cero de energía: E(n={sysm.n_manifold}, l>=3) + KRb(N=0) = "
          f"{sysm.E_manifold:.12e} E_h")

    # ------------------------------------------------------ 2. barrido
    if args.reuse and os.path.exists(npz_path):
        d = np.load(npz_path)
        R, ctx, wman = d["R"], d["context"], d["w_manifold"]
        print(f"\n[2] BARRIDO — reutilizando {npz_path} ({len(R)} puntos)")
    else:
        R = np.arange(R_min, R_max, args.step)
        R = np.append(R, R_max)
        print(f"\n[2] BARRIDO ADIABÁTICO — {len(R)} puntos, paso {args.step} a0")
        t0 = time.perf_counter()
        ctx = np.empty((len(R), N_CONTEXT))
        wman = np.empty((len(R), N_CONTEXT))
        mask0 = sysm.is_manifold(0)
        for i, r in enumerate(R):
            w, V = sysm.solve(float(r))
            ctx[i] = w[:N_CONTEXT]
            wman[i] = (V[:, :N_CONTEXT] ** 2)[mask0].sum(axis=0)
            if i % 25 == 0:
                print(f"    R = {r:8.2f}  ({i+1}/{len(R)})", flush=True)
        dt = time.perf_counter() - t0
        print(f"    {len(R)} diagonalizaciones en {dt:.1f} s "
              f"({dt/len(R):.2f} s por punto)")
        np.savez(npz_path, R=R, context=ctx, w_manifold=wman)
        print(f"    datos en {npz_path}")

    Y = (ctx - sysm.E_manifold) * GHZ_PER_HARTREE
    k_man = sysm.lowest_manifold_index(float(R[-1]))

    # Dos formas de decir "la curva del manifold", y con la base completa NO
    # coinciden. El índice fijo sigue UNA curva adiabática (correcto por la
    # regla de no cruce), pero esa curva cambia de carácter al atravesar los
    # cruces evitados con los estados de (n+1)d y (n+2)p, que ahora caen justo
    # en el rango de energía relevante. La curva comparable con la figura del
    # paper es la de CARÁCTER: la más baja con peso de manifold > 50% en cada R.
    k_char = np.array([int(np.argmax(wman[i] > 0.5)) if (wman[i] > 0.5).any()
                       else -1 for i in range(len(R))])
    y_fix = Y[:, k_man]
    y = np.array([Y[i, k_char[i]] if k_char[i] >= 0 else np.nan
                  for i in range(len(R))])
    print(f"    índice fijo identificado en R_max: k = {k_man}")
    print(f"    índice por CARÁCTER: recorre k = {k_char.min()}..{k_char.max()}, "
          f"{int(np.sum(np.diff(k_char) != 0))} cambios")
    n_bad = int(np.sum(wman[:, k_man] <= 0.5))
    print(f"    la curva de índice fijo NO tiene carácter de manifold en "
          f"{n_bad} de {len(R)} puntos")
    if n_bad:
        Rb = R[wman[:, k_man] <= 0.5]
        print(f"      (R in [{Rb.min():.1f}, {Rb.max():.1f}] a0)  -> se analiza "
              "la curva de CARÁCTER, y se dibujan las dos")

    # -------------------------------------------- 3. resonancia de onda p
    win = sysm.p_resonance_window(margin_factor=args.margin_factor)
    HARTREE_MEV = 27211.386245988
    print("\n" + "-" * 84)
    print("[3] RESONANCIA DE FORMA p EN LA COORDENADA DE ESTE MANIFOLD")
    print("-" * 84)
    print("    Mismo criterio que en docs/analysis_ventana_exclusion_resonancia.md:")
    print("    centro = cero del interpolante de 1/A_p en eps=k^2/2 (el polo que ve")
    print("    el modo 'inverse'); anchura = FWHM de |A_p|; margen = f x FWHM.")
    print(f"\n    eps_polo = {win['eps_pole']:.8f} E_h = "
          f"{win['eps_pole']*HARTREE_MEV:.3f} meV   (misma energía que n=24: la")
    print("    resonancia es una propiedad de e-+Rb(5S), no del manifold)")
    print(f"    R_polo   = {win['R_pole']:.3f} a0")
    print(f"    FWHM     = {win['fwhm']:.3f} a0   "
          f"(R in [{win['R_fwhm_lo']:.3f}, {win['R_fwhm_hi']:.3f}])")
    print(f"    margen   = +/- {args.margin_factor:g} x FWHM = "
          f"+/- {args.margin_factor*win['fwhm']:.3f} a0")
    in_domain = R_min <= win["R_pole"] <= R_max
    print(f"\n    VENTANA EXCLUIDA: R in [{win['R_lo']:.2f}, {win['R_hi']:.2f}] a0"
          f"   ({'DENTRO' if in_domain else 'FUERA'} del dominio barrido)")
    keep = (R < win["R_lo"]) | (R > win["R_hi"])
    print(f"    {(~keep).sum()} de {len(R)} puntos caen dentro")
    print(f"    huecos de malla de la tabla: {sysm.fermi.large_gaps()}")
    gap_pts = [r for r in R if sysm.fermi.bridges_gap(float(r))]
    print(f"    puntos con A_p interpolado A TRAVÉS del polo: {len(gap_pts)}"
          + (f"  -> R = {', '.join(f'{r:.1f}' for r in gap_pts)}" if gap_pts else ""))

    # banda de cola: dónde vuelve la curva a caber en la ventana de la Fig. 2
    off = (R > win["R_hi"]) & (np.abs(y) > abs(args.ymin))
    R_tail = float(R[off].max()) if off.any() else float(win["R_hi"])
    tail = (R > win["R_hi"]) & (R <= R_tail)
    keep_strict = keep & ~tail
    if tail.any():
        print(f"\n    BANDA DE COLA butterfly: R in ({win['R_hi']:.2f}, "
              f"{R_tail:.2f}] a0, {tail.sum()} puntos")
        print(f"    (la curva de manifold no vuelve a |E| <= {abs(args.ymin):.0f} GHz "
              "hasta ahí)")
    else:
        print("\n    sin banda de cola: la curva vuelve a la ventana de energía "
              "inmediatamente")

    # ------------------------------------ 4. umbrales asintóticos y residuo
    print("\n" + "-" * 84)
    print(f"[4] UMBRALES ASINTÓTICOS {sysm.n_s}s + KRb(N)  Y CHEQUEO EN EL BORDE")
    print("-" * 84)
    dE = sysm.delta_E_ns_ghz()
    print(f"    ΔE({sysm.n_s}s) = E({sysm.n_s}s) - E(n={sysm.n_manifold}) = "
          f"{dE:.4f} GHz     <- umbral N=0")
    print(f"    B(KRb) = {B_KRB_GHZ} GHz")
    thr = sysm.thresholds_ns_ghz(range(0, 7))
    print("\n     N | B·N(N+1) [GHz] | umbral ΔE+B·N(N+1) [GHz]")
    print("    ---|----------------|-------------------------")
    for N in range(7):
        print(f"    {N:2d} | {B_KRB_GHZ*N*(N+1):14.4f} | {thr[N]:23.4f}")
    print(f"\n    pedidos explícitamente:  N=5 -> ΔE + 30B = {thr[5]:.4f} GHz"
          f"   ;  N=6 -> ΔE + 42B = {thr[6]:.4f} GHz")

    w_e, V_e = sysm.solve(float(R[-1]))
    mask = sysm.is_manifold(0)
    N_of_state = np.array([st[2] for st in blk.states], dtype=float)
    l_of_state = np.array([st[0] for st in blk.states])
    # Con la base completa los estados más bajos son los del (n+2)p, muy por
    # debajo del ns: hay que buscar los del ns, no coger los 12 primeros.
    ns_rows = [k for k in range(len(w_e))
               if float((V_e[:, k] ** 2)[mask].sum()) <= 0.5
               and int(l_of_state[np.argmax(V_e[:, k] ** 2)]) == 0]
    rows_show = sorted(set(ns_rows[:8]) | set(range(4)))
    print(f"\n    En el borde del dominio, R = {R[-1]:.2f} a0:")
    print(f"    (se listan los primeros estados de carácter {sysm.n_s}s, más los "
          "4 más bajos del bloque)\n")
    print("     k | E-E_man [GHz] | peso manifold | l dom. | <N> | asignación   | "
          "umbral [GHz] | residuo [GHz]")
    print("    ---|---------------|---------------|--------|-----|--------------|"
          "--------------|--------------")
    residuals = {}
    for k in rows_show:
        c2 = V_e[:, k] ** 2
        wm = float(c2[mask].sum())
        Nexp = float((c2 * N_of_state).sum())
        l_dom = int(l_of_state[np.argmax(c2)])
        Ek = (w_e[k] - sysm.E_manifold) * GHZ_PER_HARTREE
        if wm > 0.5:
            print(f"    {k:2d} | {Ek:13.4f} | {wm:13.4f} | {l_dom:6d} | {Nexp:3.1f} | "
                  f"{'MANIFOLD':12s} | {'(desplazada)':>12s} | {'—':>13s}")
            continue
        N = int(round(Nexp))
        lvl = f"{sysm.levels[l_dom]}{L_NAME[l_dom]}"
        t = sysm.thresholds_level_ghz(l_dom, [N])[N]
        res = Ek - t
        if l_dom == 0:
            residuals[N] = res
        print(f"    {k:2d} | {Ek:13.4f} | {wm:13.4f} | {l_dom:6d} | {Nexp:3.1f} | "
              f"{lvl+', N='+str(N):12s} | {t:12.4f} | {res:13.4f}")
    if len(residuals) >= 2:
        rv = np.array(list(residuals.values()))
        spread = rv.max() - rv.min()
        print(f"\n    residuos: media {rv.mean():+.4f} GHz, sigma {rv.std():.4f} GHz, "
              f"dispersión (max-min) {spread:.4f} GHz")
        rel = spread / abs(rv.mean()) if rv.mean() else float("inf")
        if rel < 0.25:
            print(f"    dispersión / |media| = {rel:.3f}  ->  residuo "
                  "APROXIMADAMENTE CONSTANTE, como se esperaba a R finito")
        else:
            print(f"    ⚠️ dispersión / |media| = {rel:.3f}  ->  el residuo NO es "
                  "constante entre umbrales N.")
            print("    Eso NO es el efecto trivial de R finito: es una "
                  "discrepancia que hay que explicar.")

    # ----------------------------- 5. características fuera de la ventana
    print("\n" + "-" * 84)
    print("[5] CARACTERÍSTICAS DE LA CURVA DE CARÁCTER MANIFOLD, "
          "EXCLUYENDO LA VENTANA")
    print("-" * 84)
    left, right = R < win["R_lo"], R > win["R_hi"]
    rows = (segment_features(R[left], y[left], f"izq  R<{win['R_lo']:.1f}")
            + segment_features(R[right], y[right], f"der  R>{win['R_hi']:.1f}"))
    print("    † = el mínimo cae en la banda de cola butterfly, no comparable\n")
    print("     #  | segmento              |  R_min [a0] |  E-E_man [GHz] | "
          "prof. vs máx. der. [GHz] | ΔR al anterior")
    print("    ----|-----------------------|-------------|----------------|"
          "--------------------------|---------------")
    prev = prev_seg = None
    seps = []
    for n, r in enumerate(rows):
        same = prev is not None and r["segmento"] == prev_seg
        dR = (r["R"] - prev) if prev is not None else float("nan")
        if same:
            seps.append(dR)
        flag = "†" if win["R_hi"] < r["R"] <= R_tail else " "
        prev, prev_seg = r["R"], r["segmento"]
        dep = f"{r['depth']:24.4f}" if np.isfinite(r["depth"]) else " " * 20 + "---"
        sep_txt = (f"{dR:13.2f}{'*' if not same else ' '}"
                   if np.isfinite(dR) else "           ---")
        print(f"    {n:3d}{flag}| {r['segmento']:21s} | {r['R']:11.2f} | "
              f"{r['E']:14.4f} | {dep} | " + sep_txt)
    print("    (* = salto que atraviesa la ventana excluida; no es separación física)")
    if len(seps) >= 2:
        s = np.array(seps)
        cv = s.std() / s.mean()
        print(f"\n    separación entre mínimos CONTIGUOS: media {s.mean():.2f} a0, "
              f"sigma {s.std():.2f} a0, min {s.min():.2f}, max {s.max():.2f}")
        print(f"    coeficiente de variación = {cv:.3f}  "
              f"({'aprox. regular' if cv < 0.25 else 'NO regular'})")

    print(f"\n    rango fuera de la ventana: [{y[keep].min():.3f}, "
          f"{y[keep].max():.3f}] GHz")
    print(f"    mínimo fuera de ventana Y cola: {y[keep_strict].min():.4f} GHz "
          f"en R = {R[keep_strict][np.argmin(y[keep_strict])]:.2f} a0")
    if (~keep).any():
        print(f"    (dentro de la ventana la curva llega a {y[~keep].min():.1f} GHz, "
              "valor SIN significado físico)")

    print("\n    Texto del paper: 'shifted more than 20 GHz for R <~ 1200 a0'")
    print(f"    (nuestro dominio llega a R = {R[-1]:.1f} a0)")
    f_all = float(np.mean(np.abs(y) > 20.0))
    f_out = float(np.mean(np.abs(y[keep]) > 20.0))
    f_str = float(np.mean(np.abs(y[keep_strict]) > 20.0))
    print(f"      |E-E_man| > 20 GHz en {f_all*100:6.2f}% sin excluir nada")
    print(f"                            {f_out*100:6.2f}% excluida la ventana")
    print(f"                            {f_str*100:6.2f}% excluidas ventana + cola")
    below = R[keep_strict][np.abs(y[keep_strict]) <= 20.0]
    if len(below):
        print(f"      los {len(below)} puntos por debajo de 20 GHz están en "
              f"R in [{below.min():.1f}, {below.max():.1f}] a0")

    # ---- pozo "uno antes del más externo" (Fig. 4 del paper) ----------
    print("\n" + "-" * 84)
    print("[6] POZO 'UNO ANTES DEL MÁS EXTERNO'  (punto de comparación de la Fig. 4)")
    print("-" * 84)
    wells = [r for r in rows if not (win["R_hi"] < r["R"] <= R_tail)]
    R_class = 2.0 * sysm.n_manifold**2
    print(f"    punto de retorno clásico del manifold: 2n² = {R_class:.0f} a0")
    print(f"    dominio accesible del remapeo: hasta {R[-1]:.1f} a0 "
          f"({100*R[-1]/R_class:.1f} % de 2n²)")
    if len(wells) >= 2:
        outer, penult = wells[-1], wells[-2]
        print(f"    pozo MÁS EXTERNO:      R = {outer['R']:8.2f} a0   "
              f"E-E_man = {outer['E']:10.4f} GHz   R/2n² = {outer['R']/R_class:.3f}")
        print(f"    pozo UNO ANTES:        R = {penult['R']:8.2f} a0   "
              f"E-E_man = {penult['E']:10.4f} GHz   R/2n² = {penult['R']/R_class:.3f}")
        print(f"    separación entre los dos: {outer['R']-penult['R']:.2f} a0")
    else:
        print(f"    sólo hay {len(wells)} mínimos utilizables fuera de la ventana: "
              "no se puede identificar un 'uno antes del más externo'")
    print("    ⚠️ El pozo más externo de NUESTRA curva está limitado por el borde")
    print("    del dominio del remapeo, no por la física: si el paper mide el pozo")
    print("    más externo cerca de 2n², nuestro 'más externo' puede no ser el suyo.")

    # ------------------------------------------------------------ 7. figura
    ymask = np.where(keep[:, None], Y, np.nan)
    fig, axes = plt.subplots(2, 1, figsize=(9.5, 9.5), sharex=True)
    band_label = ("región excluida: divergencia conocida del pseudopotencial de\n"
                  "rango cero en la resonancia de forma p — ver\n" + DOC)
    for panel, (ax, ylim, ttl, masked) in enumerate((
        (axes[0], (float(min(np.nanmin(y), np.nanmin(y_fix))) * 1.05, 100.0),
         "Rango completo — se MUESTRA la divergencia que luego se excluye", False),
        (axes[1], (args.ymin, args.ymax),
         "Ventana de energía comparable a la Fig. 2 — región excluida en gris", True),
    )):
        ax.axvspan(win["R_lo"], win["R_hi"], color="0.45", alpha=0.32, zorder=0,
                   label=band_label if panel == 0 else None)
        if tail.any():
            ax.axvspan(win["R_hi"], R_tail, color="0.45", alpha=0.13, hatch="///",
                       ec="0.6", lw=0.0, zorder=0,
                       label=(f"cola butterfly, no comparable "
                              f"($R\\leq{R_tail:.0f}\\,a_0$)" if panel == 0 else None))
        ax.axvline(win["R_pole"], color="0.35", lw=1.0, ls=":", zorder=1,
                   label=(f"polo de $A_p$, $R={win['R_pole']:.1f}\\,a_0$"
                          if panel == 0 else None))
        src = ymask if masked else Y
        for c in range(Y.shape[1]):
            ax.plot(R, src[:, c], color="0.78", lw=0.6, zorder=1)
        y_show = np.where(keep, y, np.nan) if masked else y
        y_fix_show = np.where(keep, y_fix, np.nan) if masked else y_fix
        ax.plot(R, y_show, color="crimson", lw=3.0, zorder=4,
                label=f"más baja de carácter manifold $n={sysm.n_manifold}$ "
                      "(>50 % de peso)")
        if n_bad:
            ax.plot(R, y_fix_show, color="darkblue", lw=1.3, ls="--", zorder=3,
                    label=f"índice fijo k={k_man} (cambia de carácter)")
        ax.axhline(0, color="k", lw=0.8, ls=":", zorder=2)
        ax.axhline(-20, color="seagreen", lw=1.2, ls="--", zorder=2,
                   label="$-20$ GHz (texto del paper)")
        ax.set_ylim(*ylim)
        ax.set_xlim(R[0], max(R[-1], 2.0 * sysm.n_manifold**2) * 1.01)
        ax.grid(alpha=0.25)
        ax.set_ylabel(rf"$E - E_{{n={sysm.n_manifold}}}$  [GHz]")
        ax.set_title(ttl, fontsize=10)
        ax.legend(loc="lower right", fontsize=7.5)
    for N in range(7):
        axes[1].axhline(thr[N], color="tab:orange", lw=0.7, alpha=0.6)
        axes[1].axhline(B_KRB_GHZ * N * (N + 1), color="tab:blue", lw=0.7, alpha=0.6)
    for ax in axes:
        ax.axvline(2.0 * sysm.n_manifold**2, color="tab:purple", lw=1.2, ls="-.",
                   alpha=0.9, zorder=2)
        ax.axvspan(R[-1], 2.0 * sysm.n_manifold**2, color="tab:purple", alpha=0.07,
                   zorder=0)
    axes[0].plot([], [], color="tab:purple", lw=1.2, ls="-.",
                 label=f"$2n^2 = {2*sysm.n_manifold**2:.0f}\\,a_0$ (retorno "
                       "clásico; violeta: fuera del dominio del remapeo)")
    axes[0].legend(loc="lower right", fontsize=7.5)
    axes[1].annotate(
        "REGIÓN EXCLUIDA\ndivergencia conocida del pseudopotencial\n"
        "de rango cero en la resonancia de forma p\n"
        + ("(rayado: cola butterfly, tampoco comparable)\n" if tail.any() else "")
        + f"ver {DOC}",
        xy=(win["R_pole"], 0.80 * args.ymin), xycoords="data",
        xytext=(R[0] + 25, 0.78 * args.ymin), textcoords="data",
        fontsize=7.5, ha="left", va="center", color="0.10",
        bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="0.5", alpha=0.92),
        arrowprops=dict(arrowstyle="->", color="0.25", lw=1.1))
    axes[1].set_xlabel(r"$R$  [$a_0$]")
    fig.suptitle(
        f"Rb*-KRb, curvas BOP $M_J=0$, manifold $n={sysm.n_manifold}$ + "
        f"{vecinos}   "
        r"$H=H_a+B\mathbf{N}^2-\mathbf{d}\cdot\mathbf{F}_{ryd}+V_{Fermi}$"
        f"\ninterpolación de $1/A_p$ · ventana excluida "
        f"$R\\in[{win['R_lo']:.1f},\\,{win['R_hi']:.1f}]\\,a_0$ "
        f"= $R_{{polo}}\\pm{args.margin_factor:g}\\times$FWHM$({win['fwhm']:.1f}\\,a_0)$"
        f"\nnaranja: umbrales ${sysm.n_s}s$+KRb(N) · azul: manifold+KRb(N)",
        fontsize=10.5)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    fig.savefig(out_path, dpi=150)
    print(f"\n    PNG guardado en {out_path}")
    print(RULE)


if __name__ == "__main__":
    main()
