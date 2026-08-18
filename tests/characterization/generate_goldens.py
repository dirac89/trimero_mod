"""
Genera los golden files G1 y G2 del código legado.

IMPORTANTE: este script debe ejecutarse UNA sola vez, contra el código tal
como estaba ANTES de mover nada (paso 0 del refactor). Los ficheros que
produce son la referencia contra la que se verifica cada paso posterior.
Volver a ejecutarlo después de refactorizar destruiría la red de seguridad:
regeneraría la referencia a partir del código que se quiere verificar.

Uso:
    python tests/characterization/generate_goldens.py
"""
import os
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src"))
GOLDENS = Path(__file__).resolve().parent / "goldens"

from trimero.mathlib.special import Spherical, DRnl, DOlm, DPhilm, hydrogenicR  # noqa: E402
from trimero.hamiltonians.fermi import FermiPotentials  # noqa: E402


# --- G1: funciones especiales puras -----------------------------------
# (n, l) cubre los manifolds que el proyecto usa realmente: n1=5 de las
# pruebas cortas, n=24/27 de la ronda carga-dipolo, y 35-38 del trimero.
G1_NL = [(5, 0), (5, 1), (5, 2), (5, 3),
         (24, 0), (24, 3), (24, 10), (24, 23),
         (27, 0), (27, 1),
         (35, 0), (35, 2), (35, 4), (35, 34),
         (36, 2), (37, 1), (38, 0)]
G1_R = [1.0, 10.0, 100.0, 500.0, 1000.0, 2000.0]
G1_LM = [(0, 0), (1, 0), (1, 1), (2, 0), (2, 1), (2, 2),
         (3, 0), (3, 2), (4, 1), (4, 4)]
G1_THETA = [0.0, np.pi / 4, np.pi / 2, 3 * np.pi / 4, np.pi]


def generate_g1():
    radial_h, radial_d = [], []
    for (n, l) in G1_NL:
        for r in G1_R:
            radial_h.append(hydrogenicR(n, l, 1.0, r))
            radial_d.append(DRnl(n, l, r))

    ang_sph, ang_dolm, ang_dphi = [], [], []
    for (l, m) in G1_LM:
        for th in G1_THETA:
            ang_sph.append(Spherical(l, m, np.cos(th)))
            ang_dolm.append(DOlm(l, m, th))
            ang_dphi.append(DPhilm(l, m, th))

    np.savez(
        GOLDENS / "g1_special.npz",
        nl=np.array(G1_NL), r=np.array(G1_R),
        lm=np.array(G1_LM), theta=np.array(G1_THETA),
        hydrogenicR=np.array(radial_h, dtype=np.float64),
        DRnl=np.array(radial_d, dtype=np.float64),
        Spherical=np.array(ang_sph, dtype=np.float64),
        DOlm=np.array(ang_dolm, dtype=np.float64),
        DPhilm=np.array(ang_dphi, dtype=np.float64),
    )
    return len(radial_h) * 2 + len(ang_sph) * 3


# --- G2: elementos del potencial de Fermi ------------------------------
# Los cuatro pares (li, lj) cubren las cuatro ramas de Vs y de Vp:
#   (4,5) li>2  lj>2   -> hidrogénica en ambos lados
#   (1,4) li<=2 lj>2   -> tabulada a la izquierda
#   (4,1) li>2  lj<=2  -> tabulada a la derecha
#   (1,2) li<=2 lj<=2  -> tabulada en ambos
G2_LILJ = [(4, 5), (1, 4), (4, 1), (1, 2)]
# (mi, mj) incluye (0,0) a propósito: con m != 0 los armónicos se anulan en
# theta = 0 y theta = pi, que son EXACTAMENTE los dos ángulos que usa
# trimer.py (theta=0, theta1=pi). Sin el caso m=0, Vs valdría cero en la
# mayoría de los casos y el golden vigilaría ceros triviales.
G2_MIMJ = [(0, 0), (1, 1), (0, 1)]
G2_CASES = []
for (li, lj) in G2_LILJ:
    for (mi, mj) in G2_MIMJ:
        for s in (0, 1):
            for (r1, theta1) in [(1000.0, 0.0), (1500.0, np.pi), (800.0, np.pi / 3)]:
                G2_CASES.append(dict(
                    s=s, n=35, n2=35, li=li, lj=lj, mi=mi, mj=mj,
                    r1=r1, theta1=theta1,
                    As1=-16.1, Ap1=-24.0,        # long. de dispersión típicas Rb
                    wave1=3.1e-5, wave2=-2.7e-5, # valores tabulados no nulos
                    Dwave1=1.4e-7, Dwave2=-9.2e-8,
                ))


def generate_g2():
    vs, vp, vsp = [], [], []
    for c in G2_CASES:
        f = FermiPotentials(
            c["s"], c["n"], c["n2"], c["li"], c["lj"], c["mi"], c["mj"],
            c["r1"], c["theta1"], c["As1"], c["wave1"], c["wave2"],
            c["Ap1"], c["Dwave1"], c["Dwave2"],
        )
        vs.append(f.Vs())
        vp.append(f.Vp())
        vsp.append(f.Vsp())

    keys = sorted(G2_CASES[0].keys())
    params = np.array([[c[k] for k in keys] for c in G2_CASES], dtype=np.float64)
    np.savez(
        GOLDENS / "g2_fermi.npz",
        param_names=np.array(keys),
        params=params,
        Vs=np.array(vs, dtype=np.float64),
        Vp=np.array(vp, dtype=np.float64),
        Vsp=np.array(vsp, dtype=np.float64),
    )
    return len(vs) * 3


def _refuse_accidental_regeneration():
    """Los goldens son la referencia del refactor: regenerarlos desde el
    código ya refactorizado destruiría la red de seguridad sin avisar.

    Sobrescribir exige intención explícita:
        TRIMERO_REGENERATE_GOLDENS=1 python <este script>
    """
    import os
    existing = sorted(GOLDENS.glob("*.npz"))
    if existing and os.environ.get("TRIMERO_REGENERATE_GOLDENS") != "1":
        names = ", ".join(f.name for f in existing)
        raise SystemExit(
            f"Ya existen goldens ({names}).\n"
            "Regenerarlos invalida la referencia contra la que se verifica el "
            "refactor.\nSi de verdad quieres sobrescribirlos:\n"
            "    TRIMERO_REGENERATE_GOLDENS=1 python " + __file__
        )


if __name__ == "__main__":
    GOLDENS.mkdir(parents=True, exist_ok=True)
    _refuse_accidental_regeneration()
    n1 = generate_g1()
    print(f"G1 escrito: {n1} valores -> goldens/g1_special.npz")
    n2 = generate_g2()
    print(f"G2 escrito: {n2} valores -> goldens/g2_fermi.npz")
