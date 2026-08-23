"""
Factores de Franck-Condon, Fase 5 de `docs/PLAN_nonadiabatic_dynamics.md`.

Réplica de la Ec. (17) de Mellado-Alcedo, Guttridge, Cornish, Sadeghpour &
González-Férez, PRA 110, 013314 (2024) — en el plan maestro se la llamaba
"Ec. A9" por el nombre provisional del apéndice; el número real, verificado
contra el texto (arxiv.org/html/2401.09618), es Ec. (17):

    F_{i,N} = Σ_{k=d,u} ∫ (χᵢᵏ(R))* Ω_ns C^{J,k}_{n,l=0,N}(R) ψ_scat(R) R dR

con χᵢᵏ(R) la función de onda vibracional del estado final (canal k=d,u),
Ω_ns la frecuencia de Rabi de la excitación Rydberg de dos fotones (=1 en
el paper), C^{J,k}_{n,l=0,N}(R) el coeficiente electrónico de peso de la
onda parcial l=0 (Ec. 9 del paper — física electrónica del sistema
concreto, fuera del alcance genérico de este módulo, se trata aquí como un
peso `coupling_weight(R)` de entrada), y ψ_scat(R) la "función de onda de
dispersión del canal abierto inicial" (texto del paper, sin más detalle).

CONVENCIÓN: reducida en todo el módulo
-----------------------------------------
El factor extra "R dR" de la Ec. (17) (en vez del "dR" de las Fases 1-4)
indica que ψ_scat(R) está escrita en la forma NO reducida (3D, normalizada
vía ∫|ψ_scat|²R²dR=1), mientras que χᵢᵏ(R) es la forma REDUCIDA habitual
(u=R·R_l, ∫|χ|²dR=1) que usa el resto del proyecto. Con u_scat(R)≡R·ψ_scat(R)
(la forma reducida de ψ_scat), la Ec. (17) por canal se reescribe

    F = Ω_ns ∫ χᵢᵏ(R)* · C(R) · u_scat(R) dR                         (*)

que es la que implementa `franck_condon_factor` — exactamente la misma
convención reducida (∫|χ|²dR=1) que `coupling.py`, `coupled_channels.py`,
`energy_normalization.py` y `decay_rates.py`. La suma sobre canales k=d,u
de la Ec. (17) no se hace dentro de esta función (cada canal es una llamada
independiente); sumarlas es responsabilidad de quien la use.

⚠️ QUÉ ESTADO INICIAL SE USA PARA Rb*-RbCs (decisión explícita de esta
ronda, no una copia del paper)
--------------------------------------------------------------------------
El paper de 2024 define ψ_scat para SU protocolo (Cs-RbCs): dos escenarios,
(a) pinzas ópticas FUSIONADAS (Cs y RbCs comparten una trampa combinada,
Cs(42s), separación <200 nm) — el estado inicial natural ahí es el estado
FUNDAMENTAL del movimiento relativo en la trampa fusionada, un oscilador
armónico centrado en R≈0; o (b) pinzas SEPARADAS (Cs(74s), ~500 nm) — el
paper menciona este escenario pero no da su ψ_scat explícitamente tampoco.

Para Rb*-RbCs el sistema real YA DEMOSTRADO es un TERCER caso, distinto de
ambos: Ruttley, Guttridge, Baldock, González-Férez, Sadeghpour, Adams &
Cornish, "Observation of Rydberg blockade due to the charge-dipole
interaction between an atom and a polar molecule", PRL 131, 013401 (2023),
arXiv:2303.06126 — Rb y RbCs en pinzas ÓPTICAS ESPECÍFICAS POR ESPECIE,
INDEPENDIENTES (NO fusionadas), con la separación R controlada
explícitamente moviendo las pinzas ("Species-specific tweezers are used to
control the separation between the atom and molecule"). La separación
átomo-molécula demostrada para el bloqueo Rydberg (Rb→52s) fue
R_am=310(40) nm, y el propio paper reporta la dispersión experimental de
la alineación relativa entre pinzas tiro a tiro: "an estimated standard
deviation of 50 nm in each coordinate".

Esto es físicamente DISTINTO del escenario de pinzas fusionadas: no hay una
trampa combinada con un oscilador armónico de movimiento relativo centrado
en R=0. En vez de eso, el átomo y la molécula ocupan trampas INDEPENDIENTES
separadas por una distancia R0 FIJADA experimentalmente (el parámetro de
control del experimento), con una incertidumbre de posicionamiento relativo
σ medida directamente (50 nm). El modelo natural, adaptado a este sistema
—NO al del paper— es entonces un PAQUETE DE ONDAS GAUSSIANO centrado en la
separación fijada R0, de anchura σ dada por esa incertidumbre de
posicionamiento medida, no un oscilador armónico centrado en R=0:

    u_scat(R) = (πσ²)^{-1/4} exp[-(R-R0)²/(2σ²)]                     (**)

`gaussian_wavepacket` implementa (**), normalizada NUMÉRICAMENTE sobre la
malla real (no sólo la fórmula analítica, válida para R0≫σ — se verifica
en los tests que la aproximación es buena para los R0,σ de interés). R0 y σ
son parámetros de entrada de este módulo, NO fijados a Rb*-RbCs aquí: fijar
R0 al cruce evitado real y σ al valor de Ruttley et al. (2023) es tarea de
la Fase 6, que además exigiría antes decidir a qué n de Rydberg (y por
tanto qué escala de R0) corresponde el análisis — no se hace aquí.
"""
from typing import Optional

import numpy as np

__all__ = ["gaussian_wavepacket", "franck_condon_factor"]


def gaussian_wavepacket(R_grid: np.ndarray, R0: float, sigma: float) -> np.ndarray:
    """
    u_scat(R), Ec. (**) del docstring del módulo, normalizada NUMÉRICAMENTE
    en la malla dada (Σu²h=1 exacto en la malla, no sólo en el límite
    analítico R0≫σ).

    Args:
        R_grid: malla uniforme de R (se verifica), R≥0.
        R0: centro del paquete (separación fijada, p.ej. de pinzas).
        sigma: anchura (incertidumbre de posicionamiento).

    Returns:
        u_scat en los puntos de R_grid, forma reducida, ∫u²dR=1.
    """
    R_grid = np.asarray(R_grid, dtype=float)
    steps = np.diff(R_grid)
    if not np.allclose(steps, steps[0], rtol=1e-9, atol=1e-12):
        raise ValueError(
            "gaussian_wavepacket requiere una malla de R uniforme "
            f"(pasos observados: {steps})")
    h = float(steps[0])

    u = np.exp(-((R_grid - R0) ** 2) / (2.0 * sigma ** 2))
    norm = np.sqrt(np.sum(u ** 2) * h)
    return u / norm


def franck_condon_factor(
    R_grid: np.ndarray,
    v_final: np.ndarray,
    u_scat: np.ndarray,
    coupling_weight: Optional[np.ndarray] = None,
    omega_ns: float = 1.0,
) -> float:
    """
    F, Ec. (*) del docstring del módulo, para UN canal.

    Args:
        R_grid: malla uniforme de R (se verifica).
        v_final: autovector PURO (Σv²=1, convención de `eigh`) del estado
            final χᵢᵏ, forma reducida.
        u_scat: función de onda inicial YA física y normalizada en la
            misma malla (∫u²dR=1), p.ej. de `gaussian_wavepacket`.
        coupling_weight: C(R) de la Ec. (17), o None (=1: solapamiento de
            Franck-Condon puro, sin peso electrónico).
        omega_ns: Ω_ns, por defecto 1.0 (como en el paper).

    Returns:
        F (real; complejo si `v_final`/`u_scat`/`coupling_weight` lo son).
    """
    R_grid = np.asarray(R_grid, dtype=float)
    steps = np.diff(R_grid)
    if not np.allclose(steps, steps[0], rtol=1e-9, atol=1e-12):
        raise ValueError(
            "franck_condon_factor requiere una malla de R uniforme "
            f"(pasos observados: {steps})")
    h = float(steps[0])

    chi_final = np.asarray(v_final) / np.sqrt(h)
    weight = np.ones_like(R_grid) if coupling_weight is None else np.asarray(coupling_weight)

    return omega_ns * np.sum(np.conj(chi_final) * weight * np.asarray(u_scat)) * h
