"""
PUNTO 1 de la ronda híbrida: transformación θ=0→θ=π de UN perturbador de
V_Fermi, verificado por fuerza bruta contra ψ y ∇ψ numéricos en (0,0,-R).

DERIVACIÓN FORMAL (operador paridad 𝒫)
--------------------------------------
El perturbador individual en θ=π está en R⃗_π = -Rẑ. Con el operador paridad
𝒫|r̄⟩ = |-r̄⟩, 𝒫² = 1, que actúa sobre el gradiente como 𝒫∇𝒫† = -∇:

    𝒫 δ³(r⃗-Rẑ) 𝒫†     = δ³(r⃗+Rẑ)
    𝒫 ∇⃖δ³(r⃗-Rẑ)∇⃗ 𝒫† = (-∇⃖) δ³(r⃗+Rẑ) (-∇⃗) = ∇⃖δ³(r⃗+Rẑ)∇⃗

Los DOS signos del término bilineal p-wave se cancelan, luego

    V(R⃗_π) = 𝒫 V(R⃗_0) 𝒫† .

Sobre estados hidrogenoides reales ⟨l m|𝒫 = (-1)^l ⟨l m|, así que

    ⟨l₁m₁|V(θ=π)|l₂m₂⟩ = (-1)^{l₁+l₂} · ⟨l₁m₁|V(θ=0)|l₂m₂⟩        (*)

con las MISMAS reglas de selección que en θ=0 (Δm=0, |m| ≤ 1; V_s sólo m=0).
No hay ningún signo extra asociado al GRADIENTE más allá del producto
(-1)^{l₁+l₂}: los dos −1 de 𝒫∇𝒫† ya se cancelaron entre sí.

CONSISTENCIA CON linear_trimer: sumando el elemento de θ=0 y el de θ=π,

    V_total_simétrico = [1 + (-1)^{l₁+l₂}] · V(θ=0),

que es EXACTAMENTE el `parity_factor` de `SymmetricLinearTrimer`. La fórmula
(*) es pues el ingrediente que aquél asumía sin verificar de forma aislada;
aquí se usa SOLA, porque en el sistema híbrido el Rb neutro está en θ=π SIN
compañero en θ=0.

REFERENCIA INDEPENDIENTE
------------------------
`brute_force_fermi_element_at` evalúa ψ y ∇ψ numéricamente en el punto
cartesiano (0,0,z0) —z0=-R para θ=π— por diferencias finitas centradas de la
función de onda completa (hidrogenoide × Y_lm desde `lpmv`). No usa ninguna
forma cerrada del módulo de producción, igual que el test F2 de
`test_fermi_krb.py` hacía para θ=0.

Ejecutar:  poetry run pytest tests/systems/hybrid_neutral_polar/test_transformacion_theta_pi.py -s
"""

import numpy as np
import pytest
from scipy.special import eval_genlaguerre, gammaln, lpmv

from trimero.basis.radial import RadialBasis
from trimero.systems.rb_neutral_perturber.fermi_krb import (
    FermiPseudopotential,
    ScatteringLengths,
)

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")

N_MANIFOLD = 35


# ---------------------------------------------------------------- referencia
def _Ylm(l, m, theta, phi):
    """Y_lm complejo con fase de Condon-Shortley, construido desde lpmv."""
    am = abs(m)
    c = np.sqrt(
        (2 * l + 1) / (4 * np.pi) * np.exp(gammaln(l - am + 1) - gammaln(l + am + 1))
    )
    val = c * lpmv(am, l, np.cos(theta)) * np.exp(1j * am * phi)
    if m < 0:
        val = (-1) ** am * np.conj(val)
    return val


def hydrogenic_R_ref(n, l, r):
    """Hidrogenoide RECONSTRUIDA aquí, independiente del módulo de producción."""
    rho = 2.0 * r / n
    log_norm = 0.5 * (
        3.0 * np.log(2.0 / n) + gammaln(n - l) - np.log(2.0 * n) - gammaln(n + l + 1)
    )
    return float(
        np.exp(log_norm + l * np.log(rho) - rho / 2.0) * eval_genlaguerre(n - l - 1, 2 * l + 1, rho)
    )


def _psi_at(n, l, m, x, y, z):
    r = np.sqrt(x * x + y * y + z * z)
    return hydrogenic_R_ref(n, l, r) * _Ylm(l, m, np.arccos(z / r), np.arctan2(y, x))


def _grad_psi_at(n, l, m, point, h):
    """∇ψ en un punto cartesiano arbitrario por diferencias finitas centradas."""
    g = []
    for axis in range(3):
        p = list(point)
        q = list(point)
        p[axis] += h
        q[axis] -= h
        g.append((_psi_at(n, l, m, *p) - _psi_at(n, l, m, *q)) / (2.0 * h))
    return np.array(g)


def brute_force_fermi_element_at(n_i, l_i, m_i, n_j, l_j, m_j, point, A_s, A_p,
                                 h=0.05):
    """
    ⟨i| 2π A_s δ³(r-point) + 6π A_p δ³(r-point) ∇·∇ |j⟩ evaluado numéricamente
    en `point` CUALQUIERA. Sin formas cerradas de producción.
    """
    psi_i = _psi_at(n_i, l_i, m_i, *point)
    psi_j = _psi_at(n_j, l_j, m_j, *point)
    vs = 2.0 * np.pi * A_s * np.conj(psi_i) * psi_j
    gi = _grad_psi_at(n_i, l_i, m_i, point, h)
    gj = _grad_psi_at(n_j, l_j, m_j, point, h)
    vp = 6.0 * np.pi * A_p * np.dot(np.conj(gi), gj)
    return complex(vs + vp)


# ---------------------------------------------------------------- fixtures
RADIAL = RadialBasis(n_manifold=N_MANIFOLD)
SCAT = ScatteringLengths(n_star_table=float(N_MANIFOLD))
FERMI = FermiPseudopotential(
    RADIAL, SCAT, n_manifold=N_MANIFOLD, uniform_n_star=float(N_MANIFOLD)
)


def phase_pi(l1: int, l2: int) -> float:
    """(-1)^{l1+l2}, el peso geométrico del perturbador INDIVIDUAL en θ=π."""
    return (-1.0) ** (l1 + l2)


# ================================================================ TEST P1a
def test_p1a_theta_cero_sanidad():
    """Primero reproducir F2 en n=35: las cerradas funcionan en θ=0."""
    R = 900.0
    A_s, A_p = SCAT.scattering(R, float(N_MANIFOLD))
    casos = [(4, 0, 7, 0), (3, 0, 3, 0), (0, 0, 5, 0), (5, 1, 8, 1),
             (3, -1, 6, -1), (2, 1, 5, 1)]
    print(f"\n  SANIDAD θ=0 (punto (0,0,+{R:g})), A_s={A_s:.6f}, A_p={A_p:.4f}")
    print("   l1 m1 l2 m2 |   cerrada        |  fuerza bruta    |   rel")
    print("   ------------|------------------|------------------|--------")
    vals = []
    for (l1, m1, l2, m2) in casos:
        n1, n2 = RADIAL.n_of_l(l1), RADIAL.n_of_l(l2)
        closed = FERMI.electron_element(l1, m1, l2, m2, R)
        bf = brute_force_fermi_element_at(
            n1, l1, m1, n2, l2, m2, (0.0, 0.0, R), A_s, A_p)
        vals.append((l1, m1, l2, m2, closed, bf.real))
    escala = max(abs(c) for *_, c, _ in vals)
    floor = 1e-12 * escala
    worst = max(abs(c - b) / max(abs(c), abs(b), floor) for *_, c, b in vals)
    for (l1, m1, l2, m2, closed, bf) in vals:
        print(f"   {l1:2d}{m1:3d}{l2:3d}{m2:3d} | {closed:+.10e} | "
              f"{bf:+.10e} | {abs(closed-bf)/max(abs(closed),abs(bf),floor):.1e}")
    print(f"   peor diferencia relativa = {worst:.2e}")
    assert worst < 1e-5, worst


# ================================================================ TEST P1b
def test_p1b_transformacion_pi_fuerza_bruta():
    """
    El corazón del punto 1: cerrada(θ=π) := (-1)^{l1+l2}·cerrada(θ=0)
    contra fuerza bruta en (0,0,-R).
    """
    R = 900.0
    A_s, A_p = SCAT.scattering(R, float(N_MANIFOLD))
    casos = [
        (4, 0, 7, 0),    # manifold-manifold, m=0 (V_s y V_p radial)
        (3, 0, 3, 0),    # diagonal m=0
        (0, 0, 5, 0),    # vecino s-manifold, m=0
        (9, 0, 10, 0),   # par (l1+l2 impar) -> fase -1
        (5, 1, 8, 1),    # Π, transverso
        (3, -1, 6, -1),  # Π negativo
        (2, 1, 5, 1),    # vecino d-manifold, Π
        (5, 2, 7, 2),    # |m|=2: cero exacto en AMBOS polos
        (4, 0, 7, 1),    # Δm≠0: cero exacto en ambos polos
    ]
    crudos = []
    for (l1, m1, l2, m2) in casos:
        n1, n2 = RADIAL.n_of_l(l1), RADIAL.n_of_l(l2)
        f = phase_pi(l1, l2)
        closed_pi = f * FERMI.electron_element(l1, m1, l2, m2, R)
        bf_pi = brute_force_fermi_element_at(
            n1, l1, m1, n2, l2, m2, (0.0, 0.0, -R), A_s, A_p)
        crudos.append((l1, m1, l2, m2, f, closed_pi, bf_pi))

    escala = max(abs(c) for *_, c, _ in crudos)
    floor = 1e-12 * escala

    print(f"\n  TRANSFORMACIÓN θ=0→π, punto sur (0,0,-{R:g})")
    print("   l1 m1 l2 m2 | fase | cerrada(π)       | fuerza bruta(π)  |   rel")
    print("   ------------|------|------------------|------------------|--------")
    worst = 0.0
    for (l1, m1, l2, m2, f, closed_pi, bf_pi) in crudos:
        # Parte imaginaria: sólo ruido numérico de las diferencias finitas.
        assert abs(bf_pi.imag) <= 1e-8 * escala + 1e-30, (l1, m1, l2, m2, bf_pi)
        rel = abs(closed_pi - bf_pi.real) / max(abs(closed_pi), abs(bf_pi.real),
                                                floor)
        worst = max(worst, rel)
        nota = "  <- cero exacto" if abs(closed_pi) < floor else ""
        print(f"   {l1:2d}{m1:3d}{l2:3d}{m2:3d} | {f:+.0f}  | {closed_pi:+.10e} | "
              f"{bf_pi.real:+.10e} | {rel:.1e}{nota}")
    print(f"   escala del mayor elemento = {escala:.2e}")
    print(f"   peor diferencia relativa = {worst:.2e}")
    assert escala > 0.0, "comparación vacua"
    assert worst < 1e-5, worst


# ================================================================ TEST P1c
def test_p1c_relacion_entre_polos():
    """
    Verificación DIRECTA de (*) sin pasar por las cerradas:
    bruta(θ=π) / bruta(θ=0) debe dar (-1)^{l1+l2}.
    """
    R = 700.0
    A_s, A_p = SCAT.scattering(R, float(N_MANIFOLD))
    casos = [(4, 0, 9, 0), (3, 0, 3, 0), (5, 1, 6, 1), (3, -1, 4, -1)]
    print(f"\n  COCIENTE ENTRE POLOS en R={R:g}:  bruta(0,0,-R) / bruta(0,0,+R)")
    print("   l1 m1 l2 m2 | (-1)^{l1+l2} esperado | cociente medido")
    print("   ------------|-----------------------|------------------")
    for (l1, m1, l2, m2) in casos:
        n1, n2 = RADIAL.n_of_l(l1), RADIAL.n_of_l(l2)
        bf0 = brute_force_fermi_element_at(
            n1, l1, m1, n2, l2, m2, (0.0, 0.0, R), A_s, A_p).real
        bfpi = brute_force_fermi_element_at(
            n1, l1, m1, n2, l2, m2, (0.0, 0.0, -R), A_s, A_p).real
        ratio = bfpi / bf0
        esperado = phase_pi(l1, l2)
        print(f"   {l1:2d}{m1:3d}{l2:3d}{m2:3d} |           {esperado:+.0f}"
              f"           | {ratio:+.12f}")
        assert abs(ratio - esperado) < 1e-7, (l1, m1, l2, m2, ratio)


# ================================================================ TEST P1d
def test_p1d_consistencia_con_parity_factor_trimer():
    """
    Elemento(θ=0) + Elemento(θ=π) = [1+(-1)^{l1+l2}]·Elemento(θ=0):
    exactamente el `parity_factor` que `SymmetricLinearTrimer` tiene
    implementado para los DOS perturbadores (la suma de ambos individuales).
    """
    from trimero.systems.rb_neutral_perturber.linear_trimer import (
        SymmetricLinearTrimer,
    )

    tri = SymmetricLinearTrimer(n_manifold=N_MANIFOLD, radial_source="hydrogenic",
                                n_perturbers=2)
    R = 800.0
    ls = sorted(set(list(range(3, 12)) + [0, 1, 2]))
    print(f"\n  Identidad con el trímero (radiales hidrogenoides), R={R:g}:")
    print("   suma elem(0)+elem(π)  vs  parity_factor·elem(0)")
    worst = 0.0
    for l1 in ls:
        for l2 in ls:
            e0 = tri.pseudo.electron_element(l1, 0, l2, 0, R)
            suma = e0 + phase_pi(l1, l2) * e0
            con_factor = tri.parity_factor(l1, l2) * e0
            worst = max(worst, abs(suma - con_factor))
    print(f"   max|suma - parity_factor·elem(0)| sobre {len(ls)}² pares = {worst:.1e}")
    assert worst == 0.0
