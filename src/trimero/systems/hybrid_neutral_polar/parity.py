"""
Transformación θ=0→θ=π del elemento individual de V_Fermi.

Derivación formal con el operador paridad 𝒫 (𝒫|r̄⟩=|−r̄⟩, 𝒫∇𝒫†=−∇):

    V(−R⃗) = 𝒫 V(R⃗) 𝒫†

porque δ³ se traslada bajo conjugación y los DOS signos del término
bilineal p-wave se cancelan. Sobre estados hidrogenoides reales
⟨l m|𝒫 = (−1)^l ⟨l m|, luego

    ⟨l₁m₁|V(θ=π)|l₂m₂⟩ = (−1)^{l₁+l₂} ⟨l₁m₁|V(θ=0)|l₂m₂⟩

con las mismas reglas de selección que en θ=0 (Δm=0, |m|≤1; V_s sólo m=0).

Sumando ambos perturbadores del trímero simétrico se recupera exactamente
[1+(−1)^{l₁+l₂}]·V(θ=0), el `parity_factor` de `SymmetricLinearTrimer`.
Aquí la fase se usa SOLA: en el sistema híbrido el Rb neutro está en θ=π
sin compañero en θ=0.

Verificación por fuerza bruta contra ψ y ∇ψ numéricos en (0,0,−R):
tests/systems/hybrid_neutral_polar/test_transformacion_theta_pi.py
(peor desviación relativa 1.2e−6, ruido de diferencias finitas; cociente
bruta(θ=π)/bruta(θ=0) = (−1)^{l₁+l₂} a 12 dígitos).
"""

__all__ = ["phase_pi"]


def phase_pi(l1: int, l2: int) -> float:
    """
    Peso geométrico del par (l₁,l₂) para UN perturbador individual en θ=π:
    (−1)^{l₁+l₂}. ±1.0 exacto en coma flotante (exponente entero).
    """
    return (-1.0) ** (l1 + l2)
