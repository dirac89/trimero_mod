import numpy as np
from scipy.special import lpmv, sph_harm, factorial, genlaguerre, gammaln

# Spherical harmonics (solo parte asociada de Legendre, normalización)
def Spherical(l, m, x1):
    m = int(m)
    l = int(l)
    if m < 0:
        # lpmv requiere |m|, y el signo se ajusta manualmente
        return (
            lpmv(abs(m), l, x1)
            * (-1) ** (-m)
            * factorial(l + m) / factorial(l - m)
        )
    else:
        return lpmv(m, l, x1)

# Derivada radial de la función de onda hidrogenoide
# Necesita la función hydrogenicR definida abajo
def DRnl(n, l, r):
    n = int(n)
    l = int(l)
    if l < n - 1:
        xn = 2.0 * r / n
        N1 = 2.0 * np.sqrt(factorial(n - l - 1))
        N2 = n ** 2 * np.sqrt(factorial(n + l))
        Nnl = N1 / N2
        term1 = ((l / r) - (1.0 / n)) * hydrogenicR(n, l, 1.0, r)
        laguerre = genlaguerre(n - l - 2, 2 * l + 2)
        term2 = (xn ** l) * np.exp(-r / n) * laguerre(xn)
        return term1 - Nnl * term2 * 2.0 / n
    else:
        return 0.0

# Derivada angular DOlm
def DOlm(l, m, theta):
    x = np.cos(theta)
    term1 = 0.0
    term2 = 0.0
    l = int(l)
    m = int(m)
    if l > 0:
        C1 = np.sqrt((2.0 * l + 1.0) * factorial(l - m))
        C2 = np.sqrt((4.0 * np.pi) * factorial(l + m))
        Clm = C1 / C2 * 0.5
        # term1
        if m + 1 >= 0 and m + 1 <= l:
            term1 = lpmv(m + 1, l, x)
        elif m + 1 < 0 and -m - 1 <= l:
            term1 = (-1) ** (1 + m) * factorial(l + 1 + m) / factorial(l - 1 - m) * lpmv(-m - 1, l, x)
        # term2
        if m - 1 < 0 and -m + 1 <= l:
            term2 = (l + m) * (l - m + 1) * lpmv(1 - m, l, x) * (-1) ** (1 - m) * factorial(l - 1 + m) / factorial(l + 1 - m)
        elif m - 1 >= 0 and m - 1 <= l:
            term2 = (l + m) * (l - m + 1) * lpmv(m - 1, l, x)
        return Clm * (term1 - term2)
    else:
        return 0.0

# Derivada angular DPhilm
def DPhilm(l, m, theta):
    x = np.cos(theta)
    term1 = 0.0
    term2 = 0.0
    l = int(l)
    m = int(m)
    if l > 0:
        C1 = np.sqrt((2.0 * l + 1.0) * factorial(l - m))
        C2 = np.sqrt((4.0 * np.pi) * factorial(l + m))
        Clm = -C1 / C2 * 0.5
        # term1
        if m + 1 >= 0 and m + 1 <= l + 1:
            term1 = lpmv(m + 1, l + 1, x)
        elif m + 1 < 0 and -m - 1 <= l + 1:
            term1 = lpmv(-m - 1, l + 1, x) * (-1) ** (1 + m) * factorial(l + 2 + m) / factorial(l - m)
        # term2
        if m - 1 < 0 and -m + 1 <= l + 1:
            term2 = (l - m + 1) * (l - m + 2) * lpmv(1 - m, l + 1, x) * (-1) ** (1 - m) * factorial(l + m) / factorial(l + 2 - m)
        elif m - 1 >= 0 and m - 1 <= l + 1:
            term2 = (l - m + 1) * (l - m + 2) * lpmv(m - 1, l + 1, x)
        return Clm * (term1 + term2)
    else:
        return 0.0

# Función radial hidrogenoide (Z=1 por defecto)
def hydrogenicR(n, l, Z, r):
    n = int(n)
    l = int(l)
    # Normalización
    rho = 2.0 * Z * r / n
    norm = np.sqrt((2.0 * Z / n) ** 3 * factorial(n - l - 1) / (2.0 * n * factorial(n + l)))
    laguerre = genlaguerre(n - l - 1, 2 * l + 1)
    return norm * np.exp(-rho / 2) * rho ** l * laguerre(rho) 