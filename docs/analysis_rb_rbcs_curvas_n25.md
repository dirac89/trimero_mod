# Curvas completas Rb\*+RbCs, n=25

**Fecha:** 2026-08-23

**Modelo:** carga–dipolo puro, `N_max=6`, sin pseudopotencial de Fermi.

## Curvas BOP

Malla completa `R=400..1800 a0`, paso `5 a0`, base
`(25,l>=3)+26d+27p+28s`.

| bloque | mínimo global | E(1800 a0) | mínimos locales |
|---|---:|---:|---:|
| M_J=0 | -56.1414 GHz en 400 a0 | -1.8377 GHz | 8 |
| M_J=1 | -51.8942 GHz en 400 a0 | -1.5816 GHz | 7 |

En el extremo exterior ambos estados tienen peso de manifold `1.0000`. La
figura usa una ventana hasta `-60 GHz`, necesaria para no recortar el pozo de
RbCs. Que el mínimo global caiga en el borde inferior significa que este
barrido no demuestra un mínimo cerrado por debajo de `400 a0`; no se extrapola.

## Orientación

Para `M_J=0` se usó malla fina `100..800 a0` con paso `5 a0` y cola hasta
`1800 a0` con paso `50 a0`.

- 161 puntos, sin NaN.
- Rango: `<cos(theta_d)> = [-0.888294, 0.946571]`.
- La cota matemática `[-1,1]` se satisface.
- Se detectan 29 saltos con `|Delta cos|>0.02`; la figura los representa con
  trazo fino. La cola de paso 50 a0 no debe interpretarse como espectroscopía
  fina ni como seguimiento inequívoco de cruces.

## Archivos

- Datos: `plots/rb_rbcs_polar/data/`.
- Figuras: `plots/rb_rbcs_polar/figures/`.
- Convención global: `plots/README.md`.
