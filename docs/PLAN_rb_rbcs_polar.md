# Implementación de Rb\*+RbCs polar

**Fecha:** 2026-08-23

**Estado:** implementación estructural y convergencia inicial completadas

## Contrato físico

El sistema es Rb Rydberg + RbCs tratado como rotor rígido polar:

```text
H_ad(R) = H_A + B N² - d·[F_ion(R) + F_elec(R)]
```

No contiene pseudopotencial de Fermi. Usa la misma base electrónica que el
sistema Rb\*+KRb: `(n,l>=3)+(n+1)d+(n+2)p+(n+3)s`, bloques de `M_J`, y selección
de la curva por peso de manifold. Los parámetros actualmente adoptados son
`B=490.17 MHz` y `d=1.225 D`; antes de una predicción cuantitativa deben
asociarse explícitamente a isotopólogo y estado vibracional.

## Arquitectura aplicada

- `systems/polar_molecule.py`: catálogo único de KRb y RbCs con conversiones.
- `systems/polar_rydberg/`: motor carga–dipolo independiente de la especie y
  desacoplado de Fermi.
- `systems/rb_rbcs_polar/`: configuración pública de Rb\*+RbCs.
- `systems/rb_krb_polar/BOPSystem`: compatibilidad histórica; conserva el
  camino Fermi requerido por tests antiguos, pero deja de ser el motor de los
  scripts polares de producción.
- `HybridNeutralPolar`: comparte el catálogo RbCs y su límite `fermi=False`
  coincide con el nuevo sistema polar puro.

## Validación

1. Conversiones moleculares y rechazo de parámetros no físicos.
2. Igualdad bit a bit del motor genérico KRb con `BOPSystem(fermi=False)`.
3. Igualdad bit a bit de Rb\*+RbCs con el límite polar del híbrido.
4. Regresión de la Fig. 1 de KRb sin cambios.
5. Hermiticidad, conservación de `M_J` y límite asintótico cubiertos por los
   tests compartidos del operador carga–dipolo.

## Convergencia obtenida

Antes de aceptar curvas RbCs se ejecutará:

```bash
poetry run python scripts/check_polar_convergence.py \
  --molecule rbcs --n-manifold 25 --n-max 4 6 8 --r 500 1000 1500
```

Resultado para `n=25`, `M_J=0`:

| N_max | E(500) [GHz] | E(1000) [GHz] | E(1500) [GHz] |
|---:|---:|---:|---:|
| 4 | -49.09514246 | -45.27324783 | -7.33979301 |
| 6 | -49.26531658 | -45.29208507 | -7.34298235 |
| 8 | -49.26756684 | -45.29235428 | -7.34298493 |

El peor cambio `N_max=6→8` es `0.00225027 GHz`: `N_max=6` satisface tanto
el criterio de producción de `0.1 GHz` como el criterio fino de `0.01 GHz` en
estos tres puntos. Queda por comprobar `M_J=1` al producir su curva completa.

Producción:

```bash
poetry run python scripts/compute_bop_curve.py --molecule rbcs --n-manifold 25 --n-max 6 --mj 0 1
poetry run python scripts/compute_orientation_curve.py --molecule rbcs --n-manifold 25 --n-max 6
```

Los resultados se guardan en `plots/rb_rbcs_polar/` con especie, `B`, `d`,
`n`, `N_max`, `M_J` y criterio de carácter dentro del `.npz`.

## Criterios de aceptación física

- No usar `domain_bounds()`, remapeo `k(R)` ni ventana de resonancia p.
- Verificar `R→∞`: umbrales `E_Rb+B_RbCs N(N+1)`.
- Mantener `|<cos theta_d>|<=1` y registrar regiones de seguimiento ambiguo.
- Comparar KRb/RbCs con la misma base y malla, sin imponer de antemano mayor
  profundidad: el mayor dipolo y menor `B` sólo son tendencias a comprobar.
- No alimentar dinámica no adiabática hasta fijar `N_max` y validar la curva.
