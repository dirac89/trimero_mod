# Organización de resultados

```text
plots/
├── rb_krb_polar/
│   ├── data/       # matrices y curvas serializadas (.npz)
│   └── figures/    # figuras de producción (.png)
├── rb_rbcs_polar/
│   ├── data/
│   └── figures/
├── rb_neutral_perturber/
│   ├── data/
│   └── figures/
├── hybrid_neutral_polar/
│   ├── data/
│   └── figures/
├── geometry/       # esquemas geométricos compartidos
└── archive/        # resultados históricos; no reorganizar ni usar como producción
```

Los scripts de producción escriben en esta estructura. Los `.npz` polares
nuevos incorporan especie molecular, constantes, base y criterio de carácter
para impedir reutilizaciones incompatibles.

## Rb\*+RbCs disponible

- `data/fig1_ad_MJ0_n25.npz`: BOP completa, R=400–1800 a0, paso 5 a0.
- `data/fig1_ad_MJ1_n25.npz`: mismo barrido para M_J=1.
- `data/orientation_MJ0_n25.npz`: orientación, malla fina 100–800 a0 y cola
  hasta 1800 a0.
- `figures/bop_MJ0_MJ1_n25_Nmax6.png`.
- `figures/orientation_MJ0_n25_Nmax6.png`.
