# El autovalor más bajo del híbrido NO es una curva de ligadura: es el umbral 37p

**Fecha**: 2026-08-22
**Autor**: Javier Aguilera
**Relevancia**: Revoca la interpretación física de las primeras curvas del híbrido (`analysis_hibrido_fases12_golden_regen.md` §4). Las curvas E0(R2) producidas son el umbral atómico del vecino 37p vestido por H_mol, no una curva BOP del manifold. Establece el criterio de selección correcto antes de cualquier producción física del híbrido.
**Tipo**: analysis

## Resumen

Las curvas E0(R2) a R1∈{600,900,1100} a₀ (profundidades −102…−111 GHz, sin
mínimos locales) se generaron seleccionando **el autovalor más bajo sin más**.
El diagnóstico pedido muestra que ese estado es **99.7 % vecino 37p** (l=1,
m=0), cuyo nivel desnudo está a **−102.360 GHz** bajo el manifold n=35,l≥3.
La forma de la curva es ese umbral constante más un vestido de H_mol que
decae con R2 (−8.6 GHz en R2=500, −0.4 GHz en R2=1500); V_Fermi aporta
sólo −0.18 GHz. La coincidencia numérica con la "cota polar pura" de L-b
(−102.54 GHz) era la misma cosa en ambos sistemas: el umbral p vestido,
no ligadura de manifold.

Consecuencia: hay que seleccionar la curva del híbrido **por carácter**
(peso >50 % en el manifold), como en el sistema polar. Existe y es
seleccionable: queda a −79.5 / −71.7 / −71.5 GHz en R2=500/1000/1500
(peso de manifold 95–99 %). Además, N_max=2 no está convergido para
cuantitativa en R2 pequeño (Δ≈0.6 GHz al subir a 4), aunque sí para
tendencias.

## Palabras Clave

- carácter de autoestado
- umbral atómico (vecino 37p)
- criterio de selección de curva
- convergencia en N_max
- sistema híbrido

---

## 0. Nota metodológica: dos bugs nuestros de diagnóstico (no de producción)

Durante el diagnóstico dos errores propios distorsionaron resultados
intermedios; quedan registrados porque explican idas y venidas:

1. **Constante de unidades**: usamos Eh→GHz = 657.97 en los scripts de
   diagnóstico en vez de `GHZ_PER_HARTREE` = 6 579 683.92 (×10⁴). Producción
   siempre fue correcta (`spectrum_ghz` usa la constante del módulo). Detectado
   comparando matriz a matriz: `max|H_manual − H_metodo| = 0.0`.
2. **Etiqueta de estado**: la tupla de estado del bloque es `(l, m_l, N, M_N)`
   — **no incluye n**. Leer `st[0]` como `n` produjo etiquetas absurdas
   ("1s=0.99"); lo correcto: `l≥3` ⇒ manifold n=35; `l=0,1,2` ⇒ vecinos
   38s/37p/36d.

## 1. T1 — Carácter del estado que da E0

Referencias atómicas desnudas respecto al manifold (con defectos reales del
átomo modelo del código):

```
36d : E-E_man =   -54.019 GHz
37p : E-E_man =  -102.360 GHz      <- el sospechoso
38s : E-E_man =   -20.267 GHz
```

Descomposición del autovector de E0 (R1=900, H completo):

| R2 (a₀) | MANIFOLD | 37p | 36d | 38s | rotor N=0/N=1/N=2 |
|---|---|---|---|---|---|
| 500 | 0.0027 | **0.9971** | 0.0001 | 0.0000 | 0.349 / 0.498 / 0.154 |
| 1000 | 0.0028 | **0.9971** | 0.0000 | 0.0000 | 0.542 / 0.411 / 0.048 |
| 1500 | 0.0033 | **0.9966** | 0.0000 | 0.0000 | 0.761 / 0.232 / 0.007 |

Top componentes en R2=500: `(l=1,m=0,N=1):0.495`, `(1,0,N=0):0.348`,
`(1,0,N=2):0.151`. **Dominio absoluto: 37p (>99 %), m=0**, con rotor
sustancialmente mezclado por el campo del electrón (N=1 domina en R2=500;
N=0 en R2≥1000). La hipótesis "E0 = energía atómica desnuda del vecino" se
**confirma**: E0 en H_A solo es −102.36043 GHz = exactamente el nivel 37p.

## 2. T2 — Descomposición de la profundidad

Mínimo del espectro en cada pieza (R1=900):

| pieza | R2=500 | R2=1000 | R2=1500 |
|---|---|---|---|
| (a) H_A solo | −102.36043 | −102.36043 | −102.36043 |
| (b) H_A+H_mol | −110.95794 | −103.88410 | −102.78290 |
| (c) H_A+V_Fermi | −102.54232 | −102.54232 | −102.54232 |
| (d) completo | −111.10797 | −104.04045 | −102.95174 |

Lectura:

- **Los −102.4 GHz base son H_A por sí solo** (el nivel 37p desnudo).
- **H_mol aporta toda la dependencia con R2**: −8.60 GHz en R2=500,
  −1.52 en 1000, −0.42 en 1500 (vestido del campo eléctrico del electrón
  Rydberg sobre el rotor; decae al abrir R2).
- **V_Fermi aporta −0.18 GHz casi constante** (y explica las diferencias
  ~0.1–0.2 GHz entre las tres curvas R1 de producción).
- La curva E0(R2) "sin mínimos locales" es entonces **umbral constante +
  cola monótona de vestido Stark**: no hay pozo que interpretar.

## 3. T3 — Convergencia en N_max (R1=900)

| | E0 (R2=500) | E0 (R2=1000) | E_carácter (R2=500) | E_carácter (R2=1000) |
|---|---|---|---|---|
| N_max=2 (dim 307) | −111.10797 | −104.04045 | −79.47540 | −71.73432 |
| N_max=4 (dim 835) | −111.69137 | −104.05634 | −80.10430 | −71.78772 |
| Δ | **−0.583** | −0.016 | **−0.629** | −0.053 |

En R2=1000 ambos objetos están convergidos a <0.06 GHz; **en R2=500 el error
de truncamiento de rotor es ~0.6 GHz (0.5 %)** para cualquier criterio. Para
afirmaciones cuantitativas cerca del mínimo R2 hay que subir N_max (y
estudiar la convergencia sistemáticamente); para tendencias y estructura,
N_max=2 cualitativamente suficiente.

## 4. T4 — Valoración

**La curva E0 producida NO es una curva BOP de ligadura molecular: es el
umbral atómico del vecino 37p vestido por H_mol.** Confirmado por tres vías
independientes: carácter >99 % 37p, descomposición por piezas (base −102.36 =
nivel desnudo), y ausencia de estructura (monótona).

Además, la lectura previa "converge hacia la cota polar pura ~−102.8 GHz"
era una coincidencia numérica mal razonada: el test L-b comparaba el mismo
objeto (umbral p vestido) en dos módulos distintos; consistencia válida como
test de código, etiqueta física incorrecta.

**Criterio correcto**: seguir la curva de **CARÁCTER de manifold** (>50 %
peso en n=35, l≥3, elegido por R), como en el sistema polar. Ese objeto
existe, es seleccionable y estable en los tres puntos:

| R2 (a₀) | E_carácter | peso manifold | mezcla principal |
|---|---|---|---|
| 500 | −79.47540 | 0.954 | 36d 3.9 % |
| 1000 | −71.73432 | 0.954 | 36d 4.0 % |
| 1500 | −71.54387 | 0.989 | 36d 0.9 % |

Su física es distinta a la de E0: es el manifold hidrogenoide degenerado
partido linealmente por los campos de la molécula (dipolo permanente del
KRb + campo del electrón) — análogo invertido del fig1 polar. Ojo: su
interpretación como "ligadura" exige además seguimiento adiabático por
carácter punto a punto (aquí idx 7 en los tres puntos, pero el índice no
es garantía) y convergencia en N_max.

## Conclusiones y plan

1. Descartar como física las tres curvas E0(R2) de
   `plots/hybrid_neutral_polar/data/hybrid_curves_R1*.npz`: trazaban el umbral
   37p. No borrarlas (son consistentes como datos), pero reetiquetar o
   regenerar con el criterio de carácter.
2. Modificar `scripts/compute_hybrid_curves.py`: seleccionar por peso de
   manifold >50 % (guardar también el peso y el desglose por vecino), y
   añadir barrido de convergencia N_max∈{2,4,6} antes de producir.
3. Revisar la etiqueta del test de límite L-b (documentar que compara
   umbrales-p vestidos entre módulos).
4. Los goldens MJ0/MJ1 y la regresión del sistema polar **no se tocan**
   (ronda cerrada; este análisis no afecta al sistema polar: allí el
   criterio de carácter ya era parte del método).

## Referencias

- [`analysis_hibrido_fases12_golden_regen.md`](analysis_hibrido_fases12_golden_regen.md) — §4 es el objeto corregido aquí.
- [`analysis_base_correcta_3_vecinos.md`](analysis_base_correcta_3_vecinos.md) — composición de vecinos y su papel en el sistema polar.
- Scripts de diagnóstico (temporal):
  `/var/folders/.../opencode/wigner_fix/diag_caracter_e0_v2.py`.
