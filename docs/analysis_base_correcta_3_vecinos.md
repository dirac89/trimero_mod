# La base electrónica correcta: manifold + (n+1)d + (n+2)p + (n+3)s

> ⚠️ **AVISO PARCIAL (2026-08-20).** La corrección central de este documento —la
> base es manifold + (n+1)d + (n+2)p + (n+3)s— **sigue siendo válida y se usa**
> en el cálculo definitivo. Lo que queda invalidado es todo lo que depende del
> **pseudopotencial de Fermi**, que no forma parte del Hamiltoniano de
> Aguilera-Fernández et al. 2015: la §2.1 (acortamiento del dominio del remapeo
> por el (n+2)p) y la §5 (curva BOP con V_Fermi, ventana de resonancia, cola
> butterfly) **no aplican a esa comparación**. El §5.1 —que hace falta
> identificar la curva por CARÁCTER y no por índice fijo— sí sigue aplicando,
> porque el problema lo causan los cruces evitados con (n+1)d y (n+2)p, que
> están en H_mol. Ver `docs/analysis_fig1_carga_dipolo_sin_fermi.md`.

**Fecha**: 2026-08-19
**Autor**: Javier Aguilera
**Relevancia**: Corrige la composición de la base electrónica usada en TODA la
sesión anterior. El texto de Aguilera-Fernández et al. (2015) especifica **tres**
niveles vecinos individuales, no uno; hasta ahora sólo incluíamos (n+3)s.
Se mide qué cambia y qué no.
**Tipo**: analysis

## Resumen

El paper que da las Fig. 1 y Fig. 4 que queremos reproducir —Aguilera-Fernández,
Sadeghpour, Schmelcher & González-Férez, *J. Phys.: Conf. Ser.* **635**, 012023
(2015), arXiv:1507.07972— define la base como

> «the (n, l≥3) degenerate manifold, and the energetically neighboring levels
> (n+1)d, (n+2)p, and (n+3)s»

Nuestra base tenía manifold + (n+3)s solamente: **faltaban (n+1)d y (n+2)p**.
Corregido, con el manifold como parámetro y los vecinos derivados de n.

Tres conclusiones, y la tercera es la que importa:

1. **El coste crece poco**: dim(M_J=0) 1016 → 1064 (+4.7 %), 0.888 → 0.946 s por
   punto (+6.5 %).
2. **Los dos observables ya calculados apenas cambian**: ⟨cos θ_d⟩ se mueve
   ≤2.5×10⁻³ y el cruce de 0.78 sigue entre 380 y 400 a₀; la convergencia en N
   del estado de manifold es 10⁻¹³–10⁻¹¹ relativa en ambas bases. **La base
   incompleta era, para esto, una aproximación razonable.**
3. **Pero el dominio del remapeo k(R) se acorta 48.8 a₀** (n=24) y **50.8 a₀**
   (n=25), porque el nivel más ligado de la base pasa a ser el (n+2)p. Eso sí es
   un cambio real, y afecta a hasta dónde llegan nuestras curvas.

## Palabras clave

- base electrónica, manifold cuasi-degenerado, niveles vecinos
- (n+1)d, (n+2)p, (n+3)s
- convergencia en N, orientación ⟨cos θ_d⟩
- dominio del remapeo semiclásico k(R)

## 1. Qué se generalizó

Ya no hay ningún «24, 27» ni «25, 28» escrito a mano. La composición de la base
vive en **una sola función**:

```python
# src/trimero/systems/rb_defects.py
NEIGHBOR_DN = {0: 3, 1: 2, 2: 1}          # Δn = 3 - l

def neighbor_levels(n_manifold: int) -> dict:
    """{l: n} de los vecinos: {0: n+3, 1: n+2, 2: n+1}."""
    return {l: n_manifold + dn for l, dn in NEIGHBOR_DN.items()}
```

y todo lo demás la consume:

| pieza | qué depende del manifold | cómo se resuelve ahora |
|---|---|---|
| `CoupledBasis` | l máximo del manifold, y qué vecinos hay | `manifold_l_max=n-1`, `neighbor_l=(0,1,2)` |
| `RadialBasis` | n de cada radial | `neighbor_n_eff` = round(n\*) de cada vecino |
| `rydberg_diagonal` | energías diagonales de H_a | `neighbor_levels(n_manifold)` |
| `FermiPseudopotential` / `n_star_of_l` | n\* del remapeo k(R) | `neighbor_levels` + `n_star_nl` |
| `ChargeDipoleHamiltonian` / `RydbergElectronField` | — | ya era genérico en l: recorre `radial.l_values` |
| `BOPSystem` | todo lo anterior a la vez | `BOPSystem(n_manifold=25)` |

**Detalle que hace que esto sea barato**: como los tres vecinos tienen l = 0, 1,
2 y el manifold empieza en l = 3, **l sigue identificando unívocamente el
nivel**. El estado sigue siendo `(l, m_l, N, M_N)`; no hubo que meter n en la
tupla ni tocar los Hamiltonianos, que ya eran genéricos en l.

`neighbors=(0,)` reproduce exactamente la base incompleta anterior, y se usa en
este documento para medir la diferencia.

## 2. Tamaño de la base y coste (n=24, N_max=6)

Composición electrónica:

| bloque | l | estados m_l | n\* | radial hidrogenoide |
|---|---|---|---|---|
| manifold n=24 | 3..23 | 567 | 24.000000 | n = 24 |
| 25d | 2 | 5 | 23.653793 | n = 24 |
| 26p | 1 | 3 | 23.351184 | n = 23 |
| 27s | 0 | 1 | 23.867894 | n = 24 |
| **total** | | **576 = 24²** | | |

Que el total sea exactamente **24²** es una comprobación bonita: la base son
todos los estados hidrogenoides de n=24 en número, aunque tres de ellos sean en
realidad niveles de n distinto con defecto cuántico.

| | base incompleta | base correcta | cambio |
|---|---|---|---|
| electrónicos | 568 | 576 | +8 |
| dim total (×49) | 27 832 | 28 224 | +392 |
| **dim(M_J=0), N_max=6** | **1016** | **1064** | **+48 (+4.7 %)** |
| dim(M_J=0), N_max=8 | 1640 | 1704 | +64 |
| construir H | 0.799 s | 0.846 s | +5.9 % |
| `eigvalsh` | 0.089 s | 0.100 s | +12 % |
| **por punto de barrido** | **0.888 s** | **0.946 s** | **+6.5 %** |

**Valoración sin adornos: el crecimiento es modesto.** +4.7 % en dimensión y
+6.5 % en tiempo no cambia la viabilidad de nada. El coste sigue dominado por
construir H (89 % del tiempo), no por diagonalizar; añadir 48 filas al bloque no
mueve esa aguja. Un barrido completo pasa de ~3.3 a ~3.5 minutos.

### 2.1 Lo que sí cambia de verdad: el dominio del remapeo

Los puntos de retorno clásicos `2n*²` de cada nivel de la base (n=24):

| nivel | n\* | 2n\*² [a₀] |
|---|---|---|
| **26p** | 23.351184 | **1090.6** ← el más corto |
| 25d | 23.653793 | 1119.0 |
| 27s | 23.867894 | 1139.4 |
| manifold 24f | 24.000000 | 1152.0 |

El dominio del remapeo lo fija el par más ligado presente en la matriz, que
**deja de ser (n+3)s y pasa a ser (n+2)p**:

| | n=24 | n=25 |
|---|---|---|
| base incompleta | [105.61, **1138.92**] a₀ | [106.37, **1236.32**] a₀ |
| base correcta | [105.61, **1090.16**] a₀ | [106.37, **1185.50**] a₀ |
| acortamiento | **−48.8 a₀** | **−50.8 a₀** |

⚠️ **Lectura correcta de este límite.** No es «donde se acaba la física». Es
donde deja de estar definido el **remapeo semiclásico k(R)** con el que leemos
las tablas de longitudes de dispersión: pasado el punto de retorno clásico del
par más ligado, k² < 0 y no hay A_s, A_p que leer. `H_mol` (carga-dipolo) **no**
tiene esa restricción — se evalúa con la función de onda hidrogenoide real, que
decae exponencialmente pero no se anula ahí. Resolverlo queda fuera de esta
ronda por decisión explícita; lo que sí se hace es dejar de describirlo mal.

Consecuencia práctica inmediata: un test que evaluaba el Fermi en R=1100 a₀
(`test_f3d_magnitude_vs_charge_dipole`) dejó de ser válido y se movió a
R=1080 a₀.

## 3. Tabla I: qué es cada fila

La Tabla I son **energías atómicas puras**: no dependen de la base ni requieren
diagonalizar nada, así que sus valores no cambian (siguen siendo los de
`docs/analysis_verificacion_tabla_I.md` §9, con δ₀(ns)=3.13180 y peor
desviación 0.0040 GHz). Lo que sí había que confirmar es **qué papel juega cada
fila en la base**, y ahora está confirmado:

| n | l | E_nl [E_h] | n\* | ΔE vs manifold [GHz] | ¿en la base? |
|---|---|---|---|---|---|
| 25 | 3 | −8.000000000000000e−04 | 25.0000000 | 447.78404 | no — contraste |
| 28 | 0 | −8.085207285236110e−04 | 24.8679178 | 391.72034 | no — contraste |
| 26 | 2 | −8.226318540067425e−04 | 24.6537076 | 298.87360 | no — contraste |
| 27 | 1 | −8.431955228646972e−04 | 24.3512274 | 163.57116 | no — contraste |
| **24** | **3** | **−8.680555555555555e−04** | **24.0000000** | **0.00000** | **SÍ — manifold** |
| **27** | **0** | **−8.776913467604867e−04** | **23.8678936** | **63.40046** | **SÍ — vecino l=0** |
| **25** | **2** | **−8.936519942688221e−04** | **23.6537928** | **168.41648** | **SÍ — vecino l=2** |
| **26** | **1** | **−9.169637842916073e−04** | **23.3511842** | **321.80069** | **SÍ — vecino l=1** |

**Los tres vecinos del paper —25d, 26p, 27s— están los tres en la Tabla I: 3 de
3.** Los otros cuatro niveles son el **siguiente de cada serie** (25f, 26d, 27p,
28s), es decir puntos de contraste que acotan el manifold por arriba y **no
entran en la base**. La tabla está construida simétricamente alrededor del
manifold: cuatro niveles por encima y tres por debajo, y son precisamente los
tres de debajo los que se incluyen.

## 4. Los dos chequeos que dependían de la base incompleta

### 4.1 Convergencia en N (N_max=6 vs N_max=8)

⚠️ **El criterio importa más de lo que parecía.** El paper habla de una
diferencia relativa < 2×10⁻⁶, y hasta ahora la medíamos sobre «los 10
autovalores más bajos». Con la base completa **eso deja de medir lo que
queremos**:

```
composición de los 10 más bajos en R=900 a0 (l dominante : E-E_man [GHz])
  base incompleta: l=0:-63.7  l=0:-61.3  l=0:-56.9  l=0:-50.2  l=18:-41.7
                   l=0:-41.3  l=15:-39.9 l=17:-35.6 l=0:-30.3  l=0:-30.0
  base correcta:   l=1:-322.6 l=1:-320.0 l=1:-319.6 l=1:-319.5 l=1:-315.6
                   l=1:-315.1 l=1:-315.1 l=1:-308.9 l=1:-308.4 l=1:-308.4
```

Los diez más bajos pasan a ser **todos del 26p** (−322 a −308 GHz), que está
321.8 GHz por debajo del manifold y apenas se acopla al rotor. Por eso su
convergencia en N es trivialmente perfecta: no es una mejora del cálculo, es
otra pregunta.

Midiendo las dos cosas por separado:

| base | criterio | R [a₀] | peor \|ΔE\| [MHz] | peor \|ΔE\|/\|E\| | ¿< 2×10⁻⁶? |
|---|---|---|---|---|---|
| incompleta | 10 más bajos | 600 | 18.974 | 3.072×10⁻⁶ | **NO** |
| correcta | 10 más bajos | 600 | 19.649 | 3.178×10⁻⁶ | **NO** |
| incompleta | 10 más bajos | 900 | 0.771 | 1.342×10⁻⁷ | sí |
| correcta | 10 más bajos | 900 | 0.000 | 3.98×10⁻¹⁴ | sí (espurio, ver arriba) |
| incompleta | **manifold (k=0)** | 600 | 0.000 | 7.60×10⁻¹⁴ | **sí** |
| correcta | **manifold (k=0)** | 600 | 0.000 | 1.18×10⁻¹³ | **sí** |
| incompleta | **manifold (k=4)** | 900 | 0.0001 | 1.73×10⁻¹¹ | **sí** |
| correcta | **manifold (k=52)** | 900 | 0.0002 | 3.61×10⁻¹¹ | **sí** |

**Veredicto: completar la base NO cambia la conclusión de convergencia.** El
estado físicamente relevante —el más bajo con carácter de manifold— está
convergido a 10⁻¹³–10⁻¹¹ relativo en las dos bases. El incumplimiento a R=600
con el criterio «10 más bajos» sigue ahí y sigue viniendo de un estado casi
degenerado suelto, no de un fallo de la base.

(El índice del estado de manifold sube de k=4 a k=52 en R=900 simplemente porque
ahora hay 48 estados de 26p y 25d por debajo.)

### 4.2 Orientación ⟨cos θ_d⟩

Autoestado más bajo con carácter de manifold, bloque M_J=0. `‖C‖_F` pasa de
17.1957 a 17.7331 sólo por haber más estados.

| R [a₀] | ⟨cos⟩ incompleta | ⟨cos⟩ correcta | Δ⟨cos⟩ |
|---|---|---|---|
| 110 | +0.930172 | +0.930173 | +6.8×10⁻⁷ |
| 170 | +0.900148 | +0.900150 | +1.8×10⁻⁶ |
| 230 | +0.868343 | +0.868363 | +2.0×10⁻⁵ |
| 290 | +0.842665 | +0.842667 | +2.3×10⁻⁶ |
| 320 | +0.826686 | +0.826801 | +1.2×10⁻⁴ |
| 380 | +0.807021 | +0.807119 | +9.8×10⁻⁵ |
| 400 | +0.769873 | +0.769989 | +1.2×10⁻⁴ |
| 440 | +0.725637 | +0.725682 | +4.5×10⁻⁵ |
| 480 | +0.712773 | +0.712722 | −5.1×10⁻⁵ |
| 520 | +0.684758 | +0.684694 | −6.3×10⁻⁵ |
| **560** | +0.672983 | +0.670482 | **−2.5×10⁻³** |

- Peor \|Δ⟨cos⟩\| en la rama suave: **2.5×10⁻³**, y sólo en el último punto,
  R=560 a₀, ya cerca de la resonancia de onda p.
- **El cruce de 0.78 sigue entre R = 380 y 400 a₀** en las dos bases.
- El valor de borde (R=110 a₀) es 0.930172 → 0.930173.

**Veredicto: el efecto de añadir 25d y 26p sobre la orientación es pequeño**,
del orden de 10⁻⁵–10⁻⁴ en casi todo el rango. La conclusión de
`analysis_verificacion_tabla_I.md` §11 —que no se puede confirmar ni refutar el
0.78 del paper porque nuestra curva no tiene máximo interior— **se mantiene sin
cambios**.

### 4.3 Entonces, ¿era grave la base incompleta?

Para **estos dos observables, no**: los cambios están entre 10⁻¹³ y 10⁻³
relativos. Para el **dominio en R, sí**: se perdían/ganaban 49-51 a₀ de alcance
y, sobre todo, se atribuía el límite al nivel equivocado. Y en general no es una
respuesta que se pueda extrapolar: los estados 26p y 25d están 322 y 168 GHz por
debajo del manifold, así que dominan el fondo del espectro y cambian qué es «el
k-ésimo autovalor» en cualquier análisis que ordene por energía — como se acaba
de ver en §4.1.

## 5. Curva BOP de M_J=0 para n=25, con la base correcta

Base: manifold (25, l≥3) + 26d + 27p + 28s. `dim(M_J=0) = 1113`, 217 puntos,
paso 5 a₀, **245.8 s** (1.13 s por punto, con `eigh` porque hacen falta los
autovectores).

```
dominio del remapeo: R ∈ [106.3734, 1185.4957] a0
  par más restrictivo 27p-27p, retorno clásico 2n*² = 1185.99 a0
  R_max queda 0.49 a0 por debajo: lo fija el borde de la tabla (R'=2448 a0)
cero de energía: E(n=25, l≥3) + KRb(N=0) = -8.000000000000e-04 E_h
```

### 5.1 ⚠️ El índice fijo deja de identificar la curva del manifold

Éste es el hallazgo no anticipado de esta ronda, y obliga a cambiar el método.

Hasta ahora identificábamos la curva del manifold por su **índice ordenado**,
fijado en el borde del dominio: por la regla de no cruce, el k-ésimo autovalor
ordenado ES la k-ésima curva adiabática, así que basta identificarlo una vez
(`analysis_curva_bop_MJ0.md` §2). Eso sigue siendo cierto **como afirmación
sobre curvas adiabáticas**, pero con la base completa deja de servir para lo que
lo usábamos: los estados de 26d (−148.9 GHz) y 27p (−284.2 GHz) caen justo en el
rango de energía relevante, y la curva de índice fijo **cambia de carácter** al
atravesar sus cruces evitados.

Con `k = 53` (el índice de carácter manifold en R_max = 1185.5 a₀):

| R [a₀] | E(k=53) [GHz] | peso de manifold | E de la más baja con carácter manifold |
|---|---|---|---|
| 200 | −31.44 | 1.00 | −33.66 (k=52) |
| 300 | −21.17 | **0.00** | −16.60 (k=54) |
| 400 | −22.22 | **0.00** | −11.01 (k=54) |
| 500 | −22.56 | **0.00** | −13.60 (k=54) |
| 900 | −37.55 | 0.99 | −41.29 (k=52) |
| 1000 | −33.98 | **0.01** | −34.50 (k=52) |

**La curva de índice fijo no tiene carácter de manifold en 111 de los 217 puntos
del barrido** (R ∈ [110, 1040] a₀). El primer barrido, hecho con ese criterio,
daba «99 % del dominio desplazado más de 20 GHz»: estaba midiendo un estado
26d/27p, no el manifold.

**Método corregido**: se guarda el peso de manifold de las 60 curvas más bajas
en cada R y se toma la **curva de CARÁCTER** — la más baja con peso > 50 % en
cada R. Recorre k = 0..55 con 27 cambios de índice, y su peso de manifold es
1.000 en prácticamente todo el dominio. La figura dibuja las dos para que se vea
dónde se separan.

### 5.2 Y con esto, el efecto de completar la base vuelve a ser pequeño

Comparando **curva de carácter contra curva de carácter** entre las dos bases:

| R [a₀] | base correcta [GHz] | base incompleta [GHz] | Δ |
|---|---|---|---|
| 255 | −21.8274 | −21.8295 | +0.0021 |
| 265 | −20.3060 | −20.3145 | +0.0084 |
| 400 | −11.0080 | −11.0101 | +0.0021 |
| 405 | −10.9995 | −11.0021 | +0.0027 |
| 1015 | −34.7675 | −34.9090 | +0.1415 |
| 1115 | −32.2370 | −32.3592 | +0.1222 |
| 1180 | −29.2464 | −29.4724 | +0.2259 |

**Diferencias de 0.002 a 0.23 GHz**, creciendo con R. Es decir: igual que en
n=24 (§4), añadir 26d y 27p **apenas mueve la energía de la curva**. Lo que
rompía era el criterio de selección, no la física.

(Nota: la tabla de mínimos del primer barrido de n=25 —mínimos en 265, 990 y
1180 a₀— era del índice fijo y queda descartada. Incluso con la base incompleta,
el índice fijo k=6 se desviaba ~0.5 GHz de la curva de carácter en R≈265 a₀.)

### 5.3 Resonancia de onda p, recalculada para n=25

```
eps_polo = 0.00090887 E_h = 24.732 meV     <- IDÉNTICO al de n=24, bit a bit
R_polo   = 585.183 a0                       (era 562.771 en n=24)
FWHM     = 14.258 a0                        (era 13.187)
VENTANA EXCLUIDA: R ∈ [556.67, 613.70] a0   = R_polo ± 2×FWHM
  11 de 217 puntos dentro; 4 con A_p interpolado a través del polo
BANDA DE COLA butterfly: R ∈ (613.70, 705.00] a0, 19 puntos
```

La energía del polo es la misma porque la resonancia es una propiedad de
e⁻+Rb(5S), no del manifold; lo que cambia es su posición en R, porque
ε = E_manifold + 1/R. Dentro de la ventana la curva llega a **−6459 GHz**, sin
significado físico (misma limitación de rango cero de siempre, no corregida por
decisión: ver `analysis_ventana_exclusion_resonancia.md` §6).

### 5.4 Umbrales asintóticos y residuo en el borde

`ΔE(28s) = E(28s) − E(n=25) = −56.0637 GHz` (umbral N=0), `B(KRb) = 1.114 GHz`:

| N | B·N(N+1) [GHz] | umbral ΔE+B·N(N+1) [GHz] |
|---|---|---|
| 0 | 0.0000 | −56.0637 |
| 1 | 2.2280 | −53.8357 |
| 2 | 6.6840 | −49.3797 |
| 3 | 13.3680 | −42.6957 |
| 4 | 22.2800 | −33.7837 |
| **5** | **33.4200** | **−22.6437** |
| **6** | **46.7880** | **−9.2757** |

En el borde del dominio, R = 1185.50 a₀:

| k | E−E_man [GHz] | peso manifold | ⟨N⟩ | asignación | umbral [GHz] | residuo [GHz] |
|---|---|---|---|---|---|---|
| 0 | −284.4965 | 0.0001 | 0.0 | 27p, N=0 | −284.2129 | −0.2837 |
| 1 | −282.1669 | 0.0001 | 1.0 | 27p, N=1 | −281.9849 | −0.1820 |
| **48** | **−56.2534** | 0.0028 | 0.0 | **28s, N=0** | −56.0637 | **−0.1897** |
| **49** | **−54.0183** | 0.0034 | 1.0 | **28s, N=1** | −53.8357 | **−0.1826** |
| **50** | **−49.5649** | 0.0033 | 2.0 | **28s, N=2** | −49.3797 | **−0.1852** |
| **51** | **−42.8827** | 0.0034 | 3.0 | **28s, N=3** | −42.6957 | **−0.1870** |
| **52** | **−33.9728** | 0.0037 | 4.0 | **28s, N=4** | −33.7837 | **−0.1891** |
| **54** | **−22.8355** | 0.0042 | 5.0 | **28s, N=5** | −22.6437 | **−0.1918** |
| **57** | **−9.4391** | 0.0043 | 6.0 | **28s, N=6** | −9.2757 | **−0.1634** |

```
residuos (serie 28s): media -0.1841 GHz, sigma 0.0089 GHz
                      dispersión (max-min) 0.0284 GHz
                      dispersión/|media| = 0.154  ->  APROXIMADAMENTE CONSTANTE
```

**Los siete umbrales del 28s se reproducen con un residuo constante de
−0.184 ± 0.009 GHz**, dispersión del 15 % de su valor. Es el comportamiento
esperado a R finito: la curva no ha llegado a su asíntota. **No hay
discrepancia que reportar aquí.** Para comparar, en n=24 el residuo era
≈ −0.19 GHz con la misma estructura; que salga prácticamente el mismo número
para dos manifolds distintos refuerza que es efecto de R finito y no un error de
los umbrales.

(La asignación se hace por peso de autovector y ⟨N⟩, no por orden: con la base
completa los 48 estados más bajos son de 27p y 26d, y los del 28s empiezan en
k=48. Los residuos del 27p —−0.284, −0.182, …— no son constantes entre sí, pero
eso es de esperar: l=1 tiene tres componentes m_l que el campo separa, así que
no hay un solo umbral por N.)

### 5.5 Características de la curva y criterio de los 20 GHz

Curva de carácter, excluyendo la ventana. Extremos locales por segmento
contiguo:

| # | segmento | R_min [a₀] | E−E_man [GHz] | prof. vs máx. derecho [GHz] | ΔR |
|---|---|---|---|---|---|
| 0 | der (R > 613.7) | 1015.00 | −34.7675 | 2.5561 | — |
| 1 | der (R > 613.7) | 1115.00 | −32.2370 | — | 100.00 |

**En el segmento izquierdo (R < 556.7 a₀) no hay ningún mínimo local estricto**:
la curva sube monótonamente de −110.6 GHz (R=110) a un **máximo en R ≈ 405 a₀**
(−11.0 GHz) y a partir de ahí baja monótonamente hasta la ventana de resonancia
(−15.3 GHz en R=555). Es una forma distinta de la de n=24, que sí tenía mínimos
en 235, 400 y 455 a₀.

| criterio | \|E−E_man\| > 20 GHz |
|---|---|
| sin excluir nada | 70.51 % |
| excluida la ventana | 71.84 % |
| excluidas ventana + cola butterfly | **68.98 %** |

Los 58 puntos por debajo de 20 GHz están confinados a **R ∈ [270, 555] a₀**.
Es el mismo patrón que en n=24 (73 % / R ∈ [275, 535] a₀): la afirmación del
paper «shifted more than 20 GHz for R ≲ 1200 a₀» **se sostiene**, y como allí,
la conclusión es robusta frente a la exclusión (68.98–71.84 %, 2.9 puntos de
margen). Nuestro dominio llega a 1185.5 a₀, así que cubre casi todo el rango que
enuncia el texto.

## 6. Figura

`plots/bop_curve_MJ0_n25_excluded.png`, mismo estilo que
`bop_curve_MJ0_excluded.png`, con el cero de energía del paper
(E_{n=25,l≥3} + KRb(N=0)):

- **Panel superior**: rango completo, se muestra la divergencia de −6459 GHz.
- **Panel inferior**: ventana de energía comparable a la Fig. 2.
- Banda **gris** = ventana excluida R ∈ [556.7, 613.7] a₀, rotulada con el
  motivo; banda **rayada** = cola butterfly hasta 705 a₀; línea de puntos = polo
  de A_p en 585.2 a₀.
- **Rojo grueso**: curva de carácter manifold. **Azul discontinuo**: curva de
  índice fijo k=53, dibujada a propósito para que se vea que en R ≈ 270–1040 a₀
  **no es la misma curva**.
- Líneas horizontales: naranja = umbrales 28s+KRb(N), azul = manifold+KRb(N).
- **Línea vertical violeta en 2n² = 1250 a₀** con la franja sombreada que va
  desde el borde del dominio (1185.5 a₀) hasta ahí: es la parte del retorno
  clásico que **no** alcanzamos con el remapeo actual.

## 7. Pozo «uno antes del más externo»

```
punto de retorno clásico del manifold: 2n² = 1250 a0
dominio accesible del remapeo:         hasta 1185.5 a0  (94.8 % de 2n²)

pozo MÁS EXTERNO:   R = 1115.00 a0   E-E_man = -32.2370 GHz   R/2n² = 0.892
pozo UNO ANTES:     R = 1015.00 a0   E-E_man = -34.7675 GHz   R/2n² = 0.812
separación:         100.00 a0
```

**El dato pedido es R ≈ 1015 a₀ (0.812 · 2n²)**, con el más externo en 1115 a₀
(0.892 · 2n²).

⚠️ **Con una salvedad que no conviene esconder**: nuestro pozo «más externo» lo
es dentro del dominio accesible, que se corta en 1185.5 a₀ — un 5.2 % antes de
2n². Si la Fig. 4 del paper mide el pozo más externo cerca de 2n², puede que su
«más externo» sea uno que a nosotros nos queda fuera, y entonces nuestro «uno
antes» tampoco sería el suyo. La comparación completa de la Fig. 4 (varios n)
queda fuera de esta ronda, como se pidió.

## 8. Pendiente

1. **Rehacer n=24 con la base correcta.** `plots/bop_curve_MJ0.npz`,
   `bop_curve_MJ0.png` y `bop_curve_MJ0_excluded.png` se generaron con la base
   incompleta **y con el criterio de índice fijo**, y el código ya no los
   reproduce. Por §4 y §5.2 se espera que la curva cambie ≲0.25 GHz, pero el
   dominio se acorta a 1090.16 a₀ y hay que rehacer la identificación por
   carácter. No se ha hecho aquí porque el encargo acotaba la curva a n=25.
2. **M_J=1**: sigue sin abordarse, como se pidió.
3. **Extender el dominio de H_mol más allá del remapeo**: sabemos que el tope
   en R es una limitación de cómo leemos las tablas de dispersión, no física.
   No se toca en esta ronda por decisión explícita.
4. La forma radial de los vecinos sigue siendo hidrogenoide de n entero
   (27s→24, 26p→23, 25d→24). Con tres vecinos en vez de uno, esa aproximación
   ahora afecta a tres niveles; el error en la extensión radial es del 1-3 %.

## Archivos tocados

- `src/trimero/systems/rb_defects.py` — `neighbor_levels()`, `n_star_nl()` (nuevos).
- `src/trimero/basis/quantum.py` — `CoupledBasis(manifold_l_max, neighbor_l)`.
- `src/trimero/basis/radial.py` — `neighbor_l` / `neighbor_n_eff` en vez de `n_s_eff`.
- `src/trimero/hamiltonians/charge_dipole.py` — `rydberg_diagonal` usa `neighbor_levels`.
- `src/trimero/hamiltonians/fermi_krb.py` — `n_star_of_l` generalizado a los tres vecinos.
- `src/trimero/simulation/bop_system.py` — `BOPSystem(n_manifold, neighbors)`.
- `scripts/run_basis_correction_check.py` — nuevo: §2, §3, §4 de este documento.
- `scripts/run_bop_manifold.py` — barrido adiabático parametrizado por manifold.
- `scripts/run_validation_paper.py` — pasa a usar `BOPSystem`.
- `tests/basis/test_basis_enumeration.py`, `tests/hamiltonians/test_fermi_krb.py`,
  `tests/simulation/test_bop_system.py` — actualizados.

## Referencias

- Aguilera-Fernández, J., Sadeghpour, H. R., Schmelcher, P. & González-Férez, R.,
  *Ultralong-Range Rb-KRb Rydberg Molecules: Selected Aspects of Electronic
  Structure, Orientation and Alignment*, J. Phys.: Conf. Ser. **635**, 012023
  (2015), arXiv:1507.07972 — **define la base**.
- González-Férez, R., Sadeghpour, H. R. & Schmelcher, P., *New J. Phys.* **17**,
  013021 (2015) — Tabla I, Fig. 2.
- `docs/analysis_verificacion_tabla_I.md`, `docs/analysis_ventana_exclusion_resonancia.md`.
