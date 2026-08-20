# Fig. 1 de Aguilera-Fernández 2015: H_ad = H_A + H_mol, sin pseudopotencial de Fermi

**Fecha**: 2026-08-20
**Autor**: Javier Aguilera
**Relevancia**: Corrige la premisa de fondo de varias rondas de esta sesión. El
Hamiltoniano de ese paper **no lleva pseudopotencial de Fermi**. Con eso
desaparece toda la maquinaria de limitación de dominio y el rango completo
R ∈ [400, 1800] a₀ se calcula directo. Incluye por fin **M_J = 1**.
**Tipo**: analysis

## Resumen

`H_ad = H_A + H_mol` (Ec. 1 del paper): Rydberg libre más carga-dipolo, con KRb
tratado como **dipolo puntual** que siente el campo del Rydberg. Nada de
dispersión de contacto.

Consecuencia inmediata y muy buena: **no hay remapeo k(R), ni ventana de
exclusión de la resonancia de onda p, ni cola butterfly, ni tope en 2n²a₀**.
H_A y H_mol se evalúan con la función de onda hidrogenoide, que está definida en
todo R. Las curvas salen enteras y continuas de 400 a 1800 a₀, sin extensiones
aproximadas ni recortes.

Resultados para n=25, base completa (25,l≥3 + 26d + 27p + 28s):

| | M_J = 0 | M_J = 1 |
|---|---|---|
| dim(bloque) | 1113 | 1106 |
| pozo más profundo | **−23.100 GHz** (R = 400 a₀) | **−19.015 GHz** (R = 400 a₀) |
| mínimos locales | **8** | **7** |
| separación media entre mínimos | 100.7 a₀ | 97.5 a₀ |
| E(R = 1800 a₀) | **−0.338 GHz** | **−0.276 GHz** |
| cruza −10 GHz subiendo | R = 1295 a₀ | R = 1240 a₀ |

**Los tres puntos de comparación que pediste caen dentro de lo que describes de
la figura real**: profundidad de 20–25 GHz (tenemos 23.1), subida hacia cero en
R ≈ 1200–1400 a₀ (tenemos el tramo empinado entre 1150 y 1400, cruzando −10 GHz
en 1295), y del orden de media docena de oscilaciones antes de esa subida
(tenemos 8 y 7).

⚠️ **No he podido ver la imagen adjunta** —no llegó ningún adjunto—, así que la
comparación de §5 es contra tu descripción en texto, no contra la figura.

## Palabras clave

- H_ad = H_A + H_mol, Ec. 1
- carga-dipolo, dipolo puntual, KRb polar
- M_J = 0 y M_J = 1
- sin pseudopotencial de Fermi

## 1. La premisa que estaba mal, y qué queda invalidado

### 1.1 Qué pasó

Hace varias rondas se integró `hamiltonians/fermi_krb.py`, el pseudopotencial de
Fermi de onda s y p, como tercer término del Hamiltoniano de Rb*-KRb. **Eso no
corresponde a este sistema.** El pseudopotencial de contacto modela un
perturbador **NEUTRO** que dispersa al electrón Rydberg (la línea de
Aguilera-Fernández 2016). Aquí el perturbador es **polar**, y el paper lo trata
como un dipolo puntual en el campo eléctrico del Rydberg: `-d·F_ryd`, que es
justo lo que ya hace `hamiltonians/charge_dipole.py`.

### 1.2 Qué documentos quedan afectados, y en qué grado

**Su contenido técnico sigue siendo válido**: la física del pseudopotencial, la
resonancia de forma p, la interpolación consciente del polo y la ventana de
exclusión están bien hechas y bien medidas. Lo que estaba mal era **aplicarlas a
la comparación con Aguilera-Fernández 2015 / González-Férez 2015**. Siguen
siendo el material de referencia si algún día se trabaja la línea de perturbador
neutro.

| documento | qué sigue valiendo | qué NO aplica a esta comparación |
|---|---|---|
| `analysis_pseudopotencial_fermi_krb.md` | la derivación de V_s/V_p, el remapeo k(R) y sus aproximaciones | **todo su uso como término del Hamiltoniano de Rb*-KRb** |
| `analysis_resonancia_onda_p.md` | la resonancia ³P^o a 24.8 meV es real y está bien localizada en la tabla | el pozo butterfly de −6.5 THz no existe en este modelo |
| `analysis_interpolacion_polo_Ap.md` | la interpolación de 1/A_p a través del polo es correcta y está bien justificada | su efecto sobre «nuestras curvas BOP» |
| `analysis_ventana_exclusion_resonancia.md` | el criterio de ventana (polo ± 2×FWHM) está bien construido y testado | **no hay nada que excluir**: sin V_Fermi no hay resonancia en la curva |
| `analysis_extension_dominio_fig1.md` | la medida de que V_Fermi NO es despreciable en el borde es correcta | **la pregunta entera desaparece**: sin V_Fermi no hay borde de dominio que extender |
| `analysis_base_correcta_3_vecinos.md` | **la base (manifold + (n+1)d + (n+2)p + (n+3)s) y el criterio de CARÁCTER: se usan aquí** | §2.1 (dominio) y §5 (curva con V_Fermi) |
| `analysis_curva_bop_MJ0.md` | la metodología adiabático/diabático | los números, que llevan V_Fermi |

Cada uno de esos ficheros lleva ya un aviso al principio apuntando aquí.

### 1.3 Lo que se conserva de todo aquello

No fue tiempo perdido: tres cosas de esas rondas **se usan tal cual** en este
cálculo.

1. **La base correcta**: manifold (n,l≥3) + (n+1)d + (n+2)p + (n+3)s, con el
   manifold como parámetro (`neighbor_levels`, `BOPSystem`).
2. **El criterio de identificación por CARÁCTER** en vez de índice fijo. Sigue
   siendo necesario **aunque no haya V_Fermi**: los estados de 26d y 27p siguen
   produciendo cruces evitados con las curvas del manifold, y un índice fijo
   cambiaría de objeto por el camino. Medido aquí: el índice de la curva de
   carácter va de k=55 (M_J=0) en R=1800 a valores distintos a R menor.
3. **δ₀(ns) = 3.13180** y toda la verificación de la Tabla I, que es energía
   atómica pura y no dependía de nada de esto.

`fermi_krb.py` y sus tests **no se han tocado**: siguen ahí, verdes, para el otro
problema.

## 2. El cálculo

`scripts/run_fig1_charge_dipole.py`. El Hamiltoniano se arma a mano para que se
vea de un vistazo qué términos hay — **no** se llama a
`BOPSystem.hamiltonian(fermi=False)`, y el script **ni siquiera importa**
`fermi_krb`:

```python
def hamiltonian_ad(sysm, M_J, R):
    """
    H_ad(R) = H_A + H_mol, escrito término a término.

    H_A   : diagonal, energías Rydberg del manifold y de los tres vecinos.
    H_mol : B·N² - d·F_ion(R) - d·F_elec(R), construido por
            ChargeDipoleHamiltonian.

    No hay un tercer término. Esto no es `fermi=False`: es que aquí el
    pseudopotencial no forma parte del modelo.
    """
    blk = sysm.block(M_J)
    return np.diag(sysm.rydberg_diagonal(M_J)) + sysm.hmol.build(blk, R)
```

Parámetros: n=25, δ₀(ns)=3.13180, N_max=6, cero de energía
E(n=25,l≥3) + E_KRb(N=0) = −8.0×10⁻⁴ E_h. Malla R = 400…1800 a₀ en pasos de
5 a₀ (281 puntos), **sin un solo hueco**. Coste: 306.7 s (M_J=0) y 107.6 s
(M_J=1).

Umbrales asintóticos, de rondas anteriores y sin cambios (son energías
atómicas):

```
ΔE(28s) + 30B  (N=5) = -22.6437 GHz
ΔE(28s) + 42B  (N=6) =  -9.2757 GHz
```

## 3. Verificación en el borde superior: la curva tiende a cero

Sin dominio truncado, la curva del manifold debe ir al cero de energía
(manifold + KRb(N=0)) cuando R → ∞. Comprobado:

| M_J | E(1800 a₀) [GHz] | k | peso de manifold | E(1700 a₀) [GHz] | pendiente [GHz/100 a₀] |
|---|---|---|---|---|---|
| 0 | **−0.33761** | 55 | **1.0000** | −0.59177 | +0.25415 |
| 1 | **−0.27628** | 52 | **1.0000** | −0.46239 | +0.18611 |

- Las dos curvas están a **menos de 0.34 GHz de cero** en el borde, viniendo de
  −23 y −19 GHz: han recorrido el 98.5 % del camino.
- El peso de manifold es **1.0000**: en el borde el estado es manifold puro, no
  hay mezcla residual con 26d/27p/28s.
- La pendiente es positiva y pequeña, consistente con una aproximación suave a
  la asíntota. Sin ningún escalón ni artefacto: es el contraste directo con las
  rondas anteriores, donde el borde del dominio metía un salto de 12–14 GHz.

## 4. La figura

`plots/fig1_ad_MJ0_MJ1_n25.png`, dos paneles (M_J=0 y M_J=1), ejes R ∈ [400,
1800] a₀ y V(R) ∈ [−25, +1] GHz:

- **136 curvas finas grises**: el resto del bloque, sin filtrar — todas las que
  entran en la ventana de energía.
- **Trazo grueso negro**: la más baja con carácter de manifold (>50 % de peso).
- **Naranja discontinua**: ΔE_28s + 30B (N=5) = −22.64 GHz.
  **Verde discontinua**: ΔE_28s + 42B (N=6) = −9.28 GHz.
- **Vertical violeta**: 2n² = 1250 a₀, sólo como referencia de escala — ya no es
  ningún límite de nada.

## 5. Números para comparar con la figura real

### 5.1 Profundidad del pozo

| | M_J = 0 | M_J = 1 |
|---|---|---|
| valor más profundo de la curva | **−23.100 GHz** en R = 400 a₀ | **−19.015 GHz** en R = 400 a₀ |
| mínimo local más profundo (interior) | −21.597 GHz en R = 450 a₀ | −17.881 GHz en R = 460 a₀ |
| rango completo de la curva | [−23.100, −0.338] GHz | [−19.015, −0.276] GHz |

**Comparación**: describes ~20–25 GHz en la figura. Nuestro M_J=0 da **23.1
GHz**, dentro de ese intervalo. Salvedad: el valor más profundo está en R=400
a₀, que es el **borde izquierdo del rango**, no un mínimo interior — la curva
sigue bajando hacia R menor. El mínimo local más profundo dentro del rango es
21.6 GHz. Las dos cifras caen en 20–25 GHz, así que la comparación aguanta con
cualquiera de los dos criterios.

También encaja que la curva del manifold **quede por encima del umbral N=5**
(−22.64 GHz) en casi todo el rango y sólo lo roce en R ≲ 420 a₀.

### 5.2 Oscilaciones antes de la subida

**M_J = 0: 8 mínimos locales**, entre R = 450 y 1155 a₀:

| R_mín [a₀] | E [GHz] | máx. siguiente [a₀] | amplitud pico-valle [GHz] | ΔR al mín. anterior |
|---|---|---|---|---|
| 450 | −21.597 | 495 | 1.417 | — |
| 515 | −20.444 | 560 | 1.192 | 65 |
| 590 | −19.655 | 635 | 1.074 | 75 |
| 670 | −19.122 | 725 | 1.029 | 80 |
| 765 | −18.809 | 825 | 1.103 | 95 |
| 875 | −18.735 | 950 | 1.377 | 110 |
| 1015 | −19.292 | 1130 | 2.169 | 140 |
| 1155 | −17.165 | — | — | 140 |

**M_J = 1: 7 mínimos locales**, entre R = 460 y 1045 a₀, con la misma
estructura y amplitudes ~30 % menores (0.98 → 1.00 GHz, mínimo 0.68 en R=685).

Dos cosas medidas que conviene tener aunque no las hayas pedido:

- **La separación entre mínimos crece monótonamente**: 65 → 75 → 80 → 95 → 110
  → 140 → 140 a₀ (media 100.7 a₀). Es el comportamiento esperado de una
  modulación ligada a los nodos de la función de onda Rydberg.
- **La amplitud pico-valle NO decrece monótonamente**: baja de 1.42 a 1.03 GHz
  hasta R≈670 a₀ y luego vuelve a subir hasta 2.17 GHz en el último pozo antes
  de la subida. Si la figura real muestra amplitud decreciente en todo el rango,
  ahí habría una discrepancia que mirar; con la descripción en texto que tengo
  no puedo decidirlo.

### 5.3 Forma de la subida final hacia cero

| nivel cruzado (subiendo) | M_J = 0 | M_J = 1 |
|---|---|---|
| −10 GHz | R = **1295** a₀ | R = **1240** a₀ |
| −5 GHz | R = 1385 a₀ | R = 1325 a₀ |
| −2 GHz | R = 1510 a₀ | R = 1450 a₀ |
| −1 GHz | R = 1615 a₀ | R = 1565 a₀ |

**Comparación**: describes la subida hacia cero «alrededor de R ~ 1200–1400 a₀».
El tramo empinado de nuestra curva va de **R ≈ 1155 a₀** (último mínimo) a
**R ≈ 1400 a₀** (donde ya está por encima de −4 GHz), con el punto medio
—el cruce de −10 GHz, que es la mitad de la caída total desde el último
mínimo— en **1295 a₀**. Cae dentro del intervalo que describes.

La forma es una sigmoide limpia: pendiente máxima cerca de 1250–1300 a₀ y cola
suave hacia cero, sin estructura. Nótese que el tramo empinado **atraviesa
2n² = 1250 a₀**, lo que tiene sentido físico: es donde la molécula deja de estar
dentro de la nube electrónica.

### 5.4 M_J = 1 frente a M_J = 0

M_J=1 es sistemáticamente **menos ligada** (−19.0 vs −23.1 GHz de máxima
profundidad, −4.1 GHz de diferencia), tiene **un mínimo menos** (7 vs 8),
oscilaciones **~30 % más pequeñas** y **sube hacia cero antes** (55 a₀ antes en
el cruce de −10 GHz). Las dos curvas tienen la misma forma general.

## 6. Archivos tocados

- `scripts/run_fig1_charge_dipole.py` — nuevo. No importa `fermi_krb`.
- `plots/fig1_ad_MJ0_n25.npz`, `plots/fig1_ad_MJ1_n25.npz`,
  `plots/fig1_ad_MJ0_MJ1_n25.png`.
- Avisos de premisa añadidos a: `analysis_pseudopotencial_fermi_krb.md`,
  `analysis_resonancia_onda_p.md`, `analysis_interpolacion_polo_Ap.md`,
  `analysis_ventana_exclusion_resonancia.md`,
  `analysis_extension_dominio_fig1.md`, `analysis_base_correcta_3_vecinos.md`,
  `analysis_curva_bop_MJ0.md`.
- **Sin tocar**: `src/trimero/hamiltonians/fermi_krb.py` y sus tests.

`poetry run pytest -m "not slow"` → 44 passed.

## 7. Pendiente

1. **Ver la figura real** y afinar §5: sin ella la comparación es contra una
   descripción en texto.
2. **Otros n** (24, 26, 27…) para la Fig. 4 y para las tendencias con n. Ahora
   es barato: no hay dominio que gestionar y el rango es el mismo para todos.
3. **R < 400 a₀**: la curva sigue bajando en el borde izquierdo; el paper
   empieza en 400 a₀, pero conviene saber dónde está el mínimo real.
4. La amplitud pico-valle no decrece monótonamente (§5.2) — contrastar contra la
   figura cuando esté disponible.

## Referencias

- Aguilera-Fernández, J., Sadeghpour, H. R., Schmelcher, P. & González-Férez, R.,
  *Ultralong-Range Rb-KRb Rydberg Molecules: Selected Aspects of Electronic
  Structure, Orientation and Alignment*, J. Phys.: Conf. Ser. **635**, 012023
  (2015), arXiv:1507.07972 — **Ec. 1 y Fig. 1**.
- González-Férez, R., Sadeghpour, H. R. & Schmelcher, P., *New J. Phys.* **17**,
  013021 (2015) — Tabla I y Fig. 2.
- Marinescu, M., Sadeghpour, H. R. & Dalgarno, A., *Phys. Rev. A* **49**, 982
  (1994) — potencial modelo del que salen los defectos cuánticos de H_A.
- `docs/analysis_base_correcta_3_vecinos.md` — base y criterio de carácter, que
  sí se usan aquí.
- `docs/analysis_verificacion_tabla_I.md` — δ₀(ns) = 3.13180.
