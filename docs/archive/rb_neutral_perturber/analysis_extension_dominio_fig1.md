# Fig. 1(a): curvas M_J=0 de n=24 y n=25, y por qué NO se extiende el dominio

> ⚠️ **AVISO DE PREMISA (2026-08-20).** Este documento trata el
> **pseudopotencial de Fermi**, que **NO forma parte del Hamiltoniano de
> Aguilera-Fernández et al. 2015 / González-Férez et al. 2015** para Rb*-KRb.
> Su Ec. 1 es `H_ad = H_A + H_mol`, con KRb como **dipolo puntual**, no como
> centro de dispersión de contacto. Su contenido técnico **sigue siendo
> válido** —y sigue aplicando a la línea de perturbador NEUTRO
> (Aguilera-Fernández 2016)— pero **su uso para comparar con esos dos papers
> partía de una premisa equivocada**. Ver
> `docs/analysis_fig1_carga_dipolo_sin_fermi.md` §1.

**Fecha**: 2026-08-20
**Autor**: Javier Aguilera
**Relevancia**: Reproducción de la Fig. 1(a) de Aguilera-Fernández et al. (2015),
que va de 400 a 1800 a₀. Nuestro dominio del remapeo k(R) llega a 1090 a₀ (n=24)
y 1185 a₀ (n=25). Se evaluó la opción de extender con H_a+H_mol y **se
descartó**, con medida.
**Tipo**: analysis

## Resumen

**El paso previo de verificación falla, y la extensión queda descartada.** La
idea era: más allá del dominio del remapeo, poner V_Fermi ≡ 0 y seguir con
H_a+H_mol, asumiendo que ahí el pseudopotencial ya es despreciable. Medido en el
empalme, **no lo es**: quitarlo desplaza la curva **14.28 GHz (n=24) y 11.86 GHz
(n=25)**, un 43 % y un 41 % de su energía de ligadura. La ventana vertical de la
figura que queremos comparar es de decenas de GHz, así que la discontinuidad
sería del orden de la profundidad de los pozos que se comparan.

La razón física es que el borde del dominio **no cae en la cola de la función de
onda sino casi encima de su lóbulo externo** (|ψ(R)|² al 89 % de su máximo
externo), que es justo donde el pseudopotencial de contacto es grande.

Lo que sí se entrega: las dos curvas (n=24 y n=25) con H completo dentro de su
dominio, con la base correcta y el criterio de carácter, y la comparación
cualitativa con las tres tendencias del texto del paper.

## Palabras clave

- extensión del dominio, remapeo semiclásico k(R)
- pseudopotencial de Fermi, onda s
- lóbulo externo de la función de onda Rydberg
- Fig. 1(a), M_J=0, n=24 y n=25

## 1. El paso previo: criterio, medida y veredicto

### 1.1 Criterio, fijado antes de mirar los números

Escrito en `scripts/run_fig1_curves.py` como constantes de módulo, no ajustado a
posteriori:

```python
# Criterio de aceptación de la extensión, fijado ANTES de mirar los números:
# el salto en el empalme al poner V_Fermi = 0 debe ser pequeño frente a la
# escala vertical de la figura que queremos comparar (~decenas de GHz).
JUMP_TOL_GHZ = 1.0
JUMP_TOL_REL = 0.05          # ... y frente a la propia ligadura de la curva
```

Se miden tres cosas en los últimos 150 a₀ del dominio de cada manifold:

| medida | qué dice |
|---|---|
| ‖V_Fermi‖_F | tamaño absoluto del pseudopotencial en el bloque M_J=0 |
| ‖H_mol − diag(B·N(N+1))‖_F | la parte de H_mol que depende de R, para comparar |
| **salto** = E_carácter(con V_Fermi) − E_carácter(sin V_Fermi) | **el error que la extrapolación mete justo en el empalme** |

La tercera es la que decide; las dos primeras dan contexto.

### 1.2 Resultado (n=24, dominio hasta 1090.16 a₀)

| R [a₀] | ‖V_Fermi‖ | ‖H_mol−d‖ | V_F/H_mol | E con Fermi | E sin Fermi | **salto** |
|---|---|---|---|---|---|---|
| 940.2 | 1.400×10⁻⁵ | 1.547×10⁻⁵ | 0.905 | −40.398 | −21.969 | **−18.429** |
| 990.2 | 1.194×10⁻⁵ | 1.429×10⁻⁵ | 0.836 | −37.552 | −20.166 | −17.386 |
| 1040.2 | 9.773×10⁻⁶ | 1.243×10⁻⁵ | 0.786 | −37.310 | −19.544 | −17.766 |
| 1080.2 | 7.552×10⁻⁶ | 1.152×10⁻⁵ | 0.656 | −34.734 | −19.470 | −15.264 |
| **1089.7** | **6.984×10⁻⁶** | **1.126×10⁻⁵** | **0.621** | **−33.536** | **−19.261** | **−14.275** |

### 1.3 Resultado (n=25, dominio hasta 1185.50 a₀)

| R [a₀] | ‖V_Fermi‖ | ‖H_mol−d‖ | V_F/H_mol | E con Fermi | E sin Fermi | **salto** |
|---|---|---|---|---|---|---|
| 1035.5 | 1.095×10⁻⁵ | 1.353×10⁻⁵ | 0.810 | −34.369 | −18.987 | −15.381 |
| 1085.5 | 9.378×10⁻⁶ | 1.185×10⁻⁵ | 0.791 | −32.314 | −17.525 | −14.790 |
| 1135.5 | 8.105×10⁻⁶ | 1.096×10⁻⁵ | 0.739 | −32.073 | −17.129 | −14.944 |
| 1175.5 | 6.253×10⁻⁶ | 1.020×10⁻⁵ | 0.613 | −29.704 | −17.030 | −12.675 |
| **1185.0** | **5.789×10⁻⁶** | **9.979×10⁻⁶** | **0.580** | **−28.701** | **−16.841** | **−11.860** |

### 1.4 Veredicto

```
criterio: |salto| < 1.0 GHz  Y  < 5 % de |E|

n=24:  |salto| = 14.275 GHz = 42.6 % de |E| = 33.536 GHz   ->  NO PASA
n=25:  |salto| = 11.860 GHz = 41.3 % de |E| = 28.701 GHz   ->  NO PASA
```

**Hay que ser preciso sobre qué falla y qué no.** La condición literal que se
pidió comprobar —«tendencia decreciente clara»— **sí se cumple**: ‖V_Fermi‖ cae
un factor 2 en los últimos 150 a₀ y su cociente con H_mol baja de 0.90 a 0.62.
Lo que falla es la **magnitud**: una cosa puede estar decreciendo y seguir siendo
grande. En el empalme V_Fermi sigue valiendo el 58-62 % de la parte de H_mol que
depende de R, y el 41-43 % de la ligadura de la curva.

## 2. Por qué falla: el borde no está en la cola, está en el lóbulo externo

La intuición «a R grande el electrón ya no llega, así que el contacto es
despreciable» es correcta **para R ≫ 2n²**, pero el borde de nuestro dominio no
está ahí. El dominio lo corta el punto de retorno clásico del par más ligado de
la base, el (n+2)p, que está a un 5 % de 2n² del manifold:

| n | R_max del remapeo | 2n² del manifold | R_max / 2n² |
|---|---|---|---|
| 24 | 1090.16 a₀ | 1152 a₀ | 0.946 |
| 25 | 1185.50 a₀ | 1250 a₀ | 0.948 |

Y ahí la densidad electrónica está prácticamente en su máximo externo. Medido
como |ψ(R)|² del manifold normalizado al máximo del lóbulo más externo:

| n | R [a₀] | \|ψ(R)\|² / máx. lóbulo externo |
|---|---|---|
| 24 | 940 | 0.021 |
| 24 | 990 | 0.351 |
| 24 | 1040 | **0.943** |
| 24 | 1089.7 (borde) | **0.894** |
| 25 | 1035 | 0.000 |
| 25 | 1085 | 0.435 |
| 25 | 1135 | **0.960** |
| 25 | 1185.0 (borde) | **0.890** |

El borde del dominio cae **sobre el lóbulo externo**, no después de él. El
pseudopotencial de contacto es proporcional a esa densidad, así que ahí está
cerca de su máximo en la región externa — exactamente lo contrario de
despreciable.

### 2.1 Y lo que se estaría tirando es la onda s, que es la que liga

Descomponiendo el pseudopotencial en sus dos términos (n=24):

| R [a₀] | E completo | E sin Fermi | E sólo onda s | E sólo onda p |
|---|---|---|---|---|
| 940.2 | −40.398 | −21.969 | −35.504 | −30.981 |
| 1040.2 | −37.310 | −19.544 | −34.250 | −23.137 |
| 1089.7 | −33.536 | −19.261 | **−29.607** | −23.133 |

y para n=25:

| R [a₀] | E completo | E sin Fermi | E sólo onda s | E sólo onda p |
|---|---|---|---|---|
| 1135.5 | −32.073 | −17.129 | −29.661 | −19.872 |
| 1185.0 | −28.701 | −16.841 | **−25.653** | −19.838 |

De los ~14 GHz que se pierden en n=24, **~10.3 GHz los aporta la onda s sola**.
No es una corrección fina que se pueda descartar: es el término de contacto que
crea el pozo de la molécula Rydberg de rango ultralargo. Poner V_Fermi ≡ 0 no
deja «la misma curva un poco desplazada», deja una curva de naturaleza distinta
—sólo carga-dipolo— que no es comparable con la del paper.

## 3. Decisión y qué se entrega en su lugar

**No se extiende.** La región R > R_max se calcula igualmente con V_Fermi ≡ 0 y
se **dibuja punteada y rotulada como RECHAZADA**, con el salto marcado
explícitamente en la figura (línea vertical + triángulo), porque es la evidencia
visual de por qué la opción no vale. No se usa para ninguna conclusión
cuantitativa, y no se presenta como continuación de la curva.

La implementación del interruptor sí se queda, porque es lo que permite medir:

```python
    def hamiltonian(self, R, M_J=0, fermi=True):
        """
        `fermi=False` pone V_Fermi ≡ 0 EXACTAMENTE (matriz nula, no un valor
        devuelto por error ni una extrapolación de las tablas).
        """
        blk = self.block(M_J)
        H = np.diag(self.rydberg_diagonal(M_J)) + self.hmol.build(blk, R)
        if fermi:
            H = H + self.fermi.build(blk, R)
        return H
```

### 3.1 Vías alternativas para llegar a 1800 a₀, sin implementar

Se listan porque la pregunta sigue abierta, no porque se hayan probado:

1. **Regenerar las tablas de dispersión para n=24/25** en vez de remapear las de
   n=35 (la «opción (ii)» de `analysis_pseudopotencial_fermi_krb.md`). Elimina
   el remapeo y con él su límite; es la vía limpia y la más cara.
2. **Empalmar con la forma asintótica de la onda s** en vez de con cero:
   A_s(k) → a_s constante cuando k→0, de modo que
   V_s ∝ 2π a_s |ψ(R)|², que sí se puede evaluar más allá del punto de retorno
   con la función de onda hidrogenoide. Requiere justificar que k es
   suficientemente pequeño ahí y decidir qué hacer con la onda p.
3. **Aceptar el recorte** y comparar con la Fig. 1(a) sólo en R ≲ 1185 a₀,
   diciendo explícitamente qué parte de la figura no se cubre.

La 2 parece la más razonable en relación coste/beneficio, pero no se ha
implementado ni verificado: queda como propuesta.

## 4. Curvas de n=24 y n=25 con H completo

Ambas con la base correcta (manifold + (n+1)d + (n+2)p + (n+3)s) y el criterio
de **carácter** (la más baja con peso de manifold > 50 % en cada R), no el
índice fijo. Esto cierra de paso el punto que quedaba pendiente de la ronda
anterior: **n=24 rehecho con la base completa**.

| | n=24 | n=25 |
|---|---|---|
| dim(M_J=0) | 1064 | 1113 |
| dominio del remapeo | [400, **1090.16**] a₀ | [400, **1185.50**] a₀ |
| 2n² | 1152 a₀ | 1250 a₀ |
| ventana p excluida | [536.4, 589.1] a₀ | [556.7, 613.7] a₀ |
| cola butterfly (no comparable) | hasta 685.0 a₀ | hasta 705.0 a₀ |
| puntos con V_Fermi | 140 | 159 |
| puntos de la extensión rechazada | 71 | 62 |

`plots/fig1_MJ0_n24_n25.npz` guarda para cada n: `R_in/E_in/K_in/W_in` (región
con V_Fermi) y `R_ex/E_ex/K_ex/W_ex` (región sin él), separados a propósito para
que no se puedan concatenar por descuido.

### 4.1 La figura

`plots/fig1_MJ0_n24_n25.png`, dos paneles sobre R ∈ [400, 1800] a₀:

- **Trazo grueso continuo**: curva con H completo, dentro del dominio del
  remapeo. Rojo n=24, azul n=25.
- **Punteado fino**: extensión con V_Fermi ≡ 0, **rotulada RECHAZADA** en la
  leyenda. Está dibujada para que se vea el problema, no como continuación.
- **Línea vertical continua + triángulo** en el borde del dominio: el salto de
  14.3 GHz (n=24) y 11.9 GHz (n=25). Es lo que hay que mirar.
- **Banda llena** de color: ventana excluida de la resonancia p.
  **Banda rayada**: cola butterfly, donde los mínimos no son comparables.
- **Punteado vertical**: fin del dominio del remapeo. **Raya-punto**: 2n²a₀.

Nada está recortado en silencio: las tres regiones —válida, cola, extendida—
están marcadas y etiquetadas.

## 5. Comparación con las tendencias del texto del paper

Pozos de cada curva, excluyendo ventana de resonancia **y** cola butterfly:

**n=24** (2n² = 1152 a₀)

| # | R_pozo [a₀] | E−E_man [GHz] | profundidad [GHz] | R/2n² |
|---|---|---|---|---|
| 0 | 690.00 | −119.3661 | 78.7735 | 0.599 |
| 1 | 925.00 | −40.7704 | 3.5099 | 0.803 |
| 2 | 1030.00 | −37.3691 | no medible | 0.894 |

**n=25** (2n² = 1250 a₀)

| # | R_pozo [a₀] | E−E_man [GHz] | profundidad [GHz] | R/2n² |
|---|---|---|---|---|
| 0 | 710.00 | −115.2670 | 80.8392 | 0.568 |
| 1 | 1015.00 | −34.7675 | 2.5561 | 0.812 |
| 2 | 1115.00 | −32.2370 | no medible | 0.892 |

«No medible» = el pozo más externo no tiene un máximo local a su derecha dentro
del dominio, así que su profundidad no está definida con estos datos.

⚠️ **Cuántos pozos hay de verdad para comparar: dos por curva, y sólo uno con
profundidad medible.** El pozo #0 (690 / 710 a₀) está pegado al borde del
criterio de cola butterfly, así que su profundidad de ~80 GHz sigue
contaminada por la resonancia y no es un pozo trilobite limpio. Con eso, la
comparación descansa sobre muy poca base estadística, y así hay que leerla.

### (1) ¿Los pozos se desplazan a mayor R al aumentar n? — **SÍ**

| pozo (desde fuera) | n=24 | n=25 | Δ |
|---|---|---|---|
| más externo | 1030.0 a₀ | 1115.0 a₀ | **+85.0 a₀** |
| segundo | 925.0 a₀ | 1015.0 a₀ | **+90.0 a₀** |
| tercero (contaminado) | 690.0 a₀ | 710.0 a₀ | +20.0 a₀ |

**Consistente, y con un detalle que refuerza el resultado**: en unidades de
2n²a₀ las posiciones son casi idénticas —0.894 vs 0.892 el más externo, 0.803 vs
0.812 el segundo—, es decir los pozos **escalan con n²**, que es exactamente lo
que se espera de una estructura fijada por la función de onda Rydberg. El
desplazamiento absoluto de +85/+90 a₀ es del orden de
2·(25²−24²) = 98 a₀, como debe ser.

### (2) ¿La profundidad decrece con n? — **SÍ, pero con un solo pozo medido**

| pozo | n=24 | n=25 | |
|---|---|---|---|
| segundo desde fuera | 3.5099 GHz | 2.5561 GHz | **decrece** |
| más externo | no medible | no medible | — |
| tercero (contaminado) | 78.7735 GHz | 80.8392 GHz | crece (no fiable) |

La única comparación limpia da **2.556 / 3.510 = 0.728**, frente al 
(24/25)⁶ = 0.783 que daría un escalado n⁻⁶ (el esperado si la profundidad va
como |ψ(R)|² ∝ n⁻³ al cuadrado, evaluado en el pozo correspondiente). Van en la
misma dirección y en el mismo orden, pero **con un solo pozo, malla de 5 a₀ y
profundidades de pocos GHz, esto no da para afirmar el exponente** — sólo para
decir que el signo de la tendencia es el correcto.

El pozo #0 va al revés (78.8 → 80.8 GHz), pero es el contaminado por la cola de
la resonancia: no cuenta como contraejemplo ni como confirmación.

### (3) ¿La amplitud de la oscilación decrece con R? — **SÍ**

| n | profundidades de dentro a fuera [GHz] | |
|---|---|---|
| 24 | 78.773 → 3.510 → (no medible) | decreciente |
| 25 | 80.839 → 2.556 → (no medible) | decreciente |

Las dos curvas caen más de un orden de magnitud en amplitud entre el pozo
interno y el siguiente, y la estructura se aplana hacia el borde del dominio.
Visualmente (panel superior de la figura) es lo más claro de las tres
tendencias: entre 800 y 1090 a₀ las dos curvas suben suavemente con ondulaciones
de pocos GHz sobre un fondo de ~−40 a −32 GHz.

### Veredicto conjunto

**Las tres tendencias que describe el texto son consistentes con lo que dan
nuestras dos curvas.** Con las salvedades, que son grandes y conviene repetirlas:

1. Sólo son **dos manifolds y dos pozos limpios por curva**; el texto describe
   una familia de curvas sobre varios n.
2. **La comparación no cubre la mitad exterior de la Fig. 1(a)**: nuestro
   dominio acaba en 1090 / 1185 a₀ y la figura llega a 1800 a₀. La región
   R > R_max de nuestra figura es la extensión rechazada, no un resultado.
3. La tendencia (2) se apoya en **un único par de profundidades**.

Lo que sí se puede afirmar sin reservas es lo del escalado n²: los pozos caen
en la misma posición relativa a 2n²a₀ en los dos manifolds, con tres cifras
significativas de acuerdo en el más externo.

## 6. Archivos tocados

- `src/trimero/simulation/bop_system.py` — `hamiltonian/eigvals/solve(..., fermi=)`
  con V_Fermi ≡ 0 exacto, y `character_curve()` (nuevos).
- `scripts/run_fig1_curves.py` — nuevo: comprobación previa, barridos, tendencias
  y figura.
- `tests/simulation/test_bop_system.py` — `test_s5_fermi_off_es_cero_exacto`
  (nuevo): comprueba que `fermi=False` resta exactamente V_Fermi y que el salto
  en el borde NO es despreciable.
- `plots/fig1_MJ0_n24.npz`, `plots/fig1_MJ0_n25.npz`,
  `plots/fig1_MJ0_n24_n25.png`.

`poetry run pytest -m "not slow"` → 44 passed.

## 7. Pendiente

1. **Llegar a 1800 a₀ de verdad**: elegir entre las tres vías de §3.1. Sin eso,
   la comparación con la Fig. 1(a) se queda en su mitad interior.
2. **Más manifolds** (n=26, 27, …) para que las tendencias de §5 tengan base
   estadística, y para la Fig. 4.
3. **M_J=1**: sigue sin abordarse.
4. Las curvas y figuras de n=24 de rondas anteriores
   (`plots/bop_curve_MJ0*.png`) siguen siendo de la base incompleta y del
   criterio de índice fijo. La curva buena de n=24 es ahora la de
   `plots/fig1_MJ0_n24.npz`; las viejas quedan obsoletas.

## Referencias

- Aguilera-Fernández, J., Sadeghpour, H. R., Schmelcher, P. & González-Férez, R.,
  *J. Phys.: Conf. Ser.* **635**, 012023 (2015), arXiv:1507.07972 — Fig. 1(a).
- González-Férez, R., Sadeghpour, H. R. & Schmelcher, P., *New J. Phys.* **17**,
  013021 (2015) — Fig. 2 y el texto de las tendencias.
- `docs/analysis_base_correcta_3_vecinos.md` — base correcta y criterio de carácter.
- `docs/analysis_ventana_exclusion_resonancia.md` — ventana de la resonancia p.
- `docs/analysis_pseudopotencial_fermi_krb.md` — el remapeo k(R) y sus opciones.
