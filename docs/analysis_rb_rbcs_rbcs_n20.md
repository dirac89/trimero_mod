# Rb*+RbCs+RbCs en el manifold n=20

## Modelo

El sistema contiene dos rotores rígidos RbCs etiquetados, con
`B=490.17 MHz` y `d=1.225 D`, acoplados al campo del ion y del electrón
Rydberg, entre sí mediante la interacción dipolo–dipolo y a un campo DC
paralelo al eje global `+Z`. La base es

`(l, m_l, N1, M1, N2, M2)`, con `M_J=m_l+M1+M2`.

No se incluyen potenciales químicos de corto alcance, dispersión RbCs–RbCs,
movimiento nuclear ni simetrización bosónica.

## Geometrías

- `symmetric`: `z1=-R`, `z2=+R`, de modo que `R12=2R` y `Vdd` cambia durante
  el barrido.
- `unilateral`: `z1=R`, `z2=R+D`; `D=300 a0` por defecto, configurable, y
  `Vdd` es constante a lo largo del corte.

Por ello la comparación no aísla únicamente un efecto angular: también cambia
la distancia intermolecular. Este hecho debe conservarse en títulos y
metadatos de cualquier figura.

## Producción prevista

```bash
poetry run python scripts/compute_double_rbcs_curves.py \
  --geometry both --separation 300 --n-manifold 20 --n-max 3 \
  --mj 0 1 --fields 0 100 300 500 \
  --rmin 200 --rmax 1800 --step 10
```

La base tiene dimensiones aproximadas 4664/4620 para `N_max=3` y
11084/10999 para `N_max=4`, en `M_J=0/1`. Se usa `eigsh` con desplazamiento,
selección por carácter del manifold y seguimiento desde `R_max` hacia dentro
por máximo solapamiento.

Antes de considerar definitivas las curvas debe ejecutarse la convergencia:

```bash
poetry run python scripts/check_double_rbcs_convergence.py \
  --n-manifold 20 --n-max 3 4 --geometry symmetric \
  --r 400 800 1400 --mj 0 1 --fields 0 500
```

El criterio es `|ΔE|<0.01 GHz` y diferencias inferiores a `0.01` en ambas
orientaciones. Los puntos que fallen deben repetirse con `N_max=5`.

## Benchmark y convergencia inicial

En `R=800 a0`, geometría simétrica, `M_J=0`, `F=0`, la construcción y
diagonalización selectiva tardan aproximadamente 9 s para `N_max=3`, 25 s
para `N_max=4` y 97 s para `N_max=5` en la máquina de desarrollo. Los
resultados son:

| N_max | dimensión | E (GHz) | cos(theta1) | cos(theta2) |
|---:|---:|---:|---:|---:|
| 3 | 4664 | -64.606344 | +0.094729 | -0.094729 |
| 4 | 11084 | -66.060188 | +0.112822 | -0.112822 |
| 5 | 22364 | -66.531816 | +0.122274 | -0.122274 |
| 6 | 40292 | -66.650543 | +0.125940 | -0.125940 |

El paso `3→4` cambia `1.453844 GHz`, `4→5` cambia `0.471628 GHz` y `5→6`
cambia `0.118727 GHz`. Por tanto este punto **no está convergido en energía**
con `N_max=6`. La orientación sí cumple el umbral en `5→6`, con un cambio de
`0.003667`.

Las razones sucesivas de los incrementos energéticos son `0.324` y `0.252`.
Una extrapolación de Aitken con los tres últimos puntos estima
`E(N_max→∞)≈-66.6905 GHz`, pero esta extrapolación es sólo un diagnóstico y no
sustituye una diagonalización convergida. Las curvas con `N_max=3` son de
exploración y no deben presentarse como resultados definitivos.

### Límite computacional observado

`N_max=7` produce 66800 estados para `M_J=0`. Se probaron dos algoritmos:

1. `eigsh` shift-invert: la factorización LU de `H-sigma I` no terminó en un
   tiempo razonable.
2. Búsqueda matricial libre sobre `[(H-E_n)GHz-sigma]^2`: evita la LU y fue
   validada contra diagonalización exacta en una base pequeña, pero tampoco
   convergió con rapidez suficiente en el espectro denso de `N_max=7`.

La producción completa queda condicionada a introducir una base rotacional
contraída (estados pendulares locales) o una reducción de simetría. Aumentar
directamente a `N_max=8` daría 103908 estados y no es una ruta de producción
viable con el solver actual.

## Evaluación de la base pendular contraída

Se implementó una base local que diagonaliza cada rotor en el campo axial del
ion más el campo externo y conserva un número configurable de autoestados por
`M_i`. Si se conservan todos, el espectro coincide con el Hamiltoniano
primitivo dentro de `2e-12 Eh` en las pruebas pequeñas.

Para `n=20`, corte primitivo `N=8`, geometría simétrica, `R=800 a0`, `M_J=0`
y `F=0`:

| modos por M | dimensión | E (GHz) | cos(theta1) | tiempo |
|---:|---:|---:|---:|---:|
| 1 | 4148 | -8.600613 | -0.353842 | 9.7 s |
| 2 | 15008 | -29.055632 | -0.021976 | 25.9 s |
| 3 | 30244 | -46.160829 | -0.080725 | 113.8 s |
| 4 | 47624 | -60.432113 | +0.053376 | 437.0 s |

La energía se acerca monótonamente a la referencia primitiva, pero cuatro
modos todavía quedan a más de 6 GHz y ya cuestan más que `N_max=6`. El campo
electrónico transversal mezcla demasiados estados pendulares para que esta
contracción local sea una ruta de producción suficiente por sí sola.

El benchmark es reproducible con:

```bash
poetry run python scripts/check_double_rbcs_contracted.py \
  --n-manifold 20 --primitive-n-max 8 --rotor-keep 1 2 3 4 \
  --geometry symmetric --r 800 --mj 0 --field 0
```

### Contracción mediante estados naturales

Se implementó la alternativa basada en las matrices de densidad reducida de
cada rotor. Usando como fuente el autoestado primitivo `N_max=6` y proyectando
una base objetivo con corte `N_max=8`, en el mismo punto de control:

| estados naturales por M | dimensión | E (GHz) | cos(theta1) | peso manifold |
|---:|---:|---:|---:|---:|
| 1 | 4148 | -51.734279 | +0.831425 | -- |
| 2 | 15008 | -66.618616 | +0.125855 | 0.971402 |
| 3 | 30244 | -66.644196 | +0.126791 | 0.971117 |
| 4 | 47624 | -66.650354 | +0.125704 | 0.971060 |

Entre tres y cuatro estados naturales, `dE=0.006158 GHz` y
`dcos=0.001087`; ambos satisfacen el criterio de `0.01`. El resultado de cuatro
estados difiere de la fuente primitiva `N_max=6` en sólo `0.000189 GHz`.
La simetría esperada también se conserva:
`cos(theta2)=-cos(theta1)` dentro de la precisión numérica.

Este resultado demuestra convergencia **de la proyección respecto al número de
estados naturales en este punto**, pero no demuestra todavía convergencia
absoluta respecto a `N_max`: la base natural se entrenó con una fuente
`N_max=6`, cuya energía aún cambia apreciablemente frente a la extrapolación.
Tampoco valida por sí solo la transferencia de la base a otros radios, campos,
geometrías o valores de `M_J`.

El benchmark se reproduce (reutilizando opcionalmente una fuente almacenada)
con:

```bash
poetry run python scripts/check_double_rbcs_natural.py \
  --source-n-max 6 --target-n-max 8 --natural-keep 2 3 4 \
  --geometry symmetric --r 800 --mj 0 --field 0
```

El siguiente paso numérico es validar bases naturales entrenadas en radios
ancla y transferidas a puntos vecinos. Sólo después de establecer el espaciado
necesario de esos anclajes debe iniciarse el barrido completo de producción.

## Barrido exploratorio completo

Se completó el barrido primitivo con `N_max=3` para las dos geometrías,
`M_J=0,1`, campos `0,100,300,500 V/m` y 161 radios entre `200` y `1800 a0`.
Los 16 archivos NPZ están en `plots/rb_rbcs_rbcs_polar/data/`; no contienen
`NaN`, todas las orientaciones pertenecen a `[-1,1]` y todos los estados
seleccionados tienen peso de manifold superior a `0.5`.

Los puntos cuyo solapamiento con el radio anterior es inferior a `0.7` se
conservan en `warning_R` y aparecen como cruces rojas en las figuras. El número
de avisos por campo es:

| geometría | M_J | F=0 | F=100 | F=300 | F=500 |
|---|---:|---:|---:|---:|---:|
| simétrica | 0 | 9 | 10 | 9 | 19 |
| simétrica | 1 | 23 | 24 | 27 | 17 |
| unilateral | 0 | 22 | 54 | 42 | 49 |
| unilateral | 1 | 31 | 43 | 32 | 19 |

La geometría unilateral presenta claramente más regiones de mezcla y
seguimiento ambiguo. Estos datos sirven para localizar pozos, anticruces y
regiones de orientación, pero continúan etiquetados como **exploratorios** por
dos motivos independientes: `N_max=3` no está convergido en energía y un
número apreciable de radios no supera el umbral de continuidad. Las cruces
rojas no deben interpretarse como una BOP adiabática inequívoca.

Las figuras resultantes son:

- `plots/rb_rbcs_rbcs_polar/figures/comparison_symmetric_n20_Nmax3.png`
- `plots/rb_rbcs_rbcs_polar/figures/comparison_unilateral_n20_Nmax3.png`

La prueba de transferencia de una base natural fija desde `R=800 a0` hasta
`R=700 a0` falló cuantitativamente: produjo `-76.0861 GHz`, mientras que la
referencia primitiva `N_max=6` es `-82.1401 GHz`, y además cambió la
orientación. Por tanto no es correcto interpolar anclas separadas `200 a0`.
La ruta de alta precisión requiere actualizar adaptativamente las matrices de
densidad reducida o introducir una contracción conjunta de los dos rotores.

## Continuación adaptativa en cruces densos

El seguimiento de `compute_double_rbcs_curves.py` incorpora sub-stepping
adaptativo para los puntos que siguen por debajo de `--overlap` después de
agotar la duplicación de `k`. Primero se conserva sin cambios el intento
normal con `k, 2k, ... --max-k`. Si todavía falla, el código biseca el hueco
desde el último radio con continuidad aceptable, utiliza los autoestados de
los puntos medios como semillas y vuelve a intentar el radio del grid.

Los puntos auxiliares no se escriben en los NPZ: `R`, energías, orientaciones
y espectros mantienen exactamente las 161 posiciones originales. La bisección
se detiene con `--min-substep` (por defecto `0.5 a0`) o
`--bisect-max-depth` (por defecto 6). Si no consigue recuperar el umbral, se
mantiene el comportamiento anterior: el punto se acepta y queda registrado en
`warning_R`.

Los datasets de dos rotores que contienen puntos de solapamiento bajo deben
regenerarse **sin `--reuse`**: geometrías `symmetric` y `unilateral`,
`M_J=0,1` y los cuatro campos `0,100,300,500 V/m`. Por ejemplo:

```bash
poetry run python scripts/compute_double_rbcs_curves.py \
  --geometry both --n-manifold 20 --n-max 3 --mj 0 1 \
  --fields 0 100 300 500 --rmin 200 --rmax 1800 --step 10 \
  --min-substep 0.5 --bisect-max-depth 6 --workers 4
```

No debe añadirse `--reuse` a ese comando, porque los NPZ existentes fueron
calculados antes del sub-stepping. Los restantes datasets del proyecto no
cambian y pueden seguir conservándose o reutilizándose con `--reuse`.

## Regeneración con sub-stepping adaptativo (25-ago-2026)

Se regeneraron los 16 datasets que tenían al menos un punto con
`overlap < 0.7`. Se usaron `--bisect-max-depth 3` y `--min-substep 2.0` para
limitar el coste recursivo. El canario `symmetric, M_J=1, F=300 V/m` se
ejecutó con los defaults `k=80`, `max-k=320`; tardó `1919.5 s` y redujo los
avisos de 27 a 9. El resto se ejecutó con la configuración original
`k=40`, `max-k=160`, generalmente en parejas de procesos.

| geometría | M_J | F (V/m) | bajos antes | bajos después | tiempo (s) |
|---|---:|---:|---:|---:|---:|
| simétrica | 0 | 0 | 9 | 10 | 1471.4 |
| simétrica | 0 | 100 | 10 | 11 | 1563.8 |
| simétrica | 0 | 300 | 9 | 11 | 1885.3 |
| simétrica | 0 | 500 | 19 | 13 | 2040.5 |
| simétrica | 1 | 0 | 23 | 14 | 1385.2 |
| simétrica | 1 | 100 | 24 | 21 | 1675.1 |
| simétrica | 1 | 300 | 27 | 9 | 1919.5 |
| simétrica | 1 | 500 | 17 | 9 | 981.5 |
| unilateral | 0 | 0 | 22 | 15 | 1469.4 |
| unilateral | 0 | 100 | 54 | 21 | 1851.8 |
| unilateral | 0 | 300 | 42 | 38 | 2888.5 |
| unilateral | 0 | 500 | 49 | 40 | 3061.0 |
| unilateral | 1 | 0 | 31 | 15 | 1396.2 |
| unilateral | 1 | 100 | 43 | 30 | 1965.6 |
| unilateral | 1 | 300 | 32 | 25 | 3449.8 |
| unilateral | 1 | 500 | 19 | 17 | 3020.7 |

El total pasa de 430 a 299 puntos bajos: 131 puntos menos (`30.5 %`). La suma
de `elapsed_seconds` de los 16 NPZ es `32025.1 s` (`8.90 h`) de cómputo por
caso. Esta cifra suma procesos que se ejecutaron en paralelo y no equivale al
tiempo de pared. Hubo además un intento de cuatro procesos abortado tras
`42.7 min` de pared por contención; no produjo datasets y no está incluido en
la suma anterior.

El sub-stepping es especialmente efectivo en los casos con discontinuidades
más densas, pero no garantiza una reducción monótona: tres casos simétricos
con `M_J=0` aumentaron ligeramente su número de avisos. La geometría
unilateral con `F=300,500 V/m` conserva muchos puntos ambiguos incluso tras la
regeneración conservadora. Las figuras se regeneraron desde los nuevos NPZ
usando `--reuse`; `N_max=3` continúa siendo exploratorio y no convergido.

Se hizo una segunda pasada sobre los seis casos menos concluyentes usando los
defaults `k=80`, `max-k=320`, y manteniendo `--bisect-max-depth 3` y
`--min-substep 2.0`. Los valores «antes» de la tabla siguiente son los
resultados de la primera pasada con `k=40`, `max-k=160`; no son los valores de
los datasets originales previos al sub-stepping.

| geometría | M_J | F (V/m) | bajos antes | bajos después | tiempo (s) |
|---|---:|---:|---:|---:|---:|
| simétrica | 0 | 0 | 10 | 10 | 1728.4 |
| simétrica | 0 | 100 | 11 | 11 | 2721.1 |
| simétrica | 0 | 300 | 11 | 11 | 2618.1 |
| unilateral | 0 | 300 | 38 | 32 | 4972.6 |
| unilateral | 0 | 500 | 40 | 40 | 5734.5 |
| unilateral | 1 | 300 | 25 | 25 | 3175.1 |

En conjunto, esta segunda pasada reduce los puntos bajos de 135 a 129 y suma
`20949.8 s` (`5.82 h`) de cómputo por caso. Sólo `unilateral, M_J=0,
F=300 V/m` mejora, de 38 a 32 puntos. Por tanto, los valores reducidos
`k=40`, `max-k=160` **no eran la causa general** de la falta de mejora: los
tres casos simétricos y dos de los tres unilaterales permanecen sin cambios al
duplicar la ventana espectral. En el caso unilateral restante sí ocultaban una
mejora moderada, a costa de un tiempo apreciablemente mayor. Los otros diez
datasets no se regeneraron.
