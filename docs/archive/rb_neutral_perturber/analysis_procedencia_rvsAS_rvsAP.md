# Procedencia de los datos `rvsAS.dat`, `rvsAP.dat` y `rvsR*.dat`

**Fecha**: 2026-08-18
**Autor**: Javier Aguilera
**Relevancia**: Estos archivos alimentan directamente la construcción del Hamiltoniano
(`FermiPotentials`). Sin conocer su procedencia y sus unidades no se puede validar
la física ni reproducir los resultados.

**Status**: ⚠️ **PUNTO ABIERTO** — el script generador NO se ha encontrado.
No se adopta ninguna interpretación definitiva hasta resolverlo (contacto con
autores originales / revisión de la tesis en papel).

## Resumen

Se realizó una búsqueda exhaustiva de arqueología en el repositorio (todas las ramas,
todos los commits, todos los blobs, objetos inalcanzables, reflog, stash) y en el disco
local. **No existe ni ha existido nunca en este repositorio ningún script generador**
(Fortran, Python, MATLAB, Mathematica u otro) que produzca estos `.dat`. Los archivos
entraron ya construidos en el commit raíz.

La única documentación sobre el contenido de la columna 2 es **un solo comentario de
una línea** en el C++ original, presente desde el commit raíz.

## 1. Búsqueda del script generador — NEGATIVO

### Alcance de la búsqueda

| Comprobación | Resultado |
|---|---|
| `git branch -a` / `git for-each-ref` | Solo `master` y `migrate-python` (+ sus remotos). No hay más refs. |
| `git log --all` | 10 commits en total (2018-03-05 → 2025-07-24). |
| Listado de **todos** los paths que existieron en algún commit | Ni un `.f`, `.f90`, `.for`, `.m`, `.nb`, ni ningún `.py` generador. |
| Nombres tipo `generate`, `scattering`, `phase`, `quantum_defect`, `bahrim`, `thumm`, `effective_range` | Cero coincidencias en nombres de archivo en todo el historial. |
| `git fsck --unreachable --dangling` | Solo un *stash* WIP del 2026-08-18 (`54eeb37`, `e6dc3bd`) y sus blobs: `poetry.lock`, `pyproject.toml`, `src/trimer.py`, blob vacío. Nada nuevo. |
| `git reflog --all` / `git stash list` / `git notes` | Nada relevante. El repo se clonó de `github.com:dirac89/trimero_mod.git`. |
| Búsqueda en disco (`~/Documents`, `~/Desktop`, `~/Downloads`) de `*.f`, `*.f90`, `*.for`, `*.nb`, `*.m`, `rvs*`, y `grep -r "rvsAS"` | Cero coincidencias fuera del repositorio. |

### Contenido real de los commits

El commit raíz `ea14a75` ("trimero modificado", 2018-03-05) **ya contiene** los
`.dat`, junto con código C++ (GSL) y artefactos de build de Eclipse (`.cproject`,
`Debug/`, `.o`). El repositorio nace como un *snapshot* de un proyecto Eclipse ya
en marcha, no desde cero.

Nota: en `ea14a75` y `4912a4c` los archivos vivían en `src/rvsAS.dat`, `src/rvsAP.dat`
y `src/Wavefunction/rvsR*.dat`; se movieron a `data/Wavefunction/` en `117a7ce`.

### Los datos nunca se han modificado

El hash del blob de `rvsAS.dat` es **`d5e2ae8b3998f0b54e91dec002f6832dc86576d2` en los
10 commits**, y el de `rvsAP.dat` es `8ebbf7c2503019fc59b3ce834904176bea5997a8` en los
10. Byte-idénticos desde el primer commit hasta hoy. No hay ninguna revisión intermedia
que pudiera contener pistas del generador.

**Conclusión**: el generador vivió fuera de control de versiones (probablemente en la
máquina del autor original del código C++, anterior a este repo).

## 2. Documentación de la columna 2 — UNA SOLA LÍNEA

La búsqueda de comentarios/docstrings en todo el historial (C++ y Python, todas las
ramas) devuelve **exactamente una** línea que documenta estos archivos.

En `src/trimero.cpp` (línea 232 en `2c0a4d6`; línea 201 ya en el commit raíz `ea14a75`):

```cpp
// getting the radial possition and the s an p-wave scattering lengths
R  = gsl_matrix_get(As,row,0) ;
AS = gsl_matrix_get(As,row,1) ;
AP = gsl_matrix_get(Ap,row,1) ;
```

Es decir, la intención declarada por el autor es:
- **columna 1** → posición radial `R`
- **columna 2** de `rvsAS.dat` → longitud de dispersión de onda s
- **columna 2** de `rvsAP.dat` → magnitud de dispersión de onda p

No existe ningún otro comentario, README, docstring ni nota — ni en el C++, ni en el
Python migrado, ni en ningún commit antiguo — que amplíe esto, indique unidades,
o mencione la fuente de los datos. El `README.md` solo lista los nombres de archivo
como "necesarios"; no describe su contenido.

Otros comentarios del C++ que dan contexto (no procedencia):
- `trimero.cpp:169` — `s=1; //s=0 only s-wave, s=1 s-wave + p-wave`
- `trimero.cpp:53` — `const int lc = 2 ;//lc=2`

## 3. Contenido numérico real de los archivos (verificado)

`rvsAS.dat`: 776 filas × 2 columnas.

| Magnitud | Valor |
|---|---|
| Columna 1 (R) | 111 → 2448, paso constante de 3 (u.a.) |
| Columna 2, primer valor (R=111) | `+3.1289135977` |
| Columna 2, último valor (R=2448) | `-15.9929682993` |
| Monotonía | estrictamente decreciente en todo el rango |
| Cruce por cero | R = 426 |

**Corrección a la hipótesis de partida**: el archivo **no** contiene "valores +3.1
decrecientes". El valor `+3.1289` es únicamente el **primer punto** de la tabla. La
columna atraviesa cero en R=426 y termina en **−15.99**, que sí cae dentro del rango
de literatura para la longitud de dispersión triplete e⁻–Rb (−16 a −20 a₀), en el
extremo de R grande. La coincidencia numérica de `3.1289` con la constante base
`muns_Rb ≈ 3.1312` afecta solo al primer punto de 776 y no se sostiene para el resto
de la tabla.

`rvsAP.dat`: 776 filas × 2 columnas, mismo eje R.

| Magnitud | Valor |
|---|---|
| Columna 2, rango | −680411.61 a +410651.63 |
| Primer valor (R=111) | `+324.9025` |
| Último valor (R=2448) | `−28433.0231` |
| Cruce por cero | R = 750 |

Los valores de `rvsAP.dat` presentan una divergencia entre R≈561 (+9229) y R≈798
(−82737), con órdenes de magnitud enormes — comportamiento compatible con una cantidad
cúbica (volumen de dispersión, ~a₀³) y con el paso por una resonancia, pero **esto no
está documentado en ninguna parte del repositorio**.

## Conclusiones y aplicación al proyecto

1. **El script generador no existe en este repositorio ni en el disco local.** No es
   una cuestión de buscar mejor: el historial completo son 10 commits y el listado
   exhaustivo de paths no contiene ningún candidato.
2. **La única evidencia documental interna** es el comentario
   `// getting the radial possition and the s an p-wave scattering lengths`,
   que identifica la columna 2 como longitud/magnitud de dispersión, no como defecto
   cuántico.
3. **La evidencia numérica** (asíntota en −15.99 a₀ a R grande) es consistente con esa
   lectura, no con la de un defecto cuántico efectivo.
4. **Punto abierto**: sigue sin estar documentado el método de cálculo, la referencia
   bibliográfica de los desplazamientos de fase usados, ni las unidades exactas.
   Pendiente de resolver por vía externa (autores originales / tesis en papel).
   **No se modifica `fermi_potentials.py` ni ningún otro código hasta entonces.**

## Referencias

- Repositorio: `git@github.com:dirac89/trimero_mod.git`
- Commit raíz: `ea14a75b1d9e1ff697fe93e3889cb34d9939a85e` (2018-03-05)
- Blob invariante `rvsAS.dat`: `d5e2ae8b3998f0b54e91dec002f6832dc86576d2`
- Blob invariante `rvsAP.dat`: `8ebbf7c2503019fc59b3ce834904176bea5997a8`
- Bahrim & Thumm — longitud de dispersión triplete e⁻–Rb (referencia externa citada por
  el usuario; **no** presente ni citada en el repositorio)
