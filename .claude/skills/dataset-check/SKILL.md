---
name: dataset-check
description: Audita un .npz de resultados antes de darlo por bueno o commitearlo — metadatos de molécula, NaN, carácter de manifold, solapamiento y continuidad de la rama. Usa esto tras cada barrido y antes de commitear plots/.
---

# Skill: dataset-check

## Propósito

Un barrido puede terminar sin error y estar mal: la rama se pierde en un cruce
evitado, el peso de manifold cae y la curva pasa a ser otro objeto físico. Pasó,
y la reparación son los commits `a0730b1` y `b85295d`.

## Trigger

`/dataset-check [ruta.npz]`, o: «audita el resultado», «¿esta curva es buena?»,
y automáticamente como paso 4 de `/sweep`.

## Flujo

1. **Elegir el fichero**. Con argumento, ése. Sin argumento, el más reciente:
   ```bash
   ls -t plots/*/data/*.npz | head -1
   ```
2. **Delegar en el subagente `dataset-auditor`**, pasándole la ruta. Ahí está el
   detalle de las comprobaciones por esquema (polar, polar doble, neutro).
3. **Presentar el veredicto** tal cual, sin suavizarlo, y actuar:

| veredicto | qué hacer |
|---|---|
| apto | seguir; ya se puede commitear |
| apto con reservas | mirar la figura en los `R` señalados antes de usar la curva en una publicación |
| no apto | recalcular el tramo con `--min-substep` menor (o `--step` menor). **No** subir `--weight`/`--overlap` para que encaje |

## Lo que se comprueba (resumen)

- **Identidad**: `molecule`, `B_hz`, `d_debye`, `n_manifold`, `N_max`, `M_J`,
  `schema_version`. Un `.npz` polar sin `molecule` **se rechaza**; y las
  constantes deben cuadrar con la especie (KRb 1.114 GHz / 0.566 D; RbCs
  490.17 MHz / 1.225 D). También que el nombre y el directorio digan lo mismo
  que el contenido.
- **Sanidad numérica**: sin NaN/Inf, `R` creciente, `E` real y finita, saltos
  anómalos localizados con su `R`.
- **Carácter**: `W ≥ character_weight` en todo el rango y `W → ~1` en el borde
  superior de `R`.
- **Seguimiento**: `overlap ≥ overlap_threshold`; puntos por debajo del umbral
  catastrófico (rama resembrada); todos los avisos ya grabados en
  `warning_R`/`warning_overlap`/`warning_weight`; cambios del índice `K`.

## Nota

Los `.npz` son resultados verificados y **se commitean** (`plots/` nunca va al
`.gitignore`). Precisamente por eso se auditan antes: un dataset malo commiteado
es caro de detectar después.
