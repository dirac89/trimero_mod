---
name: run-simulation
description: Alias histórico de /sweep. Ejecutar un barrido de curvas — usa /sweep, que declara sistema y molécula, estima el coste antes de gastar CPU y audita el resultado.
---

# Skill: run-simulation (alias)

**Usa [`/sweep`](../sweep/SKILL.md).** Este comando se conserva sólo porque el
nombre está en documentos y transcripciones antiguas.

Invoca `/sweep` y sigue su protocolo: declarar la identidad del cálculo →
estimar y anunciar el coste → ejecutar → auditar el `.npz`.

> La versión anterior de esta skill pedía `n1` y `dc_field_au` y hablaba de dos
> sistemas. Esos parámetros son del camino **legado congelado**
> (`Trimer_energies_field`), y hoy el repositorio tiene siete paquetes de
> sistemas: seguirla llevaba a lanzar el cálculo equivocado. Se deja constancia
> del error en vez de borrarlo (regla §7 de `.claude/CLAUDE.md`).
