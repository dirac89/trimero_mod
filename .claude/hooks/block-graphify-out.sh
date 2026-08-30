#!/usr/bin/env bash
# PreToolUse/Bash — bloquea un `git add`/`git stage` que arrastre graphify-out/.
#
# graphify-out/ es un artefacto derivado y regenerable con /graphify: está en
# .gitignore y fue purgado del histórico el 2026-08-21 (regla §10 de
# .claude/CLAUDE.md). Un `git add -f` o un `git add .` desde dentro del
# directorio lo devolvería al repositorio sin que se note.
#
# Dos precisiones, ambas aprendidas de falsos positivos reales al instalarlo:
#
#  1. Sólo se vigilan `add`/`stage`, no `commit`: son los que meten ficheros en
#     el índice, y un commit sólo puede incluir lo que ya esté allí. Vigilar el
#     commit bloqueaba mensajes que se limitaban a MENCIONAR graphify-out.
#  2. Se exige que `git add` abra un segmento de comando y que el artefacto sea
#     argumento suyo, en vez de buscar las cadenas sueltas por el comando
#     entero — si no, un comando que sólo NOMBRA la orden (documentación, o una
#     tubería hacia este mismo script, cuyo nombre ya contiene «graphify-out»)
#     quedaba bloqueado.
#  3. El nombre tiene que ser un COMPONENTE de ruta completo, no una subcadena:
#     si no, añadir este mismo fichero (block-graphify-out.sh) se bloqueaba solo.
set -euo pipefail

cmd=$(jq -r '.tool_input.command // ""')
dir="graphify-out"

# Segmento que empieza por `git add`/`git stage`: principio de línea o tras ; & | && || (
segmento=$(printf '%s' "$cmd" \
  | grep -Eo '(^|[;&|(]|&&|\|\|)[[:space:]]*git[[:space:]]+(add|stage)([[:space:]][^;&|]*)?' \
  | head -1 || true)

bloquear=0
for tok in $segmento; do
  # Normalizar a /tok/ para exigir componente completo:
  #   graphify-out/nodes.json  -> /graphify-out/... ✔ bloquea
  #   hooks/block-graphify-out.sh -> /hooks/block-graphify-out.sh/ ✘ permite
  case "/$tok/" in
    */"$dir"/*) bloquear=1 ;;
  esac
done

if [ "$bloquear" -eq 1 ]; then
  cat <<'JSON'
{
  "hookSpecificOutput": {
    "hookEventName": "PreToolUse",
    "permissionDecision": "deny",
    "permissionDecisionReason": "El artefacto derivado no se versiona (regla §10 de .claude/CLAUDE.md): es regenerable con /graphify y fue purgado del histórico el 2026-08-21. Se mantiene en local; no entra en git."
  }
}
JSON
fi

exit 0
