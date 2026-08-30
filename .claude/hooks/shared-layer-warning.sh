#!/usr/bin/env bash
# PostToolUse/Edit|Write — avisa al tocar la capa compartida.
#
# mathlib/ y basis/ son de las que dependen TODOS los sistemas y los golden
# files de caracterización. El wigner_3j canónico (586776d, 642fabd) obligó a
# regenerar goldens y a rehacer las fases 1-2 del híbrido.
set -euo pipefail

path=$(jq -r '.tool_input.file_path // .tool_response.filePath // ""')

case "$path" in
  *src/trimero/mathlib/*|*src/trimero/basis/*)
    cat <<'JSON'
{
  "systemMessage": "⚠️  Capa compartida (mathlib/ o basis/): de aquí dependen todos los sistemas y los golden files. Suite COMPLETA antes de commitear (poetry run pytest, ~11 min).",
  "hookSpecificOutput": {
    "hookEventName": "PostToolUse",
    "additionalContext": "Se ha modificado la capa compartida (src/trimero/mathlib/ o src/trimero/basis/). Los golden files de tests/systems/rb_neutral_perturber/characterization/ dependen de ella con rtol=1e-12. Antes de commitear hay que correr la suite completa (poetry run pytest, ~11 min), no sólo -m \"not slow\". Si un golden se mueve, es un cambio de física: hay que parar y reportarlo, no regenerar el golden."
  }
}
JSON
    exit 0
    ;;
esac

exit 0
