#!/usr/bin/env bash
# PreToolUse/Bash (git commit) — recuerda la suite completa. Avisa, no bloquea.
#
# Regla §5 de .claude/CLAUDE.md: `-m "not slow"` mientras se itera, suite
# completa antes de commitear. Los 16 tests `slow` son los que protegen los
# goldens del legado y las regresiones end-to-end.
set -euo pipefail

cmd=$(jq -r '.tool_input.command // ""')

case "$cmd" in
  *"git commit"*)
    cat <<'JSON'
{
  "systemMessage": "Recordatorio: la suite COMPLETA (poetry run pytest, ~11 min, incluye los 16 tests slow) va antes de commitear, no sólo -m \"not slow\".",
  "hookSpecificOutput": {
    "hookEventName": "PreToolUse",
    "additionalContext": "Antes de este commit debe haberse ejecutado la suite completa (poetry run pytest, ~11 min). Si sólo se corrió -m \"not slow\", faltan los 16 tests slow: los goldens de caracterización del legado y las regresiones end-to-end. Si no se ha hecho, dilo antes de commitear."
  }
}
JSON
    exit 0
    ;;
esac

exit 0
