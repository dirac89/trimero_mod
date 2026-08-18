"""Pone `src/` en el path para toda la suite.

Nota para el paso 2 del refactor: cuando `src/trimero/` sea un paquete
instalado, este fichero deja de ser necesario y los tests importarán
`from trimero...` directamente.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
