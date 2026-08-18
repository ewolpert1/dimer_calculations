"""Portable defaults for external optimisation programs."""

from __future__ import annotations

import os
import shutil


def _executable_path(environment_variable: str, executable: str) -> str | None:
    """Return an explicitly configured path or locate an executable."""
    configured_path = os.environ.get(environment_variable)
    if configured_path:
        return configured_path
    return shutil.which(executable)


# Paths.
CAGES_FOLDER_PATH = "cages"
TRIALDEHYDES_FILE = "trialdehydes.txt"
TRIAMINES_FILE = "triamines.txt"
DIALDEHYDES_FILE = "dialdehydes.txt"
DIAMINES_FILE = "diamines.txt"

# External optimisation programs. Set the corresponding environment variable
# when the program is not available on PATH. ``SCHRODINGER_PATH`` is the root
# of a Schrödinger installation; the other two values are executable paths.
SCHRODINGER_PATH = os.environ.get("SCHRODINGER_PATH")
GULP_PATH = _executable_path("GULP_PATH", "gulp")
XTB_PATH = _executable_path("XTB_PATH", "xtb")
