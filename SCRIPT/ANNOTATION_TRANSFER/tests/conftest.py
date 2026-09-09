from __future__ import annotations

import sys
from pathlib import Path


TEST_DIR = Path(__file__).resolve().parent
MODULE_DIR = TEST_DIR.parent
REPOSITORY_ROOT = TEST_DIR.parents[2]

sys.path.insert(0, str(MODULE_DIR))


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "external_tools: requires local BLAST+ and MAFFT executables",
    )
    config.addinivalue_line(
        "markers",
        "integration: exercises genePrediction.sh/downstream pipeline helpers",
    )
