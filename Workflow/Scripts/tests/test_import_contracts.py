"""The env-driven import rules for the five shared sibling modules, in one table.

Each module runs in a different conda env, so its import list is a contract, not a style
preference: the two config modules are imported by the Snakefile into the DRIVER env (a
third-party import there breaks every snakemake invocation, --dry-run included), and
_common/_gumbel_fit are imported by scripts running in rnahybrid.yaml, which ships neither
pandas nor numpy.
"""
from __future__ import annotations
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # Workflow/Scripts

SCRIPTS = Path(__file__).resolve().parent.parent

#----- (module, the exact top-level imports it is allowed, why it is pinned) -----#
IMPORT_CONTRACTS = [
    ("_intarna_config", set(),
     "the Snakefile imports it at DAG-build time, in the driver env"),
    ("_rnacalibrate_config", {"import math", "import bisect"},
     "same driver-env import; needs log/sqrt and bisect only"),
    ("_gumbel_fit", {"import math", "from collections import namedtuple"},
     "runs in rnahybrid.yaml, which has neither pandas nor numpy"),
    ("_common", {"import shutil", "from pathlib import Path"},
     "imported by the RNAhybrid arm's scripts, whose env ships no pandas"),
    ("_length", {"import numpy as np", "import pandas as pd"},
     "postprocess-only; may use the scientific stack"),
    ("_report", {"import pandas as pd", "from _common import ensure_parent"},
     "enhance-report scaffold; postprocess-only, so pandas is allowed"),
]


#----- Top-level imports only: an indented one is inside a function and never runs at import -----#
def _top_level_imports(module_name):
    source = (SCRIPTS / f"{module_name}.py").read_text()
    return {line.strip() for line in source.splitlines()
            if line.startswith(("import ", "from ")) and not line.startswith("from __future__")}


@pytest.mark.parametrize("module_name,allowed,reason",
                         IMPORT_CONTRACTS,
                         ids=[row[0] for row in IMPORT_CONTRACTS])
def test_module_keeps_its_import_contract(module_name, allowed, reason):
    extra = sorted(_top_level_imports(module_name) - allowed)
    assert not extra, f"{module_name} may not import {extra}: {reason}."


def test_the_table_covers_every_shared_sibling_module():
    # A new _-prefixed sibling is an env decision; it must be pinned here, not left implicit.
    on_disk = {p.stem for p in SCRIPTS.glob("_*.py") if not p.stem.startswith("__")}
    # _rnahybrid_worker is a GNU Parallel entry point, not a shared import target.
    on_disk -= {"_rnahybrid_worker"}
    assert on_disk == {row[0] for row in IMPORT_CONTRACTS}
