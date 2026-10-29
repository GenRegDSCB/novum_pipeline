"""Load/save Config/config.yaml as structured data for the UI bridge.

Pure and dependency-light (stdlib + PyYAML only) so it is unit-testable without
importing pywebview/PyQt. `launcher.py` delegates its API methods here."""
from __future__ import annotations

import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "data-prep"))
from config_writer import atomic_write  # noqa: E402  stdlib-only

# libyaml's C parser is ~12x faster on this config and is what the conda env ships; the
# pure-Python classes are only a fallback for an interpreter built without it.
try:
    from yaml import CSafeLoader as _Loader, CSafeDumper as _Dumper
except ImportError:                                  # pragma: no cover - env without libyaml
    from yaml import SafeLoader as _Loader, SafeDumper as _Dumper


#----- Parse the on-disk config into a plain dict; never raises to the UI -----#
def load_config(path):
    p = Path(path)
    try:
        text = p.read_text()
    except OSError:
        return {}
    try:
        data = yaml.load(text, Loader=_Loader)
    except yaml.YAMLError as exc:
        return {"__error__": str(exc)}
    return data if isinstance(data, dict) else {}


#----- Emit `config` as YAML with safe quoting, written atomically -----#
def save_config(path, config):
    p = Path(path)
    text = yaml.dump(config, Dumper=_Dumper, sort_keys=False,
                     default_flow_style=False, allow_unicode=True)
    p.parent.mkdir(parents=True, exist_ok=True)
    atomic_write(p, text)
    return {"ok": True, "path": str(p)}
