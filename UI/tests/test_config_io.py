from __future__ import annotations
import os
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # UI/ on path

import config_io


def test_roundtrip_preserves_structure(tmp_path):
    cfg = {
        "queries": {"ecoli": "a.fa"}, "targets": {"ecoli": "g.fna"},
        "threads": 16, "max_suboptimal_hits": 1, "seed": None,
        "rnahybrid": {"max_hybrid_energy": -18},
        "intarna": {"accessibility_variant": "on",
                    "output": {"columns": "id1,id2,E", "overlap": "B"}},
        "plots": {"type": None}, "results_dir": "Data/Results/",
    }
    p = tmp_path / "config.yaml"
    assert config_io.save_config(p, cfg)["ok"] is True
    assert config_io.load_config(p) == cfg

def test_special_characters_survive_roundtrip(tmp_path):
    # The #10 regression: ':' '#' and a leading '-' must not corrupt the YAML.
    cfg = {"results_dir": "a: b # c", "note": "-leading", "variant": "on"}
    p = tmp_path / "config.yaml"
    config_io.save_config(p, cfg)
    assert config_io.load_config(p) == cfg

def test_missing_file_returns_empty(tmp_path):
    assert config_io.load_config(tmp_path / "nope.yaml") == {}

def test_bad_yaml_returns_error(tmp_path):
    p = tmp_path / "config.yaml"
    p.write_text("queries: [unterminated\n")
    assert "__error__" in config_io.load_config(p)

def test_helix_and_energy_set_survive_roundtrip(tmp_path):
    cfg = {
        "intarna": {
            "model": "B",
            "energy_set": "Turner99",
            "helix": {
                "min_bp": 3, "max_bp": 8, "max_internal_loop": None,
                "min_unpaired_probability": None, "max_energy": -1.5,
                "full_energy": True,
            },
        },
    }
    p = tmp_path / "config.yaml"
    assert config_io.save_config(p, cfg)["ok"] is True
    assert config_io.load_config(p) == cfg

def test_save_is_atomic_no_partial(tmp_path):
    p = tmp_path / "config.yaml"
    config_io.save_config(p, {"a": 1})
    # no leftover temp files in the directory
    assert [x.name for x in tmp_path.iterdir()] == ["config.yaml"]

def test_save_keeps_the_existing_files_mode(tmp_path):
    # A temp file + os.replace used to hand config.yaml the temp file's 0600.
    p = tmp_path / "config.yaml"
    for mode in (0o644, 0o640):
        p.write_text("a: 0\n")
        p.chmod(mode)
        config_io.save_config(p, {"a": 1})
        assert p.stat().st_mode & 0o777 == mode

def test_save_gives_a_new_file_the_umask_mode(tmp_path):
    old = os.umask(0o027)
    try:
        config_io.save_config(tmp_path / "config.yaml", {"a": 1})
    finally:
        os.umask(old)
    assert (tmp_path / "config.yaml").stat().st_mode & 0o777 == 0o640
