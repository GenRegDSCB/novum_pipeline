from __future__ import annotations

import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import tempfile
from pathlib import Path

from config_writer import replace_block, update_path_blocks

SAMPLE = """\
#----- miRNAs -----#
queries:
  escherichia: Data/Raw/RNAs/validated_escherichia.fa

#----- Genomes -----#
targets:
  escherichia: Data/Raw/genomes/old.fna

threads: 16
intarna:
  seed:
    query_range: null
""".splitlines()

def test_replace_populated_block_preserves_everything_else():
    out = replace_block(SAMPLE, "queries", {
        "veillonella_parvula": "Data/Raw/RNAs/veillonella_parvula.fa",
        "fusobacterium_nucleatum": "Data/Raw/RNAs/fusobacterium_nucleatum.fa",
    })
    text = "\n".join(out)
    assert "  veillonella_parvula: Data/Raw/RNAs/veillonella_parvula.fa" in text
    assert "  fusobacterium_nucleatum: Data/Raw/RNAs/fusobacterium_nucleatum.fa" in text
    assert "validated_escherichia.fa" not in text          # old entry gone
    assert "#----- miRNAs -----#" in text                  # comment above kept
    assert "#----- Genomes -----#" in text                 # terminator kept
    assert "threads: 16" in text                           # unrelated param kept
    assert "    query_range: null" in text                 # nested key untouched

def test_does_not_match_nested_query_range():
    # 'query_range:' is nested and must NOT be treated as the 'queries' block.
    out = replace_block(SAMPLE, "targets", {"x": "y.fna"})
    text = "\n".join(out)
    assert "  x: y.fna" in text
    assert "    query_range: null" in text
    assert "old.fna" not in text

def test_empty_block_inserts_without_removing():
    lines = ["queries:", "", "#----- next -----#", "threads: 1"]
    out = replace_block(lines, "queries", {"a": "p.fa"})
    assert out == ["queries:", "  a: p.fa", "", "#----- next -----#", "threads: 1"]

INTERIOR_COMMENT = """\
queries:
  old_one: Data/Raw/RNAs/old_one.fa
# a note the user wrote at column 0
  old_two: Data/Raw/RNAs/old_two.fa
targets:
  t: g.fna
""".splitlines()


def _entries_under(lines, key):
    """The indented entries between `key:` and the next top-level key."""
    out = []
    for ln in lines[lines.index(f"{key}:") + 1:]:
        if ln.strip() == "" or ln.lstrip().startswith("#"):
            continue
        if ln[:1] not in (" ", "\t"):
            break
        out.append(ln.strip())
    return out


def test_column_zero_comment_inside_a_block_does_not_strand_stale_entries():
    # A comment at column 0 used to terminate the block, so the entries below it
    # survived alongside the new ones — a duplicate YAML key that PyYAML loads
    # without complaint (last wins), silently running the pipeline on the old path.
    out = replace_block(INTERIOR_COMMENT, "queries", {"new": "n.fa"})
    assert _entries_under(out, "queries") == ["new: n.fa"]


def test_column_zero_comment_inside_a_block_is_preserved():
    # Dropping the stale entries must not take the user's comment with them.
    out = replace_block(INTERIOR_COMMENT, "queries", {"new": "n.fa"})
    assert "# a note the user wrote at column 0" in out
    assert _entries_under(out, "targets") == ["t: g.fna"]   # next block untouched


def test_indented_comment_inside_a_block_is_preserved_too():
    lines = ["queries:", "  # human samples", "  s1: a.fa", "threads: 1"]
    out = replace_block(lines, "queries", {"new": "n.fa"})
    assert _entries_under(out, "queries") == ["new: n.fa"]
    assert "  # human samples" in out
    assert out[-1] == "threads: 1"


def test_missing_key_raises():
    try:
        replace_block(["threads: 1"], "queries", {"a": "b"})
    except ValueError as e:
        assert "queries" in str(e)
    else:
        raise AssertionError("expected ValueError")


def test_update_path_blocks_keeps_the_configs_mode():
    # A temp file + os.replace used to hand config.yaml the temp file's 0600.
    d = Path(tempfile.mkdtemp())
    config = d / "config.yaml"
    for mode in (0o644, 0o640):
        config.write_text("queries:\n  a: a.fa\ntargets:\n  a: a.fna\n")
        config.chmod(mode)
        update_path_blocks(config, {"b": "b.fa"}, {"b": "b.fna"})
        assert config.stat().st_mode & 0o777 == mode
        assert "  b: b.fna" in config.read_text().splitlines()
    assert [p.name for p in d.iterdir()] == ["config.yaml"]   # no temp left behind

if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for t in tests:
        t()
        print(f"ok  {t.__name__}")
    print(f"\n{len(tests)} passed")
