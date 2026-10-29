"""The three enhance_* reports share one sectioning rule.

A scoring column is rendered in its own block and therefore never echoed in the generic
metadata header. Pinned across all three because the reports drifted once already: the
consensus report showed IntaRNA's energies as plain metadata while the per-arm report gave
them a block, and nobody chose that.
"""
from __future__ import annotations
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # Workflow/Scripts

import enhance_consensus
import enhance_intarna
import enhance_rnahybrid
import tidy_intarna
import tidy_rnahybrid
from _report import INTARNA_DUPLEX_COLUMNS, RNAHYBRID_DUPLEX_COLUMNS


REPORTS = (
    (enhance_rnahybrid, [tidy_rnahybrid.ENERGY_COLUMNS]),
    (enhance_intarna,   [tidy_intarna.ENERGY_COLUMNS, tidy_intarna.SEED_COLUMNS]),
    (enhance_consensus, [tidy_rnahybrid.ENERGY_COLUMNS,
                         tidy_intarna.ENERGY_COLUMNS, tidy_intarna.SEED_COLUMNS]),
)


def test_every_report_takes_its_sections_from_the_producers():
    # The groups are the tidy steps' own, so a column added to a tidy schema lands in its
    # block rather than falling through to the metadata lines.
    for module, groups in REPORTS:
        assert module.SECTION_COLUMNS == groups


def _consensus_frame():
    row = {column: "x" for column in RNAHYBRID_DUPLEX_COLUMNS + INTARNA_DUPLEX_COLUMNS}
    row["subseqDP"] = "ACGU&UGCA"
    row["hybridDP"] = "((((&))))"
    for group in enhance_consensus.SECTION_COLUMNS:
        row.update({column: 1.0 for column in group})
    row["miRNA"], row["Gene"] = "m1", "g1"
    return pd.DataFrame([row])


def test_consensus_blocks_both_arms_rather_than_only_intarna(tmp_path):
    src, out = tmp_path / "consensus.csv", tmp_path / "consensus.txt"
    _consensus_frame().to_csv(src, index=False)
    enhance_consensus.enhance_results(src, out)
    text = out.read_text()

    # The metadata header runs to the first blank line; every scored column must sit past it.
    header = text.split("\n\n", 1)[0]
    for group in enhance_consensus.SECTION_COLUMNS:
        for column in group:
            assert f"{column}:" in text, f"{column} missing from the report"
            assert f"{column}:" not in header, f"{column} echoed as metadata"

    assert "miRNA: m1" in header and "Gene: g1" in header


def test_an_all_na_group_is_not_rendered_as_a_section(tmp_path):
    # The seed columns are all-NA under --noSeed, and intersect reindexes them back in even
    # when IntaRNA never wrote them. An empty block is not a section in either report.
    frame = _consensus_frame()
    for column in tidy_intarna.SEED_COLUMNS:
        frame[column] = pd.NA

    src, out = tmp_path / "c.csv", tmp_path / "c.txt"
    frame.to_csv(src, index=False)
    enhance_consensus.enhance_results(src, out)
    text = out.read_text()

    for column in tidy_intarna.SEED_COLUMNS:
        assert column not in text
    assert "E_hybrid:" in text          # a populated group still renders
    assert "Energy:" in text


def test_a_partly_populated_group_still_renders_in_full(tmp_path):
    # Only a wholly empty group drops; one NA column inside a live block stays, so the
    # block's column list does not silently vary per file.
    frame = _consensus_frame()
    frame["ED1"] = pd.NA

    src, out = tmp_path / "c.csv", tmp_path / "c.txt"
    frame.to_csv(src, index=False)
    enhance_consensus.enhance_results(src, out)
    assert "ED1: NA" in out.read_text()
