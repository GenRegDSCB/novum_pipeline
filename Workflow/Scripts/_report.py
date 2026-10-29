"""Shared driver for the per-record enhance_* reports: a metadata block, one block per scoring
group, then each arm's duplex. Imports pandas, so only postprocess-env scripts pull it in."""

import pandas as pd

from _common import ensure_parent


LABEL_TARGET = "target 5' "
LABEL_MIRNA  = "miRNA  3' "
LABEL_INDENT = "          "
MARK_3PRIME  = "   3'"
MARK_5PRIME  = "   5'"
SEPARATOR    = "-" * 39


#----- Renders one cell as text, mapping pandas NA to the literal "NA" so the report never shows "nan" -----#
def fmt(value):
    return "NA" if pd.isna(value) else str(value)


#----- Raises unless every column in `columns` is present, naming all the missing ones at once -----#
def require_columns(annotated, columns, label):
    missing = [column for column in columns if column not in annotated.columns]
    if missing:
        raise ValueError(f"{label} missing required column(s): {', '.join(missing)}")


#----- Column name → itertuples position. `row` is a positional tuple, so every lookup below
#      resolves names once here rather than per record -----#
def _index_of(annotated):
    return {col: idx for idx, col in enumerate(annotated.columns)}


#----- Resolves `columns` to (name, itertuples index) pairs, skipping any that are absent -----#
def column_pairs(annotated, columns):
    index_of = _index_of(annotated)
    return [(column, index_of[column]) for column in columns if column in index_of]


#----- Writes one `name: value` line per pair — the report's metadata-block format, shared by the driver below and the arm-specific energy/seed blocks -----#
def write_pairs(fh, pairs, row):
    for column, idx in pairs:
        fh.write(f"{column}: {fmt(row[idx])}\n")


#----- Resolves each column group to its (name, index) pairs. One report section per group, in
#      the order given, dropping any group that is absent or holds no value at all: the seed
#      columns are all-NA under --noSeed, and consensus reindexes them back in even when
#      IntaRNA never wrote them, so "present" alone would section an empty block -----#
def column_blocks(annotated, groups):
    blocks = []
    for group in groups:
        pairs = column_pairs(annotated, group)
        if pairs and any(annotated[column].notna().any() for column, _idx in pairs):
            blocks.append(pairs)
    return blocks


#----- Writes each resolved group as its own block, a blank line ahead of each -----#
def write_blocks(fh, blocks, row):
    for block in blocks:
        fh.write("\n")
        write_pairs(fh, block, row)


#----- Per record: metadata (every column not in `skip_cols`), `render_record`, the separator -----#
def stream_records(annotated, output_path, skip_cols, render_record):
    index_of = _index_of(annotated)
    metadata_pairs = [(col, idx) for col, idx in index_of.items() if col not in skip_cols]

    with ensure_parent(output_path).open("w", encoding="utf-8") as fh:
        for row in annotated.itertuples(index=False, name=None):
            write_pairs(fh, metadata_pairs, row)
            render_record(fh, row, index_of)
            fh.write(f"\n{SEPARATOR}\n\n")


# Duplex-string columns each renderer consumes.
RNAHYBRID_DUPLEX_COLUMNS = ["Target_unmatches", "Target_matches", "miRNA_matches", "miRNA_unmatches"]
INTARNA_DUPLEX_COLUMNS   = ["subseqDP", "hybridDP"]


#----- Writes RNAhybrid's 4-line ASCII duplex (target 5'->3' on top, miRNA 3'->5' below).
#      Takes the whole row and reads its own columns, so no caller respells them in order -----#
def render_rnahybrid_duplex(fh, row, index_of):
    target_unmatches, target_matches, mirna_matches, mirna_unmatches = (
        row[index_of[column]] for column in RNAHYBRID_DUPLEX_COLUMNS
    )
    fh.write(f"{LABEL_TARGET}{fmt(target_unmatches)}{MARK_3PRIME}\n")
    fh.write(f"{LABEL_INDENT}{fmt(target_matches)}\n")
    fh.write(f"{LABEL_INDENT}{fmt(mirna_matches)}\n")
    fh.write(f"{LABEL_MIRNA}{fmt(mirna_unmatches)}{MARK_5PRIME}\n")


#----- Writes IntaRNA's 3-line dot-bracket duplex from subseqDP/hybridDP (miRNA shown reversed) -----#
def render_intarna_duplex(fh, row, index_of):
    subseq_dp, hybrid_dp = (row[index_of[column]] for column in INTARNA_DUPLEX_COLUMNS)
    # IntaRNA writes both fields as "target&query"; partition yields ("text", "", "") when
    # the separator is absent, so a malformed field degrades to an empty query strand.
    target_seq, _, query_seq = str(subseq_dp).partition("&")
    target_dp, _, _          = str(hybrid_dp).partition("&")
    n = max(len(target_seq), len(target_dp))
    indicator = "".join(
        "|" if i < len(target_dp) and target_dp[i] == "(" else " "
        for i in range(n)
    )
    query_display = query_seq[::-1]
    fh.write(f"{LABEL_TARGET}{target_seq}{MARK_3PRIME}\n")
    fh.write(f"{LABEL_INDENT}{indicator}\n")
    fh.write(f"{LABEL_MIRNA}{query_display}{MARK_5PRIME}\n")


_DUPLEX_COLUMNS = {render_rnahybrid_duplex: RNAHYBRID_DUPLEX_COLUMNS,
                   render_intarna_duplex:   INTARNA_DUPLEX_COLUMNS}


#----- One enhance report: metadata, the `sections` blocks, then each duplex after a blank line.
#      The metadata skip set and the column guard both derive from `sections` and `renderers` -----#
def write_report(annotated_csv_path, output_path, sections, renderers, label):
    annotated = pd.read_csv(annotated_csv_path)
    duplex_columns = [column for render in renderers for column in _DUPLEX_COLUMNS[render]]
    require_columns(annotated, duplex_columns, label)
    blocks = column_blocks(annotated, sections)

    def render_record(fh, row, index_of):
        write_blocks(fh, blocks, row)
        for render in renderers:
            fh.write("\n")
            render(fh, row, index_of)

    stream_records(annotated, output_path, set(duplex_columns).union(*sections), render_record)
