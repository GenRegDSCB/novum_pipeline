import sys
import pandas as pd

from _common import ensure_parent
from annotate import ANNOTATION_COLUMNS
from tidy_intarna import OUTPUT_COLUMNS as INTARNA_TIDY_COLUMNS
from tidy_rnahybrid import OUTPUT_COLUMNS as RNAHYBRID_TIDY_COLUMNS


# Columns both arms already carry, taken once from the RNAhybrid side and dropped from IntaRNA
# before the merge. The four annotation columns come from annotate.py on both arms (imported,
# so adding one there cannot leave a _x/_y pair here); Gene_length does NOT — IntaRNA's comes
# from the target FASTA, RNAhybrid's from its own output column 2.
SHARED_ANNOTATION = [*ANNOTATION_COLUMNS, "Gene_length"]


#----- One arm's tidy columns as the consensus carries them: the pair keys and the shared
#      annotation drop out (kept once, above) and Position disambiguates to the arm -----#
def _arm_columns(tidy_columns, position_name, after_position=()):
    carried_once = {"miRNA", "Gene", *SHARED_ANNOTATION}
    columns = []
    for column in tidy_columns:
        if column in carried_once:
            continue
        columns.append(position_name if column == "Position" else column)
        if column == "Position":
            columns.extend(after_position)
    return columns


# Final consensus column order (minimal-rename: only Position is disambiguated). Derived from
# the two producers, so a column added to either tidy schema reaches the consensus instead of
# being dropped by the reindex below without a word.
CONSENSUS_COLUMNS = [
    "miRNA", "Gene",
    *SHARED_ANNOTATION,
    *_arm_columns(RNAHYBRID_TIDY_COLUMNS, "Position_rnahybrid"),
    *_arm_columns(INTARNA_TIDY_COLUMNS, "Position_intarna", ["Site_offset_nt"]),
]


#----- Keeps one row per `keys` group: the minimum `rank_col` (NaNs dropped first), tie-broken
#      by `tiebreak_col`. The dropna is load-bearing — an all-NaN group must not survive via keep="first" -----#
def select_best_hit(df, keys, rank_col, tiebreak_col):
    ranked = df.dropna(subset=[rank_col]).sort_values([rank_col, tiebreak_col], kind="mergesort")
    return ranked.drop_duplicates(subset=keys, keep="first")


#----- Inner-joins the two arms' representative hits per (miRNA, Gene), each chosen by its own gated energy; shared annotation kept from RNAhybrid, Position disambiguated -----#
def intersect_annotations(rnahybrid_csv, intarna_csv, output_csv):
    keys = ["miRNA", "Gene"]
    rnahybrid = pd.read_csv(rnahybrid_csv)
    intarna = pd.read_csv(intarna_csv)

    rh_best = select_best_hit(rnahybrid, keys, "Energy", "P_value")
    in_best = select_best_hit(intarna, keys, "E_hybrid", "E")

    rh_best = rh_best.rename(columns={"Position": "Position_rnahybrid"})
    in_best = in_best.rename(columns={"Position": "Position_intarna"})

    in_best = in_best.drop(columns=SHARED_ANNOTATION, errors="ignore")

    merged = rh_best.merge(in_best, on=keys, how="inner")

    # The join is pair-granular, so measure the site disagreement rather than gate on it.
    # Signed: positive means IntaRNA's site lies downstream.
    merged["Site_offset_nt"] = (
        (merged["Position_intarna"] - merged["Position_rnahybrid"]) * merged["Gene_length"]
    ).round(1)

    absent = [c for c in CONSENSUS_COLUMNS if c not in merged.columns]
    if absent:
        print(
            "intersect: consensus columns absent from inputs (filled NA): "
            + ", ".join(absent)
            + " — expected when intarna.output.columns is trimmed, or when the top-level "
              "seed is null and the seed columns are withheld.",
            file=sys.stderr,
        )
    merged = merged.reindex(columns=CONSENSUS_COLUMNS)

    merged.to_csv(ensure_parent(output_csv), index=False)
    print(
        f"intersect: {len(merged)} consensus pair(s) "
        f"from {len(rh_best)} RNAhybrid × {len(in_best)} IntaRNA best-hits",
        file=sys.stderr,
    )


#----- Snakemake entry point: unpacks the injected `snakemake` object and calls the pure logic above -----#
def run_from_snakemake(snakemake):
    intersect_annotations(
        rnahybrid_csv=snakemake.input.rnahybrid,
        intarna_csv=snakemake.input.intarna,
        output_csv=snakemake.output.annotated,
    )

if "snakemake" in globals():
    run_from_snakemake(snakemake)
