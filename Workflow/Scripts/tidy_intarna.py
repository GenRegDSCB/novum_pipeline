import pandas as pd

from _common import iter_fasta_records, parse_header_id, ensure_parent
from _intarna_config import INTARNA_SEED_COLUMNS
from _length import add_length_corrected, corrected_name


# The two IntaRNA column groups, in report order. Spelled once here because the enhance
# reports render each as its own block and OUTPUT_COLUMNS orders the CSV from them.
ENERGY_COLUMNS = ["E", "E_hybrid", corrected_name("E_hybrid"), "ED1", "ED2", "Pu1", "Pu2"]
SEED_COLUMNS   = list(INTARNA_SEED_COLUMNS)

OUTPUT_COLUMNS = [
    "miRNA", "Gene", *ENERGY_COLUMNS,
    "Gene_length", "Start1", "End1", "Start2", "End2", "Position",
    "subseqDP", "hybridDP",
    *SEED_COLUMNS,
]

# Columns that must be present in IntaRNA CSV output for this script to work.
# E_hybrid is the gate and sort key, so a trimmed intarna.output.columns must fail here
# rather than silently produce an ungated file.
REQUIRED_COLS = {"id1", "id2", "start1", "end1", "start2", "end2", "E", "E_hybrid"}

# Rows per read batch: ~100 MB, already below the surviving set's own footprint (the one term
# chunking cannot remove), and throughput is flat from 100k to 1M rows — a pure memory knob.
READ_CHUNK_ROWS = 250_000

# One pair. Spelled once because the cap is applied twice — once per batch, once globally —
# and the two must group identically or the pre-cap stops being a no-op on the result.
CAP_KEYS = ["miRNA", "Gene"]

# IntaRNA writes uppercase NAN in the seed columns under --noSeed (top-level seed: null),
# and that spelling is outside pandas' default NA set — the five columns would land as
# Python strings and cost 34% of the frame. Declaring it types them as float, so they also
# reach the enhance report as "NA" rather than the literal string "NAN".
INTARNA_NA_VALUES = ["NAN"]


#----- Streams the target FASTA once and returns a {Gene → sequence length} map for downstream position normalization -----#
def _parse_gene_lengths(fasta_path):
    return {parse_header_id(header): len(seq) for header, seq in iter_fasta_records(fasta_path)}


#----- Column guard, run once against the header rather than per batch. Pu1/Pu2 stay out
#      of REQUIRED_COLS: trimming intarna.output.columns is legitimate. But a floor
#      against a column that is not there would silently admit every row -----#
def _require_intarna_columns(columns, floors):
    missing = REQUIRED_COLS - set(columns)
    if missing:
        raise ValueError(f"IntaRNA CSV is missing required columns: {sorted(missing)}")

    for column, floor in floors:
        if floor is not None and column not in columns:
            raise ValueError(
                f"IntaRNA CSV is missing {column}, required by the accessibility floor "
                f"({floor}). Keep {column} in intarna.output.columns or clear the floor."
            )


#----- Everything that is row-independent and cheap on a doomed row — rename, validate, gate,
#      floor — so it runs per batch and lets all but the survivors go. The sort, the per-pair
#      cap and the annotation need the whole surviving set and stay in tidy_intarna below -----#
def _prepare_chunk(chunk, gene_lengths, max_hybrid_energy, floors):
    chunk = chunk.rename(columns={"id1": "Gene", "id2": "miRNA",
                                  "start1": "Start1", "end1": "End1",
                                  "start2": "Start2", "end2": "End2"})

    # Validated on the DISTINCT genes, before any filtering, so a target FASTA that does
    # not cover the output still fails even when the offending rows would have been gated
    # away. Row counts are only computed on the error path.
    genes = chunk["Gene"].unique()
    unknown = [gene for gene in genes if gene not in gene_lengths]
    if unknown:
        count = int(chunk["Gene"].isin(unknown).sum())
        raise ValueError(f"Gene_length missing for {count} record(s): {unknown[:5]}")

    nonpositive = [gene for gene in genes if gene_lengths[gene] <= 0]
    if nonpositive:
        raise ValueError(f"Non-positive Gene_length in target FASTA for: {nonpositive}")

    # The authoritative E_hybrid gate, then the Pu floors, which must run before the cap.
    # One mask, one take; NaN fails every comparison, so unknown accessibility never passes.
    mask = None
    if max_hybrid_energy is not None:
        mask = chunk["E_hybrid"] <= max_hybrid_energy
    for column, floor in floors:
        if floor is None:
            continue
        floored = chunk[column] >= floor
        mask = floored if mask is None else mask & floored

    return chunk if mask is None else chunk[mask]


#----- The per-pair cap, on a frame already sorted by E_hybrid. Applied per batch and again
#      globally, so both go through here and can never disagree about how they group -----#
def _cap_sorted(df, max_suboptimal_hits):
    if max_suboptimal_hits is None:
        return df
    return df.groupby(CAP_KEYS, sort=False).head(int(max_suboptimal_hits))


#----- The per-pair cap, hoisted into the read loop -----#
def _precap_chunk(chunk, max_suboptimal_hits):
    if max_suboptimal_hits is None:
        return chunk
    return _cap_sorted(chunk.sort_values("E_hybrid", kind="stable"), max_suboptimal_hits)


#----- Reads IntaRNA's CSV in batches, keeping only gated/floored rows, then ranks and caps on E_hybrid and emits the CSV -----#
def tidy_intarna(input_path, target_fasta_path, output_path,
                 max_hybrid_energy=None, max_suboptimal_hits=None,
                 min_target_unpaired_probability=None,
                 min_query_unpaired_probability=None):
    # Cutoffs resolved to float once here rather than per batch, where they are loop-invariant.
    gate = None if max_hybrid_energy is None else float(max_hybrid_energy)
    floors = tuple(
        (column, None if value is None else float(value))
        for column, value in (("Pu1", min_target_unpaired_probability),
                              ("Pu2", min_query_unpaired_probability))
    )
    # E/E_hybrid typed by the C parser rather than cast afterwards: an astype() would run on
    # the whole batch, most of which the gate and cap below are about to discard.
    read_options = dict(sep=";", dtype={"id1": str, "id2": str, "E": float, "E_hybrid": float},
                        na_values=INTARNA_NA_VALUES)

    header = pd.read_csv(input_path, nrows=0, **read_options)
    _require_intarna_columns(header.columns, floors)

    gene_lengths = _parse_gene_lengths(target_fasta_path)

    # Batched, and capped per batch: how much of a batch dies in the filters alone depends entirely on the config
    prepared = [
        _precap_chunk(_prepare_chunk(chunk, gene_lengths, gate, floors), max_suboptimal_hits)
        for chunk in pd.read_csv(input_path, chunksize=READ_CHUNK_ROWS, **read_options)
    ]
    df = pd.concat(prepared, ignore_index=True)

    # Rank by the gated quantity, then cap. Capping first could discard a qualifying row
    # in favour of a better-total-E one that fails the gate.
    df = _cap_sorted(df.sort_values("E_hybrid", kind="stable"), max_suboptimal_hits)

    # Annotated only once the cap has run: nothing above reads either column, and the gate
    # alone lets through far more rows than survive. Position is a 0-1 fraction for cross-arm
    # comparability with tidy_rnahybrid, which anchors the same way.
    df = df.assign(
        Gene_length=lambda d: d["Gene"].map(gene_lengths).astype(int),
        Position=lambda d: d["Start1"].astype(float) / d["Gene_length"],
    )

    # After the cap, not before: the correction ranks a pair's genes against each other, so a
    # gene should weigh as close to once as the cap allows. Column only, never a gate — see
    # _length.py.
    df = add_length_corrected(df, "E_hybrid")

    df = df[[c for c in OUTPUT_COLUMNS if c in df.columns]]

    df.to_csv(ensure_parent(output_path), index=False)


#----- Snakemake entry point: unpacks the injected `snakemake` object and calls the pure logic above -----#
def run_from_snakemake(snakemake):
    tidy_intarna(
        input_path=snakemake.input.csv,
        target_fasta_path=snakemake.input.target,
        output_path=snakemake.output.tidy,
        max_hybrid_energy=snakemake.params.max_hybrid_energy,
        max_suboptimal_hits=snakemake.params.max_suboptimal_hits,
        min_target_unpaired_probability=snakemake.params.min_target_unpaired_probability,
        min_query_unpaired_probability=snakemake.params.min_query_unpaired_probability,
    )

if "snakemake" in globals():
    run_from_snakemake(snakemake)
