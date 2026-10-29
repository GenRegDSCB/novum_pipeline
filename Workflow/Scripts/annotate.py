import re
import sys
import pandas as pd

from _common import iter_fasta_headers, parse_header_id, ensure_parent


# Extracts key/value pairs from FASTA headers like:
# [gene=hscA] [locus_tag=b2526] [protein=[Fe-S] cluster biosynthesis chaperone HscA]
# A value may hold balanced brackets nested two deep; an unbalanced one fails its own match only.
_PLAIN = r"[^\[\]]"
_NESTED = rf"\[(?:{_PLAIN}|\[{_PLAIN}*\])*\]"
FIELD_RE = re.compile(rf"\[([^\s=\[\]]+)=((?:{_PLAIN}|{_NESTED})*)\]")

# Output column → the FASTA attribute keys feeding it, first non-empty wins. intersect.py
# imports ANNOTATION_COLUMNS so the consensus schema cannot drift from what this rule emits.
ANNOTATION_FIELDS = (
    ("gene_name",    ("gene",)),
    ("locus_tag",    ("locus_tag",)),
    ("protein_name", ("protein", "product")),
    ("protein_id",   ("protein_id",)),
)
ANNOTATION_COLUMNS = [column for column, _keys in ANNOTATION_FIELDS]


#----- First key holding a non-empty value, falling back to the last key's raw value (so an
#      empty `[gene=]` still reads as "" rather than becoming None) -----#
def _first_set(fields, keys):
    return next((fields[key] for key in keys if fields.get(key)), fields.get(keys[-1]))


#----- Pulls gene/locus/protein metadata out of FASTA headers and returns one deduped row per Gene -----#
def parse_fasta_annotations(fasta_path):
    records = []

    # Headers only: the sequences are never read here, so streaming them would build and
    # discard a multi-MB string per genome on every annotate job.
    for header in iter_fasta_headers(fasta_path):
        fields = dict(FIELD_RE.findall(header))
        records.append({
            "Gene": parse_header_id(header),
            **{column: _first_set(fields, keys) for column, keys in ANNOTATION_FIELDS},
        })

    annotations = pd.DataFrame(records)
    duplicates = annotations.duplicated(subset=["Gene"]).sum()
    if duplicates:
        print(
            f"annotate: dropped {duplicates} duplicate Gene row(s) from FASTA annotations.",
            file=sys.stderr,
        )
    annotations = annotations.drop_duplicates(subset=["Gene"])

    return annotations


#----- Left-joins the FASTA annotations onto the tidy table, inserting them after `insert_after` column -----#
def annotate_results(tidy_csv_path, fasta_path, output_path, insert_after):
    tidy = pd.read_csv(tidy_csv_path)

    # Checked before the FASTA is streamed: the anchor column depends on neither, so a
    # misconfigured insert_after should not first cost a pass over a genome-scale file.
    if insert_after not in tidy.columns:
        raise ValueError(f"tidy CSV is missing required column {insert_after!r}")

    annotations = parse_fasta_annotations(fasta_path)
    annotated = tidy.merge(annotations, on="Gene", how="left")

    cut = tidy.columns.get_loc(insert_after) + 1
    annotated = annotated[[*tidy.columns[:cut], *ANNOTATION_COLUMNS, *tidy.columns[cut:]]]

    annotated.to_csv(ensure_parent(output_path), index=False)


#----- Snakemake entry point: unpacks the injected `snakemake` object and calls the pure logic above -----#
def run_from_snakemake(snakemake):
    annotate_results(
        tidy_csv_path=snakemake.input.tidy,
        fasta_path=snakemake.input.target,
        output_path=snakemake.output.annotated,
        insert_after=snakemake.params.insert_after,
    )

if "snakemake" in globals():
    run_from_snakemake(snakemake)