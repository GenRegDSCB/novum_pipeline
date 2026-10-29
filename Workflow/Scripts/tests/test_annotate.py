from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # Workflow/Scripts on path

import pandas as pd
import pytest
import annotate  # importable because the trailing call is guarded


def _fields(header):
    return dict(annotate.FIELD_RE.findall(header))


# Verbatim E. coli K-12 headers (GCF_000005845.2): the value opens with its own bracket pair.
HSCA = ("lcl|NC_000913.3_cds_NP_417021.1_2506 [gene=hscA] [locus_tag=b2526] "
        "[db_xref=UniProtKB/Swiss-Prot:P0A6Z1] "
        "[protein=[Fe-S] cluster biosynthesis chaperone HscA] [protein_id=NP_417021.1] "
        "[location=complement(2657085..2658935)] [gbkey=CDS]")


def test_plain_fields_parse_unchanged():
    assert _fields("g1 [gene=thrL] [locus_tag=b0001] [protein=thr operon leader peptide]") == {
        "gene": "thrL", "locus_tag": "b0001", "protein": "thr operon leader peptide",
    }


def test_a_value_opening_with_brackets_is_kept_whole():
    fields = _fields(HSCA)
    assert fields["protein"] == "[Fe-S] cluster biosynthesis chaperone HscA"
    # ...without swallowing the field that follows it.
    assert fields["protein_id"] == "NP_417021.1"
    assert fields["location"] == "complement(2657085..2658935)"


@pytest.mark.parametrize("value", [
    "methylenetetrahydrofolate reductase [NAD(P)H]",       # brackets close the value
    "biotin--[acetyl-CoA-carboxylase] ligase",              # brackets mid-value
    "ditrans,polycis-undecaprenyl-diphosphate synthase [(2E,6E)-farnesyl-diphosphate specific]",
    "outer [middle [inner] middle] outer",                  # deeper than any NCBI header
])
def test_bracketed_values_survive_in_every_position(value):
    header = f"g1 [gene=x] [protein={value}] [protein_id=WP_1.1] [gbkey=CDS]"
    assert _fields(header) == {
        "gene": "x", "protein": value, "protein_id": "WP_1.1", "gbkey": "CDS",
    }


def test_a_bracketed_value_last_in_the_header_is_kept_whole():
    assert _fields("g1 [gene=x] [protein=fumarate reductase [Fe-S] protein]") == {
        "gene": "x", "protein": "fumarate reductase [Fe-S] protein",
    }


def test_empty_values_and_unspaced_fields_still_parse():
    assert _fields("g1 [gene=][locus_tag=b1]") == {"gene": "", "locus_tag": "b1"}


def test_free_text_between_or_after_fields_is_not_a_value():
    assert _fields("g1 [gene=abc] putative kinase") == {"gene": "abc"}
    assert _fields("g1 [gene=abc] putative [Fe-S] kinase [protein=k]") == {
        "gene": "abc", "protein": "k",
    }


def test_an_unbalanced_value_never_swallows_the_next_field():
    assert _fields("g1 [protein=foo [bar] [gene=x] [locus_tag=b1]")["gene"] == "x"


def _annotate(tmp_path, headers, genes):
    fasta = tmp_path / "t.fna"
    fasta.write_text("".join(f">{header}\nAUGC\n" for header in headers))
    tidy = tmp_path / "tidy.csv"
    pd.DataFrame({"miRNA": ["m1"] * len(genes), "Gene": genes,
                  "P_value": [0.01] * len(genes)}).to_csv(tidy, index=False)
    out = tmp_path / "out.csv"
    annotate.annotate_results(tidy, fasta, out, insert_after="P_value")
    return pd.read_csv(out, keep_default_na=False)


def test_annotate_results_writes_the_whole_bracketed_protein_name(tmp_path):
    annotated = _annotate(tmp_path, [HSCA], ["lcl|NC_000913.3_cds_NP_417021.1_2506"])
    assert annotated.loc[0, "protein_name"] == "[Fe-S] cluster biosynthesis chaperone HscA"
    assert annotated.loc[0, "gene_name"] == "hscA"


def test_annotate_results_joins_on_the_token_the_tools_report(tmp_path):
    # RNAhybrid 2.1.2 and IntaRNA 3.4.1 both report the first whitespace token, so a
    # free-text description before the attributes must not leak into the join key.
    annotated = _annotate(tmp_path, ["gene1 putative kinase [gene=abc] [protein=kinase]"],
                          ["gene1"])
    assert annotated.loc[0, "gene_name"] == "abc"
    assert annotated.loc[0, "protein_name"] == "kinase"
