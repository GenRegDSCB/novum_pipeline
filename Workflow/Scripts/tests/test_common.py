from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # Workflow/Scripts on path

import pytest
import _common


@pytest.mark.parametrize("header, expected", [
    # NCBI cds_from_genomic: the id, then bracketed attributes.
    ("lcl|NC_000913.3_cds_NP_417021.1_2506 [gene=hscA] [protein=[Fe-S] cluster]",
     "lcl|NC_000913.3_cds_NP_417021.1_2506"),
    ("lcl|X_cds_1", "lcl|X_cds_1"),
    # No bracket attributes: RNAhybrid and IntaRNA both report "gene1" for this record.
    ("gene1 hypothetical protein", "gene1"),
    ("gene1\thypothetical protein", "gene1"),
    ("gene1 putative kinase [gene=abc]", "gene1"),
])
def test_parse_header_id_is_the_token_the_tools_report(header, expected):
    assert _common.parse_header_id(header) == expected


@pytest.mark.parametrize("header", [
    "lcl|X_cds_1 [gene=a]", "gene1 hypothetical protein", "hsa-miR-21-5p MIMAT0000076",
])
def test_targets_and_queries_share_one_key(header):
    assert _common.parse_header_id(header) == _common.query_key(header)


def test_tidy_intarna_gene_lengths_key_on_the_reported_token(tmp_path):
    # tidy_intarna joins IntaRNA's id1 against these keys; a description in the header
    # would otherwise raise "Gene_length missing" for every record.
    import tidy_intarna
    fasta = tmp_path / "t.fna"
    fasta.write_text(">gene1 hypothetical protein\nAUGCAUGC\n>lcl|X_cds_2 [gene=a]\nAUG\n")
    assert tidy_intarna._parse_gene_lengths(fasta) == {"gene1": 8, "lcl|X_cds_2": 3}
