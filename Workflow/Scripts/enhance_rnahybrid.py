from _report import write_report, render_rnahybrid_duplex
from tidy_rnahybrid import ENERGY_COLUMNS


SECTION_COLUMNS  = [ENERGY_COLUMNS]
DUPLEX_RENDERERS = [render_rnahybrid_duplex]


#----- Per-record report: metadata, the energy block, then the RNAhybrid duplex -----#
def enhance_results(annotated_csv_path, output_path):
    write_report(annotated_csv_path, output_path, SECTION_COLUMNS, DUPLEX_RENDERERS, "annotated CSV")


#----- Snakemake entry point: unpacks the injected `snakemake` object and calls the pure logic above -----#
def run_from_snakemake(snakemake):
    enhance_results(
        annotated_csv_path=snakemake.input.annotated,
        output_path=snakemake.output.enhanced,
    )

if "snakemake" in globals():
    run_from_snakemake(snakemake)
