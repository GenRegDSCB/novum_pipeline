from _report import write_report, render_intarna_duplex, render_rnahybrid_duplex
from tidy_intarna import ENERGY_COLUMNS as INTARNA_ENERGY_COLUMNS, SEED_COLUMNS
from tidy_rnahybrid import ENERGY_COLUMNS as RNAHYBRID_ENERGY_COLUMNS


SECTION_COLUMNS  = [RNAHYBRID_ENERGY_COLUMNS, INTARNA_ENERGY_COLUMNS, SEED_COLUMNS]
DUPLEX_RENDERERS = [render_rnahybrid_duplex, render_intarna_duplex]


#----- Per-record report: merged metadata, both arms' scoring blocks, then both duplexes -----#
def enhance_results(annotated_csv_path, output_path):
    write_report(annotated_csv_path, output_path, SECTION_COLUMNS, DUPLEX_RENDERERS, "consensus CSV")


#----- Snakemake entry point: unpacks the injected `snakemake` object and calls the pure logic above -----#
def run_from_snakemake(snakemake):
    enhance_results(
        annotated_csv_path=snakemake.input.annotated,
        output_path=snakemake.output.enhanced,
    )

if "snakemake" in globals():
    run_from_snakemake(snakemake)
