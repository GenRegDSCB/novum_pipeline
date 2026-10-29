import os
import sys

#----- Handle to the config file -----#
configfile: "Config/config.yaml"

#----- Stdlib-only config validators, imported into the driver env so a bad config fails at
#      DAG-build time. Never import a script that ends in run_from_snakemake() here -----#
sys.path.insert(0, os.path.join(workflow.basedir, "Workflow", "Scripts"))
from _intarna_config import (
    ACCESSIBILITY_PU_FLOOR_KEYS,
    INTARNA_DEFAULT_ACCESSIBILITY_SEARCH_DEPTH,
    opt,
    reject_removed_keys,
    tool_config,
    validate_intarna_config,
)
from _rnacalibrate_config import validate_rnacalibrate_config

#----- Config values: opt() wherever a default exists (an explicit YAML null would flow on as
#      None), plain .get() where null is meaningful -----#
targets             = config["targets"]
rnacalibrate_config = config.get("rnacalibrate", {}) or {}
rnahybrid_config    = config.get("rnahybrid", {}) or {}
intarna_config      = config.get("intarna", {}) or {}
results_dir         = opt(config, "results_dir", "Data/Results").rstrip("/")
max_target_length   = opt(rnacalibrate_config, "max_target_length", 50000)
shared_threads             = int(opt(config, "threads", 1))
shared_max_suboptimal_hits = config.get("max_suboptimal_hits")
shared_seed                = config.get("seed")
rnahybrid_max_hybrid_energy = rnahybrid_config.get("max_hybrid_energy")
intarna_max_hybrid_energy   = intarna_config.get("max_hybrid_energy")

SAMPLE_DIR         = results_dir + "/{sample}"
SAMPLE_VARIANT_DIR = SAMPLE_DIR + "/{variant}"

#----- `queries:` keyed like `targets:`, or a legacy single `query:` broadcast to every target -----#
def _resolve_queries(cfg, target_keys):
    queries_cfg = cfg.get("queries")
    if queries_cfg is not None:
        missing = set(target_keys) - set(queries_cfg)
        extra   = set(queries_cfg) - set(target_keys)
        if missing or extra:
            raise ValueError(
                "queries/targets keys must match. "
                f"Missing in queries: {sorted(missing)}; extra in queries: {sorted(extra)}."
            )
        return dict(queries_cfg)
    legacy = cfg.get("query")
    if legacy is None:
        raise ValueError("Config must define either `queries:` (mapping) or `query:` (single path).")
    return {sample: legacy for sample in target_keys}

queries = _resolve_queries(config, targets.keys())


#----- Per-sample miRNA query and CDS target FASTAs -----#
def sample_query(wc):
    return queries[wc.sample]

def sample_target(wc):
    return targets[wc.sample]

#----- Optional plotting configuration -----#
plots_config  = config.get("plots", {}) or {}
plots_type    = plots_config.get("type")
plots_enabled = plots_type is not None
plots_slug    = (plots_type or "none").replace(",", "-").replace(" ", "")

#----- Each arm's variant axis as a (with, without) pair of output directory names -----#
W_CALIBRATION    = "w_calibration"
WO_CALIBRATION   = "wo_calibration"
W_ACCESSIBILITY  = "w_accessibility"
WO_ACCESSIBILITY = "wo_accessibility"
CALIBRATION_LITERALS   = (W_CALIBRATION,   WO_CALIBRATION)
ACCESSIBILITY_LITERALS = (W_ACCESSIBILITY, WO_ACCESSIBILITY)
VARIANT_LABELS = {W_CALIBRATION: "calibrated", WO_CALIBRATION: "uncalibrated"}

#----- An arm's on/off/both key (bare YAML on/off arrives as a bool) to its directory literals -----#
def variants_for(cfg, key, label, literals):
    with_literal, without_literal = literals
    mode = cfg.get(key, "both")
    if isinstance(mode, bool):
        mode = "on" if mode else "off"
    modes = {"on": [with_literal], "off": [without_literal],
             "both": [with_literal, without_literal]}
    if mode not in modes:
        raise ValueError(f"{label} must be one of {list(modes)}, got {mode!r}")
    return modes[mode]

variants = variants_for(rnacalibrate_config, "calibration_variant",
                        "rnacalibrate.calibration_variant", CALIBRATION_LITERALS)
intarna_variants = variants_for(intarna_config, "accessibility_variant",
                                "intarna.accessibility_variant", ACCESSIBILITY_LITERALS)

#----- A removed config key is fatal: it is inert but looks live -----#
reject_removed_keys(config)

#----- Validate both tool blocks: fatals raise, warnings are logged -----#
for _warning in validate_intarna_config(
    intarna_config, shared_seed,
    rnahybrid_max_hybrid_energy=rnahybrid_max_hybrid_energy,
    max_suboptimal_hits=shared_max_suboptimal_hits,
    accessibility_on=W_ACCESSIBILITY in intarna_variants,
):
    logger.warning(_warning)

for _warning in validate_rnacalibrate_config(
    rnacalibrate_config,
    calibration_on=W_CALIBRATION in variants,
):
    logger.warning(_warning)

#----- --outNumber floor for the accessibility arm -----#
intarna_search_depth = opt(intarna_config, "accessibility_search_depth",
                           INTARNA_DEFAULT_ACCESSIBILITY_SEARCH_DEPTH)


RNAHYBRID_VARIANT_RE = "|".join(CALIBRATION_LITERALS)
INTARNA_VARIANT_RE   = "|".join(ACCESSIBILITY_LITERALS)

#----- Shared by rule rnacalibrate's output and calibration_input -----#
RNACALIBRATE_JSON = SAMPLE_DIR + "/" + W_CALIBRATION + "/rnacalibrate.json"


#----- Every {variant} is one of the four arm literals; each rule narrows to its own arm's two -----#
wildcard_constraints:
    variant = f"{RNAHYBRID_VARIANT_RE}|{INTARNA_VARIANT_RE}"


#----- The calibration JSON for w_calibration only. Keep `[]`, not None: rnahybrid.py reads
#      an empty input as the uncalibrated branch -----#
def calibration_input(wc):
    if wc.variant == W_CALIBRATION:
        return RNACALIBRATE_JSON.format(sample=wc.sample)
    return []


#----- species/distribution reach wo_calibration only: w_calibration ignores both, and params
#      are a rerun trigger -----#
def _uncalibrated_only(value):
    return lambda wc: None if wc.variant == W_CALIBRATION else value

rnahybrid_species      = _uncalibrated_only(rnahybrid_config.get("species", "3utr_human"))
rnahybrid_distribution = _uncalibrated_only(rnahybrid_config.get("distribution"))


#----- The predicate every IntaRNA resolver below branches on -----#
def _accessibility_on(wc):
    return wc.variant == W_ACCESSIBILITY


#----- {variant} to IntaRNA's --acc mode (C = accessibility, N = none) -----#
def intarna_acc_mode(wc):
    return "C" if _accessibility_on(wc) else "N"


#----- --outMaxE bounds total E, which equals E_hybrid only under acc=N; under acc=C
#      tidy_intarna gates on E_hybrid instead -----#
def intarna_outmaxe(wc):
    return None if _accessibility_on(wc) else intarna_max_hybrid_energy


#----- acc=C ranks by total E, not the gated E_hybrid, so report deeper there. A floor, never
#      a ceiling: null still means all hits -----#
def intarna_outnumber(wc):
    if not _accessibility_on(wc) or shared_max_suboptimal_hits is None:
        return shared_max_suboptimal_hits
    return max(intarna_search_depth, shared_max_suboptimal_hits)


#----- Per-side Pu floors, withheld under acc=N, which computes no Pu -----#
def _pu_floor_resolver(key):
    return lambda wc: intarna_config.get(key) if _accessibility_on(wc) else None

intarna_min_pu_target, intarna_min_pu_query = (
    _pu_floor_resolver(key) for key in ACCESSIBILITY_PU_FLOOR_KEYS
)


#----- Consensus compares each arm's w_* variant, or wo_* when only that ran -----#
def _preferred_variant(available, literals):
    with_literal, without_literal = literals
    return with_literal if with_literal in available else without_literal

CONSENSUS_RNAHYBRID_VARIANT = _preferred_variant(variants, CALIBRATION_LITERALS)
CONSENSUS_INTARNA_VARIANT   = _preferred_variant(intarna_variants, ACCESSIBILITY_LITERALS)


#----- output finale ‒ pipeline conclusion -----#
rule all:
    input:
        expand(
            f"{SAMPLE_VARIANT_DIR}/rnahybrid_enhanced.txt",
            sample=targets.keys(),
            variant=variants,
        ),
        expand(
            f"{SAMPLE_VARIANT_DIR}/intarna_enhanced.txt",
            sample=targets.keys(),
            variant=intarna_variants,
        ),
        expand(
            f"{SAMPLE_DIR}/consensus/consensus_enhanced.txt",
            sample=targets.keys(),
        ),
        *(expand(
            f"{SAMPLE_DIR}/plots_{plots_slug}.pdf",
            sample=targets.keys(),
        ) if plots_enabled else [])

#----- Dynamically calibrates the statistics based on the target sequence -----#
rule rnacalibrate:
    input:
        query=sample_query,
        target=sample_target
    output:
        calibration=RNACALIBRATE_JSON
    conda:
        "Workflow/Envs/rnahybrid.yaml"
    threads:
        # per-miRNA calls are independent and run concurrently
        shared_threads
    params:
        k=opt(rnacalibrate_config, "k", 10000),
        max_target_length=max_target_length,
        randomize_targets=opt(rnacalibrate_config, "randomize_targets", False),
        rng_seed=rnacalibrate_config.get("rng_seed"),
        # a ladder edit re-runs this rule, and rnahybrid through the JSON input
        length_anchors=rnacalibrate_config.get("length_anchors"),
        max_internal_loop=rnahybrid_config.get("max_internal_loop"),
        max_bulge_loop=rnahybrid_config.get("max_bulge_loop"),
        seed=shared_seed
    script:
        "Workflow/Scripts/rnacalibrate.py"

#----- Run RNAhybrid on multiple cores -----#
rule rnahybrid:
    wildcard_constraints:
        variant = RNAHYBRID_VARIANT_RE
    input:
        query=sample_query,
        target=sample_target,
        calibration=calibration_input
    output:
        compact=f"{SAMPLE_VARIANT_DIR}/rnahybrid_output.tsv"
    conda:
        "Workflow/Envs/rnahybrid.yaml"
    params:
        species=rnahybrid_species,
        max_suboptimal_hits=shared_max_suboptimal_hits,
        max_internal_loop=rnahybrid_config.get("max_internal_loop"),
        max_bulge_loop=rnahybrid_config.get("max_bulge_loop"),
        max_hybrid_energy=rnahybrid_max_hybrid_energy,
        pvalue_threshold=rnahybrid_config.get("pvalue_threshold"),
        seed=shared_seed,
        distribution=rnahybrid_distribution,
        max_target_length=max_target_length
    threads:
        shared_threads
    script:
        "Workflow/Scripts/rnahybrid.py"

#----- Filter, select, arrange the RNAhybrid's output -----#
rule tidy_rnahybrid:
    wildcard_constraints:
        variant = RNAHYBRID_VARIANT_RE
    input:
        compact=f"{SAMPLE_VARIANT_DIR}/rnahybrid_output.tsv"
    output:
        tidy=f"{SAMPLE_VARIANT_DIR}/tidy_output.csv"
    conda:
        "Workflow/Envs/postprocess.yaml"
    script:
        "Workflow/Scripts/tidy_rnahybrid.py"

#----- Parse metadata from the target file (genome) and merge them with the tidied output -----#
rule annotate_rnahybrid:
    wildcard_constraints:
        variant = RNAHYBRID_VARIANT_RE
    input:
        tidy=f"{SAMPLE_VARIANT_DIR}/tidy_output.csv",
        target=sample_target
    output:
        annotated=f"{SAMPLE_VARIANT_DIR}/rnahybrid_annotated.csv"
    conda:
        "Workflow/Envs/postprocess.yaml"
    params:
        insert_after="P_value"
    script:
        "Workflow/Scripts/annotate.py"

#----- Visualize RNAhybrid alignments + all the data -----#
rule enhance_rnahybrid:
    wildcard_constraints:
        variant = RNAHYBRID_VARIANT_RE
    input:
        annotated=f"{SAMPLE_VARIANT_DIR}/rnahybrid_annotated.csv"
    output:
        enhanced=f"{SAMPLE_VARIANT_DIR}/rnahybrid_enhanced.txt"
    conda:
        "Workflow/Envs/postprocess.yaml"
    script:
        "Workflow/Scripts/enhance_rnahybrid.py"

#----- Run IntaRNA on multiple cores (native --threads, no external chunking) -----#
rule intarna:
    wildcard_constraints:
        variant = INTARNA_VARIANT_RE
    input:
        query=sample_query,
        target=sample_target
    output:
        csv=f"{SAMPLE_VARIANT_DIR}/intarna_output.csv"
    conda:
        "Workflow/Envs/intarna.yaml"
    threads:
        shared_threads
    params:
        # tool keys only: a downstream-only key here would rerun IntaRNA on every edit
        intarna           = tool_config(intarna_config),
        acc                 = intarna_acc_mode,
        max_suboptimal_hits = intarna_outnumber,
        out_max_energy      = intarna_outmaxe,
        seed                = shared_seed
    script:
        "Workflow/Scripts/intarna.py"

#----- Filter, rename, and normalize the IntaRNA CSV output -----#
rule tidy_intarna:
    wildcard_constraints:
        variant = INTARNA_VARIANT_RE
    input:
        csv=f"{SAMPLE_VARIANT_DIR}/intarna_output.csv",
        target=sample_target
    output:
        tidy=f"{SAMPLE_VARIANT_DIR}/intarna_tidy.csv"
    resources:
        mem_mb = 1200
    conda:
        # needs pandas, never the IntaRNA binary
        "Workflow/Envs/postprocess.yaml"
    params:
        # the authoritative E_hybrid gate on both variants
        max_hybrid_energy   = intarna_max_hybrid_energy,
        max_suboptimal_hits = shared_max_suboptimal_hits,
        # withheld under acc=N; applied before the cap
        min_target_unpaired_probability = intarna_min_pu_target,
        min_query_unpaired_probability  = intarna_min_pu_query
    script:
        "Workflow/Scripts/tidy_intarna.py"

#----- Parse metadata from the target FASTA and merge with the IntaRNA tidy output -----#
rule annotate_intarna:
    wildcard_constraints:
        variant = INTARNA_VARIANT_RE
    input:
        tidy=f"{SAMPLE_VARIANT_DIR}/intarna_tidy.csv",
        target=sample_target
    output:
        annotated=f"{SAMPLE_VARIANT_DIR}/intarna_annotated.csv"
    conda:
        "Workflow/Envs/postprocess.yaml"
    params:
        insert_after="E"
    script:
        "Workflow/Scripts/annotate.py"

#----- Visualize IntaRNA alignments and energy data as a per-record human-readable report -----#
rule enhance_intarna:
    wildcard_constraints:
        variant = INTARNA_VARIANT_RE
    input:
        annotated=f"{SAMPLE_VARIANT_DIR}/intarna_annotated.csv"
    output:
        enhanced=f"{SAMPLE_VARIANT_DIR}/intarna_enhanced.txt"
    conda:
        "Workflow/Envs/postprocess.yaml"
    script:
        "Workflow/Scripts/enhance_intarna.py"

#----- Intersect the two arms' annotated tables into the consensus (both-tools) set -----#
rule intersect:
    input:
        rnahybrid=f"{SAMPLE_DIR}/{CONSENSUS_RNAHYBRID_VARIANT}/rnahybrid_annotated.csv",
        intarna=f"{SAMPLE_DIR}/{CONSENSUS_INTARNA_VARIANT}/intarna_annotated.csv"
    output:
        annotated=f"{SAMPLE_DIR}/consensus/consensus_annotated.csv"
    conda:
        "Workflow/Envs/postprocess.yaml"
    script:
        "Workflow/Scripts/intersect.py"

#----- Render the consensus set as a per-record report: merged metadata + both duplexes -----#
rule enhance_consensus:
    input:
        annotated=f"{SAMPLE_DIR}/consensus/consensus_annotated.csv"
    output:
        enhanced=f"{SAMPLE_DIR}/consensus/consensus_enhanced.txt"
    conda:
        "Workflow/Envs/postprocess.yaml"
    script:
        "Workflow/Scripts/enhance_consensus.py"

#----- Build configurable ggplot2 plots from the annotated tables -----#
rule build_plots:
    input:
        annotated=expand(
            f"{SAMPLE_VARIANT_DIR}/rnahybrid_annotated.csv",
            variant=variants,
            allow_missing=True,
        )
    output:
        pdf=f"{SAMPLE_DIR}/plots_{plots_slug}.pdf"
    conda:
        "Workflow/Envs/plots.yaml"
    params:
        # opt() because the shipped config writes `locus:`/`gene:`/`protein:` as bare nulls.
        type=plots_type,
        basesize=opt(plots_config, "basesize", 12),
        pvalue_threshold=opt(plots_config, "pvalue_threshold", []),
        locus=opt(plots_config, "locus", []),
        gene=opt(plots_config, "gene", []),
        protein=opt(plots_config, "protein", []),
        per_mirna_top_n=opt(plots_config, "per_mirna_top_n", 12),
        variant_labels=[VARIANT_LABELS[v] for v in variants]
    script:
        "Workflow/Scripts/build_plots.R"
