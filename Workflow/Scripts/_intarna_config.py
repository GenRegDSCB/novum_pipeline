"""IntaRNA config schema: defaults, seed derivation, and DAG-build-time validation.

Imported by BOTH `intarna.py` (inside the IntaRNA conda env) and the Snakefile (at
DAG-build time, in the driver env). That second consumer is why this module must stay
**dependency-free** — stdlib only: the driver env ships none of the arms' scientific
dependencies. Keeping it separate from `intarna.py` also means the driver never imports
a script that ends in `run_from_snakemake(snakemake)`.
"""

INTARNA_ENERGY_SETS = ("Turner99", "Turner04", "Andronescu07")
INTARNA_DEFAULT_ENERGY_SET = "Turner04"

# IntaRNA's own defaults, mirrored because validate_intarna_config reproduces two of its
# constraints.
INTARNA_DEFAULT_MODEL = "X"
INTARNA_DEFAULT_HELIX_MAX_BP = 10

# --outNumber accepts [0,1000] (0 = report nothing), so 1000 is how `max_suboptimal_hits:
# null` ("all hits") is spelled.
INTARNA_MAX_OUTNUMBER = 1000

# Ours, not IntaRNA's: the --outNumber floor on w_accessibility only, where acc=C ranks by
# total E while the gate is on E_hybrid.
INTARNA_DEFAULT_ACCESSIBILITY_SEARCH_DEPTH = 20

# Median of (IntaRNA E_hybrid - RNAhybrid mfe) over matched, UNFILTERED acc=N hits. It
# equalizes the gates under acc=N only; never re-derive it from w_accessibility output.
# See CLAUDE.md "Per-arm energy gates".
INTARNA_RNAHYBRID_ENERGY_OFFSET = 5.10

# Drift past this trips warning (d). A usability threshold, not a confidence bound — the
# median is known to SE ~0.03, so 1 kcal/mol is a deliberate choice, not noise.
ENERGY_OFFSET_TOLERANCE = 1.0

# Keys nothing reads any more, mapped to the message explaining the hand migration. An
# inert-but-present key is the exact failure this refactor ends, so presence is fatal.
REMOVED_TOP_LEVEL_KEYS = {
    "max_total_energy": (
        "max_total_energy has been removed. It filtered pure hybridization energy on the "
        "RNAhybrid arm but TOTAL interaction energy (E_hybrid + ED1 + ED2) on the IntaRNA "
        "arm, so one number meant two different bars — roughly a 2x stricter one for "
        "IntaRNA. There is no single correct automatic migration; replace it with the "
        "per-arm keys:\n"
        "    rnahybrid:\n"
        "      max_hybrid_energy: -18      # -e ; the literature MFE threshold\n"
        "    intarna:\n"
        "      max_hybrid_energy: -12.9    # -18 + 5.10 measured convention offset"
    ),
}


#----- The only `intarna:` keys build_command reads. Everything else in the block steers a
#      DOWNSTREAM rule (the tidy gate, the per-side Pu floors) or the DAG itself, and
#      Snakemake reruns on any change to a rule's params — so handing rule intarna the whole
#      block makes a tidy_intarna knob invalidate the hour-long arm. Filter through
#      tool_config() instead; a key added here must be one build_command actually consumes -----#
INTARNA_TOOL_KEYS = (
    "prediction_mode", "model", "energy_set", "max_interaction_length", "max_loop_size",
    "seed", "helix", "accessibility", "output",
)


#----- The tool-flag subset of the `intarna:` block, for rule intarna's params -----#
def tool_config(cfg):
    cfg = cfg or {}
    return {key: cfg[key] for key in INTARNA_TOOL_KEYS if key in cfg}


#----- Raises on a top-level key nothing reads any more. Called from the Snakefile at DAG-build time, before any rule runs -----#
def reject_removed_keys(config):
    for key, message in REMOVED_TOP_LEVEL_KEYS.items():
        if key in (config or {}):
            raise ValueError(message)


#----- cfg.get() that survives an explicit null: YAML `key:` parses to None, which a plain .get(key, default) would pass through and stringify into the flag as "None" -----#
def opt(cfg, key, default):
    value = cfg.get(key)
    return default if value is None else value


#----- Parses the shared top-level seed "x,y" into IntaRNA's seedBP+seedQRange shape; raises if the width is outside [2,20] -----#
def derive_seed_from_string(seed_str):
    if seed_str is None:
        return None
    try:
        parts = str(seed_str).split(",")
        if len(parts) != 2:
            raise ValueError
        x, y = int(parts[0].strip()), int(parts[1].strip())
    except ValueError:
        raise ValueError(f"top-level seed must be 'x,y' integers, got: {seed_str!r}")
    bp = y - x + 1
    if not (2 <= bp <= 20):
        raise ValueError(
            f"top-level seed '{seed_str}' gives seedBP={bp}, which is outside IntaRNA's [2,20] range"
        )
    return {"length": bp, "query_range": f"{x}-{y}"}


# IntaRNA's five seed columns, in its own spelling. The single source: intarna.py withholds
# them from --outCsvCols under --noSeed and tidy_intarna names its report block from them.
INTARNA_SEED_COLUMNS = ("seedStart1", "seedEnd1", "seedE", "seedStart2", "seedEnd2")


#----- The --outCsvCols list with the seed columns dropped when no seed is enforced. Under
#      --noSeed IntaRNA writes the literal NAN in all five on every row — ~12% of a
#      panel-scale CSV that is written, re-read and parsed to carry nothing -----#
def resolve_output_columns(columns, derived_seed):
    if not columns or derived_seed is not None:
        return columns
    kept = [name for name in str(columns).split(",")
            if name.strip() not in INTARNA_SEED_COLUMNS]
    return ",".join(kept)


#----- (d) The two arms' cutoffs must sit ~one convention offset apart, or the consensus silently intersects non-comparable bars -----#
def _energy_comparability_warnings(cfg, rnahybrid_max_hybrid_energy):
    intarna_cutoff = cfg.get("max_hybrid_energy")

    # Both unset: a deliberate "no gating anywhere" — nothing to warn about.
    if intarna_cutoff is None and rnahybrid_max_hybrid_energy is None:
        return []

    # Exactly one unset — a bigger mismatch than two cutoffs merely drifting apart, and
    # where a half-finished migration off max_total_energy lands.
    if intarna_cutoff is None or rnahybrid_max_hybrid_energy is None:
        null_key, null_arm, set_key, set_value = (
            ("intarna", "IntaRNA", "rnahybrid", rnahybrid_max_hybrid_energy)
            if intarna_cutoff is None
            else ("rnahybrid", "RNAhybrid", "intarna", intarna_cutoff)
        )
        return [
            f"{null_key}.max_hybrid_energy is null ({null_arm} arm ungated on "
            f"hybridization energy) while {set_key}.max_hybrid_energy is {set_value}. "
            f"The consensus intersects a gated arm with an ungated one."
        ]

    expected = rnahybrid_max_hybrid_energy + INTARNA_RNAHYBRID_ENERGY_OFFSET
    drift = abs(intarna_cutoff - expected)
    if drift <= ENERGY_OFFSET_TOLERANCE:
        return []
    return [
        f"intarna.max_hybrid_energy is {intarna_cutoff} while rnahybrid.max_hybrid_energy "
        f"is {rnahybrid_max_hybrid_energy}, which is {expected:.2f} on IntaRNA's E_hybrid "
        f"scale (measured offset {INTARNA_RNAHYBRID_ENERGY_OFFSET:+.2f} kcal/mol). The two "
        f"arms are gating at non-comparable bars, {drift:.2f} kcal/mol apart, which biases "
        f"the consensus toward whichever arm is looser."
    ]


# The accessibility window/span knobs, in Config/config.yaml order. Single source of the
# key→flag mapping: intarna.py emits from it, the conflict check below quotes from it.
# Spelled out rather than derived from a "q"/"t" rule so `grep -r qAccW` still finds it.
ACCESSIBILITY_WINDOW_FLAGS = (
    ("window", "--accW"),
    ("max_bp_span", "--accL"),
    ("query_window", "--qAccW"),
    ("query_max_bp_span", "--qAccL"),
    ("target_window", "--tAccW"),
    ("target_max_bp_span", "--tAccL"),
)

# Each shared accessibility key paired with the per-side counterparts that override it.
# Measured on 3.4.1: a shared value beside a DIFFERING per-side one aborts the run; equal
# values are fine. UI/src/data.js mirrors this table; a test asserts they agree.
ACCESSIBILITY_WINDOW_CONFLICTS = (
    ("window", "query_window"),
    ("window", "target_window"),
    ("max_bp_span", "query_max_bp_span"),
    ("max_bp_span", "target_max_bp_span"),
)


#----- (f) IntaRNA aborts on a shared accessibility window set beside a differing per-side one,
#      rather than letting either win; catch it before the hour-long RNAhybrid arm starts.
#      Deliberately NOT checked: `*_max_bp_span <= *_window` — IntaRNA already reports that
#      clearly per side, and a copy here would have to re-derive the side resolution -----#
def _reject_accessibility_conflicts(acc):
    flag_for = dict(ACCESSIBILITY_WINDOW_FLAGS)
    for shared_key, side_key in ACCESSIBILITY_WINDOW_CONFLICTS:
        shared, side = acc.get(shared_key), acc.get(side_key)
        if shared is None or side is None or shared == side:
            continue
        raise ValueError(
            f"intarna.accessibility.{shared_key} is {shared} and {side_key} is {side}. "
            f'IntaRNA rejects a shared window alongside a differing per-side one ("ERROR: '
            f'{flag_for[shared_key]} and {flag_for[side_key]} are set to different '
            f'values"). Set either {shared_key} or {side_key}, not both.'
        )


#----- (b) IntaRNA reads an unknown energy set as a parameter-FILE path, fails to load it, and exits 255 asking the user to file a bug report. Reject it here instead -----#
def _reject_unknown_energy_set(cfg):
    energy_set = opt(cfg, "energy_set", INTARNA_DEFAULT_ENERGY_SET)
    if energy_set not in INTARNA_ENERGY_SETS:
        raise ValueError(
            f"intarna.energy_set must be one of {list(INTARNA_ENERGY_SETS)}, got {energy_set!r}"
        )


#----- (a) Under model B the seed must fit inside a single helix or IntaRNA hard-errors — reachable by widening the seed or by lowering the cap beneath it. Inert under every other model -----#
def _reject_seed_wider_than_helix(model, helix, derived_seed, seed_str):
    if model != "B" or derived_seed is None:
        return
    max_bp = opt(helix, "max_bp", INTARNA_DEFAULT_HELIX_MAX_BP)
    if derived_seed["length"] > max_bp:
        raise ValueError(
            f"seed {seed_str!r} needs {derived_seed['length']} base pairs, but "
            f"intarna.helix.max_bp is {max_bp}. Under intarna.model: B the seed must "
            f"fit inside one helix. Raise intarna.helix.max_bp to at least "
            f"{derived_seed['length']}, or shorten the top-level seed."
        )


#----- (c) The --helix* flags are emitted unconditionally because IntaRNA ignores them cleanly under any model but B; report the mismatch instead of withholding them -----#
def _inert_helix_warnings(model, helix):
    if model == "B":
        return []
    # `0` is a legitimate value (max_internal_loop: 0 = pure stacks), so compare
    # against None/False rather than testing truthiness.
    inert = sorted(k for k, v in helix.items() if v is not None and v is not False)
    if not inert:
        return []
    return [
        f"intarna.helix {inert} is set but intarna.model is {model!r}: helix "
        f"parameters apply only to model B and will be ignored."
    ]


#----- (e) w_accessibility withholds --outMaxE, so the per-pair cap is all that bounds the intermediate CSV. accessibility_search_depth cannot stand in: it is a FLOOR on --outNumber, which a null cap ignores -----#
def _unbounded_output_warnings(max_suboptimal_hits, accessibility_on):
    if not accessibility_on or max_suboptimal_hits is not None:
        return []
    return [
        "intarna: w_accessibility withholds --outMaxE (its energy gate lives downstream "
        "in tidy_intarna.py), so with max_suboptimal_hits: null IntaRNA writes up to "
        f"{INTARNA_MAX_OUTNUMBER} hits/pair, bounded only by its own --outMaxE default "
        "of 0. Set max_suboptimal_hits to bound intarna_output.csv; "
        "intarna.accessibility_search_depth cannot — it is a floor on --outNumber, which "
        "a null cap ignores entirely."
    ]


# The per-side Pu floors, target side first. Single source of the key pair: the Snakefile
# builds one resolver per entry, and (g)/(h) below name them from it.
ACCESSIBILITY_PU_FLOOR_KEYS = ("min_target_unpaired_probability", "min_query_unpaired_probability")


#----- The per-side floors actually set, named for a message. Shared by (g) and (h) so the
#      two can never disagree about what counts as "set" -----#
def _floors_in_use(cfg):
    return ", ".join(
        "intarna." + key for key in ACCESSIBILITY_PU_FLOOR_KEYS if cfg.get(key) is not None
    )


#----- (g) --outMinPu and the per-side floors both bound accessibility, but on different quantities, and setting both compounds them silently -----#
def _duplicate_pu_floor_warnings(cfg):
    out_min_pu = (cfg.get("output") or {}).get("min_unpaired_probability")
    named = _floors_in_use(cfg)
    if out_min_pu is None or not named:
        return []
    return [
        f"intarna: output.min_unpaired_probability ({out_min_pu}) requires EVERY interacting "
        f"position to reach it, on both sides, at the tool; {named} floors the whole site's "
        "Pu1/Pu2 per side in tidy_intarna. Different bars, and they compound — a site must "
        "clear both."
    ]


#----- (h) The floors are withheld under acc=N, so with accessibility off they are wholly inert while still reading as live config -----#
def _inert_accessibility_floor_warnings(cfg, accessibility_on):
    named = _floors_in_use(cfg)
    if accessibility_on or not named:
        return []
    return [
        f"intarna: {named} set under accessibility_variant: off. IntaRNA computes no "
        "unpaired probabilities under acc=N, so these floors are ignored entirely. Set "
        "accessibility_variant to 'on' or 'both' to make them take effect."
    ]


#----- DAG-build-time validation: raises on a config IntaRNA would die on, returns advisory warnings.
#      Each check is a lettered helper above — (a)/(b)/(f) raise, (c)/(d)/(e)/(g)/(h) return warnings — so
#      adding one means writing a helper and a line here, never growing this body -----#
def validate_intarna_config(cfg, seed_str, rnahybrid_max_hybrid_energy=None,
                            max_suboptimal_hits=None, accessibility_on=False):
    cfg = cfg or {}
    helix = cfg.get("helix") or {}
    model = opt(cfg, "model", INTARNA_DEFAULT_MODEL)

    # Every fatal runs before any warning is collected — a warning list a raise discards
    # is wasted work.
    _reject_unknown_energy_set(cfg)
    _reject_accessibility_conflicts(cfg.get("accessibility") or {})
    # The only validation of the top-level seed string, so it must run for every model.
    derived_seed = derive_seed_from_string(seed_str)
    _reject_seed_wider_than_helix(model, helix, derived_seed, seed_str)

    return (
        _energy_comparability_warnings(cfg, rnahybrid_max_hybrid_energy)
        + _unbounded_output_warnings(max_suboptimal_hits, accessibility_on)
        + _inert_helix_warnings(model, helix)
        + _duplicate_pu_floor_warnings(cfg)
        + _inert_accessibility_floor_warnings(cfg, accessibility_on)
    )
