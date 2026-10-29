"""RNAcalibrate ladder schema: cell geometry and DAG-build-time validation.

Stdlib only: the Snakefile imports it into the driver env, which ships none of the arms'
scientific deps. `rnacalibrate.py` and `rnahybrid.py` import it too.

Mirrors `_intarna_config.py`: lettered helpers, fatals first, then warnings concatenated.
"""

import bisect
import math

# An anchor is a cell centre AND a fit point; they separate here. Above 1500 nt fits are
# unreliable even when they converge (theta non-monotone, residuals 5.2% -> 17.3%, -nan at
# 2100/3000/6000), so higher anchors are declared but extrapolated. Measured property of
# this binary, not a preference — hence a constant, not a config key.
FIT_CEILING_NT = 1500

# Seven cells, six fits. 900 is off-geometric on purpose: 735-1162 is the densest cell.
DEFAULT_LENGTH_ANCHORS = (76, 150, 300, 600, 900, 1500, 3000)

# Below three points alpha is undetermined.
MIN_FITTED_ANCHORS = 3


#----- Ladder from a comma-separated string or a YAML list; None = no ladder, feature off -----#
def parse_anchors(spec):
    if spec is None:
        return None
    if isinstance(spec, str):
        text = spec.strip()
        if not text:
            return None
        parts = [part.strip() for part in text.split(",")]
    else:
        parts = [str(part).strip() for part in spec]
        if not parts:
            return None

    anchors = []
    for part in parts:
        try:
            value = int(part)
        except (TypeError, ValueError):
            raise ValueError(
                f"rnacalibrate.length_anchors: {part!r} is not a whole number of "
                f"nucleotides (got {spec!r})."
            ) from None
        if value <= 0:
            raise ValueError(
                f"rnacalibrate.length_anchors: anchors must be positive lengths, got {value}."
            )
        anchors.append(value)

    for lower, upper in zip(anchors, anchors[1:]):
        if upper <= lower:
            raise ValueError(
                "rnacalibrate.length_anchors must be strictly ascending with no "
                f"duplicates, got {anchors}."
            )
    return tuple(anchors)


#----- Cell boundaries: geometric midpoints, since assignment is nearest in log space -----#
def cell_edges(anchors):
    return tuple(math.sqrt(lower * upper) for lower, upper in zip(anchors, anchors[1:]))


#----- The cell this length belongs to. Outermost cells are OPEN-ENDED; a length exactly on
#      an edge joins the upper anchor -----#
def assign_anchor(length, anchors, edges):
    return anchors[bisect.bisect_right(edges, length)]


#----- Anchors RNAcalibrate is actually invoked at; the rest are extrapolated -----#
def fitted_anchors(anchors):
    return tuple(anchor for anchor in anchors if anchor <= FIT_CEILING_NT)


#----- (a) Without -s there is no RNG and -l is dropped entirely (verified inert), so every
#      anchor would hold bit-identical numbers. Fatal only if calibration runs, else (d) -----#
def _reject_ladder_without_shuffling(anchors, cfg, calibration_on):
    if calibration_on and not cfg.get("randomize_targets"):
        raise ValueError(
            "rnacalibrate.length_anchors is set but randomize_targets is false. Without "
            "-s there is no RNG and RNAcalibrate drops -l entirely, so every anchor would "
            "produce identical xi/theta and the per-anchor tree would be meaningless. Set "
            "randomize_targets: true, or remove length_anchors."
        )


#----- (b) The one shape constraint that depends on the ceiling; parse_anchors covers the rest -----#
def _reject_too_few_fitted_anchors(anchors):
    fitted = fitted_anchors(anchors)
    if len(fitted) < MIN_FITTED_ANCHORS:
        raise ValueError(
            f"rnacalibrate.length_anchors needs at least {MIN_FITTED_ANCHORS} anchors at "
            f"or below the fit ceiling of {FIT_CEILING_NT} nt; got {list(fitted)} out of "
            f"{list(anchors)}. Anchors above the ceiling are never fitted, so they cannot "
            "contribute to the per-miRNA alpha that every extrapolated cell is read from."
        )


#----- (c) Expected on the shipped ladder: the top anchor gives the long tail its own cell -----#
def _above_ceiling_warnings(anchors):
    above = [anchor for anchor in anchors if anchor > FIT_CEILING_NT]
    if not above:
        return []
    return [
        f"rnacalibrate: anchors {above} sit above the fit ceiling of {FIT_CEILING_NT} nt. "
        "They will be extrapolated from each miRNA's own alpha, never fitted — fits above "
        "the ceiling are measured unreliable on this binary. This is expected on the "
        "default ladder, which declares 3000 precisely so the long tail gets its own cell."
    ]


#----- (d) The block is inert when the arm never runs. Takes the RESOLVED flag, not the raw
#      key: YAML 1.1 parses a bare `off` as False, which no string compare here would catch -----#
def _inert_ladder_warnings(calibration_on):
    if not calibration_on:
        return [
            "rnacalibrate: length_anchors is set but calibration_variant is 'off', so no "
            "calibration runs and the ladder is wholly inert."
        ]
    return []


#----- (e) Re-draws still work unpinned, but the retry COUNT varies run to run, so two runs
#      of one config can extrapolate different cells -----#
def _unpinned_ladder_warnings(cfg):
    if cfg.get("rng_seed") is None:
        return [
            "rnacalibrate: length_anchors is set with rng_seed: null. The Tier 2 re-draws "
            "still work, but which anchors need them is no longer reproducible, so two "
            "runs of this config can differ in how many cells were extrapolated."
        ]
    return []


#----- DAG-build-time validation. (a)/(b) raise, (c)/(d)/(e) return warnings; adding a check
#      means a helper plus one line here, never growing this body -----#
def validate_rnacalibrate_config(cfg, calibration_on=True):
    cfg = cfg or {}
    anchors = parse_anchors(cfg.get("length_anchors"))
    # No ladder means every check below is inert, so each helper can assume one exists.
    if not anchors:
        return []

    # Every fatal runs before any warning is collected — a warning list a raise discards
    # is wasted work.
    _reject_ladder_without_shuffling(anchors, cfg, calibration_on)
    _reject_too_few_fitted_anchors(anchors)

    return (
        _above_ceiling_warnings(anchors)
        + _inert_ladder_warnings(calibration_on)
        + _unpinned_ladder_warnings(cfg)
    )
