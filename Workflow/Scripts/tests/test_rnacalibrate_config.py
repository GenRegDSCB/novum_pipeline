from __future__ import annotations
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # Workflow/Scripts

import _rnacalibrate_config as rcc

LADDER = (76, 150, 300, 600, 900, 1500, 3000)

# This module's stdlib-only contract is pinned in test_import_contracts.py, with the other four.


# --- parse_anchors ---

def test_parse_anchors_accepts_the_comma_separated_string_form():
    assert rcc.parse_anchors("76,150,300,600,900,1500,3000") == LADDER

def test_parse_anchors_accepts_a_list_and_tolerates_whitespace():
    assert rcc.parse_anchors([76, 150, 300]) == (76, 150, 300)
    assert rcc.parse_anchors(" 76 , 150 ,300 ") == (76, 150, 300)

def test_parse_anchors_treats_absent_and_blank_as_no_ladder():
    assert rcc.parse_anchors(None) is None
    assert rcc.parse_anchors("") is None
    assert rcc.parse_anchors("   ") is None

@pytest.mark.parametrize("spec", ["76,abc,300", "76,-150,300", "76,0,300",
                                  "76,150,150,300", "76,300,150", "76,1.5,300"])
def test_parse_anchors_rejects_malformed_ladders(spec):
    with pytest.raises(ValueError):
        rcc.parse_anchors(spec)


# --- cell geometry ---

def test_cell_edges_are_geometric_midpoints():
    assert [round(e) for e in rcc.cell_edges(LADDER)] == [107, 212, 424, 735, 1162, 2121]

def test_assign_anchor_uses_the_geometric_midpoint_not_the_arithmetic_one():
    edges = rcc.cell_edges(LADDER)
    # The arithmetic midpoint of 600 and 900 is 750; the geometric one is 734.8.
    assert rcc.assign_anchor(734, LADDER, edges) == 600
    assert rcc.assign_anchor(736, LADDER, edges) == 900
    assert rcc.assign_anchor(340, LADDER, edges) == 300
    assert rcc.assign_anchor(746, LADDER, edges) == 900

def test_assign_anchor_at_the_top_cell_boundary():
    edges = rcc.cell_edges(LADDER)
    assert rcc.assign_anchor(2120, LADDER, edges) == 1500
    assert rcc.assign_anchor(2122, LADDER, edges) == 3000

def test_assign_anchor_puts_a_length_exactly_on_an_edge_in_the_upper_cell():
    # The shipped ladder's edges are all irrational, so no integer length can ever sit on
    # one — the documented tie convention is only reachable at an exact geometric midpoint.
    anchors = (100, 400)
    edges = rcc.cell_edges(anchors)
    assert edges == (200.0,)
    assert rcc.assign_anchor(199, anchors, edges) == 100
    assert rcc.assign_anchor(200, anchors, edges) == 400
    assert rcc.assign_anchor(201, anchors, edges) == 400

def test_assign_anchor_outer_cells_are_open_ended():
    edges = rcc.cell_edges(LADDER)
    assert rcc.assign_anchor(27, LADDER, edges) == 76
    assert rcc.assign_anchor(17676, LADDER, edges) == 3000

def test_cell_edges_of_a_single_anchor_ladder_is_empty():
    assert rcc.cell_edges((900,)) == ()
    assert rcc.assign_anchor(5, (900,), ()) == 900


# --- the fit ceiling ---

def test_fit_ceiling_is_fifteen_hundred():
    assert rcc.FIT_CEILING_NT == 1500

def test_fitted_anchors_stops_at_the_ceiling_inclusive():
    assert rcc.fitted_anchors(LADDER) == (76, 150, 300, 600, 900, 1500)

def test_default_ladder_declares_one_anchor_above_the_ceiling():
    assert rcc.DEFAULT_LENGTH_ANCHORS == LADDER
    assert len(rcc.fitted_anchors(rcc.DEFAULT_LENGTH_ANCHORS)) == 6


# --- validation ---

def _cfg(**overrides):
    base = {"length_anchors": "76,150,300,600,900,1500,3000",
            "randomize_targets": True, "rng_seed": 1}
    base.update(overrides)
    return base

def test_no_checks_fire_when_no_ladder_is_configured():
    assert rcc.validate_rnacalibrate_config({"randomize_targets": False}) == []

def test_a_fatal_when_the_ladder_is_set_without_randomize_targets():
    with pytest.raises(ValueError, match="randomize_targets"):
        rcc.validate_rnacalibrate_config(_cfg(randomize_targets=False))

def test_a_fatal_message_explains_that_every_anchor_would_be_identical():
    with pytest.raises(ValueError, match="identical"):
        rcc.validate_rnacalibrate_config(_cfg(randomize_targets=False))

def test_a_is_not_fatal_when_calibration_is_off():
    # No calibration runs, so the ladder is inert whatever randomize_targets says, and
    # (d) is the check that reports it.
    warnings = rcc.validate_rnacalibrate_config(_cfg(randomize_targets=False),
                                                calibration_on=False)
    assert any("wholly inert" in w for w in warnings)

def test_a_still_fires_when_calibration_is_explicitly_on():
    with pytest.raises(ValueError, match="randomize_targets"):
        rcc.validate_rnacalibrate_config(_cfg(randomize_targets=False), calibration_on=True)

def test_b_fatal_on_a_malformed_ladder():
    with pytest.raises(ValueError):
        rcc.validate_rnacalibrate_config(_cfg(length_anchors="76,150,abc"))

def test_b_fatal_when_fewer_than_three_anchors_sit_at_or_below_the_ceiling():
    # Only 76 and 150 are fitted, so Tier 3 never gets its 3 points. Also covers the
    # "every anchor above the ceiling" case.
    with pytest.raises(ValueError, match="at or below"):
        rcc.validate_rnacalibrate_config(_cfg(length_anchors="76,150,3000,6000"))
    with pytest.raises(ValueError, match="at or below"):
        rcc.validate_rnacalibrate_config(_cfg(length_anchors="2000,3000,6000"))

def test_c_warns_that_above_ceiling_anchors_are_extrapolated_never_fitted():
    warnings = rcc.validate_rnacalibrate_config(_cfg())
    assert any("3000" in w and "extrapolated" in w for w in warnings)

def test_c_is_silent_when_every_anchor_is_at_or_below_the_ceiling():
    warnings = rcc.validate_rnacalibrate_config(_cfg(length_anchors="76,150,300,600,900,1500"))
    assert not any("extrapolated" in w for w in warnings)

def test_d_warns_when_the_ladder_is_set_with_calibration_off():
    warnings = rcc.validate_rnacalibrate_config(_cfg(), calibration_on=False)
    assert any("calibration_variant" in w for w in warnings)

def test_d_is_silent_when_the_calibrated_arm_runs():
    assert not any("wholly inert" in w
                   for w in rcc.validate_rnacalibrate_config(_cfg(), calibration_on=True))

def test_d_reads_the_resolved_flag_not_the_raw_key():
    # YAML 1.1 parses a bare `off` as False, so the raw key is unusable here; the Snakefile
    # resolves it through variants_for and passes the result in.
    warnings = rcc.validate_rnacalibrate_config(_cfg(calibration_variant=False),
                                                calibration_on=True)
    assert not any("wholly inert" in w for w in warnings)

def test_e_warns_when_the_ladder_is_set_without_an_rng_seed():
    warnings = rcc.validate_rnacalibrate_config(_cfg(rng_seed=None))
    assert any("rng_seed" in w for w in warnings)

def test_e_is_silent_when_the_seed_is_zero():
    # 0 is a legitimate pinned seed; only None is unpinned.
    warnings = rcc.validate_rnacalibrate_config(_cfg(rng_seed=0))
    assert not any("rng_seed" in w for w in warnings)
