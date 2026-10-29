// Novum Pipeline — reference data

//----- Dropdown option lists. Only `value` reaches config.yaml; `label` is display-only -----//
window.SPECIES_OPTIONS = [
  { value: "3utr_human", label: "Human (3′UTR)" },
  { value: "3utr_fly",   label: "Fly (3′UTR)" },
  { value: "3utr_worm",  label: "Worm (3′UTR)" }
];

window.VARIANT_MODES = [
  { value: "on",   label: "On"   },
  { value: "off",  label: "Off"  },
  { value: "both", label: "Both" }
];

window.INTARNA_PREDICTION_MODES = [
  { value: "H", label: "H · heuristic" },
  { value: "M", label: "M · exact"     },
  { value: "S", label: "S · seed-only" }
];

window.INTARNA_MODELS = [
  { value: "X", label: "X · seed-extension MFE" },
  { value: "S", label: "S · legacy MFE"          },
  { value: "B", label: "B · helix-block"         },
  { value: "P", label: "P · ensemble"            }
];

window.INTARNA_ENERGY_SETS = [
  { value: "Turner99",     label: "Turner99"     },
  { value: "Turner04",     label: "Turner04"     },
  { value: "Andronescu07", label: "Andronescu07" }
];

window.INTARNA_OVERLAPS = [
  { value: "N", label: "N · none"          },
  { value: "T", label: "T · target only"   },
  { value: "Q", label: "Q · query only"    },
  { value: "B", label: "B · both"          }
];

window.PLOT_TYPES = [
  { value: "pvalue_distribution", label: "p-value distribution" },
  { value: "position_pvalue",     label: "position × p-value" },
  { value: "position_energy",     label: "position × energy" },
  { value: "volcano",             label: "volcano" },
  { value: "pvalue_ecdf",         label: "p-value ECDF" },
  { value: "calibration_delta",   label: "calibration Δ" },
  { value: "per_mirna",           label: "per-miRNA" },
  { value: "position_density",    label: "position density" },
  { value: "seed_class",          label: "seed class" }
];

window.PLOT_TYPE_VALUES = window.PLOT_TYPES.map((p) => p.value);

//----- React list keys. A counter, not the clock alone: two rows added inside the same
//      millisecond would otherwise share an id and React would reuse the wrong row -----//
let _idSeq = 0;
window.mintId = (prefix) => `${prefix}${(_idSeq++).toString(36)}${Date.now().toString(36)}`;

//----- Queries and Targets are one keyed picker over two directories; only the nouns and
//      the directory strings differ. `footerNote` is appended after the row count -----//
window.KEYED_FILE_SECTIONS = [
  { key: "queries", num: "1", title: "Queries", noun: "query", many: "queries",
    sub: "One query FASTA per sample (joined by key)", api: "list_queries",
    pickerHead: "Query FASTA · Data/Raw/RNAs/", pathHead: "Query FASTA path",
    placeholderPath: "Data/Raw/RNAs/validated_*.fa", emptyHint: "No files in Data/Raw/RNAs/",
    footerNote: " · keys must match target keys" },
  { key: "targets", num: "2", title: "Targets", noun: "target", many: "targets",
    sub: "One or more (CDS, RNAs, full genomes)", api: "list_targets",
    pickerHead: "Target genome · Data/Raw/genomes/", pathHead: "Target genome path",
    placeholderPath: "Data/Raw/genomes/*.fna", emptyHint: "No files in Data/Raw/genomes/",
    footerNote: "" },
];

//----- The three plots pair-lists — same (value, optional miRNA) editor, different nouns -----//
window.PLOT_PAIR_LISTS = [
  { key: "locus",   heading: "Locus tags", noun: "locus",   valueHeader: "Locus tag",
    valuePlaceholder: "locus tag (required), e.g. FE838_RS16060" },
  { key: "gene",    heading: "Genes",      noun: "gene",    valueHeader: "Gene",
    valuePlaceholder: "gene name (required), e.g. yegH" },
  { key: "protein", heading: "Proteins",   noun: "protein", valueHeader: "Protein",
    valuePlaceholder: "protein name (required)" },
];

//----- Ordered field specs, one per config section. kind: "n" nullable {set,…} field
//      (scalar or pair — the shape is read off the default), "p" plain value, "v" variant
//      string (YAML 1.1 turns bare on/off into booleans), "c" constant. Key ORDER is the
//      order written to config.yaml (safe_dump(sort_keys=False)) — keep it matching
//      Config/config.yaml. Both directions of config.js walk these, so a key listed here
//      can never be emitted without being hydrated back, or vice versa -----//
window.CONFIG_BLOCKS = {
  shared: [
    ["max_suboptimal_hits", "n"],
    ["seed", "n"],
  ],
  rnacalibrate: [
    ["calibration_variant", "v"],
    ["k", "p"],
    ["max_target_length", "p"],
    ["randomize_targets", "p"],
    ["rng_seed", "n"],
    ["length_anchors", "s"],
  ],
  rnahybrid: [
    ["species", "p"],
    ["max_hybrid_energy", "n"],
    ["max_internal_loop", "n"],
    ["max_bulge_loop", "n"],
    ["pvalue_threshold", "n"],
    ["distribution", "n"],
  ],
  // `type` and the locus/gene/protein lists have their own shapes and are read by hand.
  plots: [
    ["basesize", "p"],
    ["pvalue_threshold", "p"],
    ["per_mirna_top_n", "p"],
  ],
};

//----- The `intarna:` sub-blocks, same entry shape as CONFIG_BLOCKS above -----//
window.INTARNA_BLOCKS = {
  top: [
    ["accessibility_variant", "v"],
    ["prediction_mode", "p"],
    ["model", "p"],
    ["energy_set", "p"],
    ["max_hybrid_energy", "n"],
    ["min_target_unpaired_probability", "n"],
    ["min_query_unpaired_probability", "n"],
    ["accessibility_search_depth", "n"],
    ["max_interaction_length", "n"],
    ["max_loop_size", "n"],
  ],
  helix: [
    ["min_bp", "n"], ["max_bp", "n"], ["max_internal_loop", "n"],
    ["min_unpaired_probability", "n"], ["max_energy", "n"],
    ["full_energy", "p"],
  ],
  seed: [
    ["max_energy", "n"], ["max_hybrid_energy", "n"], ["min_unpaired_probability", "n"],
    ["forbid_gu", "p"], ["forbid_gu_at_ends", "p"],
    ["target_range", "n"], ["max_unpaired_bases", "n"],
    ["report_best_only", "p"],
  ],
  accessibility: [
    ["window", "n"], ["max_bp_span", "n"],
    ["query_window", "n"], ["query_max_bp_span", "n"],
    ["target_window", "n"], ["target_max_bp_span", "n"],
    ["forbid_lonely_pairs", "p"], ["forbid_gu_at_ends", "p"],
  ],
  output: [
    ["max_delta_energy", "n"], ["overlap", "p"], ["min_unpaired_probability", "n"],
    ["forbid_lonely_pairs", "p"], ["forbid_gu_at_ends", "p"],
    ["columns", "c"],
  ],
};

//----- The nullable keys of a spec, for the section renderers and the Hero "set params"
//      counter. Derived, never hand-listed: a key can't be in one list and not the other -----//
const nullableKeys = (spec) => spec.filter(([, kind]) => kind === "n").map(([k]) => k);

// On window: read across files (yaml.js's blank-line rule, config.js's emit loop).
window.SHARED_NULLABLE_KEYS    = nullableKeys(window.CONFIG_BLOCKS.shared);
window.RNAHYBRID_NULLABLE_KEYS = nullableKeys(window.CONFIG_BLOCKS.rnahybrid);

// Local: only countOptionalParams and helixInertKeys below read these.
const INTARNA_TOP_NULLABLE_KEYS    = nullableKeys(window.INTARNA_BLOCKS.top);
const INTARNA_HELIX_NULLABLE_KEYS  = nullableKeys(window.INTARNA_BLOCKS.helix);
const INTARNA_SEED_NULLABLE_KEYS   = nullableKeys(window.INTARNA_BLOCKS.seed);
const INTARNA_ACC_NULLABLE_KEYS    = nullableKeys(window.INTARNA_BLOCKS.accessibility);
const INTARNA_OUTPUT_NULLABLE_KEYS = nullableKeys(window.INTARNA_BLOCKS.output);

//----- `distribution` is supplied by RNAcalibrate only when calibration is the sole variant;
//      under "both" the wo_calibration arm still runs on the static distribution/species -----//
window.isDistributionForced = (cfg) => cfg.rnacalibrate.calibration_variant === "on";

//----- IntaRNA seed sub-fields only affect a run when a shared seed is set (intarna.py:
//      _add_seed_flags). Via resolveNullable, not `.set`: a toggled-on-but-blank seed saves
//      as null, and would otherwise enable sub-fields for a run getting --noSeed -----//
window.isSeedEnforced = (cfg) => window.resolveNullable(cfg.seed) != null;

//----- Helix parameters only reach IntaRNA under model B; under X/S/P they are inert.
//      Mirrors _intarna_config.py:INTARNA_DEFAULT_HELIX_MAX_BP -----//
window.INTARNA_DEFAULT_HELIX_MAX_BP = 10;
window.isHelixActive = (cfg) => cfg.intarna.model === "B";

//----- The accessibility block (and the per-side Pu floors) only bite when the variant is on -----//
window.isAccessibilityActive = (cfg) => cfg.intarna.accessibility_variant !== "off";

//----- IntaRNA's own per-side accessibility defaults, named once: INITIAL_CONFIG's
//      placeholders, the "Global" pill's reset and the AUTO (…) tags all read them here -----//
window.INTARNA_DEFAULT_ACC_WINDOW = 150;
window.INTARNA_DEFAULT_ACC_BP_SPAN = 100;

//----- Mirrors _intarna_config.py check (a): under model B the seed must fit inside one
//      helix. Returns the offending numbers, else null. Both fields go through
//      resolveNullable so a set-but-blank input collapses to null exactly as it does on save -----//
window.helixSeedConflict = (cfg) => {
  if (!window.isHelixActive(cfg)) return null;
  const seed = window.resolveNullable(cfg.seed);
  if (seed == null) return null;
  const [a, b] = String(seed).split(",").map(Number);
  const seedBP = b - a + 1;
  const maxBP = window.resolveNullable(cfg.intarna.helix.max_bp)
    ?? window.INTARNA_DEFAULT_HELIX_MAX_BP;
  return seedBP > maxBP ? { seedBP, maxBP } : null;
};

//----- Mirrors check (c). Inert values are still serialized — a model toggle must not
//      destroy the user's tuning — so the form names them instead of dropping them -----//
window.helixInertKeys = (cfg) => {
  if (window.isHelixActive(cfg)) return [];
  const keys = INTARNA_HELIX_NULLABLE_KEYS.filter((k) => cfg.intarna.helix[k].set);
  return cfg.intarna.helix.full_energy ? keys.concat("full_energy") : keys;
};

//----- The two accessibility sides — same window/span editor, different nouns and flags, in
//      the shape of PLOT_PAIR_LISTS above. Field keys derive as `${side}_window` /
//      `${side}_max_bp_span`, matching INTARNA_BLOCKS.accessibility and Config/config.yaml -----//
window.ACCESSIBILITY_SIDES = [
  { side: "query",  label: "Query (miRNA)", name: "Query",  wFlag: "--qAccW", lFlag: "--qAccL" },
  { side: "target", label: "Target (CDS)",  name: "Target", wFlag: "--tAccW", lFlag: "--tAccL" },
];

//----- Mirrors _intarna_config.py:ACCESSIBILITY_WINDOW_CONFLICTS — IntaRNA aborts if a
//      shared accessibility key is set beside a DIFFERING per-side one; equal or unset is
//      fine. test_ui_config_contract.py asserts this table equals the Python one -----//
window.ACCESSIBILITY_WINDOW_CONFLICTS = [
  ["window", "query_window"],
  ["window", "target_window"],
  ["max_bp_span", "query_max_bp_span"],
  ["max_bp_span", "target_max_bp_span"],
];

//----- The conflicting (shared, per-side) pairs in the current config, for the live note -----//
window.accessibilityConflicts = (cfg) => {
  const acc = cfg.intarna.accessibility;
  const conflicts = [];
  for (const [sharedKey, sideKey] of window.ACCESSIBILITY_WINDOW_CONFLICTS) {
    const shared = window.resolveNullable(acc[sharedKey]);
    const side = window.resolveNullable(acc[sideKey]);
    if (shared == null || side == null || shared === side) continue;
    conflicts.push({ sharedKey, sideKey, shared, side });
  }
  return conflicts;
};

//----- Count optional (nullable) params the user has set, against the number in play. Inert
//      groups (no shared seed, accessibility Off) and the locked distribution field are
//      excluded from BOTH sides, so the ratio only counts params that can affect this run -----//
window.countOptionalParams = (cfg) => {
  let set = 0, total = 0;
  const tally = (obj, keys, { active = true, locked = () => false } = {}) => {
    if (!active) return;
    for (const k of keys) {
      if (locked(k)) continue;
      total += 1;
      if (obj[k]?.set) set += 1;
    }
  };

  tally(cfg, window.SHARED_NULLABLE_KEYS);
  tally(cfg.rnahybrid, window.RNAHYBRID_NULLABLE_KEYS, {
    locked: (k) => k === "distribution" && window.isDistributionForced(cfg),
  });

  tally(cfg.intarna, INTARNA_TOP_NULLABLE_KEYS);
  tally(cfg.intarna.seed, INTARNA_SEED_NULLABLE_KEYS, {
    active: window.isSeedEnforced(cfg),
  });
  tally(cfg.intarna.helix, INTARNA_HELIX_NULLABLE_KEYS, {
    active: window.isHelixActive(cfg),
  });
  tally(cfg.intarna.accessibility, INTARNA_ACC_NULLABLE_KEYS, {
    active: window.isAccessibilityActive(cfg),
  });
  tally(cfg.intarna.output, INTARNA_OUTPUT_NULLABLE_KEYS);

  return { set, total };
};

//----- Default --outCsvCols string. The form never exposes `intarna.output.columns`, so the
//      writer emits this constant over whatever is on disk (kind "c") -----//
window.INTARNA_OUTPUT_COLUMNS_DEFAULT =
  "id1,id2,start1,end1,start2,end2,subseqDP,hybridDP,E,E_hybrid,ED1,ED2,Pu1,Pu2,seedStart1,seedEnd1,seedE,seedStart2,seedEnd2";

//----- Canonical initial state; every field absent from a loaded config keeps its value
//      here. Plots toggle on the size of `plots.types`: an empty list serializes to
//      `type: null` (skip), any selection emits the full plot block -----//
window.INITIAL_CONFIG = {
  queries: [],
  targets: [],
  threads: 16,
  max_suboptimal_hits: { set: false, value: 1 },
  seed:                { set: false, a: 2, b: 7 },
  hadMaxTotalEnergy: false,
  rnacalibrate: {
    calibration_variant: "on",
    k: 10000,
    max_target_length: 50000,
    randomize_targets: true,
    rng_seed: { set: false, value: 1 },   // null = don't pin the clock
    length_anchors: "76,150,300,600,900,1500,3000",
  },
  rnahybrid: {
    species: "3utr_human",
    max_hybrid_energy: { set: true, value: -18 },
    max_internal_loop: { set: false, value: 1 },
    max_bulge_loop:    { set: false, value: 1 },
    pvalue_threshold:  { set: false, value: 0.05 },
    distribution:      { set: false, a: 0, b: 1 }
  },
  intarna: {
    accessibility_variant: "on",
    prediction_mode: "H",
    model: "X",
    energy_set: "Turner04",
    max_hybrid_energy:          { set: true, value: -12.9 },
    // Per-side Pu floors, applied downstream in tidy_intarna (not IntaRNA flags) and
    // inert when the accessibility variant is off.
    min_target_unpaired_probability: { set: false, value: 0.001 },
    min_query_unpaired_probability:  { set: false, value: 0.001 },
    accessibility_search_depth: { set: true, value: 20 },
    // Placeholder `value`s mirror IntaRNA's own defaults, so flipping a field on is
    // initially a no-op. Except max_interaction_length, whose default (0) *means* auto.
    max_interaction_length: { set: false, value: 10 },
    max_loop_size:          { set: false, value: 10 },
    helix: {
      min_bp:                   { set: false, value: 2 },
      max_bp:                   { set: false, value: window.INTARNA_DEFAULT_HELIX_MAX_BP },
      max_internal_loop:        { set: false, value: 0 },
      min_unpaired_probability: { set: false, value: 0 },
      max_energy:               { set: false, value: 0 },
      full_energy: false
    },
    seed: {
      max_energy:               { set: false, value: 0 },
      max_hybrid_energy:        { set: false, value: 999 },
      min_unpaired_probability: { set: false, value: 0 },
      forbid_gu: false,
      forbid_gu_at_ends: false,
      target_range:             { set: false, value: "" },
      max_unpaired_bases:       { set: false, value: 0 },
      report_best_only: false
    },
    accessibility: {
      window:      { set: false, value: window.INTARNA_DEFAULT_ACC_WINDOW },
      max_bp_span: { set: false, value: window.INTARNA_DEFAULT_ACC_BP_SPAN },
      query_window:       { set: false, value: window.INTARNA_DEFAULT_ACC_WINDOW },
      query_max_bp_span:  { set: false, value: window.INTARNA_DEFAULT_ACC_BP_SPAN },
      target_window:      { set: false, value: window.INTARNA_DEFAULT_ACC_WINDOW },
      target_max_bp_span: { set: false, value: window.INTARNA_DEFAULT_ACC_BP_SPAN },
      forbid_lonely_pairs: false,
      forbid_gu_at_ends: false
    },
    output: {
      max_delta_energy:         { set: false, value: 100 },
      overlap: "B",
      min_unpaired_probability: { set: false, value: 0 },
      forbid_lonely_pairs: false,
      forbid_gu_at_ends: false
    }
  },
  plots: {
    types: window.PLOT_TYPES.map((p) => p.value),
    basesize: 12,
    pvalue_threshold: 0.01,
    per_mirna_top_n: 12,
    locus: [],
    gene: [],
    protein: []
  },
  results_dir: "Data/Results/"
};
