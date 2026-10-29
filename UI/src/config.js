// Novum Pipeline — config <-> form-state bridge.
//
// `configToObject` is the sole writer of the config schema: api.save_config dumps
// it, and yaml.js renders the preview from it. Both directions walk the same
// window.CONFIG_BLOCKS / INTARNA_BLOCKS specs, so a key cannot be written without
// being read back — the failure that once dropped `rnacalibrate.rng_seed` silently.

//----- Nullable form field -> serialized value. Pair-shaped fields (`a`/`b`) emit "a,b";
//      unset/forced/blank fields collapse to null. data.js's derived warnings route through
//      it too, so they cannot contradict what Save writes -----//
function resolveNullable(field, forced = false) {
  if (forced || !field.set) return null;
  if ("a" in field && "b" in field) {
    const { a, b } = field;
    return a === "" || b === "" || a == null || b == null ? null : `${a},${b}`;
  }
  return field.value === "" || field.value == null ? null : field.value;
}
window.resolveNullable = resolveNullable;

//----- One config block, emitted from its spec. `forced(key)` collapses a field to null
//      whatever its toggle says (RNAhybrid's distribution, once calibration supplies it).
//      "s" is a plain string whose blank is a meaningful null, not a missing value -----//
function _block(state, spec, forced = () => false) {
  return Object.fromEntries(spec.map(([key, kind]) => [
    key,
    kind === "n" ? resolveNullable(state[key], forced(key))
    : kind === "s" ? (String(state[key] ?? "").trim() || null)
    : kind === "c" ? window.INTARNA_OUTPUT_COLUMNS_DEFAULT
    : state[key],
  ]));
}

//----- Rows whose value is blank are dropped, so a lone miRNA never reaches the config -----//
function _pairListToArray(items) {
  return (items || [])
    .map((it) => ({ v: (it.value || "").trim(), m: (it.mirna || "").trim() }))
    .filter((it) => it.v)
    .map((it) => (it.m ? `${it.v},${it.m}` : it.v));
}

//----- form-state cfg -> plain config object (what launcher.save_config/run_pipeline receive) -----//
window.configToObject = function configToObject(cfg) {
  const obj = {};
  const B = window.CONFIG_BLOCKS;

  obj.queries = {};
  cfg.queries.forEach((q) => { obj.queries[(q.key || "").trim() || "unnamed_query"] = q.path || ""; });
  obj.targets = {};
  cfg.targets.forEach((t) => { obj.targets[(t.key || "").trim() || "unnamed_target"] = t.path || ""; });

  // The shared nullables ride at the top level, beside `threads`.
  obj.threads = cfg.threads;
  Object.assign(obj, _block(cfg, B.shared));

  obj.rnacalibrate = _block(cfg.rnacalibrate, B.rnacalibrate);

  const distForced = window.isDistributionForced(cfg);
  obj.rnahybrid = _block(cfg.rnahybrid, B.rnahybrid,
                         (key) => key === "distribution" && distForced);

  // Helix values are written even while the form greys them out (model !== "B"), so
  // toggling the model never destroys the user's tuning; the backend warns they are inert.
  const it = cfg.intarna;
  const BI = window.INTARNA_BLOCKS;
  obj.intarna = {
    ..._block(it, BI.top),
    helix:         _block(it.helix,         BI.helix),
    seed:          _block(it.seed,          BI.seed),
    accessibility: _block(it.accessibility, BI.accessibility),
    output:        _block(it.output,        BI.output),
  };

  // `type` alone enables plots (Snakefile: `plots_enabled = plots_type is not None`);
  // every other key is read only inside rule build_plots. The rest of the block is
  // therefore inert under `type: null`, so always write it — disabling plots must not
  // destroy the user's basesize/thresholds/locus lists.
  const pl = cfg.plots;
  const isAll = pl.types.length > 0 && window.PLOT_TYPE_VALUES.every((t) => pl.types.includes(t));
  obj.plots = {
    type: pl.types.length === 0 ? null : (isAll ? "all" : pl.types.join(",")),
    ..._block(pl, B.plots),
    locus: _pairListToArray(pl.locus),
    gene: _pairListToArray(pl.gene),
    protein: _pairListToArray(pl.protein),
  };

  obj.results_dir = cfg.results_dir;
  return obj;
};

//----- inverse: parsed config dict -> form state (start from INITIAL_CONFIG defaults) -----//

//----- Inverse of resolveNullable. Dispatches on the DEFAULT's shape, not the raw value: a
//      pair arrives from YAML as the string "2,7", which cannot tell you which it is -----//
function _toNullable(v, dflt) {
  if (v === null || v === undefined) return { ...dflt, set: false };
  if (!("a" in dflt && "b" in dflt)) return { set: true, value: v };
  const [a, b] = String(v).split(",");
  return { set: true, a: Number(a), b: Number(b) };
}

//----- PyYAML (YAML 1.1) turns bare on/off into booleans; coerce variants back to strings -----//
function _variantStr(v, dflt) {
  if (v === true) return "on";
  if (v === false) return "off";
  return typeof v === "string" ? v : dflt;
}

//----- Reads one block back into form state, in place. "c" is a constant the form never
//      edits, so it is deliberately not hydrated. "s" keys on PRESENCE, not on null: an
//      explicit null blanks the field, where an absent key keeps the default -----//
function _hydrateBlock(rawBlock, state, spec) {
  if (!rawBlock) return;
  spec.forEach(([key, kind]) => {
    if (kind === "n") state[key] = _toNullable(rawBlock[key], state[key]);
    else if (kind === "v") state[key] = _variantStr(rawBlock[key], state[key]);
    else if (kind === "s") state[key] = key in rawBlock ? String(rawBlock[key] ?? "") : state[key];
    else if (kind === "p" && rawBlock[key] != null) state[key] = rawBlock[key];
  });
}

function _arrayToPairList(items, prefix) {
  if (!Array.isArray(items)) return [];
  return items.map((s) => {
    const [value, mirna] = String(s).split(",");
    return { id: window.mintId(prefix), value: value || "", mirna: mirna || "" };
  });
}

//----- parsed config dict -> form state; anything absent keeps its INITIAL_CONFIG default -----//
window.hydrateConfig = function hydrateConfig(raw) {
  const cfg = JSON.parse(JSON.stringify(window.INITIAL_CONFIG));
  if (!raw || typeof raw !== "object") return cfg;
  const B = window.CONFIG_BLOCKS;

  ["queries", "targets"].forEach((side) => {
    if (raw[side] && typeof raw[side] === "object") {
      cfg[side] = Object.entries(raw[side]).map(([key, path]) =>
        ({ id: window.mintId(side[0]), key, path: path || "" }));
    }
  });

  if (typeof raw.threads === "number") cfg.threads = raw.threads;
  _hydrateBlock(raw, cfg, B.shared);
  // configToObject rebuilds the file from form state, so opening and saving an old config
  // migrates it. Flag the removed key so the substitution is visible, not merely correct.
  // UI-only: configToObject builds its output explicitly, so this never reaches the YAML.
  cfg.hadMaxTotalEnergy = raw.max_total_energy !== undefined;

  _hydrateBlock(raw.rnacalibrate, cfg.rnacalibrate, B.rnacalibrate);
  _hydrateBlock(raw.rnahybrid, cfg.rnahybrid, B.rnahybrid);

  if (raw.intarna) {
    const ri = raw.intarna, ci = cfg.intarna;
    const BI = window.INTARNA_BLOCKS;
    _hydrateBlock(ri,               ci,               BI.top);
    _hydrateBlock(ri.helix,         ci.helix,         BI.helix);
    _hydrateBlock(ri.seed,          ci.seed,          BI.seed);
    _hydrateBlock(ri.accessibility, ci.accessibility, BI.accessibility);
    _hydrateBlock(ri.output,        ci.output,        BI.output);
  }

  if (raw.plots) {
    const rp = raw.plots;
    // Read the surrounding settings whether or not `type` is null — see the writer.
    cfg.plots.types = rp.type == null ? []
      : rp.type === "all" ? window.PLOT_TYPE_VALUES.slice()
      : String(rp.type).split(",").map((s) => s.trim()).filter(Boolean);
    _hydrateBlock(rp, cfg.plots, B.plots);
    window.PLOT_PAIR_LISTS.forEach(({ key }) => {
      if (rp[key] != null) cfg.plots[key] = _arrayToPairList(rp[key], key[0]);
    });
  }

  if (raw.results_dir != null) cfg.results_dir = raw.results_dir;
  return cfg;
};
