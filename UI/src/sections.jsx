// Novum Pipeline — section UIs
//
// One numbered card per config section, in the order they appear on screen. Every card
// takes `{cfg, setCfg}` and edits its own slice; the schema itself lives in data.js.

//----- Numbered card shell shared by all eight sections -----//
function SectionCard({ num, title, sub, anchor, children, right }) {
  return (
    <section className="card" data-anchor={anchor}>
      <div className="card-head">
        <div className="num">{num}</div>
        <div className="title">{title}</div>
        <div className="sub">{sub}</div>
        <div className="spacer" />
        {right}
      </div>
      <div className="card-body">{children}</div>
    </section>
  );
}

//----- 1./2. Queries and Targets: one keyed picker driven by KEYED_FILE_SECTIONS. `rowHint`
//      renders the per-row note only Targets wants (does this key name a query?) -----//
function KeyedFileSection({ cfg, setCfg, spec, rowHint }) {
  const items = useApiList(spec.api);
  const hasApi = items !== undefined;
  const rows = cfg[spec.key];
  const { add, update, remove } = makeListEditor(
    rows,
    (next) => setCfg({ ...cfg, [spec.key]: next }),
    () => ({ id: window.mintId(spec.key[0]), key: "", path: "" }),
    1,
  );
  return (
    <SectionCard
      num={spec.num}
      title={spec.title}
      sub={spec.sub}
      anchor={spec.key}
      right={<button className="btn sm" onClick={add}><PlusIcon /> Add {spec.noun}</button>}
    >
      <div className="stack sm">
        <div className={`eyebrow head ${keyedRowClass(hasApi)}`}>
          <span>Sample key</span>
          <span>{hasApi ? spec.pickerHead : spec.pathHead}</span>
          <span></span>
          {!hasApi && <span></span>}
        </div>
        {rows.length === 0 && (
          <div className="hint italic" style={{ padding: "8px 2px" }}>
            No {spec.many} yet — click <b>Add {spec.noun}</b>.
          </div>
        )}
        {rows.map((row, i) => (
          <div key={row.id} className="stack xs">
            <KeyedFileRow
              keyVal={row.key}
              pathVal={row.path}
              items={items}
              onKey={(v) => update(i, { key: v })}
              onPath={(v) => update(i, { path: v })}
              onRemove={() => remove(i)}
              canRemove={rows.length > 1}
              placeholderKey="e.g. escherichia"
              placeholderPath={spec.placeholderPath}
              emptyHint={spec.emptyHint}
            />
            {rowHint && rowHint(row)}
          </div>
        ))}
        <div className="hint" style={{ paddingTop: 4 }}>
          {countOf(rows.length, spec.noun, spec.many)}{spec.footerNote}
        </div>
      </div>
    </SectionCard>
  );
}

//----- The two specs, destructured in table order -----//
const [QUERIES_SPEC, TARGETS_SPEC] = window.KEYED_FILE_SECTIONS;

function QueriesSection(props) {
  return <KeyedFileSection {...props} spec={QUERIES_SPEC} />;
}

function TargetsSection({ cfg, setCfg }) {
  const queryKeys = cfg.queries.map((q) => q.key).filter(Boolean);
  const rowHint = (t) => (
    <div className={`hint ${t.key && !queryKeys.includes(t.key) ? "warn" : ""}`}
         style={{ padding: "2px 4px" }}>
      {!t.key ? <span>Enter a key shared with a query</span>
        : queryKeys.includes(t.key) ? <span>✓ matches query <code>{t.key}</code></span>
        : <span>No query named <code>{t.key}</code></span>}
    </div>
  );
  return <KeyedFileSection cfg={cfg} setCfg={setCfg} spec={TARGETS_SPEC} rowHint={rowHint} />;
}

//----- 3. Shared parameters — apply to both arms -----//
function SharedParamsSection({ cfg, setCfg }) {
  const set = (patch) => setCfg({ ...cfg, ...patch });
  const setNullable = (key, next) => set({ [key]: next });
  return (
    <SectionCard
      num="3"
      title="Shared parameters"
      sub="Apply to both RNAhybrid and IntaRNA"
      anchor="shared"
    >
      {cfg.hadMaxTotalEnergy && (
        <InlineNote tone="warn">
          · max_total_energy was removed — it meant hybridization energy on RNAhybrid but
          total energy on IntaRNA. The per-arm cutoffs now live on each tool's card; saving
          rewrites the file.
        </InlineNote>
      )}
      <div className="grid-2">
        <Subcard label="Threads">
          <input
            className="input mono"
            type="number"
            min={1}
            step="1"
            value={cfg.threads}
            onChange={(e) => set({ threads: Number(e.target.value) })}
          />
        </Subcard>
        {/* "auto", not "loose": a null here means each arm falls back to its own default. */}
        <NullableField name="Max suboptimal hits" kind="integer"
          min={1} nullHint="null · auto"
          value={cfg.max_suboptimal_hits}
          onChange={(v) => setNullable("max_suboptimal_hits", v)} />
        <NullableField name="Seed" doc="nucleotides start,end" kind="pair"
          min={1} nullHint="null · auto"
          pairLabels={["start nt", "end nt"]}
          value={cfg.seed}
          onChange={(v) => setNullable("seed", v)} />
      </div>
    </SectionCard>
  );
}

//----- 4. RNAcalibrate — variant-based; its whole body is inert when the variant is off -----//
function RNACalibrateSection({ cfg, setCfg }) {
  const r = cfg.rnacalibrate;
  const set = (patch) => setCfg({ ...cfg, rnacalibrate: { ...r, ...patch } });
  const inert = r.calibration_variant === "off";
  return (
    <SectionCard
      num="4"
      title="RNAcalibrate"
      sub="Calibrate hybridization against random targets"
      anchor="rnacalibrate"
      right={
        inert ? <span className="eyebrow inert-badge">Parameters inactive</span> : null
      }
    >
      <div className="stack lg">
        <Subcard label="Calibration variant">
          <SegmentedControl
            options={window.VARIANT_MODES}
            value={r.calibration_variant}
            onChange={(v) => set({ calibration_variant: v })}
            ariaLabel="Calibration variant"
          />
        </Subcard>
        {inert && (
          <div className="callout">
            <span className="badge" aria-hidden>i</span>
            <div>
              Calibration is off. The parameters below won’t affect the run — set
              RNAhybrid’s <code>distribution</code> instead.
            </div>
          </div>
        )}
        <DisableGroup disabled={inert} className="stack lg">
          <div className="grid-2">
            <Subcard label="k (sample size)">
              <input
                className="input mono"
                type="number"
                step="1"
                min={0}
                value={r.k}
                onChange={(e) => set({ k: Number(e.target.value) })}
              />
            </Subcard>
            <Subcard label="Maximum target length" hint="nt">
              <input
                className="input mono"
                type="number"
                step="1"
                min={0}
                value={r.max_target_length}
                onChange={(e) => set({ max_target_length: Number(e.target.value) })}
              />
            </Subcard>
          </div>
          <Subcard
            label="Randomize targets"
            hint="shuffle target sequences before calibration"
            right={
              <Toggle
                on={r.randomize_targets}
                onChange={(v) => set({ randomize_targets: v })}
                ariaLabel="Randomize targets"
              />
            }
          />
          <Subcard
            label="Target-length anchors"
            hint="nt, comma-separated — one Gumbel null per cell"
          >
            <input
              className="input mono"
              type="text"
              value={r.length_anchors}
              onChange={(e) => set({ length_anchors: e.target.value })}
            />
          </Subcard>
        </DisableGroup>
      </div>
    </SectionCard>
  );
}

//----- 5. RNAhybrid -----//
function RNAHybridSection({ cfg, setCfg }) {
  const r = cfg.rnahybrid;
  const set = (patch) => setCfg({ ...cfg, rnahybrid: { ...r, ...patch } });
  const setF = (key, next) => set({ [key]: next });
  const distLocked = window.isDistributionForced(cfg);
  // Under "both" species/distribution still drive the wo_calibration arm, and only it.
  const armNote = cfg.rnacalibrate.calibration_variant === "both" ? " · wo_calibration arm only" : "";

  return (
    <SectionCard
      num="5"
      title="RNAhybrid"
      sub="Target prediction parameters"
      anchor="rnahybrid"
    >
      <div className="stack lg">
        {distLocked ? (
          <Subcard label="Species model" locked
                   right={<span className="auto-tag locked">LOCKED · overridden</span>}>
            <div className="ghost-input">
              overridden · distribution is provided by RNAcalibrate
            </div>
          </Subcard>
        ) : (
          <Subcard label="Species model" hint={`3′UTR distribution${armNote}`}>
            <select
              className="select"
              value={r.species}
              onChange={(e) => set({ species: e.target.value })}
            >
              {window.SPECIES_OPTIONS.map((o) => (
                <option key={o.value} value={o.value}>{o.label}</option>
              ))}
            </select>
          </Subcard>
        )}

        <div className="grid-2">
          <NullableField name="Max hybridization energy" doc="kcal/mol · MFE"
            max={0} nullHint="null · no cutoff"
            value={r.max_hybrid_energy} onChange={(v) => setF("max_hybrid_energy", v)} />
          <NullableField name="Max internal loop" doc="nt" kind="integer"
            min={0}
            value={r.max_internal_loop} onChange={(v) => setF("max_internal_loop", v)} />
          <NullableField name="Max bulge loop" doc="nt" kind="integer"
            min={0}
            value={r.max_bulge_loop} onChange={(v) => setF("max_bulge_loop", v)} />
          <NullableField name="p-value threshold"
            min={0}
            value={r.pvalue_threshold} onChange={(v) => setF("pvalue_threshold", v)} />
          <NullableField name="Distribution" doc={`Gumbel xi,theta${armNote}`} kind="pair"
            pairNumber="float" min={0}
            pairLabels={["xi", "theta"]}
            value={r.distribution}
            onChange={(v) => setF("distribution", v)}
            forcedNull={distLocked}
            forcedNullHint="null · distribution is provided by RNAcalibrate"
          />
        </div>
      </div>
    </SectionCard>
  );
}

//----- 6. IntaRNA -----//
function IntaRNASection({ cfg, setCfg }) {
  const it = cfg.intarna;
  const set      = (patch) => setCfg({ ...cfg, intarna: { ...it, ...patch } });
  const setSeed  = (patch) => set({ seed:          { ...it.seed,          ...patch } });
  const setHelix = (patch) => set({ helix:         { ...it.helix,         ...patch } });
  const setAcc   = (patch) => set({ accessibility: { ...it.accessibility, ...patch } });
  const setOut   = (patch) => set({ output:        { ...it.output,        ...patch } });

  const accInert     = !window.isAccessibilityActive(cfg);
  const seedEnforced = window.isSeedEnforced(cfg);                         // enforcement follows the shared seed
  const helixActive = window.isHelixActive(cfg);                   // helix reaches IntaRNA only under model B
  const helixInert   = window.helixInertKeys(cfg);                  // mirrors backend validator (c)
  const seedConflict = window.helixSeedConflict(cfg); // mirrors backend validator (a)
  const accConflicts = window.accessibilityConflicts(cfg); // mirrors backend validator (f)

  //----- "Global" is not stored — it IS window === 0 && max_bp_span === 0 on that side, so
  //      the pill can never disagree with the fields it summarizes. Untoggling unsets both,
  //      handing that side back to IntaRNA's own defaults -----//
  const accSideGlobal = (side) => {
    const w = it.accessibility[`${side}_window`];
    const s = it.accessibility[`${side}_max_bp_span`];
    return w.set && w.value === 0 && s.set && s.value === 0;
  };
  const setSideGlobal = (side, on) => setAcc(
    on
      ? { [`${side}_window`]:      { set: true, value: 0 },
          [`${side}_max_bp_span`]: { set: true, value: 0 } }
      : { [`${side}_window`]:      { set: false, value: window.INTARNA_DEFAULT_ACC_WINDOW },
          [`${side}_max_bp_span`]: { set: false, value: window.INTARNA_DEFAULT_ACC_BP_SPAN } }
  );

  return (
    <SectionCard
      num="6"
      title="IntaRNA"
      sub="Accessibility-aware RNA–RNA interaction prediction"
      anchor="intarna"
    >
      <div className="stack lg">
        <Subcard label="Accessibility variant">
          <SegmentedControl
            options={window.VARIANT_MODES}
            value={it.accessibility_variant}
            onChange={(v) => set({ accessibility_variant: v })}
            ariaLabel="Accessibility variant"
          />
        </Subcard>

        <div className="grid-2">
          <Subcard label="Prediction mode">
            <select
              className="select"
              value={it.prediction_mode}
              onChange={(e) => set({ prediction_mode: e.target.value })}
            >
              {window.INTARNA_PREDICTION_MODES.map((o) => (
                <option key={o.value} value={o.value}>{o.label}</option>
              ))}
            </select>
          </Subcard>
          <Subcard label="Interaction model">
            <select
              className="select"
              value={it.model}
              onChange={(e) => set({ model: e.target.value })}
            >
              {window.INTARNA_MODELS.map((o) => (
                <option key={o.value} value={o.value}>{o.label}</option>
              ))}
            </select>
          </Subcard>
          <Subcard label="Energy set" hint="Nearest-Neighbor table">
            <select
              className="select"
              value={it.energy_set}
              onChange={(e) => set({ energy_set: e.target.value })}
            >
              {window.INTARNA_ENERGY_SETS.map((o) => (
                <option key={o.value} value={o.value}>{o.label}</option>
              ))}
            </select>
          </Subcard>
          <NullableField name="Max hybridization energy" doc="kcal/mol · E_hybrid"
            max={0} nullHint="null · tool default (E ≤ 0)"
            value={it.max_hybrid_energy}
            onChange={(v) => set({ max_hybrid_energy: v })} />
          <NullableField name="Min target unpaired probability" doc="Pu1 floor · gene side"
            min={0} max={1} nullHint="null · no accessibility floor"
            value={it.min_target_unpaired_probability}
            onChange={(v) => set({ min_target_unpaired_probability: v })} />
          <NullableField name="Min query unpaired probability" doc="Pu2 floor · miRNA side"
            min={0} max={1} nullHint="null · no accessibility floor"
            value={it.min_query_unpaired_probability}
            onChange={(v) => set({ min_query_unpaired_probability: v })} />
          <NullableField name="Accessibility search depth" doc="hits reported" kind="integer"
            min={1} nullHint="null · use shared cap"
            value={it.accessibility_search_depth}
            onChange={(v) => set({ accessibility_search_depth: v })} />
          <NullableField name="Max interaction length" doc="nt" kind="integer"
            min={0} autoTag="AUTO (0)"
            value={it.max_interaction_length}
            onChange={(v) => set({ max_interaction_length: v })} />
          <NullableField name="Max loop size" doc="unpaired nt" kind="integer"
            min={0} autoTag="AUTO (10)"
            value={it.max_loop_size}
            onChange={(v) => set({ max_loop_size: v })} />
        </div>

        <div>
          <div className="row-head group">
            <div className="eyebrow strong">Helix</div>
            {!helixActive && (helixInert.length > 0 ? (
              <InlineNote tone="warn">
                · {helixInert.join(", ")} set but ignored under model {it.model}
              </InlineNote>
            ) : (
              <InlineNote>· Applies only to model B</InlineNote>
            ))}
            {seedConflict && (
              <InlineNote tone="warn">
                · seed needs {seedConflict.seedBP} bp but max bp is {seedConflict.maxBP} — IntaRNA will fail
              </InlineNote>
            )}
          </div>
          <DisableGroup disabled={!helixActive} className="stack sm">
            <div className="grid-2">
              <NullableField name="Min bp" doc="[2..4]" kind="integer"
                min={2} max={4}
                value={it.helix.min_bp} onChange={(v) => setHelix({ min_bp: v })} />
              <NullableField name="Max bp" doc="[2..20]" kind="integer"
                min={2} max={20}
                value={it.helix.max_bp} onChange={(v) => setHelix({ max_bp: v })} />
              <NullableField name="Max internal loop" doc="[0..2], 0 = pure stacks" kind="integer"
                min={0} max={2}
                value={it.helix.max_internal_loop} onChange={(v) => setHelix({ max_internal_loop: v })} />
              <NullableField name="Max energy" doc="kcal/mol"
                value={it.helix.max_energy} onChange={(v) => setHelix({ max_energy: v })} />
              <DisableGroup disabled={accInert}>
                <NullableField name="Min unpaired probability" doc="[0..1]"
                  min={0} max={1}
                  value={it.helix.min_unpaired_probability} onChange={(v) => setHelix({ min_unpaired_probability: v })} />
              </DisableGroup>
            </div>
            <PillGroup label="Constraints">
              <PillToggle on={it.helix.full_energy} onChange={(v) => setHelix({ full_energy: v })}
                label="Full helix energy (max energy counts E_init, ED, dangles)" />
            </PillGroup>
          </DisableGroup>
        </div>

        <div>
          <div className="row-head group">
            <div className="eyebrow strong">Seed</div>
            <div className="spacer" />
            <span className={`status-text ${seedEnforced ? "on" : ""}`}>
              {seedEnforced ? "Enforced via shared seed" : "Disabled — set a shared seed"}
            </span>
          </div>
          <DisableGroup disabled={!seedEnforced} className="stack sm">
            <div className="grid-2">
              <NullableField name="Max unpaired bases" kind="integer"
                min={0}
                value={it.seed.max_unpaired_bases} onChange={(v) => setSeed({ max_unpaired_bases: v })} />
              <NullableField name="Target range" doc="'from-to,...'"
                value={it.seed.target_range} onChange={(v) => setSeed({ target_range: v })} />
              <NullableField name="Max energy"
                max={999}
                value={it.seed.max_energy} onChange={(v) => setSeed({ max_energy: v })} />
              <NullableField name="Max hybrid energy"
                max={999}
                value={it.seed.max_hybrid_energy} onChange={(v) => setSeed({ max_hybrid_energy: v })} />
              <NullableField name="Min unpaired probability" doc="[0..1]"
                min={0} max={1}
                value={it.seed.min_unpaired_probability} onChange={(v) => setSeed({ min_unpaired_probability: v })} />
            </div>
            <PillGroup label="Constraints">
              <PillToggle on={it.seed.forbid_gu}         onChange={(v) => setSeed({ forbid_gu: v })}         label="Forbid G:U pairs" />
              <PillToggle on={it.seed.forbid_gu_at_ends} onChange={(v) => setSeed({ forbid_gu_at_ends: v })} label="Forbid G:U at ends" />
              <PillToggle on={it.seed.report_best_only}  onChange={(v) => setSeed({ report_best_only: v })}  label="Report best only" />
            </PillGroup>
          </DisableGroup>
        </div>

        <div>
          <div className="row-head group">
            <div className="eyebrow strong">Accessibility</div>
            {accInert && (
              <InlineNote>· variant is Off — flags below are sent but ignored</InlineNote>
            )}
            {accConflicts.length > 0 && (
              <InlineNote tone="warn">
                · {accConflicts.map(({ sharedKey, sideKey }) => `${sharedKey} vs ${sideKey}`).join(", ")}
                {" "}differ — IntaRNA will abort on a shared window set beside a differing per-side one
              </InlineNote>
            )}
          </div>
          <DisableGroup disabled={accInert} className="stack sm">
            <div className="grid-2">
              <NullableField name="Window" doc="0 = whole sequence" kind="integer"
                min={0} autoTag={`AUTO (${window.INTARNA_DEFAULT_ACC_WINDOW})`}
                value={it.accessibility.window} onChange={(v) => setAcc({ window: v })} />
              <NullableField name="Max bp span" doc="0 = use window" kind="integer"
                min={0} autoTag={`AUTO (${window.INTARNA_DEFAULT_ACC_BP_SPAN})`}
                value={it.accessibility.max_bp_span} onChange={(v) => setAcc({ max_bp_span: v })} />
            </div>

            <Subcard label="Per-side overrides"
                     hint="leave the shared pair above unset when using these — IntaRNA aborts on a mismatch">
              <div className="stack sm">
                {window.ACCESSIBILITY_SIDES.map((spec) => (
                  <div key={spec.side}>
                    <PillGroup label={spec.label}>
                      <PillToggle on={accSideGlobal(spec.side)}
                                  onChange={(v) => setSideGlobal(spec.side, v)}
                                  label="Global" />
                    </PillGroup>
                    <DisableGroup disabled={accSideGlobal(spec.side)} className="grid-2">
                      <NullableField name={`${spec.name} window`} doc={spec.wFlag} kind="integer"
                        min={0} autoTag={`AUTO (${window.INTARNA_DEFAULT_ACC_WINDOW})`}
                        value={it.accessibility[`${spec.side}_window`]}
                        onChange={(v) => setAcc({ [`${spec.side}_window`]: v })} />
                      <NullableField name={`${spec.name} max bp span`} doc={spec.lFlag} kind="integer"
                        min={0} autoTag={`AUTO (${window.INTARNA_DEFAULT_ACC_BP_SPAN})`}
                        value={it.accessibility[`${spec.side}_max_bp_span`]}
                        onChange={(v) => setAcc({ [`${spec.side}_max_bp_span`]: v })} />
                    </DisableGroup>
                  </div>
                ))}
              </div>
            </Subcard>

            <PillGroup label="Constraints">
              <PillToggle on={it.accessibility.forbid_lonely_pairs} onChange={(v) => setAcc({ forbid_lonely_pairs: v })} label="Forbid lonely pairs" />
              <PillToggle on={it.accessibility.forbid_gu_at_ends}   onChange={(v) => setAcc({ forbid_gu_at_ends: v })}   label="Forbid G:U at ends" />
            </PillGroup>
          </DisableGroup>
        </div>

        <div>
          <div className="row-head group"><div className="eyebrow strong">Output</div></div>
          <div className="grid-2">
            <NullableField name="Max ΔE" doc="cap vs MFE"
              min={0}
              value={it.output.max_delta_energy} onChange={(v) => setOut({ max_delta_energy: v })} />
            <Subcard label="Overlap">
              <select
                className="select"
                value={it.output.overlap}
                onChange={(e) => setOut({ overlap: e.target.value })}
              >
                {window.INTARNA_OVERLAPS.map((o) => (
                  <option key={o.value} value={o.value}>{o.label}</option>
                ))}
              </select>
            </Subcard>
            <NullableField name="Min unpaired probability" doc="[0..1]"
              min={0} max={1}
              value={it.output.min_unpaired_probability} onChange={(v) => setOut({ min_unpaired_probability: v })} />
            <PillGroup label="Constraints">
              <PillToggle on={it.output.forbid_lonely_pairs} onChange={(v) => setOut({ forbid_lonely_pairs: v })} label="Forbid lonely pairs" />
              <PillToggle on={it.output.forbid_gu_at_ends}   onChange={(v) => setOut({ forbid_gu_at_ends: v })}   label="Forbid G:U at ends" />
            </PillGroup>
          </div>
        </div>
      </div>
    </SectionCard>
  );
}

//----- 7. Plots — the whole block is optional; an empty type list disables it -----//
function PlotsSection({ cfg, setCfg }) {
  const pl = cfg.plots;
  const set = (patch) => setCfg({ ...cfg, plots: { ...pl, ...patch } });
  const enabled = pl.types.length > 0;
  const setEnabled = (v) => set({ types: v ? window.PLOT_TYPES.map((p) => p.value) : [] });
  const toggleType = (val) => {
    set({ types: pl.types.includes(val) ? pl.types.filter((t) => t !== val) : [...pl.types, val] });
  };

  const editors = Object.fromEntries(
    window.PLOT_PAIR_LISTS.map(({ key }) => [
      key,
      makeListEditor(
        pl[key],
        (next) => set({ [key]: next }),
        () => ({ id: window.mintId(key[0]), value: "", mirna: "" }),
      ),
    ])
  );

  return (
    <SectionCard
      num="7"
      title="Plots"
      sub="Optional · omit or disable to skip"
      anchor="plots"
      right={
        <div className="tag-row">
          <span className={`status-text ${enabled ? "on" : ""}`}>
            {enabled ? "Enabled" : "Disabled"}
          </span>
          <Toggle on={enabled} onChange={setEnabled} ariaLabel="Enable plots" />
        </div>
      }
    >
      <DisableGroup disabled={!enabled} className="stack lg">
        <PillGroup label="Plot types" hint={`${pl.types.length}/${window.PLOT_TYPES.length} on`}>
          {window.PLOT_TYPES.map((p) => (
            <PillToggle
              key={p.value}
              on={pl.types.includes(p.value)}
              onChange={() => toggleType(p.value)}
              label={p.label}
            />
          ))}
        </PillGroup>
        <div className="grid-3">
          <Subcard label="Base font size">
            <input
              className="input mono"
              type="number"
              step="1"
              min={0}
              value={pl.basesize}
              onChange={(e) => set({ basesize: Number(e.target.value) })}
            />
          </Subcard>
          <Subcard label="p-value threshold" hint="0 – 1">
            <input
              className="input mono"
              type="number"
              step="0.001"
              min={0}
              value={pl.pvalue_threshold}
              onChange={(e) => set({ pvalue_threshold: Number(e.target.value) })}
            />
          </Subcard>
          <Subcard label="Top-N miRNAs" hint="per-miRNA plot cap">
            <input
              className="input mono"
              type="number"
              min={1}
              step="1"
              value={pl.per_mirna_top_n}
              onChange={(e) => set({ per_mirna_top_n: Number(e.target.value) })}
            />
          </Subcard>
        </div>

        {window.PLOT_PAIR_LISTS.map((spec) => (
          <PairListEditor
            key={spec.key}
            heading={spec.heading}
            noun={spec.noun}
            valueHeader={spec.valueHeader}
            valuePlaceholder={spec.valuePlaceholder}
            items={pl[spec.key]}
            editor={editors[spec.key]}
          />
        ))}
      </DisableGroup>
    </SectionCard>
  );
}

//----- 8. Output -----//
function OutputSection({ cfg, setCfg }) {
  return (
    <SectionCard num="8" title="Output" sub="Destination folder" anchor="output">
      <Subcard label="Results directory" hint="created if missing">
        <PathInput
          value={cfg.results_dir}
          onChange={(v) => setCfg({ ...cfg, results_dir: v })}
          placeholder="Data/Results/"
        />
      </Subcard>
    </SectionCard>
  );
}

//----- Publish to window; see the note at the head of components.jsx -----//
Object.assign(window, {
  SectionCard, KeyedFileSection,
  QueriesSection, TargetsSection, SharedParamsSection,
  RNACalibrateSection, RNAHybridSection, IntaRNASection,
  PlotsSection, OutputSection
});
