// Novum Pipeline — main app (single-scroll layout)

const { useState, useEffect, useRef } = React;

//----- Pipeline lifecycle: idle → running → succeeded | failed | cancelled. A terminal state
//      stays until the overlay is dismissed, then acknowledge_pipeline resets the backend,
//      so a fast finish is never missed between two polls -----//

//----- Why a load_config result cannot seed the form, else null. config_io answers {} for a
//      file it cannot read or that holds no mapping, and {__error__} for unparseable YAML -----//
function loadFailure(raw) {
  if (raw && raw.__error__) return `parse error · ${raw.__error__}`;
  if (!raw || typeof raw !== "object" || !Object.keys(raw).length) {
    return "file missing, unreadable, empty or not a mapping";
  }
  return null;
}

//----- Root component: owns the config, the theme, and the pipeline state machine -----//
function App() {
  const [cfg, setCfg] = useState(() => JSON.parse(JSON.stringify(window.INITIAL_CONFIG)));
  // Save and Run rewrite config.yaml whole from `cfg`, so neither may write until a load has
  // succeeded: { ok: false, error: null } is a load in flight, a set `error` one that failed.
  const [load, setLoad] = useState({ ok: false, error: null });
  const [theme, setTheme] = useState(() => localStorage.getItem("np:theme") || "light");
  const [yamlOpen, setYamlOpen] = useState(false);
  const [toast, setToast] = useState(null);
  const toastTimer = useRef(null);

  const [saving, setSaving] = useState(false);
  const [pipelineStatus, setPipelineStatus] = useState(null);

  const pipelineState = pipelineStatus?.state ?? "idle";
  const overlayOpen = pipelineStatus !== null;

  useEffect(() => {
    document.documentElement.setAttribute("data-theme", theme);
    localStorage.setItem("np:theme", theme);
  }, [theme]);

  useEffect(() => () => clearTimeout(toastTimer.current), []);

  //----- Read config.yaml into the form. A failure leaves the form as it was and keeps
  //      writes blocked until a later reload succeeds -----//
  const loadConfig = async () => {
    const api = window.pywebview?.api;
    if (!api?.load_config) return;
    try {
      const raw = await api.load_config();
      const error = loadFailure(raw);
      if (error) return setLoad({ ok: false, error });
      setCfg(window.hydrateConfig(raw));
      setLoad({ ok: true, error: null });
    } catch (exc) {
      setLoad({ ok: false, error: String(exc?.message || exc) });
    }
  };

  // Keyed on the method, not `api`: pywebview creates `api: {}` before filling it.
  useEffect(() => {
    if (typeof window.pywebview?.api?.load_config === "function") loadConfig();
    else window.addEventListener("pywebviewready", loadConfig, { once: true });
    return () => window.removeEventListener("pywebviewready", loadConfig);
  }, []);

  //----- One toast at a time; a new one replaces whatever is on screen -----//
  const showToast = (msg, tone = "neutral") => {
    clearTimeout(toastTimer.current);
    setToast({ msg, tone });
    toastTimer.current = setTimeout(() => setToast(null), 3200);
  };

  //----- Refuses (with a toast) while no load has succeeded; true means "do not write" -----//
  const writeBlocked = () => {
    if (load.ok) return false;
    showToast(load.error ? "Config not loaded · fix the file, then Reload" : "Config still loading", "danger");
    return true;
  };

  //----- Write config.yaml without launching anything -----//
  const save = async () => {
    if (saving) return;
    const api = window.pywebview?.api;
    setSaving(true);
    try {
      if (!api) {
        showToast("Preview mode · nothing written");
        return;
      }
      if (writeBlocked()) return;
      const res = await api.save_config(window.configToObject(cfg));
      if (res?.ok) showToast(`Saved to ${res.path}`, "success");
      else showToast(`Save failed · ${res?.error || "unknown error"}`, "danger");
    } finally {
      setSaving(false);
    }
  };

  //----- Save, then launch snakemake and switch the overlay on -----//
  const run = async () => {
    if (pipelineState === "running") return;
    const api = window.pywebview?.api;
    if (!api?.run_pipeline) {
      showToast("Preview mode · cannot launch");
      return;
    }
    if (writeBlocked()) return;               // run_pipeline saves the config before launching
    // Clear any stale terminal state so the polling effect can't snap the
    // overlay back to a prior 'succeeded'/'failed' as the new run starts.
    // If the backend refuses (a previous run is still draining), surface
    // that instead of plowing into run_pipeline which would also reject.
    if (api.acknowledge_pipeline) {
      const ack = await api.acknowledge_pipeline();
      if (ack && ack.ok === false) {
        showToast("Previous run still finishing — try again", "danger");
        return;
      }
    }
    setPipelineStatus({ state: "running", elapsed: 0 });
    const res = await api.run_pipeline(window.configToObject(cfg));
    if (!res?.ok) {
      setPipelineStatus({ state: "failed", returncode: -1, elapsed: 0 });
      showToast(`Run failed · ${res?.error || "unknown error"}`, "danger");
    }
  };

  //----- SIGTERM the run; the poll below picks up the terminal state -----//
  const cancel = async () => {
    const api = window.pywebview?.api;
    if (api?.cancel_pipeline) await api.cancel_pipeline();
  };

  //----- Acknowledge the terminal state so the backend returns to idle -----//
  const dismissOverlay = async () => {
    const api = window.pywebview?.api;
    if (api?.acknowledge_pipeline) await api.acknowledge_pipeline();
    setPipelineStatus(null);
  };

  useEffect(() => {
    if (pipelineState !== "running") return;
    const api = window.pywebview?.api;
    if (!api?.pipeline_status) return;
    let cancelled = false;

    const poll = async () => {
      if (cancelled) return;
      try {
        const status = await api.pipeline_status();
        // Race: between optimistic running and run_pipeline returning, the
        // backend reports idle. Ignore so the overlay doesn't flicker back.
        if (cancelled || status.state === "idle") return;
        setPipelineStatus(status);
        if (status.state === "succeeded") {
          showToast(`Pipeline complete · ${window.formatElapsed(status.elapsed)}`, "success");
        } else if (status.state === "failed") {
          showToast(`Pipeline failed · exit ${status.returncode}`, "danger");
        } else if (status.state === "cancelled") {
          showToast("Pipeline cancelled");
        }
      } catch (_) { /* transient bridge errors retry on the next tick */ }
    };

    poll();
    const id = setInterval(poll, 1000);
    return () => { cancelled = true; clearInterval(id); };
  }, [pipelineState]);

  return (
    <div className="desk">
      <div className="window">
        <div className="window-body">
          {load.error && (
            <div className="callout" role="alert"
                 style={{ margin: "12px 24px 0", borderColor: "var(--danger)", color: "var(--danger)" }}>
              <span className="badge" aria-hidden>!</span>
              <div style={{ flex: 1 }}>
                Config/config.yaml was not loaded ({load.error}). Save and Run are blocked so
                defaults cannot overwrite it; fix the file, then reload.
              </div>
              <button className="btn sm" onClick={loadConfig}>Reload</button>
            </div>
          )}
          <div className="no-scrollbar" style={{ flex: 1, overflowY: "auto" }}>
            <div className="card-stack">
              <ErrorBoundary>
                <Hero cfg={cfg} />
                <QueriesSection cfg={cfg} setCfg={setCfg} />
                <TargetsSection cfg={cfg} setCfg={setCfg} />
                <SharedParamsSection cfg={cfg} setCfg={setCfg} />
                <RNACalibrateSection cfg={cfg} setCfg={setCfg} />
                <RNAHybridSection cfg={cfg} setCfg={setCfg} />
                <IntaRNASection cfg={cfg} setCfg={setCfg} />
                <PlotsSection cfg={cfg} setCfg={setCfg} />
                <OutputSection cfg={cfg} setCfg={setCfg} />
              </ErrorBoundary>
            </div>
          </div>
          <ActionsBar
            cfg={cfg}
            onSave={save}
            onRun={run}
            onToggleYAML={() => setYamlOpen((o) => !o)}
            theme={theme}
            onTheme={setTheme}
            pipelineState={pipelineState}
            elapsed={pipelineStatus?.elapsed}
            saving={saving}
          />
          <YAMLDrawer open={yamlOpen} onClose={() => setYamlOpen(false)} cfg={cfg} onSave={save} />
          <RunOverlay
            open={overlayOpen}
            status={pipelineStatus}
            threads={cfg.threads}
            queries={cfg.queries.length}
            targets={cfg.targets.length}
            onCancel={cancel}
            onClose={dismissOverlay}
          />
        </div>
      </div>

      {toast && <Toast tone={toast.tone}>{toast.msg}</Toast>}
    </div>
  );
}

//----- Mount -----//
const root = ReactDOM.createRoot(document.getElementById("root"));
root.render(<App />);
