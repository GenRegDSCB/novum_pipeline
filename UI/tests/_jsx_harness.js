// Drives the SHIPPED UI/src/*.jsx under Node, transpiled by the vendored Babel, so
// test_ui_config_contract.py exercises real components rather than a re-implementation.
// There is no DOM: function components are called directly and their element trees walked,
// and App runs on a minimal hooks runtime (state, refs, effects) against a mocked bridge.
//
// Usage:   node _jsx_harness.js <path/to/UI>  < input.json
// stdin:   {"scenario": "field", "section": <window name>, "raw": <config>, "field": <name>,
//           "typed": [[<"a"|"b"|"value">, <string>], ...]}
//       |  {"scenario": "app", "loads": [<load_config result>, ...], "actions": [<name>, ...],
//           "bridge": "late"?}
//          A load is returned as-is, except {"__throw__": msg} (rejects) and
//          {"__pending__": true} (never settles). Actions: "save", "run", "reload".
//          "late" mounts App on pywebview's empty `api: {}` and only then fills it and
//          fires pywebviewready, as pywebview's two injection steps can.
// stdout:  the scenario's JSON result (see each runner below).
// A JS exception is left uncaught on purpose: a non-zero exit with a stack trace is a test
// failure, not a value to compare.
"use strict";
const fs = require("fs");
const path = require("path");
const vm = require("vm");

const [, , uiDir] = process.argv;
if (!uiDir) {
  process.stderr.write("usage: node _jsx_harness.js <UI dir> < input.json\n");
  process.exit(2);
}
const input = JSON.parse(fs.readFileSync(0, "utf8"));

// The options @babel/standalone gives a <script type="text/babel"> that names none of its own.
const PAGE_BABEL_OPTIONS = {
  presets: ["react", "env"],
  plugins: ["transform-class-properties", "transform-object-rest-spread",
            "transform-flow-strip-types"],
};

//----- One page's worth of globals: the vendored libraries, then the scripts in the page's
//      own order. Each .jsx runs in its own function scope, as the comment heading
//      components.jsx says the page's text/babel tags do -----//
function loadPage({ withApp, prepare }) {
  const sandbox = { console, setImmediate };
  sandbox.window = sandbox.self = sandbox;
  vm.createContext(sandbox);
  const run = (rel) => {
    const file = path.join(uiDir, rel);
    vm.runInContext(fs.readFileSync(file, "utf8"), sandbox, { filename: file });
  };
  const runJsx = (rel) => {
    const file = path.join(uiDir, rel);
    const { code } = sandbox.Babel.transform(fs.readFileSync(file, "utf8"),
                                             { ...PAGE_BABEL_OPTIONS, filename: file });
    vm.runInContext(`(function () {\n${code}\n})();`, sandbox, { filename: file });
  };
  run("vendor/react.production.min.js");
  run("vendor/babel.min.js");
  if (prepare) prepare(sandbox);
  ["src/data.js", "src/yaml.js", "src/config.js"].forEach(run);
  ["src/components.jsx", "src/sections.jsx"].concat(withApp ? ["src/app.jsx"] : []).forEach(runJsx);
  return sandbox;
}

//----- Every element in a tree, depth first. Descends into `children` and into any prop
//      holding elements (Subcard's `right`), but never calls a function component -----//
function elements(node, out = []) {
  if (Array.isArray(node)) { node.forEach((n) => elements(n, out)); return out; }
  if (!node || typeof node !== "object" || !node.props) return out;
  out.push(node);
  Object.values(node.props).forEach((v) => elements(v, out));
  return out;
}

//----- The visible strings of a tree: string/number children only, never attribute values -----//
function text(node) {
  if (Array.isArray(node)) return node.map(text).join("");
  if (typeof node === "string" || typeof node === "number") return String(node);
  if (!node || typeof node !== "object" || !node.props) return "";
  return Object.entries(node.props)
    .filter(([k, v]) => k === "children" || (v && typeof v === "object"))
    .map(([, v]) => text(v)).join("");
}

//----- "field": render one section, type into one NullableField's inputs (re-rendering
//      between keystrokes, as React would), and report what it rendered and what Save writes -----//
function runField(sandbox) {
  let cfg = sandbox.hydrateConfig(input.raw);
  const render = () => sandbox[input.section]({ cfg, setCfg: (next) => { cfg = next; } });
  const field = () => {
    const el = elements(render()).find((e) => e.type === sandbox.NullableField
                                             && e.props.name === input.field);
    if (!el) throw new Error(`no NullableField named ${input.field} in ${input.section}`);
    return el;
  };
  const inputs = () => elements(field().type(field().props)).filter((e) => e.type === "input");

  for (const [side, typed] of input.typed || []) {
    const box = inputs()[side === "b" ? 1 : 0];
    box.props.onChange({ target: { value: typed } });
  }
  const { name, doc, kind, pairLabels, forcedNull } = field().props;
  return {
    props: { name, doc, kind, pairLabels: pairLabels || null, forcedNull: !!forcedNull },
    inputs: inputs().map((e) => ({ step: e.props.step, inputMode: e.props.inputMode })),
    sectionText: text(render()),
    hasSelect: elements(render()).some((e) => e.type === "select"),
    obj: sandbox.configToObject(cfg),
  };
}

//----- Just enough of React's hooks for App: one slot per call site, effects run after
//      each render when their deps change. Re-rendering is explicit (see runApp) -----//
function hooksRuntime(sandbox) {
  const slots = [];
  let cursor = 0;
  let effects = [];
  const changed = (prev, deps) => !prev || !deps || deps.some((d, i) => !Object.is(d, prev[i]));
  sandbox.React = Object.assign({}, sandbox.React, {
    useState(init) {
      const s = slots[cursor] || (slots[cursor] = { v: typeof init === "function" ? init() : init });
      cursor += 1;
      return [s.v, (next) => { s.v = typeof next === "function" ? next(s.v) : next; }];
    },
    useRef(init) {
      const s = slots[cursor] || (slots[cursor] = { current: init });
      cursor += 1;
      return s;
    },
    useMemo(fn, deps) {
      const s = slots[cursor] || (slots[cursor] = {});
      cursor += 1;
      if (changed(s.deps, deps)) { s.deps = deps; s.v = fn(); }
      return s.v;
    },
    useEffect(fn, deps) {
      const s = slots[cursor] || (slots[cursor] = {});
      cursor += 1;
      if (s.ran && !changed(s.deps, deps)) return;
      effects.push(() => {
        if (typeof s.cleanup === "function") s.cleanup();
        Object.assign(s, { ran: true, deps, cleanup: fn() });
      });
    },
  });
  return (Component) => {
    cursor = 0;
    effects = [];
    const tree = Component({});
    effects.forEach((e) => e());
    return tree;
  };
}

//----- "app": mount App over a mocked pywebview bridge, play the actions, and report every
//      write the bridge saw plus the final screen text -----//
async function runApp() {
  const calls = { load: 0, save: [], run: [] };
  const loads = (input.loads || []).slice();
  const api = {
    load_config: () => {
      calls.load += 1;
      const next = loads.length ? loads.shift() : {};
      if (next && next.__pending__) return new Promise(() => {});
      if (next && next.__throw__) return Promise.reject(new Error(next.__throw__));
      return Promise.resolve(JSON.parse(JSON.stringify(next)));
    },
    save_config: (obj) => { calls.save.push(obj); return Promise.resolve({ ok: true, path: "Config/config.yaml" }); },
    run_pipeline: (obj) => { calls.run.push(obj); return Promise.resolve({ ok: true, pid: 1 }); },
  };

  let mounted = null;
  let render = null;
  const late = input.bridge === "late";
  const bridge = { api: late ? {} : api };
  const listeners = [];
  const sandbox = loadPage({
    withApp: true,
    prepare: (sb) => {
      render = hooksRuntime(sb);
      const noop = () => {};
      Object.assign(sb, {
        pywebview: bridge,
        localStorage: { getItem: () => null, setItem: noop },
        document: { documentElement: { setAttribute: noop }, getElementById: () => ({}) },
        ReactDOM: { createRoot: () => ({ render: (el) => { mounted = el; } }) },
        addEventListener: (type, fn) => { if (type === "pywebviewready") listeners.push(fn); },
        removeEventListener: (type, fn) => {
          const i = listeners.indexOf(fn);
          if (i >= 0) listeners.splice(i, 1);
        },
        setTimeout: () => 0, clearTimeout: noop, setInterval: () => 0, clearInterval: noop,
      });
    },
  });
  const App = mounted.type;
  const settle = () => new Promise((resolve) => setImmediate(resolve));
  const step = async () => { render(App); await settle(); return render(App); };

  let tree = await step();
  if (late) {
    Object.assign(bridge.api, api);
    listeners.splice(0).forEach((fn) => fn());
    tree = await step();
  }
  for (const action of input.actions || []) {
    const bar = elements(tree).find((e) => e.type === sandbox.ActionsBar);
    if (action === "save") await bar.props.onSave();
    else if (action === "run") await bar.props.onRun();
    else if (action === "reload") {
      const button = elements(tree).find((e) => e.type === "button" && /Reload/.test(text(e)));
      if (!button) throw new Error("no Reload button on screen");
      await button.props.onClick();
    } else throw new Error(`unknown action ${action}`);
    tree = await step();
  }
  return { loadCalls: calls.load, saveCalls: calls.save, runCalls: calls.run, text: text(tree) };
}

(async () => {
  let result;
  if (input.scenario === "field") result = runField(loadPage({ withApp: false }));
  else if (input.scenario === "app") result = await runApp();
  else throw new Error(`unknown scenario ${input.scenario}`);
  process.stdout.write(JSON.stringify(result));
})().catch((err) => {
  process.stderr.write(`${err && err.stack || err}\n`);
  process.exit(1);
});
