// Round-trips a raw config through the SHIPPED UI/src/data.js + UI/src/config.js, run in a
// Node `vm` context, so test_ui_config_contract.py exercises the real serializer rather
// than a Python re-implementation that would only ever test itself.
//
// Usage:   node _js_roundtrip.js <path/to/data.js> <path/to/config.js>  < input.json
// stdin:   {"raw": <object to pass to hydrateConfig>}
// stdout:  {"hydrated": <hydrateConfig(raw)>, "obj": <configToObject(hydrated)>,
//           "optional": <countOptionalParams(hydrated)>}
// A JS exception is left uncaught on purpose — a non-zero exit with a stack trace is a
// test failure, not a value to compare.
"use strict";
const fs = require("fs");
const vm = require("vm");

const [, , dataJsPath, configJsPath] = process.argv;
if (!dataJsPath || !configJsPath) {
  process.stderr.write("usage: node _js_roundtrip.js <data.js> <config.js> < input.json\n");
  process.exit(2);
}

// Real browsers have `window === globalThis`; data.js/config.js only ever write through
// `window.X = ...`, so a self-referential sandbox reproduces that without a DOM.
const sandbox = {};
sandbox.window = sandbox;
vm.createContext(sandbox);
vm.runInContext(fs.readFileSync(dataJsPath, "utf8"), sandbox, { filename: dataJsPath });
vm.runInContext(fs.readFileSync(configJsPath, "utf8"), sandbox, { filename: configJsPath });

const { raw } = JSON.parse(fs.readFileSync(0, "utf8"));
const hydrated = sandbox.hydrateConfig(raw);
const obj = sandbox.configToObject(hydrated);
const optional = sandbox.countOptionalParams(hydrated);
process.stdout.write(JSON.stringify({ hydrated, obj, optional }));
