// Novum Pipeline — YAML preview renderer.
//
// Renders `configToObject(cfg)` — the same object api.save_config dumps — so the
// preview cannot drift from what Save writes. This file knows no schema.

//----- YAML scalar class, for both the highlighter and the plain-text writer -----//
const kindOf = (v) =>
  v === null || v === undefined ? "null"
  : typeof v === "boolean" ? "bool"
  : typeof v === "number" ? "num"
  : "str";

const isMap = (v) => v !== null && typeof v === "object" && !Array.isArray(v);

//----- Appends the line records for one key. `indent === 0` is top level, where scalars are
//      `kv` and maps are `k`; deeper, everything is a `nested` line -----//
function dump(lines, key, value, indent) {
  const top = indent === 0;

  if (Array.isArray(value)) {
    lines.push({ t: "nested", indent, k: key, v: "", kind: "bare" });
    if (value.length === 0) {
      lines.push({ t: "empty_list", indent: indent + 2 });
    } else {
      value.forEach((v) => lines.push({ t: "listitem", indent: indent + 2, v, kind: kindOf(v) }));
    }
    return;
  }

  if (isMap(value)) {
    lines.push(top ? { t: "k", k: key } : { t: "nested", indent, k: key, v: "", kind: "bare" });
    const entries = Object.entries(value);
    if (entries.length === 0) {
      lines.push({ t: "empty_map", indent: indent + 2 });
      return;
    }
    entries.forEach(([k, v]) => {
      if (isMap(v)) lines.push({ t: "blank" });
      dump(lines, k, v, indent + 2);
    });
    return;
  }

  lines.push(top
    ? { t: "kv", k: key, v: value, kind: kindOf(value) }
    : { t: "nested", indent, k: key, v: value, kind: kindOf(value) });
}

//----- Every top-level key opens a blank-line-separated block, except the shared scalars,
//      which ride along with `threads`. (`queries` needs no exception: it is written first,
//      so `lines.length` is still 0 when it is reached.) -----//
window.buildYAML = function buildYAML(cfg) {
  const lines = [];
  const shared = new Set(window.SHARED_NULLABLE_KEYS);
  Object.entries(window.configToObject(cfg)).forEach(([key, value]) => {
    if (lines.length && !shared.has(key)) lines.push({ t: "blank" });
    dump(lines, key, value, 0);
  });
  return lines;
};

//----- rendering — the drawer wants highlighted spans, the clipboard wants plain text -----//

const esc = (s) => String(s).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
const indent = (n) => " ".repeat(n);

//----- Turns line records into text. `key`/`val` decide the markup, so the highlighted and
//      the plain forms below share this one layout pass and cannot disagree -----//
function renderLines(lines, key, val) {
  return lines.map((L) => {
    if (L.t === "blank") return "";
    if (L.t === "k")  return `${key(L.k)}:`;
    if (L.t === "kv") return `${key(L.k)}: ${val(L.kind, L.v)}`;
    if (L.t === "empty_map")  return `${indent(L.indent)}{}`;
    if (L.t === "empty_list") return `${indent(L.indent)}[]`;
    if (L.t === "nested") {
      const head = `${indent(L.indent)}${key(L.k)}:`;
      return L.kind === "bare" ? head : `${head} ${val(L.kind, L.v)}`;
    }
    if (L.t === "listitem") return `${indent(L.indent)}- ${val(L.kind, L.v)}`;
    return "";
  }).join("\n");
}

//----- Syntax-highlighted HTML for the drawer panel -----//
window.renderYAML = function renderYAML(lines) {
  return renderLines(
    lines,
    (k) => `<span class="y-key">${esc(k)}</span>`,
    (kind, v) => {
      if (kind === "bare") return "";
      if (kind === "null") return `<span class="y-null">null</span>`;
      if (kind === "bool") return `<span class="y-bool">${v ? "true" : "false"}</span>`;
      if (kind === "num")  return `<span class="y-num">${esc(v)}</span>`;
      return `<span class="y-str">${esc(v)}</span>`;
    },
  );
};

//----- Plain text for the clipboard -----//
window.renderYAMLPlain = function renderYAMLPlain(lines) {
  return renderLines(
    lines,
    (k) => k,
    (kind, v) => {
      if (kind === "bare") return "";
      if (kind === "null") return "null";
      if (kind === "bool") return v ? "true" : "false";
      return String(v);
    },
  );
};
