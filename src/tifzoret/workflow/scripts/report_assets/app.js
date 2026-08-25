/* Tifzoret report — self-contained front-end. Vanilla JS, no dependencies.
   Reads the inlined JSON payload and renders three destinations:
     · Explore     — analysis "plates": each group's figures and its own tables
                     side by side, with a live volcano/MA drawn from the DE table
                     and a selected-gene link that highlights the gene everywhere.
     · Figures     — the assembled publication deliverables + full panel catalog.
     · Provenance  — run metadata, contrasts, warnings.
   Everything is offline; the file itself is the deliverable. */
(function () {
  "use strict";

  var DATA = JSON.parse(document.getElementById("tifzoret-data").textContent);
  var ROOT = document.documentElement;
  var app = document.getElementById("app");
  var SVGNS = "http://www.w3.org/2000/svg";

  /* ---------------------------------------------------------------- utils */
  function el(tag, props, children) {
    var node = document.createElement(tag);
    applyProps(node, props);
    if (children != null) append(node, children);
    return node;
  }
  function svg(tag, props, children) {
    var node = document.createElementNS(SVGNS, tag);
    applyProps(node, props);
    if (children != null) append(node, children);
    return node;
  }
  function applyProps(node, props) {
    if (!props) return;
    for (var k in props) {
      if (k === "class") node.setAttribute("class", props[k]);
      else if (k === "text") node.textContent = props[k];
      else if (k === "html") node.innerHTML = props[k];
      else if (k.slice(0, 2) === "on") node.addEventListener(k.slice(2).toLowerCase(), props[k]);
      else if (props[k] === true) node.setAttribute(k, "");
      else if (props[k] !== false && props[k] != null) node.setAttribute(k, props[k]);
    }
  }
  function append(node, children) {
    if (Array.isArray(children)) children.forEach(function (c) { append(node, c); });
    else if (children instanceof Node) node.appendChild(children);
    else if (children != null) node.appendChild(document.createTextNode(String(children)));
  }
  function clear(node) { while (node.firstChild) node.removeChild(node.firstChild); }
  function debounce(fn, ms) {
    var t; return function () { var a = arguments, c = this; clearTimeout(t); t = setTimeout(function () { fn.apply(c, a); }, ms); };
  }
  function cssVar(name) {
    try { return getComputedStyle(ROOT).getPropertyValue(name).trim() || FALLBACK[name] || "#888"; }
    catch (e) { return FALLBACK[name] || "#888"; }
  }
  var FALLBACK = { "--up": "#d1495b", "--down": "#1b6ca8", "--warn": "#b26a00", "--ns": "#6b7168" };

  var NUMERIC = { int: 1, num: 1, sci: 1 };
  function fmtSci(v) {
    if (v == null) return "—";
    if (v === "Inf") return "∞"; if (v === "-Inf") return "−∞";
    if (typeof v !== "number") return String(v);
    if (v === 0) return "0";
    var a = Math.abs(v);
    if (a < 1e-3 || a >= 1e5) return v.toExponential(1).replace("+", "");
    return parseFloat(v.toPrecision(3)).toString();
  }
  function fmtNum(v) {
    if (v == null) return "—";
    if (v === "Inf") return "∞"; if (v === "-Inf") return "−∞";
    if (typeof v !== "number") return String(v);
    var a = Math.abs(v);
    if (a >= 1000) return v.toLocaleString(undefined, { maximumFractionDigits: 0 });
    if (a >= 1) return v.toFixed(2);
    return v.toFixed(3);
  }
  function fmtInt(v) { return v == null ? "—" : Number(v).toLocaleString(); }
  function num(v) {
    if (v === "Inf") return Infinity; if (v === "-Inf") return -Infinity;
    return typeof v === "number" ? v : null;
  }
  function sortVal(v) {
    if (v == null) return [2, 0];
    if (v === "Inf") return [1, Infinity]; if (v === "-Inf") return [1, -Infinity];
    if (typeof v === "number") return [1, v];
    return [1, String(v).toLowerCase()];
  }
  function cmp(a, b) {
    var x = sortVal(a), y = sortVal(b);
    if (x[0] !== y[0]) return x[0] - y[0];
    return x[1] < y[1] ? -1 : x[1] > y[1] ? 1 : 0;
  }
  function splitGenes(s) { return String(s || "").split(/[,;/\s]+/).filter(Boolean); }
  function dirClass(v) { var s = String(v || "").toLowerCase(); return s.indexOf("up") === 0 ? "up" : s.indexOf("down") === 0 ? "down" : "ns"; }
  function dirGlyph(c) { return c === "up" ? "▲" : c === "down" ? "▼" : "•"; }
  function baseName(p) { return String(p || "figure").split("/").pop(); }

  /* --------------------------------------------------------------- theme */
  function initTheme() {
    var saved = null;
    try { saved = localStorage.getItem("tifzoret-theme"); } catch (e) {}
    if (!saved) saved = (window.matchMedia && window.matchMedia("(prefers-color-scheme: dark)").matches) ? "dark" : "light";
    ROOT.setAttribute("data-theme", saved);
  }
  function toggleTheme() {
    var next = ROOT.getAttribute("data-theme") === "dark" ? "light" : "dark";
    ROOT.setAttribute("data-theme", next);
    try { localStorage.setItem("tifzoret-theme", next); } catch (e) {}
    themeBtn.querySelector(".tt-label").textContent = next === "dark" ? "Light" : "Dark";
    liveCharts.forEach(function (c) { c.refresh(); });
  }

  /* --------------------------------------------------------------- shell */
  var meta = DATA.meta, counts = meta.counts;
  var views = {};
  var tabs = {};
  var themeBtn;
  var selectedGene = null;        // { name, key }
  var liveCharts = [];            // charts in the current workspace to refresh on gene-select
  var currentTable = null;        // { rehighlight } for the visible Explore table

  function tab(id, label, count) {
    var btn = el("button", { class: "tab", role: "tab", "aria-selected": id === "explore", onclick: function () { activate(id); } },
      [label, count != null ? el("span", { class: "tag", text: String(count) }) : null]);
    tabs[id] = btn; return btn;
  }
  function activate(id) {
    Object.keys(views).forEach(function (k) {
      views[k].hidden = k !== id;
      tabs[k].setAttribute("aria-selected", String(k === id));
    });
    if (id === "explore") lazyExplore();
    window.scrollTo({ top: 0 });
  }

  function buildShell() {
    var eng = meta.engine.name + (meta.engine.version ? " " + meta.engine.version : "");
    var omniInput = el("input", {
      type: "search", placeholder: "Find a gene or term across every analysis…",
      "aria-label": "Search all analyses", autocomplete: "off", spellcheck: false,
      oninput: debounce(function () { runOmni(omniInput.value.trim()); }, 130),
      onfocus: function () { if (omniInput.value.trim()) runOmni(omniInput.value.trim()); },
      onkeydown: function (e) { if (e.key === "Escape") { omniPanel.hidden = true; omniInput.blur(); } }
    });
    var omniPanel = el("div", { class: "omni-panel", hidden: true });
    var omni = el("div", { class: "omni" }, [
      el("span", { class: "prompt", text: "▸" }), omniInput,
      el("span", { class: "hint", text: "search" }), omniPanel
    ]);
    document.addEventListener("click", function (e) { if (!omni.contains(e.target)) omniPanel.hidden = true; });
    window._omni = { input: omniInput, panel: omniPanel };

    themeBtn = el("button", { class: "iconbtn", onclick: toggleTheme, "aria-label": "Toggle color theme" },
      [el("span", { class: "tt-label", text: ROOT.getAttribute("data-theme") === "dark" ? "Light" : "Dark" })]);

    var bar = el("header", { class: "appbar" }, el("div", { class: "appbar-inner" }, [
      el("div", { class: "brand" }, [
        el("span", { class: "eyebrow", text: eng + " · report" }),
        el("h1", { text: meta.title }),
        el("span", { class: "sub" }, [
          el("span", { class: "mono", text: meta.project_id }), " · ",
          el("span", { class: "mono", text: meta.analysis_set }), " · " + meta.profile + " profile"
        ])
      ]),
      omni,
      el("div", { class: "spacer" }),
      themeBtn
    ]));

    var tabbar = el("nav", { class: "tabbar" }, el("div", { class: "tabs", role: "tablist" }, [
      tab("explore", "Explore"),
      tab("figures", "Figures", counts.figures),
      tab("provenance", "Provenance")
    ]));

    views.explore = el("section", { class: "view", role: "tabpanel" });
    views.figures = el("section", { class: "view", role: "tabpanel", hidden: true });
    views.provenance = el("section", { class: "view", role: "tabpanel", hidden: true });

    app.setAttribute("aria-busy", "false");
    append(app, [bar, tabbar, views.explore, views.figures, views.provenance, lightbox()]);
    renderFigures();
    renderProvenance();
    lazyExplore();
  }

  /* =====================================================================
     EXPLORE — analysis plates (figures + tables together)
     ===================================================================== */
  var exploreBuilt = false;
  var groupsModel = [];           // [{ key,label,blurb,plate,tables,panels,contrasts }]
  var groupByKey = {};
  var tableById = {};             // id -> table object
  var tableLoc = {};              // id -> { groupKey, contrast }
  var wsState = { groupKey: null, contrast: {}, tableIdx: {} };

  function lazyExplore() { if (!exploreBuilt) { buildGroupsModel(); renderExplore(); exploreBuilt = true; } }

  function buildGroupsModel() {
    var panels = DATA.figures.panels || [];
    (DATA.analyses || []).forEach(function (g, i) {
      var gp = panels.filter(function (p) { return p.group === g.key; });
      var cset = {};
      g.tables.forEach(function (t) { if (t.contrast) cset[t.contrast] = 1; });
      gp.forEach(function (p) { if (p.contrast) cset[p.contrast] = 1; });
      // preserve the payload's contrast order where possible
      var order = (DATA.contrasts || []).map(function (c) { return c.contrast_id; });
      var contrasts = Object.keys(cset).sort(function (a, b) {
        var ia = order.indexOf(a), ib = order.indexOf(b);
        return (ia < 0 ? 99 : ia) - (ib < 0 ? 99 : ib);
      });
      var model = { key: g.key, label: g.label, blurb: g.blurb, plate: i + 1, tables: g.tables, panels: gp, contrasts: contrasts };
      groupsModel.push(model);
      groupByKey[g.key] = model;
      g.tables.forEach(function (t) { tableById[t.id] = t; tableLoc[t.id] = { groupKey: g.key, contrast: t.contrast }; });
    });
  }

  function renderExplore() {
    var v = views.explore; clear(v);
    if (!groupsModel.length) { append(v, el("div", { class: "empty-state", text: "No analyses were produced for this run." })); return; }

    var rail = el("nav", { class: "rail", "aria-label": "Analyses" }, el("div", { class: "rail-title", text: "Analyses" }));
    groupsModel.forEach(function (g) {
      var n = g.tables.reduce(function (s, t) { return s + t.total; }, 0);
      rail.appendChild(el("button", {
        class: "rail-link", "data-gid": g.key, "aria-current": false, onclick: function () { selectGroup(g.key); }
      }, [
        el("span", { class: "rl-no mono", text: pad2(g.plate) }),
        el("span", { class: "rl-name", text: g.label }),
        el("span", { class: "rl-count mono", text: n ? n.toLocaleString() : (g.panels.length + "▮") })
      ]));
    });

    var ws = el("div", { class: "workspace", id: "workspace" });
    append(v, el("div", { class: "explore" }, [rail, ws]));
    selectGroup(wsState.groupKey && groupByKey[wsState.groupKey] ? wsState.groupKey : groupsModel[0].key);
  }
  function pad2(n) { return (n < 10 ? "0" : "") + n; }

  function selectGroup(key) {
    wsState.groupKey = key;
    document.querySelectorAll(".rail-link").forEach(function (n) { n.setAttribute("aria-current", String(n.getAttribute("data-gid") === key)); });
    var active = document.querySelector('.rail-link[data-gid="' + cssEsc(key) + '"]');
    if (active) active.scrollIntoView({ block: "nearest" });
    renderWorkspace();
  }

  function currentContrast(g) {
    if (!g.contrasts.length) return null;
    var c = wsState.contrast[g.key];
    if (!c || g.contrasts.indexOf(c) < 0) c = g.contrasts[0];
    wsState.contrast[g.key] = c;
    return c;
  }

  function renderWorkspace() {
    var host = document.getElementById("workspace"); if (!host) return; clear(host);
    liveCharts = []; currentTable = null;
    var g = groupByKey[wsState.groupKey]; if (!g) return;
    var contrast = currentContrast(g);

    /* plate head */
    var head = el("div", { class: "ws-head" }, [
      el("div", { class: "plate-eyebrow" }, ["Plate", el("span", { class: "mono", text: pad2(g.plate) })]),
      el("h2", { class: "plate-title", text: g.label })
    ]);
    if (g.blurb) head.appendChild(el("p", { class: "plate-blurb", text: g.blurb }));
    head.appendChild(geneBar());
    host.appendChild(head);

    /* contrast switcher */
    if (g.contrasts.length > 1) {
      var seg = el("div", { class: "seg", role: "group", "aria-label": "Contrast" });
      g.contrasts.forEach(function (c) {
        seg.appendChild(el("button", { class: "seg-btn", "aria-pressed": c === contrast, onclick: function () { wsState.contrast[g.key] = c; renderWorkspace(); } },
          [el("span", { class: "mono", text: c })]));
      });
      host.appendChild(el("div", null, [el("span", { class: "seg-label", text: "Contrast" }), seg]));
    }

    /* figures for (group, contrast) */
    var figs = g.panels.filter(function (p) { return p.contrast == null || p.contrast === contrast; });
    var deTable = g.key === "de" ? tableForContrast(g, contrast) : null;
    if (deTable || figs.length) {
      host.appendChild(el("hr", { class: "specimen-rule" }));
      host.appendChild(el("h3", { class: "ws-block-title" }, ["Figures", el("span", { class: "count", text: figs.length + (deTable ? " + live plot" : "") })]));
      if (deTable) host.appendChild(scatterCard(deTable));
      if (figs.length) {
        var grid = el("div", { class: "fig-grid compact" });
        figs.forEach(function (p) { grid.appendChild(figCard(p, "panel", (DATA.figures.panels || []).indexOf(p))); });
        host.appendChild(grid);
      }
    }

    /* tables for (group, contrast) */
    var tables = g.tables.filter(function (t) { return t.contrast == null || t.contrast === contrast; });
    if (tables.length) {
      host.appendChild(el("hr", { class: "specimen-rule" }));
      host.appendChild(el("h3", { class: "ws-block-title" }, ["Data tables", el("span", { class: "count", text: String(tables.length) })]));
      var idx = wsState.tableIdx[g.key + "|" + contrast];
      if (idx == null || idx >= tables.length) idx = 0;
      if (tables.length > 1) {
        var pick = el("div", { class: "seg", role: "group", "aria-label": "Table" });
        tables.forEach(function (t, i) {
          pick.appendChild(el("button", { class: "seg-btn", "aria-pressed": i === idx, onclick: function () { wsState.tableIdx[g.key + "|" + contrast] = i; renderWorkspace(); } },
            [t.label, el("span", { class: "mono", text: " " + t.total.toLocaleString() })]));
        });
        host.appendChild(pick);
      }
      var wrap = el("div", { class: "tablewrap", id: "ws-table" });
      host.appendChild(wrap);
      drawTable(wrap, tables[idx], "");
    } else if (!deTable && !figs.length) {
      host.appendChild(el("div", { class: "empty-state", text: "This analysis produced no browsable tables or panels." }));
    }
  }

  function tableForContrast(g, contrast) {
    return g.tables.filter(function (t) { return t.contrast == null || t.contrast === contrast; })[0] || null;
  }

  /* selected-gene link bar */
  function geneBar() {
    if (!selectedGene) return el("span", { hidden: true });
    var where = countGeneTables(selectedGene.key);
    return el("div", { class: "genebar" }, [
      el("span", { class: "gb-label", text: "Selected gene" }),
      el("span", { class: "gb-name mono", text: selectedGene.name }),
      el("span", { class: "gb-where", text: where ? "found in " + where + " result table" + (where > 1 ? "s" : "") + " — highlighted below" : "not found in the embedded tables" }),
      el("button", { class: "gb-clear", onclick: clearGene, text: "Clear" })
    ]);
  }
  function countGeneTables(key) {
    var needle = key.toLowerCase(), n = 0;
    (DATA.analyses || []).forEach(function (g) {
      g.tables.forEach(function (t) {
        var blobs = t.search_blobs;
        for (var i = 0; i < blobs.length; i++) { if (blobs[i].indexOf(needle) >= 0) { n++; break; } }
      });
    });
    return n;
  }
  function selectGene(name) {
    name = String(name || "").trim(); if (!name) return;
    selectedGene = { name: name, key: name.toUpperCase() };
    var gb = document.querySelector(".ws-head");
    if (gb) { var old = gb.querySelector(".genebar"); var nu = geneBar(); if (old) gb.replaceChild(nu, old); else gb.appendChild(nu); }
    liveCharts.forEach(function (c) { c.refresh(); });
    if (currentTable) currentTable.rehighlight();
  }
  function clearGene() {
    selectedGene = null;
    var gb = document.querySelector(".ws-head .genebar"); if (gb) gb.parentNode.replaceChild(el("span", { hidden: true }), gb);
    liveCharts.forEach(function (c) { c.refresh(); });
    if (currentTable) currentTable.rehighlight();
  }
  function isSelected(token) { return selectedGene && String(token).toUpperCase() === selectedGene.key; }

  /* ------------------------------------------------- live volcano / MA plate */
  function colClass(cls) {
    var s = String(cls || "").toLowerCase();
    if (s.indexOf("up") >= 0) return ["up", cssVar("--up")];
    if (s.indexOf("down") >= 0) return ["down", cssVar("--down")];
    if (s === "" || s.indexOf("ns") >= 0 || s.indexOf("not") >= 0 || s.indexOf("none") >= 0 || s.indexOf("n.s") >= 0) return ["ns", cssVar("--ns")];
    return ["other", cssVar("--warn")];
  }
  function computePoints(t) {
    if (t._points) return t._points;
    var ix = {}; t.columns.forEach(function (c, i) { ix[c.key] = i; });
    var iSym = ix.gene_symbol, iId = ix.gene_id, iBase = ix.base_mean,
        iL = ix.log2_fold_change, iPadj = ix.adjusted_p_value, iP = ix.p_value,
        iCls = ix.significance_class, iDir = ix.direction;
    var pts = [], classes = {};
    t.rows.forEach(function (r) {
      var l = num(iL != null ? r[iL] : null);
      if (l == null || !isFinite(l)) return;
      var pRaw = iPadj != null ? r[iPadj] : (iP != null ? r[iP] : null);
      var p = num(pRaw);
      var base = iBase != null ? num(r[iBase]) : null;
      var clsVal = iCls != null ? r[iCls] : null;
      var cc = colClass(clsVal != null && clsVal !== "" ? clsVal : (iDir != null ? r[iDir] : ""));
      var lbl = (clsVal != null && clsVal !== "") ? String(clsVal) : cc[0];
      classes[lbl] = cc[1];
      pts.push({
        gene: (iSym != null && r[iSym]) ? r[iSym] : (iId != null ? r[iId] : ""),
        id: iId != null ? r[iId] : "",
        l2fc: l,
        negLogP: (p != null && p > 0) ? -Math.log10(p) : (p === 0 ? Infinity : null),
        logBase: (base != null && base > 0) ? Math.log10(base) : null,
        color: cc[1], klass: cc[0], clsLabel: lbl
      });
    });
    t._points = { pts: pts, classes: classes };
    return t._points;
  }

  function scatterCard(t) {
    var info = computePoints(t);
    if (!info.pts.length) return el("div", { hidden: true });
    var mode = "volcano";
    var stage = el("div");
    var note = el("p", { class: "chart-note" });
    var legend = el("div", { class: "chart-legend" });
    var toggle = el("div", { class: "chart-toggle" });
    var vBtn = el("button", { class: "seg-btn", "aria-pressed": true, text: "Volcano", onclick: function () { mode = "volcano"; sync(); } });
    var mBtn = el("button", { class: "seg-btn", "aria-pressed": false, text: "MA", onclick: function () { mode = "ma"; sync(); } });
    append(toggle, [vBtn, mBtn]);
    var card = el("div", { class: "chartcard" }, [
      el("div", { class: "chart-head" }, [
        el("h4", { text: t.contrast ? t.contrast : "Differential expression" }),
        toggle
      ]),
      stage, note, legend
    ]);
    function sync() {
      vBtn.setAttribute("aria-pressed", String(mode === "volcano"));
      mBtn.setAttribute("aria-pressed", String(mode === "ma"));
      draw();
    }
    function draw() {
      clear(stage); clear(legend);
      var made = drawScatter(info, mode);
      stage.appendChild(made.svg);
      note.textContent = made.note;
      Object.keys(info.classes).forEach(function (k) {
        legend.appendChild(el("span", { class: "lg" }, [
          el("span", { class: "sw", style: "background:" + info.classes[k] }), k
        ]));
      });
    }
    liveCharts.push({ refresh: draw });
    draw();
    return card;
  }

  var MAXPTS = 9000;
  function drawScatter(info, mode) {
    var W = 760, H = 420, ml = 56, mr = 18, mt = 16, mb = 46;
    var pts = info.pts.filter(function (p) {
      return mode === "volcano" ? (p.negLogP != null) : (p.logBase != null);
    });
    var total = pts.length, note = "";
    // Cap for smooth rendering: keep every non-ns point + a deterministic sample of ns.
    if (total > MAXPTS) {
      var strong = [], weak = [];
      pts.forEach(function (p) { (p.klass === "ns" ? weak : strong).push(p); });
      var room = Math.max(0, MAXPTS - strong.length);
      var stride = weak.length > room && room > 0 ? Math.ceil(weak.length / room) : 1;
      var kept = [];
      for (var i = 0; i < weak.length; i += stride) kept.push(weak[i]);
      // always include the selected gene if present
      if (selectedGene) weak.forEach(function (p) { if (isSelected(p.gene) || isSelected(p.id)) kept.push(p); });
      pts = strong.concat(kept);
      note = pts.length.toLocaleString() + " of " + total.toLocaleString() + " genes plotted (all classified + a sample of the rest). Every gene is in the table below and the publication PNG.";
    } else {
      note = total.toLocaleString() + " genes plotted from the DE table.";
    }

    var xf = mode === "volcano" ? function (p) { return p.l2fc; } : function (p) { return p.logBase; };
    var yf = mode === "volcano" ? function (p) { return isFinite(p.negLogP) ? p.negLogP : null; } : function (p) { return p.l2fc; };
    var xs = [], ys = [];
    pts.forEach(function (p) { xs.push(xf(p)); var y = yf(p); if (y != null) ys.push(y); });
    var xlo = Math.min.apply(null, xs), xhi = Math.max.apply(null, xs);
    var ylo = mode === "volcano" ? 0 : Math.min.apply(null, ys), yhi = Math.max.apply(null, ys);
    if (mode === "volcano") { var m = Math.max(Math.abs(xlo), Math.abs(xhi)); xlo = -m; xhi = m; }
    var xp = (xhi - xlo) * 0.04 || 1, yp = (yhi - ylo) * 0.04 || 1;
    xlo -= xp; xhi += xp; if (mode === "ma") { ylo -= yp; } yhi += yp;
    var yInf = mode === "volcano" ? yhi : null;  // p==0 rows clamp to top

    function sx(v) { return ml + (v - xlo) / (xhi - xlo) * (W - ml - mr); }
    function sy(v) { return (H - mb) - (v - ylo) / (yhi - ylo) * (H - mb - mt); }

    var root = svg("svg", { class: "chart-svg", viewBox: "0 0 " + W + " " + H, role: "img",
      "aria-label": (mode === "volcano" ? "Volcano plot" : "MA plot") + " of differential expression" });

    // gridlines + ticks
    var xt = niceTicks(xlo, xhi, 7), yt = niceTicks(ylo, yhi, 6);
    xt.forEach(function (v) {
      root.appendChild(svg("line", { class: "grid-line", x1: sx(v), y1: mt, x2: sx(v), y2: H - mb }));
      root.appendChild(svg("text", { x: sx(v), y: H - mb + 16, "text-anchor": "middle", "font-size": "11" }, fmtTick(v)));
    });
    yt.forEach(function (v) {
      root.appendChild(svg("line", { class: "grid-line", x1: ml, y1: sy(v), x2: W - mr, y2: sy(v) }));
      root.appendChild(svg("text", { x: ml - 8, y: sy(v) + 3, "text-anchor": "end", "font-size": "11" }, fmtTick(v)));
    });
    // reference line (x=0 for volcano centre / y=0 for MA)
    if (mode === "volcano" && xlo < 0 && xhi > 0) root.appendChild(svg("line", { class: "thr-line", x1: sx(0), y1: mt, x2: sx(0), y2: H - mb }));
    if (mode === "ma" && ylo < 0 && yhi > 0) root.appendChild(svg("line", { class: "thr-line", x1: ml, y1: sy(0), x2: W - mr, y2: sy(0) }));
    // axes
    root.appendChild(svg("line", { class: "axis-line", x1: ml, y1: H - mb, x2: W - mr, y2: H - mb }));
    root.appendChild(svg("line", { class: "axis-line", x1: ml, y1: mt, x2: ml, y2: H - mb }));
    root.appendChild(svg("text", { class: "axis-title", x: (ml + W - mr) / 2, y: H - 8, "text-anchor": "middle", "font-size": "12" },
      mode === "volcano" ? "log₂ fold change" : "log₁₀ mean expression"));
    root.appendChild(svg("text", { class: "axis-title", x: 14, y: (mt + H - mb) / 2, "text-anchor": "middle", "font-size": "12",
      transform: "rotate(-90 14 " + ((mt + H - mb) / 2) + ")" }, mode === "volcano" ? "−log₁₀ FDR" : "log₂ fold change"));

    // points
    var selPt = null;
    var frag = document.createDocumentFragment();
    pts.forEach(function (p) {
      var y = yf(p); if (y == null) y = yInf; if (y == null) return;
      var sel = isSelected(p.gene) || (p.id && isSelected(p.id));
      var c = svg("circle", { class: "pt" + (sel ? " sel" : ""), cx: sx(xf(p)), cy: sy(y), r: sel ? 5 : 2.4,
        fill: p.color, "fill-opacity": p.klass === "ns" ? .5 : .82 });
      c.addEventListener("mousemove", function (e) { showChartTip(e, p, mode); });
      c.addEventListener("mouseleave", hideChartTip);
      c.addEventListener("click", function () { selectGene(p.gene || p.id); });
      if (sel) selPt = { x: sx(xf(p)), y: sy(y), g: p.gene || p.id }; else frag.appendChild(c);
    });
    root.appendChild(frag);
    if (selPt) {
      root.appendChild(svg("circle", { class: "pt sel", cx: selPt.x, cy: selPt.y, r: 5, fill: cssVar("--up"), "fill-opacity": .9 }));
      root.appendChild(svg("text", { class: "sel-label", x: selPt.x + 7, y: selPt.y - 6, "font-size": "11" }, selPt.g));
    }
    return { svg: root, note: note };
  }

  var chartTip = el("div", { class: "chart-tooltip", hidden: true });
  function showChartTip(e, p, mode) {
    if (chartTip.parentNode !== document.body) document.body.appendChild(chartTip);
    clear(chartTip);
    append(chartTip, [
      el("div", null, el("span", { class: "g", text: p.gene || p.id || "—" })),
      el("div", null, ["log₂FC ", el("span", { class: "n", text: fmtNum(p.l2fc) })]),
      el("div", null, ["FDR ", el("span", { class: "n", text: fmtSci(p.negLogP != null && isFinite(p.negLogP) ? Math.pow(10, -p.negLogP) : (p.negLogP === Infinity ? 0 : null)) })]),
      el("div", null, el("span", { text: p.clsLabel }))
    ]);
    chartTip.hidden = false;
    var x = e.clientX + 14, y = e.clientY + 14;
    chartTip.style.left = Math.min(x, window.innerWidth - 250) + "px";
    chartTip.style.top = Math.min(y, window.innerHeight - 90) + "px";
  }
  function hideChartTip() { chartTip.hidden = true; }

  function niceTicks(lo, hi, n) {
    if (!(hi > lo)) return [lo];
    var span = hi - lo, step = Math.pow(10, Math.floor(Math.log10(span / n)));
    var err = (n * step) / span;
    if (err <= 0.15) step *= 10; else if (err <= 0.35) step *= 5; else if (err <= 0.75) step *= 2;
    var out = [], start = Math.ceil(lo / step) * step;
    for (var v = start; v <= hi + step * 0.001; v += step) out.push(Math.round(v / step) * step);
    return out;
  }
  function fmtTick(v) { if (v === 0) return "0"; var a = Math.abs(v); if (a >= 1000 || a < 0.01) return v.toExponential(0); return (Math.round(v * 100) / 100).toString(); }

  /* =====================================================================
     Figures view (assembled deliverables + full catalog) + lightbox
     ===================================================================== */
  var figState = { filter: "all", query: "" };

  function renderFigures() {
    var v = views.figures; clear(v);
    var asm = DATA.figures.assembled || [];
    var panels = DATA.figures.panels || [];
    if (!asm.length && !panels.length) {
      append(v, el("div", { class: "empty-state", text: "No figures were produced for this run." }));
      return;
    }
    var filters = [["all", "All"], ["assembled", "Assembled figures"], ["selected", "Selected panels"], ["catalog", "Full catalog"]];
    var controls = el("div", { class: "fig-controls" });
    filters.forEach(function (pair) {
      if (pair[0] === "assembled" && !asm.length) return;
      controls.appendChild(el("button", {
        class: "chipfilter", "aria-pressed": figState.filter === pair[0], text: pair[1],
        onclick: function () { figState.filter = pair[0]; renderFigures(); }
      }));
    });
    var search = el("input", {
      type: "search", placeholder: "Filter figures…", value: figState.query, "aria-label": "Filter figures",
      oninput: debounce(function (e) { figState.query = e.target.value.trim().toLowerCase(); drawFigGrids(); }, 120)
    });
    controls.appendChild(el("div", { class: "fig-search" }, search));
    append(v, controls);
    v.appendChild(el("div", { id: "fig-sections" }));
    drawFigGrids();
  }

  function figMatches(f) {
    if (!figState.query) return true;
    var hay = [f.label, f.constructor, f.constructor_label, f.variant_label, f.description, f.set].join(" ").toLowerCase();
    return hay.indexOf(figState.query) >= 0;
  }
  function drawFigGrids() {
    var host = document.getElementById("fig-sections"); if (!host) return; clear(host);
    var asm = (DATA.figures.assembled || []).filter(figMatches);
    var panels = (DATA.figures.panels || []).filter(figMatches);
    var showAsm = figState.filter === "all" || figState.filter === "assembled";
    var showPanels = figState.filter === "all" || figState.filter === "selected" || figState.filter === "catalog";
    var panelSet = figState.filter === "selected" ? panels.filter(function (p) { return p.selected; }) : panels;

    if (showAsm && asm.length) {
      host.appendChild(el("h2", { class: "section-title", text: "Assembled figures" }));
      host.appendChild(el("p", { class: "section-blurb", text: "Publication-ready multi-panel figures. Download the vector PDF for the manuscript, or the PNG for slides." }));
      var g1 = el("div", { class: "fig-grid assembled" });
      asm.forEach(function (f) { g1.appendChild(figCard(f, "assembled", asm.indexOf(f))); });
      host.appendChild(g1);
    }
    if (showPanels && panelSet.length) {
      var title = figState.filter === "selected" ? "Selected panels" : "Panel catalog";
      host.appendChild(el("h2", { class: "section-title", text: title }));
      host.appendChild(el("p", { class: "section-blurb", text: "Every panel the engine rendered — including alternatives not placed into a figure. A badge marks the ones used." }));
      var g2 = el("div", { class: "fig-grid" });
      panelSet.forEach(function (f) { g2.appendChild(figCard(f, "panel", (DATA.figures.panels || []).indexOf(f))); });
      host.appendChild(g2);
    }
    if (!host.children.length) host.appendChild(el("div", { class: "empty-state", text: "No figures match “" + figState.query + "”." }));
  }

  function figCard(f, kind, idx) {
    var thumb = el("div", { class: "fig-thumb", role: "button", tabindex: "0", "aria-label": "Enlarge " + (f.label || "figure"),
      onclick: function () { openLightbox(kind, idx); },
      onkeydown: function (e) { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); openLightbox(kind, idx); } }
    });
    if (f.png) thumb.appendChild(el("img", { src: f.png, alt: f.label || "", loading: "lazy" }));
    else thumb.appendChild(el("div", { class: "missing", text: "image unavailable" }));
    if (kind === "panel" && f.selected) thumb.appendChild(el("span", { class: "fig-badge", text: "used" }));
    if (kind === "assembled") thumb.appendChild(el("span", { class: "fig-badge", text: "figure" }));

    var body = el("div", { class: "fig-body" });
    body.appendChild(el("h3", { text: f.label || f.constructor_label || f.constructor || "Figure" }));
    var kindText = kind === "assembled" ? (f.set || "") : (f.constructor_label || f.constructor || "");
    if (f.variant_label && f.variant_label !== f.constructor_label) kindText += " · " + f.variant_label;
    if (f.contrast) kindText += " · " + f.contrast;
    if (kindText) body.appendChild(el("span", { class: "kind", text: kindText }));
    if (f.description) body.appendChild(el("p", { text: f.description }));

    var actions = el("div", { class: "fig-actions" });
    actions.appendChild(downloadLink("PNG", f.png, f.png_href || (f.id || f.constructor || "figure"), true));
    actions.appendChild(pdfLink(f.pdf_href));
    body.appendChild(actions);
    return el("div", { class: "fig-card" }, [thumb, body]);
  }

  function downloadLink(label, dataUri, hrefName, isData) {
    if (isData && dataUri) {
      var name = baseName(hrefName).replace(/\.[^.]+$/, "") + ".png";
      return el("a", { class: "dl", href: dataUri, download: name, text: label });
    }
    return el("a", { class: "dl", "aria-disabled": "true", text: label });
  }
  function pdfLink(href) {
    if (href) return el("a", { class: "dl", href: href, download: baseName(href), text: "PDF" });
    return el("a", { class: "dl", "aria-disabled": "true", text: "PDF" });
  }

  /* lightbox */
  var lb = {}, lbList = [], lbPos = 0;
  function lightbox() {
    lb.title = el("h3");
    lb.img = el("img", { alt: "", onclick: function () { lb.stage.classList.toggle("zoom"); } });
    lb.stage = el("div", { class: "lb-stage" }, lb.img);
    lb.png = el("a", { class: "dl", text: "PNG" });
    lb.pdf = el("a", { class: "dl", text: "PDF" });
    lb.node = el("div", { class: "lightbox", hidden: true }, [
      el("div", { class: "lb-bar" }, [
        el("button", { class: "lb-nav", text: "‹", "aria-label": "Previous", onclick: function () { stepLightbox(-1); } }),
        el("button", { class: "lb-nav", text: "›", "aria-label": "Next", onclick: function () { stepLightbox(1); } }),
        lb.title, lb.png, lb.pdf,
        el("button", { class: "lb-x", text: "✕ Close", onclick: closeLightbox })
      ]),
      lb.stage
    ]);
    lb.node.addEventListener("click", function (e) { if (e.target === lb.node || e.target === lb.stage) closeLightbox(); });
    document.addEventListener("keydown", function (e) {
      if (lb.node.hidden) return;
      if (e.key === "Escape") closeLightbox();
      else if (e.key === "ArrowLeft") stepLightbox(-1);
      else if (e.key === "ArrowRight") stepLightbox(1);
    });
    return lb.node;
  }
  function visibleFigures() {
    var out = [];
    (DATA.figures.assembled || []).forEach(function (f, i) { if (figMatches(f) && (figState.filter === "all" || figState.filter === "assembled")) out.push({ kind: "assembled", i: i, f: f }); });
    var panelShow = figState.filter === "all" || figState.filter === "selected" || figState.filter === "catalog";
    (DATA.figures.panels || []).forEach(function (f, i) {
      if (!panelShow || !figMatches(f)) return;
      if (figState.filter === "selected" && !f.selected) return;
      out.push({ kind: "panel", i: i, f: f });
    });
    return out;
  }
  function openLightbox(kind, idx) {
    // when opened from Explore, the Figures grid may be filtered — fall back to a
    // single-item list so the panel always opens.
    lbList = visibleFigures();
    lbPos = lbList.findIndex(function (x) { return x.kind === kind && x.i === idx; });
    if (lbPos < 0) { var arr = kind === "assembled" ? DATA.figures.assembled : DATA.figures.panels; lbList = [{ kind: kind, i: idx, f: arr[idx] }]; lbPos = 0; }
    showLightbox();
    lb.node.hidden = false;
  }
  function showLightbox() {
    var item = lbList[lbPos]; if (!item) return;
    var f = item.f;
    lb.stage.classList.remove("zoom");
    lb.img.src = f.png || "";
    lb.img.alt = f.label || "";
    lb.title.textContent = (f.label || f.constructor_label || "Figure") + (lbList.length > 1 ? "  (" + (lbPos + 1) + "/" + lbList.length + ")" : "");
    if (f.png) { lb.png.href = f.png; lb.png.setAttribute("download", baseName(f.png_href || f.id || "figure").replace(/\.[^.]+$/, "") + ".png"); lb.png.removeAttribute("aria-disabled"); }
    else lb.png.setAttribute("aria-disabled", "true");
    if (f.pdf_href) { lb.pdf.href = f.pdf_href; lb.pdf.setAttribute("download", baseName(f.pdf_href)); lb.pdf.removeAttribute("aria-disabled"); }
    else { lb.pdf.removeAttribute("href"); lb.pdf.setAttribute("aria-disabled", "true"); }
  }
  function stepLightbox(d) { if (!lbList.length) return; lbPos = (lbPos + d + lbList.length) % lbList.length; showLightbox(); }
  function closeLightbox() { lb.node.hidden = true; lb.img.src = ""; }

  /* =====================================================================
     Virtualized result table (both-axis scroll, gene highlighting)
     ===================================================================== */
  function cssEsc(s) { return String(s).replace(/["\\]/g, "\\$&"); }

  function colWidth(c, isFirstText) {
    if (c.kind === "direction") return 64;
    if (c.kind === "int") return 82;
    if (c.kind === "num") return 100;
    if (c.kind === "sci") return 98;
    if (c.kind === "list") return 300;
    if (c.kind === "chip") return 132;
    if (c.key === "gene_id") return 168;
    return isFirstText ? 190 : 150;   // first text column (Gene / Term / Gene set) gets room
  }
  function gridMetrics(columns) {
    var firstText = true, widths = columns.map(function (c) {
      var isText = !NUMERIC[c.kind] && c.kind !== "chip" && c.kind !== "direction" && c.kind !== "list";
      var w = colWidth(c, isText && firstText); if (isText && firstText) firstText = false; return w;
    });
    return { cols: widths.map(function (w) { return w + "px"; }).join(" "), total: widths.reduce(function (a, b) { return a + b; }, 0) };
  }

  function drawTable(host, t, presetFilter) {
    clear(host);
    var ROWH = 32;
    var st = { sortCol: t.sort.col, sortDir: t.sort.dir, filter: (presetFilter || "").toLowerCase(), facet: "" };

    var titleBox = el("div", { class: "th-title" }, [
      el("h2", { text: t.label + (t.contrast ? " — " + t.contrast : "") }),
      el("span", { class: "th-note" }, t.full_tsv
        ? [el("span", { class: "mono", text: t.note }), " · ", el("a", { href: t.full_tsv, download: baseName(t.full_tsv), text: "full TSV" })]
        : [el("span", { class: "mono", text: t.note })])
    ]);
    var filterInput = el("input", {
      type: "search", placeholder: "Filter rows…", value: presetFilter || "", "aria-label": "Filter " + t.label,
      oninput: debounce(function (e) { st.filter = e.target.value.trim().toLowerCase(); recompute(); }, 110)
    });
    var head = el("div", { class: "table-head" }, [titleBox]);
    if (t.facet_col != null) {
      var vals = {}; t.rows.forEach(function (r) { var x = r[t.facet_col]; if (x) vals[x] = 1; });
      var sel = el("select", { "aria-label": "Filter by " + t.columns[t.facet_col].label,
        onchange: function (e) { st.facet = e.target.value; recompute(); } },
        [el("option", { value: "", text: "All " + t.columns[t.facet_col].label.toLowerCase() })]
          .concat(Object.keys(vals).sort().map(function (x) { return el("option", { value: x, text: x }); })));
      head.appendChild(el("div", { class: "facet" }, sel));
    }
    head.appendChild(el("div", { class: "tfilter" }, filterInput));
    var rowcount = el("div", { class: "rowcount" });
    head.appendChild(rowcount);
    host.appendChild(head);

    var metrics = gridMetrics(t.columns);
    var thead = el("div", { class: "thead" });
    t.columns.forEach(function (c, i) {
      var isNum = NUMERIC[c.kind];
      var th = el("div", { class: "th" + (isNum ? " numeric" : ""), role: "columnheader", title: "Sort by " + c.label, onclick: function () { setSort(i); } },
        [el("span", { text: c.label }), el("span", { class: "arrow" })]);
      th.dataset.col = i;
      thead.appendChild(th);
    });
    var vbody = el("div", { class: "vbody" });
    var vscroll = el("div", { class: "vscroll" }, [thead, vbody]);
    vscroll.style.setProperty("--cols", metrics.cols);
    vscroll.style.setProperty("--tw", metrics.total + "px");
    host.appendChild(vscroll);

    var view = [];
    function recompute() {
      var f = st.filter, facet = st.facet, fc = t.facet_col;
      view = [];
      for (var i = 0; i < t.rows.length; i++) {
        if (facet && t.rows[i][fc] !== facet) continue;
        if (f && t.search_blobs[i].indexOf(f) < 0) continue;
        view.push(i);
      }
      var col = st.sortCol, dir = st.sortDir === "asc" ? 1 : -1;
      view.sort(function (a, b) { return dir * cmp(t.rows[a][col], t.rows[b][col]); });
      vbody.style.height = (view.length * ROWH) + "px";
      rowcount.textContent = t.capped
        ? view.length.toLocaleString() + " shown · " + t.rows.length.toLocaleString() + " loaded of " + t.total.toLocaleString()
        : view.length.toLocaleString() + " of " + t.total.toLocaleString();
      updateArrows();
      draw();
    }
    function updateArrows() {
      thead.querySelectorAll(".th").forEach(function (th) {
        var i = +th.dataset.col;
        if (i === st.sortCol) { th.setAttribute("data-sorted", st.sortDir); th.querySelector(".arrow").textContent = st.sortDir === "asc" ? "▲" : "▼"; }
        else { th.removeAttribute("data-sorted"); th.querySelector(".arrow").textContent = ""; }
      });
    }
    function setSort(i) {
      if (st.sortCol === i) st.sortDir = st.sortDir === "asc" ? "desc" : "asc";
      else { st.sortCol = i; st.sortDir = NUMERIC[t.columns[i].kind] ? "desc" : "asc"; }
      recompute();
    }
    function draw() {
      var offset = vbody.offsetTop;
      var top = vscroll.scrollTop - offset, h = vscroll.clientHeight;
      var first = Math.max(0, Math.floor(top / ROWH) - 6);
      var last = Math.min(view.length, Math.ceil((top + h) / ROWH) + 6);
      clear(vbody);
      for (var p = first; p < last; p++) vbody.appendChild(rowNode(t, view[p], p * ROWH));
    }
    vscroll.addEventListener("scroll", function () { window.requestAnimationFrame(draw); });
    window.addEventListener("resize", debounce(draw, 100));

    recompute();
    t._apply = function (q) { filterInput.value = q; st.filter = String(q).toLowerCase(); recompute(); filterInput.focus(); };
    currentTable = { rehighlight: draw, id: t.id };
  }

  function rowNode(t, ri, top) {
    var row = t.rows[ri];
    var hit = false;
    var tr = el("div", { class: "trow" });
    tr.style.top = top + "px";
    t.columns.forEach(function (c, i) {
      var cell = cellNode(c, row[i]);
      if (cell._hit) hit = true;
      tr.appendChild(cell);
    });
    if (hit) tr.classList.add("rowsel");
    return tr;
  }

  function cellNode(c, v) {
    var node;
    if (c.kind === "direction") {
      var dc = dirClass(v);
      return el("div", { class: "cell" }, el("span", { class: "dir " + dc, title: String(v == null ? "" : v) }, dirGlyph(dc)));
    }
    if (c.kind === "chip") {
      if (v == null || v === "") return el("div", { class: "cell muted", text: "—" });
      var s = String(v).toLowerCase();
      var extra = s.indexOf("up") >= 0 ? " class-up" : s.indexOf("down") >= 0 ? " class-down" : "";
      return el("div", { class: "cell" }, el("span", { class: "chip" + extra, text: String(v) }));
    }
    if (c.kind === "list") return listCell(c, v);
    if (NUMERIC[c.kind]) {
      var txt = c.kind === "sci" ? fmtSci(v) : c.kind === "int" ? fmtInt(v) : fmtNum(v);
      return el("div", { class: "cell numeric" + (v == null ? " muted" : ""), text: txt });
    }
    // text — gene-ish columns are clickable + highlightable
    var isGene = c.key === "gene_symbol" || c.key === "gene_id" || c.key === "regulator";
    if (v == null || v === "") return el("div", { class: "cell muted", text: "—" });
    if (isGene && isSelected(v)) {
      node = el("div", { class: "cell gene", title: String(v), onclick: function () { selectGene(v); } }, el("span", { class: "ghit", text: String(v) }));
      node._hit = true; return node;
    }
    node = el("div", { class: "cell" + (isGene ? " gene" : ""), title: String(v), text: String(v) });
    if (isGene) node.addEventListener("click", function () { selectGene(v); });
    return node;
  }

  function listCell(c, v) {
    var genes = splitGenes(v);
    if (!genes.length) { var m = el("div", { class: "cell muted", text: "—" }); return m; }
    // surface the selected gene even if it would fall outside the visible slice
    var ordered = genes, hit = false;
    if (selectedGene) {
      var found = genes.filter(function (g) { return isSelected(g); });
      if (found.length) { hit = true; ordered = found.concat(genes.filter(function (g) { return !isSelected(g); })); }
    }
    var shown = ordered.slice(0, 4);
    var span = el("span", { class: "genes", title: genes.join(" ") });
    shown.forEach(function (g, i) {
      if (i) span.appendChild(document.createTextNode(", "));
      span.appendChild(isSelected(g) ? el("span", { class: "ghit", text: g }) : document.createTextNode(g));
    });
    var wrap = el("div", { class: "cell listcell" }, span);
    if (genes.length > shown.length) {
      wrap.appendChild(el("button", { class: "more", text: "+" + (genes.length - shown.length),
        onclick: function (e) { e.stopPropagation(); showPopover(e.currentTarget, c.label, genes); } }));
    }
    wrap._hit = hit;
    return wrap;
  }

  /* gene-list popover */
  var popover = el("div", { class: "popover", hidden: true });
  document.addEventListener("click", function (e) { if (!popover.contains(e.target) && !e.target.classList.contains("more")) popover.hidden = true; });
  function showPopover(anchor, label, genes) {
    clear(popover);
    popover.appendChild(el("div", { class: "pv-title", text: label + " (" + genes.length + ")" }));
    genes.forEach(function (g) {
      popover.appendChild(el("span", { class: "g" + (isSelected(g) ? " ghit" : ""), text: g, onclick: function () { selectGene(g); popover.hidden = true; } }));
    });
    if (popover.parentNode !== document.body) document.body.appendChild(popover);
    popover.hidden = false;
    var r = anchor.getBoundingClientRect(), pw = 320;
    popover.style.left = Math.max(8, Math.min(r.left, window.innerWidth - pw - 8)) + "px";
    popover.style.top = (r.bottom + 6) + "px";
  }

  /* =====================================================================
     Omni-search — cross-analysis finder + gene locator
     ===================================================================== */
  function runOmni(q) {
    var o = window._omni, panel = o.panel;
    if (!q || q.length < 2) { panel.hidden = true; return; }
    var needle = q.toLowerCase();
    clear(panel);
    var any = false;

    // gene locator — set the selected gene and light it up everywhere
    panel.appendChild(el("button", { class: "omni-hit locate", onclick: function () { panel.hidden = true; o.input.value = ""; activate("explore"); selectGene(q); } }, [
      el("span", { class: "where", text: "Highlight “" + q + "” across every analysis" }),
      el("span", { class: "count", text: "link" })
    ]));

    var figHits = allFiguresList().filter(function (x) {
      return [x.f.label, x.f.constructor, x.f.constructor_label, x.f.variant_label, x.f.description, x.f.set].join(" ").toLowerCase().indexOf(needle) >= 0;
    });
    if (figHits.length) {
      any = true;
      panel.appendChild(el("div", { class: "omni-group-title", text: "Figures" }));
      panel.appendChild(el("button", { class: "omni-hit", onclick: function () { panel.hidden = true; figState.query = needle; activate("figures"); renderFigures(); } }, [
        el("span", { class: "where", text: figHits.length + " figure" + (figHits.length > 1 ? "s" : "") + " match" }),
        el("span", { class: "count", text: "view" })
      ]));
    }

    (DATA.analyses || []).forEach(function (g) {
      var groupHits = [];
      g.tables.forEach(function (t) {
        var n = 0, blobs = t.search_blobs;
        for (var i = 0; i < blobs.length; i++) if (blobs[i].indexOf(needle) >= 0) n++;
        if (n) groupHits.push({ t: t, n: n });
      });
      if (!groupHits.length) return;
      any = true;
      panel.appendChild(el("div", { class: "omni-group-title", text: g.label }));
      groupHits.forEach(function (hit) {
        panel.appendChild(el("button", { class: "omni-hit", onclick: function () { panel.hidden = true; revealTable(hit.t.id, q); } }, [
          el("span", { class: "where", text: hit.t.label + (hit.t.contrast ? " · " + hit.t.contrast : "") }),
          el("span", { class: "count", text: hit.n.toLocaleString() + (hit.t.capped && hit.n >= hit.t.rows.length ? "+" : "") + " hit" + (hit.n > 1 ? "s" : "") })
        ]));
      });
    });

    if (!any) panel.appendChild(el("div", { class: "omni-empty", text: "No term matches “" + q + "” in the tables — the gene link above still works." }));
    panel.hidden = false;
  }
  function allFiguresList() {
    return (DATA.figures.assembled || []).map(function (f, i) { return { kind: "assembled", i: i, f: f }; })
      .concat((DATA.figures.panels || []).map(function (f, i) { return { kind: "panel", i: i, f: f }; }));
  }

  // Navigate Explore to a specific table (group + contrast + table) and filter it.
  function revealTable(id, q) {
    activate("explore"); lazyExplore();
    var loc = tableLoc[id]; if (!loc) return;
    var g = groupByKey[loc.groupKey]; if (!g) return;
    if (loc.contrast) wsState.contrast[g.key] = loc.contrast;
    var contrast = currentContrast(g);
    var tables = g.tables.filter(function (t) { return t.contrast == null || t.contrast === contrast; });
    var idx = tables.findIndex(function (t) { return t.id === id; });
    if (idx >= 0) wsState.tableIdx[g.key + "|" + contrast] = idx;
    selectGroup(g.key);
    var t = tableById[id];
    if (t && t._apply) t._apply(q);
  }

  /* =====================================================================
     Provenance
     ===================================================================== */
  function renderProvenance() {
    var v = views.provenance; clear(v);
    var p = DATA.provenance || {};
    var grid = el("div", { class: "prov-grid" });

    var contrastCard = el("div", { class: "card wide" }, [el("h2", { text: "Contrasts" })]);
    contrastCard.appendChild(el("div", { class: "callout", text: meta.contrast_semantics }));
    if ((DATA.contrasts || []).length) {
      var tbl = el("table", { class: "prov-table" }, el("thead", null, el("tr", null, [
        el("th", { text: "Contrast" }), el("th", { text: "Numerator" }), el("th", { text: "Denominator" }), el("th", { text: "Factor" })
      ])));
      var tb = el("tbody");
      DATA.contrasts.forEach(function (c) {
        tb.appendChild(el("tr", null, [
          el("td", { text: c.contrast_id }), el("td", { text: c.numerator }),
          el("td", { text: c.denominator }), el("td", { text: c.factor })
        ]));
      });
      tbl.appendChild(tb); contrastCard.appendChild(el("div", { style: "margin-top:.7rem" }, tbl));
    }
    grid.appendChild(contrastCard);

    var warns = p.warnings || [];
    var wcard = el("div", { class: "card wide" }, el("h2", { text: "Warnings & limitations" }));
    if (!warns.length) wcard.appendChild(el("p", { class: "section-blurb", text: "None recorded." }));
    else {
      var ul = el("ul", { class: "warnlist" });
      warns.forEach(function (w) { ul.appendChild(el("li", null, [el("span", { class: "wicon", text: "▲" }), el("span", { text: String(w) })])); });
      wcard.appendChild(ul);
    }
    grid.appendChild(wcard);

    grid.appendChild(kvCard("Run", [
      ["Profile", p.profile || meta.profile],
      ["Generated (UTC)", p.generated_utc],
      ["Contrasts", (DATA.contrasts || []).length],
      ["Inputs", p.inputs_count],
      ["Result artifacts", p.results_count],
      ["Random seeds", p.random_seeds ? Object.keys(p.random_seeds).map(function (k) { return k + "=" + p.random_seeds[k]; }).join(", ") : null]
    ]));

    if (p.species || p.reference) {
      var s = p.species || {}, r = p.reference || {};
      grid.appendChild(kvCard("Reference", [
        ["Species", s.scientific_name], ["Provider", s.provider], ["Taxonomy", s.taxonomy_id],
        ["Genome build", r.genome_build], ["Annotation", r.annotation_release]
      ]));
    }

    var envRows = [], tools = p.tools || {};
    ["python", "R", "snakemake", "samtools", "featureCounts"].forEach(function (k) { if (tools[k]) envRows.push([k, tools[k]]); });
    if (p.repository) envRows.push(["Repository", (p.repository.revision || "?") + (p.repository.dirty ? " (dirty)" : "")]);
    if (p.platform) envRows.push(["Platform", p.platform]);
    if (p.environment && p.environment.container) envRows.push(["Container", p.environment.container]);
    if (envRows.length) grid.appendChild(kvCard("Environment", envRows));

    if (p.modules && p.modules.length) {
      var mcard = el("div", { class: "card" }, el("h2", { text: "Modules run" }));
      var mods = el("div", { class: "modules" });
      p.modules.forEach(function (m) { mods.appendChild(el("span", { class: "chip", text: m })); });
      mcard.appendChild(mods);
      grid.appendChild(mcard);
    }

    if (p.resource_snapshots && p.resource_snapshots.length) {
      var rc = el("div", { class: "card wide" }, el("h2", { text: "Provider snapshots" }));
      var rt = el("table", { class: "prov-table" });
      var keys = Object.keys(p.resource_snapshots[0] || {}).slice(0, 5);
      rt.appendChild(el("thead", null, el("tr", null, keys.map(function (k) { return el("th", { text: k }); }))));
      var rb = el("tbody");
      p.resource_snapshots.forEach(function (row) {
        rb.appendChild(el("tr", null, keys.map(function (k) { return el("td", { text: row[k] == null ? "—" : String(row[k]) }); })));
      });
      rt.appendChild(rb); rc.appendChild(rt); grid.appendChild(rc);
    }

    append(v, grid);
    v.appendChild(el("p", { class: "about", text:
      "This report is a single self-contained file. Figures and result tables are embedded — it browses and searches fully offline. "
      + "Vector-PDF and full-resolution downloads resolve when the file sits in its results folder. Generated by " + meta.engine.name
      + (meta.engine.version ? " " + meta.engine.version : "") + "." }));
  }
  function kvCard(title, rows) {
    var dl = el("dl", { class: "kv" });
    rows.forEach(function (r) {
      if (r[1] == null || r[1] === "") return;
      dl.appendChild(el("dt", { text: r[0] }));
      dl.appendChild(el("dd", { text: String(r[1]) }));
    });
    return el("div", { class: "card" }, [el("h2", { text: title }), dl]);
  }

  /* ----------------------------------------------------------------- go */
  initTheme();
  buildShell();
})();
