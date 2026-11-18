// Small shared interface pieces: server calls, callouts, chart key, insets, toasts.
import { css, COL } from './stage.js';

export async function api(path, body) {
  const r = await fetch('/api/' + path, body
    ? { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(body) } : undefined);
  if (!r.ok) {
    let msg = r.statusText;
    try { msg = (await r.json()).detail || msg; } catch (e) { /* keep statusText */ }
    throw new Error(msg);
  }
  return r.json();
}

// --- number formatting -------------------------------------------------------
export const sci = (x, d = 2) => {
  if (x === null || x === undefined || !isFinite(x)) return 'n/a';
  if (x === 0) return '0';
  const e = Math.floor(Math.log10(Math.abs(x)));
  if (e >= -2 && e < 4) return (+x.toPrecision(d + 1)).toString();
  return (x / 10 ** e).toFixed(d) + '×10' + sup(e);
};
const SUP = { '-': '⁻', 0: '⁰', 1: '¹', 2: '²', 3: '³', 4: '⁴', 5: '⁵', 6: '⁶', 7: '⁷', 8: '⁸', 9: '⁹' };
export const sup = (n) => String(n).split('').map((c) => SUP[c]).join('');
export const fix = (x, d = 3) => (x === null || x === undefined || !isFinite(x) ? 'n/a' : x.toFixed(d));
export const km = (x) => (Math.abs(x) >= 1e4 ? Math.round(x).toLocaleString('en') : (+x.toPrecision(3)).toString()) + ' km';
export const h = (tag, props = {}, ...kids) => {
  const el = document.createElement(tag);
  for (const [k, v] of Object.entries(props)) {
    if (k === 'class') el.className = v;
    else if (k === 'html') el.innerHTML = v;
    else if (k.startsWith('on')) el.addEventListener(k.slice(2), v);
    else el.setAttribute(k, v);
  }
  for (const k of kids.flat()) if (k !== null && k !== undefined) el.append(k);
  return el;
};

// --- toast -------------------------------------------------------------------
let toastTimer;
export function toast(msg, warn = false, ms = 5200) {
  const t = document.getElementById('toast');
  t.textContent = msg;
  t.className = warn ? 'warn' : '';
  t.hidden = false;
  clearTimeout(toastTimer);
  toastTimer = setTimeout(() => { t.hidden = true; }, ms);
}
export async function busy(msg, fn) {
  toast(msg, false, 60000);
  try { const r = await fn(); document.getElementById('toast').hidden = true; return r; }
  catch (e) { toast(e.message, true, 9000); throw e; }
}

// --- callout: a note pinned beside the selected object -----------------------
const callEl = document.getElementById('callout');
const leaders = document.getElementById('leaders');
let anchor = null;
export const callout = {
  open(title, anchorFn, build) {
    callEl.replaceChildren(h('button', { class: 'close', 'aria-label': 'Close', onclick: () => callout.close() }, '×'),
      h('h2', {}, title));
    build(callEl);
    callEl.hidden = false;
    anchor = anchorFn;
    this.onclose = null;
  },
  close() { callEl.hidden = true; anchor = null; leaders.replaceChildren(); if (this.onclose) this.onclose(); },
  get el() { return callEl; },
  place(stage) {
    if (callEl.hidden || !anchor) return;
    const v = anchor();
    if (!v) return;
    const p = stage.project(v);
    const w = callEl.offsetWidth, hgt = callEl.offsetHeight;
    // prefer the right-hand side of the object; keep clear of the title block and ruler
    let x = p.x + 46;
    if (x + w > innerWidth - 16) x = p.x - w - 46;
    x = Math.max(340, Math.min(x, innerWidth - w - 16));
    let y = Math.max(16, Math.min(p.y - 60, innerHeight - hgt - 80));
    // never sit on top of an open inset chart: slide below it, or failing that to its left
    if (!insetEl.hidden) {
      const r = insetEl.getBoundingClientRect();
      if (x + w > r.left - 8 && y < r.bottom + 8) {
        if (r.bottom + 12 + hgt < innerHeight - 76) y = r.bottom + 12;
        else x = Math.max(340, r.left - w - 14);
      }
    }
    callEl.style.left = x + 'px';
    callEl.style.top = y + 'px';
    const ax = x > p.x ? x : x + w;
    leaders.innerHTML = `<line x1="${p.x}" y1="${p.y}" x2="${ax}" y2="${Math.min(Math.max(p.y, y + 14), y + hgt - 14)}"/>`;
  },
};
export const dl = (rows) => h('dl', {}, rows.filter(Boolean).flatMap(([k, v, cls]) => [h('dt', {}, k), h('dd', cls ? { class: cls } : {}, v)]));
export const act = (label, fn, cls = '') => h('button', { class: 'act ' + cls, onclick: fn }, label);
export const why = (text) => h('p', { class: 'why' }, text);
export function slider(label, min, max, step, value, fmt, oninput, onchange) {
  const out = h('output', {}, fmt(value));
  const inp = h('input', { type: 'range', min, max, step, value });
  inp.addEventListener('input', () => { out.textContent = fmt(+inp.value); if (oninput) oninput(+inp.value); });
  if (onchange) inp.addEventListener('change', () => onchange(+inp.value));
  const el = h('label', { class: 'field' }, h('span', {}, label), out, inp);
  el.input = inp; el.out = out;
  return el;
}

// --- chart key ---------------------------------------------------------------
// items: [{group}] or [{id, label, color, style: 'line'|'dash'|'fill'|'dot'|'wide', on, toggle(on)}]
export function buildLegend(items) {
  const box = document.getElementById('legend');
  box.replaceChildren();
  for (const it of items) {
    if (it.group) { box.append(h('div', { class: 'grp' }, it.group)); continue; }
    const c = css(it.color);
    const sw = {
      line: `<svg class="sw"><line x1="0" y1="6" x2="34" y2="6" stroke="${c}" stroke-width="1.3"/></svg>`,
      wide: `<svg class="sw"><line x1="0" y1="6" x2="34" y2="6" stroke="${c}" stroke-width="3"/></svg>`,
      dash: `<svg class="sw"><line x1="0" y1="6" x2="34" y2="6" stroke="${c}" stroke-width="1.6" stroke-dasharray="5 3"/></svg>`,
      fill: `<svg class="sw"><rect x="1" y="1" width="32" height="10" fill="${c}" stroke="${css(COL.landEdge)}"/></svg>`,
      dot: `<svg class="sw"><circle cx="8" cy="6" r="2.2" fill="${c}"/><circle cx="18" cy="4" r="2.2" fill="${c}"/><circle cx="27" cy="8" r="2.2" fill="${c}"/></svg>`,
      ring: `<svg class="sw"><circle cx="17" cy="6" r="5" fill="none" stroke="${c}" stroke-width="1.4"/></svg>`,
    }[it.style || 'line'];
    const b = h('button', { 'aria-pressed': String(!!it.on), title: it.hint || '', html: sw + `<span>${it.label}</span>` });
    b.addEventListener('click', async () => {
      const on = b.getAttribute('aria-pressed') !== 'true';
      b.setAttribute('aria-pressed', String(on));
      try { await it.toggle(on); } catch (e) { b.setAttribute('aria-pressed', String(!on)); }
    });
    it.button = b;
    box.append(b);
  }
}

// --- inset chart (2D canvas) --------------------------------------------------
const insetEl = document.getElementById('inset');
export const inset = {
  open(title, height = 300) {
    const cv = h('canvas');
    const cap = h('p');
    const extra = h('div');
    insetEl.replaceChildren(h('button', { class: 'close', 'aria-label': 'Close', onclick: () => inset.close() }, '×'),
      h('h3', {}, title), cv, extra, cap);
    insetEl.hidden = false;
    const w = insetEl.clientWidth - 24;
    cv.width = w * 2; cv.height = height * 2; cv.style.height = height + 'px';
    const g = cv.getContext('2d');
    g.scale(2, 2);
    this.onclose = null;
    return { g, w, hgt: height, cap, extra, cv };
  },
  close() { insetEl.hidden = true; if (this.onclose) this.onclose(); this.onclose = null; },
  get hidden() { return insetEl.hidden; },
};

// Simple axes helper for inset plots: returns functions mapping data -> pixels.
export function axes(g, w, hgt, xr, yr, opt = {}) {
  const m = { l: 44, r: 10, t: 10, b: 30 };
  const lx = opt.logx, ly = opt.logy;
  const tx = (v) => (lx ? Math.log10(v) : v), ty = (v) => (ly ? Math.log10(v) : v);
  const X = (v) => m.l + (tx(v) - tx(xr[0])) / (tx(xr[1]) - tx(xr[0])) * (w - m.l - m.r);
  const Y = (v) => hgt - m.b - (ty(v) - ty(yr[0])) / (ty(yr[1]) - ty(yr[0])) * (hgt - m.t - m.b);
  g.clearRect(0, 0, w, hgt);
  g.strokeStyle = css(COL.rule); g.fillStyle = css(COL.soft); g.lineWidth = 0.6;
  g.font = '11px "Alegreya Sans", sans-serif';
  const ticks = (r, log) => {
    if (log) { const a = []; for (let e = Math.ceil(Math.log10(r[0])); e <= Math.floor(Math.log10(r[1])); e += Math.max(1, Math.round((Math.log10(r[1]) - Math.log10(r[0])) / 5))) a.push(10 ** e); return a; }
    const span = r[1] - r[0], step = 10 ** Math.floor(Math.log10(span / 4));
    const s = [1, 2, 5, 10].map((k) => k * step).find((k) => span / k <= 6);
    const a = []; for (let v = Math.ceil(r[0] / s) * s; v <= r[1] + 1e-12; v += s) a.push(v); return a;
  };
  const tidy = (v, r) => (Math.abs(v) < 1e-9 * (r[1] - r[0]) ? '0' : String(+v.toPrecision(3)));   // no 6.9e-18 for zero
  g.textAlign = 'center';
  for (const v of ticks(xr, lx)) { g.beginPath(); g.moveTo(X(v), m.t); g.lineTo(X(v), hgt - m.b); g.stroke(); g.fillText(lx ? '10' + sup(Math.round(Math.log10(v))) : tidy(v, xr), X(v), hgt - m.b + 12); }
  g.textAlign = 'right';
  for (const v of ticks(yr, ly)) { g.beginPath(); g.moveTo(m.l, Y(v)); g.lineTo(w - m.r, Y(v)); g.stroke(); g.fillText(ly ? '10' + sup(Math.round(Math.log10(v))) : tidy(v, yr), m.l - 4, Y(v) + 3); }
  g.textAlign = 'center';
  if (opt.xlabel) g.fillText(opt.xlabel, (m.l + w - m.r) / 2, hgt - 4);
  if (opt.ylabel) { g.save(); g.translate(10, (m.t + hgt - m.b) / 2); g.rotate(-Math.PI / 2); g.fillText(opt.ylabel, 0, 0); g.restore(); }
  g.strokeStyle = css(COL.ink); g.lineWidth = 1; g.strokeRect(m.l, m.t, w - m.l - m.r, hgt - m.t - m.b);
  return { X, Y, m };
}
export function plotLine(g, ax, xs, ys, color, width = 1.5, dash = null) {
  g.save(); g.strokeStyle = color; g.lineWidth = width; if (dash) g.setLineDash(dash);
  g.beginPath();
  let pen = false;
  xs.forEach((x, i) => { const y = ys[i]; if (y === null || !isFinite(y) || !isFinite(x)) { pen = false; return; } if (pen) g.lineTo(ax.X(x), ax.Y(y)); else { g.moveTo(ax.X(x), ax.Y(y)); pen = true; } });
  g.stroke(); g.restore();
}

export function download(name, text, type = 'text/csv') {
  const a = document.createElement('a');
  a.href = text.startsWith('data:') ? text : URL.createObjectURL(new Blob([text], { type }));
  a.download = name;
  a.click();
}
