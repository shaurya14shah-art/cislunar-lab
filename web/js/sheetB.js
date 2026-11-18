// Sheet B - the N-body sandbox with the relativity and field overlays.
// SI units throughout; the scene is scaled by one length L0 (metres per scene unit).
import { COL, css, V, THREE } from './stage.js';
import { api, busy, toast, callout, dl, act, why, slider, buildLegend, inset, axes, plotLine, sci, fix, h } from './ui.js';
import { openSchwarzschild } from './closeup.js';

const BAND_COL = { newtonian: COL.rule, '1pn-good': COL.soft, '1pn-ok': 0xd48ab6, '1pn-marginal': COL.magenta, breakdown: COL.magenta };
const TINT = [0x2f5d7c, 0x8a969d, 0xb5803a, 0x5c7d4d, 0x7a5c8f, 0x9c4a3f, 0x3f7f86];
const human = (t) => { const a = Math.abs(t); if (a < 120) return `${t.toPrecision(4)} s`; if (a < 7200) return `${(t / 60).toPrecision(4)} min`; if (a < 172800) return `${(t / 3600).toPrecision(4)} h`; if (a < 3.2e7 * 2) return `${(t / 86400).toPrecision(4)} days`; return `${(t / 3.15576e7).toPrecision(4)} years`; };
const unitOf = (T) => (T < 120 ? [1, 'seconds'] : T < 7200 ? [60, 'minutes'] : T < 172800 ? [3600, 'hours'] : T < 6.4e7 ? [86400, 'days'] : [3.15576e7, 'years']);
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const norm = (a) => Math.hypot(a[0], a[1], a[2]);

export class SheetB {
  constructor(stage, meta) {
    this.stage = stage; this.meta = meta;
    this.title = 'An isolated system in SI units. Add, remove or edit bodies; magenta marks everything relativistic.';
    this.on = { trails: true, sheet: false, lines: false, tides: false, hill: false, lagr: false, rings: true, clocks: true, rays: false };
    this.method = 'dop853'; this.pn = false; this.frameKind = 'inertial';
    this.sel = -1; this.frames = []; this.cursor = 0; this.drift = 0;
    this.g = {};
  }

  // ---- scenarios ----------------------------------------------------------------
  preset(name, extra = {}) {
    const p = this.meta.presets[name];
    return { name: p.name, gm: p.gm, radius: p.radius, pos: [0, 0, 0], vel: [0, 0, 0], spin: [...p.spin], ...extra };
  }

  // put `b` on an orbit around bodies[k]: circular if f = 1, eccentric otherwise (f scales the speed)
  orbit(bodies, b, k, d, f = 1, angle = 0) {
    const P = bodies[k], mt = P.gm + b.gm, v = f * Math.sqrt(mt / d), q = b.gm / mt;
    const e = [Math.cos(angle), Math.sin(angle), 0], t = [-Math.sin(angle), Math.cos(angle), 0];
    b.pos = P.pos.map((x, i) => x + (1 - q) * d * e[i]); b.vel = P.vel.map((x, i) => x + (1 - q) * v * t[i]);
    P.pos = P.pos.map((x, i) => x - q * d * e[i]); P.vel = P.vel.map((x, i) => x - q * v * t[i]);
    bodies.push(b);
    return bodies;
  }

  scenario(key) {
    const B = [];
    if (key === 'earth-moon') { B.push(this.preset('Earth')); this.orbit(B, this.preset('Moon'), 0, 3.844e8); }
    if (key === 'inner') {
      B.push(this.preset('Sun'));
      this.orbit(B, this.preset('Mercury'), 0, 4.6e10, Math.sqrt(1.2056), 0.4);   // starts at perihelion, e = 0.2056
      this.orbit(B, this.preset('Earth'), 0, 1.496e11, 1, 2.2);
      this.orbit(B, this.preset('Jupiter'), 0, 7.78e11, 1, 4.1);
    }
    if (key === 'ns-binary') { B.push(this.preset('Neutron star', { name: 'Neutron star A' })); this.orbit(B, this.preset('Neutron star', { name: 'Neutron star B' }), 0, 2.0e6, 0.85); }
    if (key === 'bh-wd') { B.push(this.preset('Black hole (10 Msun)', { name: 'Black hole' })); this.orbit(B, this.preset('White dwarf'), 0, 9.0e8, 0.8); }
    if (key === 'empty') { /* nothing: build your own */ }
    return B;
  }

  async enter() {
    const s = this.stage;
    s.clear();
    for (const k of ['sheet', 'lines', 'bodies', 'trails', 'tides', 'hill', 'lagr', 'rays', 'ghost']) this.g[k] = s.group(k);
    this.controls();
    this.legend();
    this.ws = new WebSocket(`ws://${location.host}/ws/nbody`);
    this.ws.onmessage = (e) => this.onMessage(JSON.parse(e.data));
    await new Promise((res) => { this.ws.onopen = res; });
    await this.load(this.scenario('earth-moon'));
  }

  leave() { callout.close(); inset.close(); this.ws?.close(); for (const d of this.dials || []) d.el.remove(); this.dials = []; }

  controls() {
    const box = document.getElementById('sheet-controls');
    const sc = h('select', {},
      h('option', { value: 'earth-moon' }, 'Earth and Moon'),
      h('option', { value: 'inner' }, 'Sun, Mercury, Earth, Jupiter'),
      h('option', { value: 'ns-binary' }, 'Two neutron stars, 2,000 km apart'),
      h('option', { value: 'bh-wd' }, 'White dwarf skimming a black hole'),
      h('option', { value: 'empty' }, 'Empty space'));
    sc.addEventListener('change', () => this.load(this.scenario(sc.value)));
    const add = h('select', {}, h('option', { value: '' }, 'Add a body…'), Object.keys(this.meta.presets).map((n) => h('option', { value: n }, n)), h('option', { value: '*' }, 'Custom body'));
    add.addEventListener('change', () => { if (add.value) this.addBody(add.value); add.value = ''; });
    const meth = h('select', {}, h('option', { value: 'dop853' }, 'High-accuracy stepping (DOP853)'), h('option', { value: 'yoshida4' }, 'Long-run stepping, no energy drift (Yoshida)'));
    meth.addEventListener('change', () => { this.method = meth.value; if (this.pn && this.method === 'yoshida4') toast('The relativistic force depends on velocity, which a symplectic splitting cannot handle; 1PN runs use the adaptive integrator.'); this.restart(); });
    const pn = h('input', { type: 'checkbox' });
    pn.addEventListener('change', () => { this.pn = pn.checked; this.restart(); });
    const fr = h('select', {}, h('option', { value: 'inertial' }, 'View: standing still'), h('option', { value: 'bary' }, 'View: centred on the balance point of all masses'), h('option', { value: 'corot' }, 'View: turning with the selected pair'));
    fr.addEventListener('change', () => { this.frameKind = fr.value; this.resetTrails(); });
    box.replaceChildren(sc, add, meth, fr, h('label', { class: 'rel-note' }, pn, ' Add Einstein’s correction to gravity (1PN)'));
  }

  legend() {
    const T = (key, after) => async (on) => { this.on[key] = on; if (this.g[key]) this.g[key].visible = on; if (after) await after(on); };
    buildLegend([
      { group: 'Motion' },
      { label: 'Trails', color: COL.ink, style: 'line', on: true, toggle: T('trails') },
      { group: 'Fields' },
      { label: 'Potential as a warped sheet', color: COL.soft, style: 'line', on: false, toggle: T('sheet', () => this.refreshField()) },
      { label: 'Field lines', color: COL.ink, style: 'line', on: false, toggle: T('lines', () => this.refreshField()) },
      { label: 'Tides: stretch and squeeze on the selected body', color: COL.unstable, style: 'wide', on: false, toggle: T('tides', () => this.refreshField()) },
      { label: 'Hill sphere and Roche limit', color: COL.stable, style: 'ring', on: false, toggle: T('hill', () => this.refreshPair()) },
      { label: 'Lagrange points of the selected pair', color: COL.unstable, style: 'dot', on: false, toggle: T('lagr', () => this.refreshPair()) },
      { group: 'Relativity' },
      { label: 'Ring colour: how relativistic each body is', color: COL.magenta, style: 'ring', on: true, toggle: T('rings') },
      { label: 'Clocks ticking at their own rates', color: COL.magenta, style: 'ring', on: true, toggle: T('clocks') },
      { label: 'Light rays past the selection', color: COL.magenta, style: 'line', on: false, toggle: T('rays', () => this.refreshRays()) },
      { label: 'Newton against Einstein, side by side', color: COL.magenta, style: 'dash', on: false, toggle: async (on) => { this.g.ghost.visible = on; if (on) await this.compare(); else inset.close(); } },
    ]);
  }

  // ---- system state -------------------------------------------------------------
  async load(bodies) {
    callout.close(); inset.close();
    this.bodies = bodies.map((b, i) => ({ ...b, color: TINT[i % TINT.length] }));
    this.sel = -1; this.t0 = 0; this.collided = new Set();
    this.frames = [{ t: 0, pos: this.bodies.map((b) => [...b.pos]), vel: this.bodies.map((b) => [...b.vel]) }];
    this.cursor = 0; this.drift = 0; this.pending = false;
    const ext = Math.max(1, ...this.bodies.map((b) => norm(b.pos)));
    this.L0 = this.bodies.length ? ext * 1.15 : 1e9;
    if (this.bodies.length === 1) this.L0 = this.bodies[0].radius * 30;
    this.stage.frameOn(new THREE.Vector3(0, 0, 0), 1.05);
    this.buildBodies();
    this.send({ cmd: 'set', bodies: this.wire(), t: 0 });
    await this.analyze();
    this.refreshField(); this.refreshPair(); this.refreshRays();
    if (!this.bodies.length) toast('Empty space. Use “Add a body…” in the title block to place the first one.');
  }

  wire() { return this.bodies.map(({ name, gm, radius, pos, vel, spin }) => ({ name, gm, radius, pos, vel, spin })); }
  send(m) { if (this.ws.readyState === 1) this.ws.send(JSON.stringify(m)); }

  // Take the frame under the playhead as the new present, drop the future, tell the server.
  commit() {
    const f = this.frameAt();
    this.bodies.forEach((b, i) => { if (f.pos[i]) { b.pos = [...f.pos[i]]; b.vel = [...f.vel[i]]; } });
    return f.t;
  }

  async restart(tNow = null) {
    const t = tNow ?? this.commit();
    this.frames = [{ t, pos: this.bodies.map((b) => [...b.pos]), vel: this.bodies.map((b) => [...b.vel]) }];
    this.cursor = 0; this.pending = false; this.drift = 0;
    this.buildBodies();
    this.send({ cmd: 'set', bodies: this.wire(), t });
    await this.analyze();
    this.refreshField(); this.refreshPair(); this.refreshRays();
  }

  buildBodies() {
    const s = this.stage;
    this.g.bodies.clear(); this.g.trails.clear();
    for (const l of [...s.labels]) s.dropLabel(l);
    for (const d of this.dials || []) d.el.remove();
    s.pickables = []; this.dials = []; this.meshes = []; this.rings = []; this.trailBuf = [];
    this.bodies.forEach((b, i) => {
      const m = s.sphere(1, b.color, this.g.bodies);
      m.userData.index = i;
      s.pickables.push(m); this.meshes.push(m);
      const ring = new THREE.Sprite(new THREE.SpriteMaterial({ map: ringTexture(), color: COL.rule, depthTest: false, sizeAttenuation: false }));
      ring.scale.set(0.06, 0.06, 1);
      this.g.bodies.add(ring); this.rings.push(ring);
      s.label(b.name, () => m.position, 'body offset', () => this.select(i));
      const el = document.createElement('div');
      el.className = 'dial';
      el.innerHTML = '<svg viewBox="-16 -16 32 32"><circle r="13" fill="#f8fbfb" stroke="#b0127a" stroke-width="1.4"/><line class="ref" x1="0" y1="0" x2="0" y2="-11" stroke="#a9c1c9" stroke-width="1.6"/><line class="hand" x1="0" y1="0" x2="0" y2="-11" stroke="#b0127a" stroke-width="1.8"/></svg><span></span>';
      s.labelLayer.appendChild(el);
      this.dials.push({ el, hand: el.querySelector('.hand'), ref: el.querySelector('.ref'), txt: el.querySelector('span'), lag: 0 });
      const geo = new THREE.BufferGeometry();
      geo.setAttribute('position', new THREE.Float32BufferAttribute(new Float32Array(900 * 3), 3));
      geo.setDrawRange(0, 0);
      const line = new THREE.Line(geo, new THREE.LineBasicMaterial({ color: b.color }));
      line.frustumCulled = false;
      this.g.trails.add(line);
      this.trailBuf.push({ line, n: 0, last: -1 });
    });
    s.onPick = (obj) => { if (obj && obj.userData.index !== undefined) this.select(obj.userData.index); };
  }

  resetTrails() { for (const t of this.trailBuf || []) { t.n = 0; t.last = -1; t.line.geometry.setDrawRange(0, 0); } }

  async analyze() {
    if (!this.bodies.length) { this.an = null; return; }
    const f = this.frameAt();
    const bodies = this.wire().map((b, i) => ({ ...b, pos: f.pos[i] || b.pos, vel: f.vel[i] || b.vel }));
    try { this.an = await api('nbody/analyze', { bodies, pair: this.pair() }); } catch (e) { /* keep the previous analysis */ }
    if (this.an) this.an.regime.forEach((r, i) => { if (this.rings[i]) this.rings[i].material.color.setHex(BAND_COL[r.band]); });
  }

  // the "pair" used for the co-rotating frame, Hill/Roche and Lagrange points:
  // the selected body and whichever other body pulls on it hardest
  pair() {
    const n = this.bodies.length;
    if (n < 2) return null;
    const i = this.sel >= 0 ? this.sel : 1;
    const f = this.frameAt();
    let best = -1, bv = -1;
    for (let j = 0; j < n; j++) if (j !== i) { const d = norm(sub(f.pos[i], f.pos[j])); const a = this.bodies[j].gm / (d * d); if (a > bv) { bv = a; best = j; } }
    return [best, i];
  }

  // ---- live integration ---------------------------------------------------------
  tChar() {
    const f = this.frames[this.frames.length - 1];
    let T = Infinity;
    for (let i = 0; i < this.bodies.length; i++) for (let j = i + 1; j < this.bodies.length; j++) {
      const d = norm(sub(f.pos[i], f.pos[j])), mt = this.bodies[i].gm + this.bodies[j].gm;
      if (mt > 0) T = Math.min(T, 2 * Math.PI * Math.sqrt(d ** 3 / mt));
    }
    return isFinite(T) ? T : 1e5;
  }

  requestMore() {
    if (this.pending || this.bodies.length === 0) return;
    this.pending = true;
    const T = this.tChar();
    this.send({ cmd: 'advance', t: 0.5 * T, n_frames: 80, method: this.method, dt: T / 500, pn: this.pn });
  }

  onMessage(m) {
    if (!m.ok) { this.pending = false; toast('Integration stopped: ' + m.error, true); this.halt = true; return; }
    if (m.set !== undefined) return;
    this.pending = false;
    for (let k = 1; k < m.t.length; k++) this.frames.push({ t: m.t[k], pos: m.pos[k], vel: m.vel[k] });
    if (this.frames.length > 5000) { const cut = this.frames.length - 5000; this.frames.splice(0, cut); this.cursor = Math.max(0, this.cursor - cut); this.resetTrails(); }
    this.drift = m.energy_drift; this.driftMethod = m.method;
    for (const c of m.collisions) {
      const key = c.names.join('+');
      if (!this.collided.has(key)) { this.collided.add(key); toast(`Collision flagged: ${c.names[0]} and ${c.names[1]} overlap at t = ${sci(c.t + m.t0)} s. They are not merged; what happens next is not physical.`, true, 12000); }
    }
  }

  frameAt() { return this.frames[Math.min(this.frames.length - 1, Math.max(0, Math.round(this.cursor)))]; }

  // express a frame's positions in the chosen display frame
  display(f) {
    const n = this.bodies.length;
    if (this.frameKind === 'inertial' || n === 0) return f.pos;
    const gm = this.bodies.map((b) => b.gm);
    const com = (idx) => { let m = 0; const c = [0, 0, 0]; for (const i of idx) { m += gm[i]; for (let k = 0; k < 3; k++) c[k] += gm[i] * f.pos[i][k]; } return m ? c.map((x) => x / m) : c; };
    if (this.frameKind === 'bary') { const c = com([...Array(n).keys()]); return f.pos.map((p) => sub(p, c)); }
    const pr = this.pairFixed || this.pair();
    if (!pr) return f.pos;
    const [i, j] = pr, o = com([i, j]);
    const d = sub(f.pos[j], f.pos[i]), dv = sub(f.vel[j], f.vel[i]);
    const ex = d.map((x) => x / norm(d)); let hz = cross(d, dv); if (norm(hz) === 0) hz = [0, 0, 1];
    const ez = hz.map((x) => x / norm(hz)), ey = cross(ez, ex);
    return f.pos.map((p) => { const r = sub(p, o); return [dot(r, ex), dot(r, ey), dot(r, ez)]; });
  }

  tick(dt, playing, speed) {
    callout.place(this.stage);
    const n = this.bodies ? this.bodies.length : 0;
    if (!n || !this.meshes) return { clock: 'No bodies yet.', drift: '', frac: 0 };
    if (playing && !this.halt) {
      this.cursor = Math.max(0, Math.min(this.frames.length - 1, this.cursor + dt * speed * 24));
      if (speed > 0 && this.cursor > this.frames.length - 50) this.requestMore();
    }
    const f = this.frameAt(), P = this.display(f), cam = this.stage.camera;
    const idx = Math.round(this.cursor);
    P.forEach((p, i) => {
      const m = this.meshes[i]; if (!m) return;
      m.position.set(p[0] / this.L0, p[1] / this.L0, p[2] / this.L0);
      const minR = 0.0065 * m.position.distanceTo(cam.position);          // never smaller than a few pixels
      const r = Math.max(this.bodies[i].radius / this.L0, minR);
      m.scale.setScalar(r); m.userData.toScale = this.bodies[i].radius / this.L0 >= minR;
      this.rings[i].position.copy(m.position); this.rings[i].visible = this.on.rings;
      const tb = this.trailBuf[i];
      if (idx !== tb.last) {                    // rebuild the trail from the buffered frames behind the playhead
        const a = tb.line.geometry.attributes.position.array, k0 = Math.max(0, idx - 899);
        let q = 0;
        for (let k = k0; k <= idx; k++) { const pp = this.display(this.frames[k])[i]; a[q++] = pp[0] / this.L0; a[q++] = pp[1] / this.L0; a[q++] = pp[2] / this.L0; }
        tb.line.geometry.setDrawRange(0, idx - k0 + 1); tb.line.geometry.attributes.position.needsUpdate = true; tb.last = idx;
      }
    });
    this.updateClocks(f.t);
    this._acc = (this._acc || 0) + dt;
    if (this._acc > 0.8 && playing) { this._acc = 0; this.analyze().then(() => { if (this.on.sheet || this.on.lines || this.on.tides) this.refreshField(); if (this.on.hill || this.on.lagr) this.refreshPair(); if (!callout.el.hidden && this.sel >= 0 && !callout.el.contains(document.activeElement)) this.fillRegime(); }); }
    const T = this.tChar();
    return {
      clock: `t = ${human(f.t)}`,
      // With 1PN on, the conserved energy of the truncated theory itself wanders by ~eps^2:
      // that is physics left out, not numerical error, so it is labelled and judged differently.
      drift: this.pn ? `1PN energy holds to ${sci(Math.abs(this.drift), 1)} (theory limit ≈ ε² = ${sci(this.worst() ** 2, 0)})`
        : `Energy drift ${sci(this.drift, 1)}, ${this.method === 'yoshida4' ? 'symplectic' : 'adaptive'}`,
      bad: Math.abs(this.drift) > (this.pn ? 30 * this.worst() ** 2 + 1e-9 : 1e-6), frac: this.frames.length > 1 ? this.cursor / (this.frames.length - 1) : 0,
    };
  }

  worst() { return this.an ? Math.max(...this.an.regime.map((q) => q.pn1_fractional_size)) : 0; }

  seek(frac) { this.cursor = frac * (this.frames.length - 1); }

  // Each dial has a grey hand (a far-away reference clock) and a magenta hand (this body's
  // own clock).  Real rate differences are parts per billion, so the *difference* is
  // magnified by a stated factor until the largest one is visible.
  updateClocks(t) {
    if (!this.an || this.an.regime.length !== this.bodies.length) return;
    const rates = this.an.regime.map((r) => r.clock_rate_offset);
    const big = Math.max(...rates.map(Math.abs), 1e-30);
    const mag = big > 0.02 ? 1 : 10 ** Math.floor(Math.log10(0.2 / big));
    const Td = this.tChar() / 3;
    this.dials.forEach((d, i) => {
      const m = this.meshes[i]; if (!m || !this.an.regime[i]) return;
      const p = this.stage.project(m.position);
      d.el.style.display = this.on.clocks && p.front ? '' : 'none';
      d.el.style.left = p.x + 'px'; d.el.style.top = p.y + 'px';
      const ref = (t / Td) * 360, own = ref * (1 + rates[i] * mag);
      d.ref.setAttribute('transform', `rotate(${ref % 360})`);
      d.hand.setAttribute('transform', `rotate(${own % 360})`);
      const us = this.an.regime[i].clock_us_per_day;
      d.txt.innerHTML = `${rates[i] <= -0.999 ? 'stopped' : Math.abs(us) < 1e7 ? sci(us, 1) + ' µs/day' : (rates[i] * 100).toFixed(1) + '%'}${mag > 1 ? `<br>lag shown ×${sci(mag, 0)}` : ''}`;
    });
  }

  // ---- selection and editing ------------------------------------------------------
  select(i) {
    this.sel = i;
    this.stage.labels.forEach((l) => l.el.classList.remove('sel'));
    const b = this.bodies[i];
    this.pairFixed = null;
    callout.open(b.name, () => this.meshes[i]?.position, (el) => {
      this.regBox = h('div');
      el.append(this.regBox);
      this.fillRegime();
      const f = this.frameAt();
      const inp = (v) => h('input', { type: 'text', value: v });
      const e = { name: inp(b.name), mass: inp((b.gm / this.meta.G).toExponential(5)), radius: inp((b.radius / 1e3).toPrecision(7)),
        pos: f.pos[i].map((x) => inp((x / 1e3).toPrecision(8))), vel: f.vel[i].map((x) => inp((x / 1e3).toPrecision(8))),
        spin: inp(norm(b.spin) > 0 ? (2 * Math.PI / norm(b.spin) / 3600).toPrecision(5) : '0') };
      el.append(h('div', { class: 'edit' },
        h('span', {}, 'Name'), e.name, h('span', {}, 'Mass, kg'), e.mass, h('span', {}, 'Radius, km'), e.radius,
        h('span', {}, 'Position, km'), h('div', { class: 'tri' }, e.pos), h('span', {}, 'Velocity, km/s'), h('div', { class: 'tri' }, e.vel),
        h('span', {}, 'Spin period, h'), e.spin));
      const acts = h('div', { class: 'acts' });
      acts.append(act('Apply these values', async () => {
        const num = (x) => { const v = parseFloat(x.value); if (!isFinite(v)) throw new Error('Every field needs a number.'); return v; };
        try {
          const t = this.commit();
          b.name = e.name.value || b.name; b.gm = num(e.mass) * this.meta.G; b.radius = num(e.radius) * 1e3;
          b.pos = e.pos.map((x) => num(x) * 1e3); b.vel = e.vel.map((x) => num(x) * 1e3);
          const sp = num(e.spin); b.spin = [0, 0, sp > 0 ? 2 * Math.PI / (sp * 3600) : 0];
          await this.restart(t); this.select(i);
        } catch (err) { toast(err.message, true); }
      }));
      acts.append(act('Close-up: exact strong-field orbits and light around this body', () => openSchwarzschild(b, this.meta), 'rel-note'));
      acts.append(act('Remove this body', async () => { const t = this.commit(); this.bodies.splice(i, 1); callout.close(); this.sel = -1; await this.restart(t); }, 'warn'));
      el.append(acts);
      el.append(why('Spin is stored and shown, but it does not enter the point-mass equations of motion: this code has no frame dragging or oblateness. The break-up fraction above is only a sanity check on the value you typed.'));
    });
    this.refreshField(); this.refreshPair(); this.refreshRays();
  }

  fillRegime() {
    const i = this.sel, r = this.an?.regime[i], b = this.bodies[i];
    if (!r || !this.regBox) return;
    const x = Math.max(r.pn1_fractional_size, 1e-12), pos = Math.min(100, (Math.log10(x) + 12) / 12 * 100);
    const u = this.an.units;
    const wasOpen = this.regBox.querySelector('details')?.open;
    this.regBox.replaceChildren(
      h('p', { class: r.band === 'newtonian' || r.band === '1pn-good' ? '' : 'rel-note' }, r.message),
      r.surface_warning ? h('p', { class: 'rel-note' }, 'Near this body’s own surface the field is strong (GM/Rc² above 0.01): 1PN is not valid there. Use the close-up for exact orbits and light.') : '',
      h('div', { class: 'meter', title: 'Size of the relativistic correction, from 10⁻¹² (left) to 1 (right)' }, h('i', { style: `left:${pos}%` }),
        h('span', { style: 'left:0%' }, '10⁻¹²'), h('span', { style: 'left:33.3%' }, '10⁻⁸'), h('span', { style: 'left:75%' }, '10⁻³'), h('span', { style: 'left:97%' }, '1')),
      dl([
        ['Depth in others’ gravity (GM/rc²)', sci(r.eps_orbit, 2)],
        ['Speed, (v/c)²', sci(r.beta2, 2)],
        ['Surface clock against a far-away clock', r.clock_rate_offset <= -0.999 ? 'stopped (event horizon)' : Math.abs(r.clock_us_per_day) < 1e7 ? sci(r.clock_us_per_day, 2) + ' µs/day' : (r.clock_rate_offset * 100).toFixed(2) + '%'],
        ['Radius ÷ Schwarzschild radius', sci(1 / r.rs_over_R, 2), r.surface_warning ? 'rel-note' : ''],
      ]),
      h('details', {}, h('summary', {}, 'More readouts and units'), dl([
        ['What 1PN itself neglects', '≈ ' + sci(r.pn1_fractional_error, 1)],
        ['Lorentz factor γ − 1', sci(r.gamma - 1 || r.beta2 / 2, 2)],
        ['Kinetic energy, relativistic ÷ Newtonian', (r.ke_ratio).toFixed(Math.abs(r.ke_ratio - 1) > 1e-4 ? 5 : 12)],
        ['Own surface, GM/Rc²', sci(r.compactness, 2)],
        ['Schwarzschild radius', r.rs < 1 ? sci(r.rs * 1000, 2) + ' mm' : sci(r.rs / 1000, 2) + ' km'],
        r.breakup_fraction !== null && r.breakup_fraction > 0 ? ['Spin ÷ break-up spin', sci(r.breakup_fraction, 1)] : null,
        this.meshes[i]?.userData.toScale ? null : ['Drawn size', 'enlarged to stay visible'],
      ]), u ? h('p', { class: 'hint' }, 'Pair units: ' + u.explain.join('; ') + '.') : ''));
    if (wasOpen) this.regBox.querySelector('details').open = true;
  }

  async addBody(name) {
    const t = this.bodies.length ? this.commit() : 0;
    const b = name === '*' ? { name: 'Custom body', gm: 1e22 * this.meta.G, radius: 5e5, pos: [0, 0, 0], vel: [0, 0, 0], spin: [0, 0, 0] } : this.preset(name);
    b.color = TINT[this.bodies.length % TINT.length];
    if (this.bodies.length) {
      // circular orbit around the heaviest body, a little outside everything already there
      const k = this.bodies.reduce((m, q, i) => (q.gm > this.bodies[m].gm ? i : m), 0);
      const ext = Math.max(...this.bodies.map((q) => norm(sub(q.pos, this.bodies[k].pos))), this.bodies[k].radius * 4);
      this.orbit(this.bodies, b, k, ext * 1.6, 1, Math.random() * 6.28);
      this.L0 = Math.max(this.L0, ext * 1.9);
    } else { this.bodies.push(b); this.L0 = b.radius * 30; }
    await this.restart(t);
    this.select(this.bodies.length - 1);
    toast(`${b.name} added on a circular orbit. Edit its position and velocity in the note beside it.`);
  }

  // ---- field overlays ---------------------------------------------------------------
  current() { const f = this.frameAt(); return this.wire().map((b, i) => ({ ...b, pos: f.pos[i], vel: f.vel[i] })); }

  async refreshField() {
    const s = this.stage;
    if (!this.bodies.length || this._fbusy) return;
    if (!(this.on.sheet || this.on.lines || this.on.tides)) { this.g.sheet.clear(); this.g.lines.clear(); this.g.tides.clear(); return; }
    this._fbusy = true;
    try {
      const bodies = this.current(), n = 64, hw = this.L0 * 1.25;
      const req = { bodies, center: [0, 0, 0], half_width: hw, n, n_lines: this.on.lines ? 26 : 0 };
      if (this.on.tides && this.sel >= 0) {            // tide felt AT the selected body, from all the others
        const others = bodies.filter((_, i) => i !== this.sel);
        if (others.length) { const r2 = await api('nbody/field', { bodies: others, center: [0, 0, 0], half_width: hw, n: 4, probe: bodies[this.sel].pos }); this.drawTides(r2.tidal, bodies[this.sel]); } else this.g.tides.clear();
      } else { this.g.tides.clear(); if (this.tideNote) { this.stage.dropLabel(this.tideNote); this.tideNote = null; } }
      const r = await api('nbody/field', req);
      this.g.sheet.clear(); this.g.lines.clear();
      if (this.on.sheet) {
        const geo = new THREE.PlaneGeometry(2 * hw / this.L0, 2 * hw / this.L0, n - 1, n - 1), a = geo.attributes.position;
        const span = Math.abs(r.phi_max - r.floor) || 1;
        for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) a.setZ((n - 1 - j) * n + i, -0.55 * (r.phi_max - r.phi[j * n + i]) / span - 0.35);
        geo.computeVertexNormals();
        this.g.sheet.add(new THREE.Mesh(geo, new THREE.MeshBasicMaterial({ color: COL.paper, side: THREE.DoubleSide, transparent: true, opacity: 0.55, polygonOffset: true, polygonOffsetFactor: 1 })));
        this.g.sheet.add(new THREE.LineSegments(new THREE.WireframeGeometry(geo), new THREE.LineBasicMaterial({ color: COL.soft, transparent: true, opacity: 0.35 })));
        if (!this.sheetNote) this.sheetNote = s.label('Potential sheet: depth ∝ −Φ, wells cut off at the deepest 1%', new THREE.Vector3(-1.2, -1.2, -0.35), 'note');
        this.sheetNote.visible = true;
      } else if (this.sheetNote) this.sheetNote.visible = false;
      if (this.on.lines && r.lines) for (const l of r.lines) s.line(l.map((p) => p.map((x) => x / this.L0)), COL.ink, this.g.lines, { opacity: 0.45 });
    } catch (e) { /* overlays are best-effort */ } finally { this._fbusy = false; }
  }

  drawTides(t, body) {
    const G = this.g.tides; G.clear();
    const c = body.pos.map((x) => x / this.L0), big = Math.max(...t.values.map(Math.abs)) || 1;
    t.values.forEach((w, k) => {
      const len = 0.16 * Math.abs(w) / big, ax = t.axes[k];
      this.stage.wide([c.map((x, q) => x - len * ax[q]), c.map((x, q) => x + len * ax[q])], w > 0 ? COL.unstable : COL.stable, 3, G);
    });
    if (this.tideNote) this.stage.dropLabel(this.tideNote);
    this.tideNote = this.stage.label(`Tide here: stretch ${sci(t.values[0], 1)} s⁻² (red), squeeze ${sci(t.values[2], 1)} s⁻² (green)`, new THREE.Vector3(c[0], c[1], c[2] + 0.2), 'note');
  }

  async refreshPair() {
    const s = this.stage;
    this.g.hill.clear(); this.g.lagr.clear();
    for (const l of this.pairLabels || []) s.dropLabel(l);
    this.pairLabels = [];
    const pr = this.pair();
    if (!pr || !(this.on.hill || this.on.lagr) || !this.an?.pair) return;
    const p = this.an.pair, f = this.display(this.frameAt());
    const iSmall = this.bodies[pr[0]].gm <= this.bodies[pr[1]].gm ? pr[0] : pr[1], iBig = iSmall === pr[0] ? pr[1] : pr[0];
    if (this.on.hill) {
      const hs = s.sphere(p.hill_radius / this.L0, COL.stable, this.g.hill, { flat: true, wire: true, opacity: 0.22 });
      hs.position.set(...f[iSmall].map((x) => x / this.L0));
      this.pairLabels.push(s.label(`Hill sphere of ${p.small}: ${sci(p.hill_radius / 1e3, 2)} km`, hs.position.clone().add(new THREE.Vector3(0, 0, p.hill_radius / this.L0)), 'note'));
      if (p.roche) {
        const rs = s.sphere(p.roche.fluid / this.L0, COL.unstable, this.g.hill, { flat: true, wire: true, opacity: 0.2 });
        rs.position.set(...f[iBig].map((x) => x / this.L0));
        this.pairLabels.push(s.label(`Roche limit for ${p.small}: ${sci(p.roche.fluid / 1e3, 2)} km (fluid), ${sci(p.roche.rigid / 1e3, 2)} km (rigid)` + (p.inside_roche_fluid ? '. It is inside!' : ''), rs.position.clone().add(new THREE.Vector3(0, 0, -p.roche.fluid / this.L0)), 'note'));
      }
    }
    if (this.on.lagr && this.frameKind === 'inertial') {
      for (const [name, q] of Object.entries(p.lagrange)) {
        const v = new THREE.Vector3(...q.map((x) => x / this.L0));
        s.mark(v, ['L4', 'L5'].includes(name) && p.mu < 0.0385 ? COL.stable : COL.unstable, this.g.lagr, 0.01);
        this.pairLabels.push(s.label(name, v, 'offset note'));
      }
      if (!p.circular_ok) toast(`The ${p.big}–${p.small} orbit has eccentricity ${p.eccentricity.toFixed(2)}: Lagrange points are true equilibria only for circular orbits, so treat these marks as a guide.`, true);
    } else if (this.on.lagr) toast('Lagrange point marks are drawn in the inertial frame; switch the frame to see them.');
  }

  async refreshRays() {
    const s = this.stage;
    this.g.rays.clear();
    for (const l of this.rayLabels || []) s.dropLabel(l);
    this.rayLabels = [];
    if (!this.on.rays || !this.bodies.length) return;
    const bodies = this.current(), k = this.sel >= 0 ? this.sel : bodies.reduce((m, q, i) => (q.gm > bodies[m].gm ? i : m), 0);
    const c = bodies[k].pos, R = bodies[k].radius;
    const spread = Math.max(6 * R, 0.25 * this.L0);
    const r = await api('nbody/rays', { bodies, origin: [c[0] - 1.6 * this.L0, c[1], c[2]], direction: [1, 0, 0], spread, n_rays: 11, length: 3.2 * this.L0 });
    let worst = 0;
    r.rays.forEach((ray) => {
      const pts = ray.pts.map((p) => p.map((x) => x / this.L0));
      s.line(pts, COL.magenta, this.g.rays, ray.trustworthy ? { opacity: 0.85 } : { dash: 0.03, opacity: 0.85 });
      if (!ray.blocked) worst = Math.max(worst, ray.deflection_arcsec);
    });
    const bad = r.rays.some((q) => !q.trustworthy);
    this.rayLabels.push(s.label(`Light bent by up to ${worst > 3600 ? (worst / 3600).toFixed(1) + '°' : sci(worst, 2) + '″'}` + (bad ? '. Dashed rays leave the weak-field model’s range: use the close-up for exact rays' : ' (too small to see; the number is real)'),
      new THREE.Vector3((c[0] + 1.3 * this.L0) / this.L0, c[1] / this.L0, c[2] / this.L0 + 0.08), 'note rel'));
  }

  // ---- Newtonian against 1PN ----------------------------------------------------------
  async compare() {
    if (this.bodies.length < 2) { toast('Add a second body to compare the two models.'); return; }
    const T = this.tChar() * 6;
    const r = await busy('Integrating the same start twice: Newtonian and 1PN…', () => api('nbody/run', { bodies: this.current(), t_end: T, compare: true, pn: false, n_frames: 500 }));
    const G = this.g.ghost; G.clear();
    this.bodies.forEach((b, i) => {
      this.stage.line(r.pos.map((f) => f[i].map((x) => x / this.L0)), COL.ink, G, { opacity: 0.8 });
      this.stage.line(r.other_pos.map((f) => f[i].map((x) => x / this.L0)), COL.magenta, G, { dash: 0.02 });
    });
    const p = inset.open('How far apart the Newtonian and 1PN bodies drift', 250);
    const [u, uname] = unitOf(T);
    const days = r.t.map((t) => t / u), ymax = Math.max(1e-3, ...r.divergence_m.flat());
    const ax = axes(p.g, p.w, p.hgt, [0, days[days.length - 1] || 1], [Math.max(ymax * 1e-8, 1e-6), ymax * 2], { logy: true, xlabel: `${uname} (six periods of the tightest pair)`, ylabel: 'separation (metres)' });
    this.bodies.forEach((b, i) => plotLine(p.g, ax, days, r.divergence_m.map((d) => Math.max(d[i], 1e-9)), css(b.color), 1.6));
    const worst = this.an ? Math.max(...this.an.regime.map((q) => q.pn1_fractional_size)) : 0;
    p.cap.textContent = `One curve per body, in its own colour. Ink trails are Newtonian, dashed magenta are 1PN. Largest correction in this system is ${sci(worst, 1)} of the Newtonian force; energy drift of this run ${sci(r.max_energy_drift, 1)}.`;
  }

  csv() {
    if (!this.frames.length || !this.bodies.length) return null;
    const head = ['t_s', ...this.bodies.flatMap((b) => ['x', 'y', 'z', 'vx', 'vy', 'vz'].map((k) => `${b.name}_${k}`))].join(',');
    return head + '\n' + this.frames.map((f) => [f.t, ...f.pos.flatMap((p, i) => [...p, ...f.vel[i]])].join(',')).join('\n');
  }
}

let _ring;
function ringTexture() {
  if (_ring) return _ring;
  const c = document.createElement('canvas'); c.width = c.height = 128;
  const g = c.getContext('2d'); g.strokeStyle = '#fff'; g.lineWidth = 6;
  g.beginPath(); g.arc(64, 64, 54, 0, 7); g.stroke();
  _ring = new THREE.CanvasTexture(c);
  return _ring;
}
