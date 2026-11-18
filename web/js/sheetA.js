// Sheet A - the three-body chart (CR3BP laboratory).
// Everything drawn here is requested from the Python core; this file only draws.
import { COL, css, V, THREE } from './stage.js';
import { api, busy, toast, callout, dl, act, why, slider, buildLegend, inset, axes, plotLine, sci, fix, km, h } from './ui.js';

const OUT_COL = [COL.stable, COL.soft, COL.unstable];           // captured, escaped, impacted
const OUT_NAME = ['captured', 'escaped', 'impacted'];

export class SheetA {
  constructor(stage, meta) {
    this.stage = stage; this.meta = meta;
    this.mu = meta.mu_earth_moon; this.C = 3.16;
    this.zHalo = 0.08; this.branch = 'first'; this.eps = 1e-6;
    this.on = { hill: true, zvs: false, family: true, halo: true, rich: false, stable: false, unstable: false, section: false };
    this.g = {};
    this.track = null; this.playhead = 0;
    this.title = 'Earth–Moon rotating frame. Lengths in Earth–Moon distances; the frame turns once a month.';
  }

  get isEM() { return Math.abs(this.mu - this.meta.mu_earth_moon) < 1e-9; }
  days(t) { return t * this.meta.tunit_s / 86400; }

  async enter() {
    const s = this.stage;
    s.clear();
    for (const k of ['hill', 'zvs', 'bodies', 'points', 'family', 'halo', 'rich', 'stable', 'unstable', 'section', 'mc', 'flight']) this.g[k] = s.group(k);
    this._view = null; this.view('all');
    this.controls();
    this.legend();
    await busy('Charting the Lagrange points…', async () => {
      await this.loadSystem();
      await this.loadHill();
    });
    this.loadFamily().then(() => this.loadHalo()).catch((e) => toast('Halo family: ' + e.message, true));
  }

  leave() { callout.close(); inset.close(); }

  // ---- title-block controls -------------------------------------------------
  controls() {
    const box = document.getElementById('sheet-controls');
    const sys = h('select', {},
      h('option', { value: this.meta.mu_earth_moon }, 'Earth–Moon  (μ = 0.012151)'),
      h('option', { value: 9.5388e-4 }, 'Sun–Jupiter  (μ ≈ 0.00095, approximate)'),
      h('option', { value: 3.0035e-6 }, 'Sun–Earth  (μ ≈ 0.000003, approximate)'),
      h('option', { value: 0.1 }, 'A heavy moon  (μ = 0.1)'),
      h('option', { value: 0.5 }, 'Equal masses  (μ = 0.5)'));
    sys.addEventListener('change', async () => {
      this.mu = +sys.value;
      callout.close();
      for (const k of ['family', 'halo', 'rich', 'stable', 'unstable', 'section', 'mc', 'flight']) this.g[k].clear();
      this.track = null;
      await busy('Recomputing for the new mass ratio…', async () => { await this.loadSystem(); await this.loadHill(); });
      if (this.isEM) this.loadFamily().then(() => this.loadHalo());
      else toast('Halo, manifold and Monte Carlo tools are set up for the Earth–Moon system. Lagrange points, stability and forbidden regions work for any mass ratio.');
    });
    this.cOut = h('p', { class: 'hint' });
    this.cSlider = slider('Energy level (Jacobi constant C; bigger C = less energy)', 2.95, 3.3, 0.001, this.C, (v) => v.toFixed(3),
      (v) => { this.C = v; this.loadHillSoon(); });
    const views = h('div', { class: 'acts', style: 'grid-template-columns:1fr 1fr' }, act('View whole system', () => this.view('all')), act('View Moon and L2', () => this.view('moon')));
    box.replaceChildren(h('label', { class: 'field' }, h('span', {}, 'System'), h('span'), sys), this.cSlider, this.cOut, views);
  }

  legend() {
    const T = (key, after) => async (on) => { this.on[key] = on; this.g[key].visible = on; if (on && after) await after(); };
    buildLegend([
      { group: 'Energy' },
      { id: 'hill', label: 'Forbidden region at this energy', color: COL.land, style: 'fill', on: true, toggle: T('hill') },
      { id: 'zvs', label: 'Edge of the forbidden region, in 3D', color: COL.landEdge, style: 'line', on: false, toggle: T('zvs', () => this.loadHill()) },
      { group: 'Parking loops near L2' },
      { id: 'family', label: 'All loop sizes (halo family)', color: COL.rule, style: 'line', on: true, toggle: T('family') },
      { id: 'halo', label: 'Selected loop (halo orbit)', color: COL.ink, style: 'wide', on: true, toggle: T('halo') },
      { id: 'rich', label: "First rough guess (Richardson's formula)", color: COL.ink, style: 'dash', on: false, toggle: T('rich') },
      { group: 'Free rides' },
      { id: 'stable', label: 'Free rides in (stable manifold)', color: COL.stable, style: 'line', on: false, toggle: T('stable', () => this.loadManifold('stable')) },
      { id: 'unstable', label: 'Free rides out (unstable manifold)', color: COL.unstable, style: 'line', on: false, toggle: T('unstable', () => this.loadManifold('unstable')) },
      { id: 'section', label: 'Target ring at the Moon line (capture curve)', color: COL.magenta, style: 'dot', on: false, toggle: T('section', () => this.loadManifold('stable', true)) },
    ]);
  }

  // ---- system: bodies and Lagrange points -------------------------------------
  async loadSystem() {
    const s = this.stage, mu = this.mu;
    this.lag = await api(`cr3bp/lagrange?mu=${mu}`);
    document.getElementById('sheet-title').textContent = this.isEM ? this.title
      : 'Rotating frame of the two bodies. Lengths in units of their separation; the frame turns once per orbit.';
    this.g.bodies.clear(); this.g.points.clear();
    for (const l of [...s.labels]) s.dropLabel(l);
    const rE = this.isEM ? this.meta.r_earth : 0.02, rM = this.isEM ? this.meta.r_moon : 0.008;
    s.sphere(rE, 0x2f5d7c, this.g.bodies).position.set(-mu, 0, 0);
    s.sphere(rM, 0x8a969d, this.g.bodies).position.set(1 - mu, 0, 0);
    const names = this.isEM ? ['Earth', 'Moon'] : ['Primary', 'Secondary'];
    s.label(names[0], V([-mu, 0, 0]), 'body offset');
    s.label(names[1], V([1 - mu, 0, 0]), 'body offset');
    for (const [name, d] of Object.entries(this.lag.points)) {
      const p = V(d.position);
      s.mark(p, d.unstable ? COL.unstable : COL.stable, this.g.points);
      s.label(name, p, 'offset', () => this.openPoint(name));
    }
    const c = this.lag.critical_jacobi;
    this.cSlider.input.min = (c.L4 - 0.06).toFixed(3); this.cSlider.input.max = (c.L1 + 0.12).toFixed(3);
    if (this.C > +this.cSlider.input.max || this.C < +this.cSlider.input.min) { this.C = c.L2 - 0.005; this.cSlider.input.value = this.C; this.cSlider.out.textContent = this.C.toFixed(3); }
  }

  // two standard views: the whole system, or the neighbourhood of the Moon and L1/L2
  view(which) {
    if (which === this._view) return;
    this._view = which;
    if (which === 'moon') this.stage.frameOn(new THREE.Vector3(1.07 - this.mu, 0, 0), 0.25);
    else this.stage.frameOn(new THREE.Vector3(0.22, 0, 0), 1.12);
  }

  openPoint(name) {
    const d = this.lag.points[name], crit = this.lag.critical_jacobi[name];
    if (['L3', 'L4', 'L5'].includes(name)) this.view('all');
    callout.open(`${name}`, () => V(d.position), (el) => {
      el.append(h('p', {}, d.verdict + '.'));
      const rows = [['Position x, y', `${fix(d.position[0], 6)}, ${fix(d.position[1], 6)}`], ['Energy number C here', fix(crit, 6)]];
      if (d.lambda !== undefined) {
        rows.push(['How fast a miss grows (λ)', fix(d.lambda, 4)], ['Sideways wobble rate (ν)', fix(d.nu, 4)], ['Up-down wobble rate (ω_z)', fix(d.omega_z, 4)]);
        if (d.e_fold_days) rows.push(['A displacement multiplies by e every', fix(d.e_fold_days, 2) + ' days']);
      }
      el.append(dl(rows));
      el.append(h('p', {}, 'Growth and wobble rates (eigenvalues)'),
        h('p', { class: 'hint' }, d.eigenvalues.map(([re, im]) => `${re.toFixed(4)} ${im < 0 ? '−' : '+'} ${Math.abs(im).toFixed(4)}i`).join(',  ')));
      el.append(why(d.unstable && d.verdict.includes('complex')
        ? 'The eigenvalues have left the imaginary axis as complex pairs with positive real part: a nudged particle spirals outward. This happens at L4 and L5 when the smaller body holds more than 3.85% of the total mass (Routh’s limit).'
        : d.unstable
        ? 'A real positive eigenvalue means a tiny push grows exponentially; the imaginary pairs are the oscillations that halo and Lyapunov orbits are built from.'
        : 'Every eigenvalue is purely imaginary: a nudged particle circles the point instead of leaving. That is why Trojan asteroids collect at L4 and L5.'));
      const acts = h('div', { class: 'acts' });
      acts.append(act(`Set the energy so the gate at ${name} is just open`, () => { this.C = crit - 0.004; this.cSlider.input.value = this.C; this.cSlider.out.textContent = this.C.toFixed(3); this.loadHill(); }));
      acts.append(act('Displace a particle by 10⁻⁸ and measure how fast it leaves', () => this.divergenceTest(name)));
      el.append(acts);
    });
  }

  async divergenceTest(name) {
    const d = this.lag.points[name];
    const rate = d.lambda ?? 0;
    const T = d.unstable && rate ? Math.min(14 / rate, 60) : 60;
    const s0 = [d.position[0] + 1e-8, d.position[1], 0, 0, 0, 0];
    const r = await busy('Flying the displaced particle…', () => api('cr3bp/propagate', { state: s0, t: T, mu: this.mu }));
    this.setTrack(r, `Particle released 10⁻⁸ from ${name}`);
    const dist = r.states.map((s) => Math.hypot(s[0] - d.position[0], s[1] - d.position[1], s[2]));
    const p = inset.open(`Distance from ${name} against time`, 240);
    const ax = axes(p.g, p.w, p.hgt, [0, T], [1e-9, Math.max(1e-2, Math.max(...dist))], { logy: true, xlabel: 'time (frame time units)', ylabel: 'distance' });
    plotLine(p.g, ax, r.t, dist, css(COL.ink), 1.6);
    if (d.unstable && rate) {
      // straight-line fit of ln(distance) where 1e-7 < d < 1e-3: the slope is the growth rate
      const xs = [], ys = [];
      r.t.forEach((t, i) => { if (dist[i] > 1e-7 && dist[i] < 1e-3) { xs.push(t); ys.push(Math.log(dist[i])); } });
      const n = xs.length, mx = xs.reduce((a, b) => a + b, 0) / n, my = ys.reduce((a, b) => a + b, 0) / n;
      const slope = xs.reduce((a, x, i) => a + (x - mx) * (ys[i] - my), 0) / xs.reduce((a, x) => a + (x - mx) ** 2, 0);
      plotLine(p.g, ax, [xs[0], xs[n - 1]], [Math.exp(my + slope * (xs[0] - mx)), Math.exp(my + slope * (xs[n - 1] - mx))], css(COL.magenta), 1, [5, 4]);
      p.cap.textContent = `Measured slope ${slope.toFixed(4)} per time unit; linear theory predicts λ = ${rate.toFixed(4)}. A straight line on this logarithmic scale is exponential growth.`;
    } else p.cap.textContent = 'No exponential growth: the particle stays within a bounded distance, as linear stability predicts.';
  }

  // ---- energy: forbidden region ----------------------------------------------
  loadHillSoon() { clearTimeout(this._ht); this._ht = setTimeout(() => this.loadHill(), 90); }

  async loadHill() {
    const lim = 1.5, n = 280;
    const r = await api(`cr3bp/hill?mu=${this.mu}&C=${this.C}&n=${n}&lim=${lim}&zlevels=${this.on.zvs ? 9 : 0}`);
    this.cOut.textContent = r.regime;
    const s = this.stage, G = this.g.hill;
    G.clear(); this.g.zvs.clear();
    // tint: land where motion is forbidden, deeper water where more speed is available
    const cv = document.createElement('canvas'); cv.width = cv.height = n;
    const g = cv.getContext('2d'), img = g.createImageData(n, n);
    const land = [0xea, 0xd9, 0xac], w0 = [0xe7, 0xf0, 0xf1], w1 = [0xb9, 0xd4, 0xde];
    for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) {
      const w = r.grid[j * n + i], o = ((n - 1 - j) * n + i) * 4;
      let c;
      if (w < 0) c = land; else { const t = Math.min(1, Math.floor(w / 0.12) * 0.2); c = w0.map((a, k) => a + (w1[k] - a) * t); }
      img.data[o] = c[0]; img.data[o + 1] = c[1]; img.data[o + 2] = c[2]; img.data[o + 3] = 255;
    }
    g.putImageData(img, 0, 0);
    const big = document.createElement('canvas'); big.width = big.height = 1200;
    const b = big.getContext('2d'); b.imageSmoothingEnabled = true; b.drawImage(cv, 0, 0, 1200, 1200);
    b.strokeStyle = 'rgba(20,48,74,0.13)'; b.lineWidth = 1;      // graticule every 0.25 length units
    for (let v = -lim; v <= lim + 1e-9; v += 0.25) { const q = (v + lim) / (2 * lim) * 1200; b.beginPath(); b.moveTo(q, 0); b.lineTo(q, 1200); b.moveTo(0, q); b.lineTo(1200, q); b.stroke(); }
    const tex = new THREE.CanvasTexture(big); tex.colorSpace = THREE.SRGBColorSpace;
    const plane = new THREE.Mesh(new THREE.PlaneGeometry(2 * lim, 2 * lim), new THREE.MeshBasicMaterial({ map: tex, side: THREE.DoubleSide, toneMapped: false, depthWrite: false }));
    plane.renderOrder = -10;      // drawn first and never hides what dips below the plane
    G.add(plane);
    for (const c of r.contours) {
      const pts = c.pts.map((p) => [p[0], p[1], c.z]);
      if (c.z === 0) s.line(pts, COL.landEdge, G);
      else { s.line(pts, COL.landEdge, this.g.zvs, { opacity: 0.7 }); s.line(pts.map((p) => [p[0], p[1], -p[2]]), COL.landEdge, this.g.zvs, { opacity: 0.35 }); }
    }
    if (this.on.zvs) for (const c of r.contours_xz) s.line(c.pts.map((p) => [p[0], 0, p[1]]), COL.landEdge, this.g.zvs, { opacity: 0.7 });
    G.visible = this.on.hill; this.g.zvs.visible = this.on.zvs;
  }

  // ---- halo family and selected halo -----------------------------------------
  async loadFamily() {
    this.family = await busy('Growing the L2 halo family (Richardson seed, differential correction, continuation)…',
      () => api(`cr3bp/family?mu=${this.mu}`));
    const G = this.g.family; G.clear();
    for (const o of this.family.ghosts) this.stage.line(o, COL.rule, G, { loop: true, opacity: 0.8 });
    G.visible = this.on.family;
  }

  async loadHalo() {
    const o = await api(`cr3bp/halo?z=${this.zHalo}&mu=${this.mu}&branch=${this.branch}`);
    this.halo = o;
    const s = this.stage;
    this.g.halo.clear(); this.g.rich.clear();
    const pts = o.points.map((p) => p.slice(0, 3));
    s.wide(pts, COL.ink, 2.6, this.g.halo, { loop: true });
    if (o.richardson) s.line(o.richardson, COL.ink, this.g.rich, { dash: 0.012, loop: true });
    this.g.rich.visible = this.on.rich; this.g.halo.visible = this.on.halo;
    if (this.haloLabel) s.dropLabel(this.haloLabel);
    this.haloLabel = s.label('Halo orbit', V(o.state0), 'offset', () => this.openHalo());
    for (const k of ['stable', 'unstable', 'section', 'mc']) this.g[k].clear();
    if (this.on.stable) this.loadManifold('stable');
    if (this.on.unstable) this.loadManifold('unstable');
    if (this.on.section) this.loadManifold('stable', true);
    return o;
  }

  haloRows() {
    const o = this.halo;
    return dl([
      ['Height above the Moon’s plane', km(o.Az_km)],
      ['Time for one lap', `${fix(o.period, 4)} (${fix(o.period_days, 2)} days)`],
      ['Energy number C', fix(o.jacobi, 5)],
      ['Instability number k', o.stability_index > 10 ? fix(o.stability_index, 1) : fix(o.stability_index, 4)],
      ['A mistake grows each lap by', o.unstable_multiplier > 1.001 ? '≈ ' + sci(o.unstable_multiplier, 1) : '1 (stable)'],
      ['The loop closes to within', sci(o.residual, 1)],
    ]);
  }

  openHalo() {
    const f = this.family;
    this.view('moon');
    callout.open('Parking loop near L2 (halo orbit)', () => V(this.halo.state0), (el) => {
      const rows = h('div', {}, this.haloRows());
      const sl = slider('Size (z at the far crossing)', Math.max(0.004, f.z_min + 0.002).toFixed(3), (f.z_max - 0.001).toFixed(3), 0.001, this.zHalo,
        (v) => km(v * this.meta.lunit_km), null, async (v) => { this.zHalo = v; await this.loadHalo(); rows.replaceChildren(this.haloRows()); });
      const br = h('select', {}, h('option', { value: 'first' }, 'Classic halo part of the family'), h('option', { value: 'last' }, 'Near-rectilinear part (towards the Moon)'));
      br.value = this.branch;
      br.addEventListener('change', async () => { this.branch = br.value; try { await this.loadHalo(); rows.replaceChildren(this.haloRows()); } catch (e) { toast(e.message, true); } });
      el.append(sl, br, rows);
      el.append(why('k measures instability: after one lap a mistake has grown by about 2k. Small loops are wild (k ≈ 600); they calm down as they get bigger and reach towards the Moon.'));
      const acts = h('div', { class: 'acts' });
      acts.append(act('Fly three laps (watch the accuracy check)', async () => this.setTrack(await busy('Integrating…', () => api('cr3bp/propagate', { state: this.halo.state0, t: 3 * this.halo.period, mu: this.mu })), 'Spacecraft on the halo, no corrections')));
      acts.append(act('Test starting errors (random trials)', () => this.openMC()));
      acts.append(act('Newton against Einstein’s correction over one lap', () => this.pnCompare(), 'rel-note'));
      acts.append(act('Crossing map at this energy (Poincaré section)', () => this.poincare()));
      el.append(acts);
    });
  }

  // ---- manifolds --------------------------------------------------------------
  async loadManifold(kind, sectionOnly = false) {
    const r = await busy(`Growing the ${kind} manifold from the monodromy eigenvectors…`,
      () => api(`cr3bp/manifold?z=${this.zHalo}&kind=${kind}&eps=${this.eps}&n=44&branch=${this.branch}&mu=${this.mu}`));
    const s = this.stage, col = kind === 'stable' ? COL.stable : COL.unstable;
    if (!sectionOnly || this.on[kind]) {
      const G = this.g[kind]; G.clear();
      for (const t of r.tubes) s.line(t.pts, t.impact ? COL.unstable : col, G, { opacity: t.impact ? 0.35 : 0.75 });
      G.visible = this.on[kind];
    }
    if (kind === 'stable' && r.section) {
      const G = this.g.section; G.clear();
      const x = 1 - this.mu, st = r.section.states;
      const ys = st.map((q) => q[1]), zs = st.map((q) => q[2]);
      const y0 = Math.min(...ys) - 0.02, y1 = Math.max(...ys, 0.02) + 0.02, z0 = Math.min(...zs) - 0.03, z1 = Math.max(...zs) + 0.03;
      s.line([[x, y0, z0], [x, y1, z0], [x, y1, z1], [x, y0, z1]], COL.magenta, G, { loop: true, opacity: 0.8 });
      s.line(st.map((q) => [q[0], q[1], q[2]]), COL.magenta, G, { loop: true });
      const geo = new THREE.BufferGeometry();
      geo.setAttribute('position', new THREE.Float32BufferAttribute(st.flatMap((q) => [q[0], q[1], q[2]]), 3));
      G.add(new THREE.Points(geo, new THREE.PointsMaterial({ color: COL.magenta, size: 5, sizeAttenuation: false })));
      G.visible = this.on.section;
      this.sectionData = r;
      if (this.on.section) {
        const c = r.convergence;
        toast(c && c.ok ? `Capture curve on the Moon plane. Convergence check: repeating with ε ten times smaller moves the curve by at most ${sci(c.relative * 100, 1)}% of its size.` : 'Capture curve on the Moon plane.', false, 9000);
        this.sectionInset(r);
      }
    }
    return r;
  }

  sectionInset(r) {
    const st = r.section.states;
    const p = inset.open('Capture curve where the stable manifold meets the Moon plane', 250);
    const ys = st.map((q) => q[1]), vs = st.map((q) => q[4]);
    const pad = (a, b) => [a - 0.08 * (b - a), b + 0.08 * (b - a)];
    const ax = axes(p.g, p.w, p.hgt, pad(Math.min(...ys), Math.max(...ys)), pad(Math.min(...vs), Math.max(...vs)), { xlabel: 'y on the plane x = 1 − μ', ylabel: 'velocity vy' });
    plotLine(p.g, ax, [...ys, ys[0]], [...vs, vs[0]], css(COL.magenta), 1.4);
    p.g.fillStyle = css(COL.magenta);
    st.forEach((q) => { p.g.beginPath(); p.g.arc(ax.X(q[1]), ax.Y(q[4]), 2, 0, 7); p.g.fill(); });
    const c = r.convergence;
    p.cap.textContent = `Each dot is one manifold trajectory (ε = ${sci(r.eps)} = ${sci(r.eps_km * 1000, 1)} m). States on this loop coast into the halo with no fuel.` +
      (c && c.ok ? ` With ε/10 the loop moves by ${sci(c.max_distance, 1)} (${sci(c.relative * 100, 1)}% of its size).` : '');
  }

  // ---- flights ----------------------------------------------------------------
  setTrack(r, what) {
    const s = this.stage, G = this.g.flight;
    G.clear();
    this.track = r; this.track.what = what; this.playhead = 0;
    this.trail = s.wide(r.states.map((q) => q.slice(0, 3)), COL.magenta, 2, G);
    this.craft = s.sphere(0.006, COL.magenta, G, { flat: true });
    this.craft.position.set(...r.states[0].slice(0, 3));
    this.playing = true;
  }

  async pnCompare() {
    const o = this.halo;
    const r = await busy('Integrating the same start twice: Newtonian and 1PN…', () => api('cr3bp/propagate', { state: o.state0, t: o.period, mu: this.mu, compare_pn: true }));
    this.setTrack(r, 'Newtonian trajectory; the 1PN one is indistinguishable at this scale');
    const p = inset.open('Gap between the Newtonian and 1PN trajectories', 240);
    const td = r.t.map((t) => this.days(t)), gap = r.pn_gap_m.map((v) => Math.max(v, 1e-6));
    const ax = axes(p.g, p.w, p.hgt, [0, td[td.length - 1]], [1e-4, Math.max(10, ...gap)], { logy: true, xlabel: 'days since the common start', ylabel: 'separation (metres)' });
    plotLine(p.g, ax, td, gap, css(COL.magenta), 1.6);
    p.cap.innerHTML = `At the start the 1PN acceleration is ${sci(r.pn_ratio, 1)} of the Newtonian one. After one revolution the two paths differ by ${sci(gap[gap.length - 1], 1)} m. The growth is the halo's own instability amplifying a tiny push; it is not relativity becoming strong.`;
  }

  async poincare() {
    const C = this.halo.jacobi;
    const r = await busy('Collecting crossings of the plane y = 0…', () => api(`cr3bp/poincare?C=${C}&mu=${this.mu}&n=30&t_max=40`));
    const p = inset.open(`Poincaré section y = 0, vy > 0 at C = ${C.toFixed(4)}`, 300);
    const all = r.orbits.flatMap((o) => o.pts);
    const xr = [0.84, 1.16], yr = [-1.2, 1.2];
    const ax = axes(p.g, p.w, p.hgt, xr, yr, { xlabel: 'x', ylabel: 'vx' });
    r.orbits.forEach((o, i) => {
      p.g.fillStyle = `hsl(${205 + (i * 47) % 150}, 55%, ${28 + (i % 3) * 10}%)`;
      for (const q of o.pts) if (q[0] > xr[0] && q[0] < xr[1] && Math.abs(q[1]) < 1.2) p.g.fillRect(ax.X(q[0]) - 0.9, ax.Y(q[1]) - 0.9, 1.8, 1.8);
    });
    p.cap.textContent = `${all.length} crossings from ${r.orbits.length} planar starts, all with the selected halo's Jacobi constant. Closed loops are regular orbits around the Moon; scattered dust is chaotic motion.`;
  }

  // ---- Monte Carlo -------------------------------------------------------------
  openMC() {
    if (!this.isEM) { toast('The injection experiment is set up for the Earth–Moon system.', true); return; }
    let sr = 0, sv = -9;
    callout.open('Starting errors', () => V(this.halo.state0), (el) => {
      el.append(h('p', {}, 'The spacecraft is placed on the stable manifold at the Moon plane, with a random error, and then coasts with no correction burns.'));
      const s1 = slider('Position error (1σ per axis)', -3, 3, 0.1, sr, (v) => km(10 ** v), (v) => { sr = v; });
      const s2 = slider('Velocity error (1σ per axis)', -4, 3, 0.1, -4, (v) => (v <= -3.95 ? 'none' : sci(10 ** v, 1) + ' mm/s'), (v) => { sv = v <= -3.95 ? -9 : v; });
      const pn = h('input', { type: 'checkbox' });
      const res = h('div');
      el.append(s1, s2, h('label', { class: 'rel-note' }, pn, ' include Einstein’s correction (1PN)'), res);
      const acts = h('div', { class: 'acts' });
      acts.append(act('Run 720 trajectories', async () => {
        const r = await busy('Flying the cloud of perturbed spacecraft…', () => api('cr3bp/montecarlo', { z: this.zHalo, sigma_r_km: 10 ** sr, sigma_v_mms: sv < -8 ? 0 : 10 ** sv, n: 120, pn: pn.checked, n_draw: 48 }));
        this.g.mc.clear();
        for (const t of r.tracks) this.stage.line(t.pts, OUT_COL[t.A], this.g.mc, { opacity: 0.6 });
        const bar = (d) => h('div', { class: 'bar' }, OUT_NAME.map((n, i) => h('i', { style: `width:${d[n] * 100}%;background:${css(OUT_COL[i])}`, title: n })));
        res.replaceChildren(h('div', { class: 'bars' },
          h('span', {}, 'A: tracks the halo'), bar(r.A), h('b', {}, Math.round(r.A.captured * 100) + '%'),
          h('span', {}, 'B: stays near the Moon'), bar(r.B), h('b', {}, Math.round(r.B.captured * 100) + '%')),
          h('p', { class: 'hint' }, `Green captured, grey escaped, red hit the Moon. ${r.n} trajectories; sampling error about ±${Math.round(r.standard_error * 100)} points. A: within ${km(r.d_tol_km)} of the halo for a full revolution. B: within ${km(r.R_B_km)} of the Moon (two Hill radii) for the same length of time.`));
      }));
      acts.append(act("Show the paper's capture curves", () => this.paperCurves()));
      el.append(acts);
    });
  }

  async paperCurves() {
    let r;
    try { r = await api('cr3bp/paper_results'); } catch (e) { toast(e.message, true); return; }
    const p = inset.open('Capture probability against position error (definition A)', 280);
    const ax = axes(p.g, p.w, p.hgt, [1e-4, 1e3], [0, 1], { logx: true, xlabel: 'position error, 1σ per axis (km)', ylabel: 'fraction captured' });
    const zs = r.halos.map((q) => q.z_far);
    const near = zs.reduce((b, z, i) => (Math.abs(z - this.zHalo) < Math.abs(zs[b] - this.zHalo) ? i : b), 0);
    r.halos.forEach((q, i) => plotLine(p.g, ax, q.sigma_km, q.pos.A, i === near ? css(COL.magenta) : css(COL.ink), i === near ? 2.2 : 0.8));
    const q = r.halos[near];
    p.cap.textContent = `One curve per halo size, from the batch run (${r.n_samples} draws × ${r.n_phases} injection points per error level). Highlighted: the size nearest yours (${km(q.Az_km)}), where half the spacecraft are captured at ${sci(q.s50.pos.A.km * 1000, 1)} m.`;
  }

  // ---- frame update, status, export ---------------------------------------------
  tick(dt, playing, speed) {
    callout.place(this.stage);
    const tr = this.track;
    if (!tr) return { clock: 'No spacecraft flying. Select a Lagrange point or the halo orbit.', drift: '', frac: 0 };
    if (playing) {
      const T = tr.t[tr.t.length - 1];
      this.playhead = Math.max(0, Math.min(1, this.playhead + dt * speed * 0.9 / Math.max(Math.abs(T), 1e-9)));
    }
    const f = this.playhead * (tr.t.length - 1), i = Math.floor(f), j = Math.min(i + 1, tr.t.length - 1), w = f - i;
    const a = tr.states[i], b = tr.states[j];
    this.craft.position.set(a[0] + w * (b[0] - a[0]), a[1] + w * (b[1] - a[1]), a[2] + w * (b[2] - a[2]));
    const t = tr.t[i] + w * (tr.t[j] - tr.t[i]);
    return {
      clock: `t = ${t.toFixed(3)}` + (this.isEM ? `  (${this.days(t).toFixed(2)} days)` : ''),
      drift: `Accuracy check: energy number changed by ${sci(tr.jacobi_drift, 1)}`, bad: tr.jacobi_drift > 1e-9, frac: this.playhead,
    };
  }

  seek(frac) { this.playhead = frac; }

  csv() {
    if (!this.track) {
      if (!this.halo) return null;
      return 'x,y,z,vx,vy,vz\n' + this.halo.points.map((p) => p.join(',')).join('\n');
    }
    return 't,x,y,z,vx,vy,vz\n' + this.track.states.map((s, i) => [this.track.t[i], ...s].join(',')).join('\n');
  }
}
