// Close-up inset: exact motion of a light companion and of light around ONE
// non-rotating mass (Schwarzschild geometry).  Lengths are in units of GM/c^2,
// so the picture is the same for every mass; the caption converts to metres.
import { COL, css } from './stage.js';
import { api, inset, slider, sci, h, toast } from './ui.js';

export async function openSchwarzschild(body, meta) {
  const M = body.gm / meta.c ** 2;            // GM/c^2 in metres
  const tM = body.gm / meta.c ** 3;           // GM/c^3 in seconds
  const p = inset.open(`Close-up of ${body.name}: exact strong-field motion`, 330);
  let rp = 8, ra = 20, orbit = null, rays = null, alive = true, k = 0;
  inset.onclose = () => { alive = false; };
  const showRays = h('input', { type: 'checkbox' });
  const s1 = slider('Closest approach', 4.1, 30, 0.1, rp, (v) => `${v.toFixed(1)} GM/c²`, null, (v) => { rp = v; if (ra < rp + 0.5) { ra = rp + 0.5; } load(); });
  const s2 = slider('Farthest point', 6, 60, 0.5, ra, (v) => `${v.toFixed(1)} GM/c²`, null, (v) => { ra = Math.max(v, rp + 0.5); load(); });
  const read = h('p');
  p.extra.append(s1, s2, h('label', { class: 'rel-note' }, showRays, ' light rays arriving from the left'), read);
  showRays.addEventListener('change', async () => { if (showRays.checked && !rays) rays = await api('schwarzschild/rays?n=14&bmin=3.5&bmax=16'); });

  async function load() {
    try { orbit = await api(`schwarzschild/orbit?rp=${rp}&ra=${ra}&n_orbits=4`); k = 0; }
    catch (e) { toast(e.message, true); return; }
    const where = body.radius > 6 * M
      ? ` For ${body.name} all of this lies deep inside the body (its surface is at ${sci(body.radius / M, 1)} GM/c²), so the picture shows what the same mass would do if it were compact.`
      : '';
    p.cap.innerHTML = `Scale: GM/c² = ${sci(M, 2)} m for this mass. The orbit turns by ${orbit.precession_deg_per_orbit.toFixed(2)}° per revolution; the weak-field formula 6πGM/c²a(1−e²) would say ${orbit.weak_field_deg_per_orbit.toFixed(2)}°.` +
      (orbit.plunged ? ' <b>This orbit plunges through the horizon.</b>' : '') + where;
  }
  await load();

  const g = p.g, W = p.w, H = p.hgt;
  function draw() {
    if (!alive) return;
    requestAnimationFrame(draw);
    if (!orbit) return;
    const R = 1.12 * Math.max(ra, 17), sc = Math.min(W, H) / (2 * R), cx = W / 2, cy = H / 2;
    const X = (x) => cx + x * sc, Y = (y) => cy - y * sc;
    g.clearRect(0, 0, W, H);
    const ring = (r, dash, col, name) => {
      g.save(); g.setLineDash(dash); g.strokeStyle = col; g.lineWidth = 1; g.beginPath(); g.arc(cx, cy, r * sc, 0, 7); g.stroke(); g.restore();
      g.fillStyle = col; g.font = 'italic 12px Alegreya, serif'; g.fillText(name, cx + r * sc * 0.72 + 3, cy - r * sc * 0.72 - 2);
    };
    ring(orbit.isco, [2, 3], css(COL.ink), 'innermost stable orbit, 6');
    ring(orbit.photon_sphere, [5, 3], css(COL.magenta), 'photon sphere, 3');
    g.fillStyle = css(COL.ink); g.beginPath(); g.arc(cx, cy, 2 * sc, 0, 7); g.fill();
    if (showRays.checked && rays) for (const r of rays.rays) {
      g.strokeStyle = css(COL.magenta); g.globalAlpha = r.captured ? 0.35 : 0.8; g.lineWidth = 0.9; g.beginPath();
      r.x.forEach((x, i) => (i ? g.lineTo(X(x), Y(r.y[i])) : g.moveTo(X(x), Y(r.y[i])))); g.stroke(); g.globalAlpha = 1;
    }
    g.strokeStyle = css(COL.ink); g.lineWidth = 1.2; g.beginPath();
    for (let i = 0; i <= k; i++) (i ? g.lineTo(X(orbit.x[i]), Y(orbit.y[i])) : g.moveTo(X(orbit.x[i]), Y(orbit.y[i])));
    g.stroke();
    g.fillStyle = css(COL.unstable); g.beginPath(); g.arc(X(orbit.x[k]), Y(orbit.y[k]), 4, 0, 7); g.fill();
    const tau = orbit.tau[k], t = orbit.t[k];
    read.innerHTML = `Clock on the orbiter: <b>${(tau * tM).toPrecision(4)} s</b> · far-away clock: <b>${(t * tM).toPrecision(4)} s</b><br>The orbiter’s clock has fallen ${(100 * (1 - tau / Math.max(t, 1e-12))).toFixed(1)}% behind.` +
      (showRays.checked && rays ? `<br>Rays aimed closer than ${rays.b_crit.toFixed(2)} GM/c² are captured (faint).` : '');
    k = (k + 3) % orbit.x.length;
  }
  draw();
}
