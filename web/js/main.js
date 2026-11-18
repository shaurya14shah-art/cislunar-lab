// Entry point: owns the stage, the two sheets, the time ruler and the exports.
import { Stage, css, COL } from './stage.js';
import { api, toast, download, callout, inset } from './ui.js';
import { SheetA } from './sheetA.js';
import { SheetB } from './sheetB.js';

const stage = new Stage(document.getElementById('stage'));
const $ = (id) => document.getElementById(id);
let sheet = null, playing = false, direction = 1, speed = 1, scrubbing = false;

function setPlaying(p) { playing = p; $('play').textContent = p ? 'Pause' : 'Play'; }

async function show(key, sheets) {
  if (sheet) sheet.leave();
  for (const b of document.querySelectorAll('#sheets button')) b.setAttribute('aria-pressed', String(b.dataset.sheet === key));
  sheet = sheets[key];
  $('sheet-title').textContent = sheet.title;
  setPlaying(key === 'B');
  direction = 1;
  try { await sheet.enter(); } catch (e) { toast('Could not start this sheet: ' + e.message, true, 15000); }
}

async function boot() {
  let meta;
  try { meta = await api('meta'); } catch (e) {
    toast('The physics server is not answering. Start it with:  python -m cislab.server', true, 600000);
    return;
  }
  const sheets = { A: new SheetA(stage, meta), B: new SheetB(stage, meta) };
  for (const b of document.querySelectorAll('#sheets button')) b.addEventListener('click', () => show(b.dataset.sheet, sheets));

  stage.onFrame = (dt) => {
    if (!sheet || !sheet.tick) return;
    if (sheet.playing) { setPlaying(true); direction = 1; sheet.playing = false; }     // a sheet may ask to start playback
    const st = sheet.tick(dt, playing && !scrubbing, direction * speed);
    if (!st) return;
    $('clock').textContent = st.clock;
    $('drift').textContent = st.drift;
    $('drift').classList.toggle('bad', !!st.bad);
    if (!scrubbing) $('scrub').value = Math.round(st.frac * 1000);
    if (playing && sheet instanceof SheetA && (st.frac >= 1 && direction > 0 || st.frac <= 0 && direction < 0)) setPlaying(false);
  };

  $('play').addEventListener('click', () => {
    if (!playing && sheet instanceof SheetA && sheet.playhead >= 1) sheet.seek(0);
    direction = 1; setPlaying(!playing);
  });
  $('rev').addEventListener('click', () => { direction = -1; setPlaying(true); });
  $('speed').addEventListener('input', (e) => { speed = 10 ** (+e.target.value); $('speed-out').textContent = '×' + (speed < 1 ? speed.toFixed(2) : speed.toFixed(speed < 10 ? 1 : 0)); });
  $('scrub').addEventListener('pointerdown', () => { scrubbing = true; });
  $('scrub').addEventListener('input', (e) => sheet?.seek(+e.target.value / 1000));
  addEventListener('pointerup', () => { scrubbing = false; });
  addEventListener('keydown', (e) => {
    if (e.target.matches('input, select, textarea')) return;
    if (e.code === 'Space') { e.preventDefault(); $('play').click(); }
    if (e.key === 'Escape') { callout.close(); inset.close(); }
  });

  $('png').addEventListener('click', () => {
    // WebGL picture plus the label layer, composited so the saved image matches the screen
    const src = stage.canvas, c = document.createElement('canvas');
    c.width = src.width; c.height = src.height;
    const g = c.getContext('2d'), k = src.width / innerWidth;
    stage.renderer.render(stage.scene, stage.camera);
    g.drawImage(src, 0, 0);
    g.textBaseline = 'middle';
    for (const l of stage.labels) {
      if (l.el.style.display === 'none') continue;
      const r = l.el.getBoundingClientRect(), cs = getComputedStyle(l.el);
      g.font = `${cs.fontStyle} ${cs.fontWeight} ${parseFloat(cs.fontSize) * k}px ${cs.fontFamily}`;
      g.lineWidth = 4 * k; g.strokeStyle = css(COL.water); g.fillStyle = cs.color;
      const x = (r.left + 5) * k, y = (r.top + r.height / 2) * k;
      g.strokeText(l.el.textContent, x, y); g.fillText(l.el.textContent, x, y);
    }
    g.font = `${12 * k}px "Alegreya Sans", sans-serif`; g.fillStyle = css(COL.soft);
    g.fillText(`Cislunar Lab · ${$('clock').textContent} · ${$('drift').textContent}`, 12 * k, c.height - 12 * k);
    download('cislunar-lab.png', c.toDataURL('image/png'));
  });
  $('csv').addEventListener('click', () => {
    const text = sheet?.csv();
    if (!text) { toast('Nothing to save yet: fly a trajectory or run the sandbox first.'); return; }
    download(sheet instanceof SheetA ? 'cr3bp-trajectory.csv' : 'nbody-history.csv', text);
  });

  const start = new URLSearchParams(location.search).get('sheet') === 'B' ? 'B' : 'A';
  await show(start, sheets);
  window.lab = { stage, sheets, get sheet() { return sheet; } };      // handy in the browser console
}
boot();
