// The stage: one Three.js scene that fills the window, plus an HTML label layer.
// Physics coordinates are used directly (z is "up"), so what you see is what
// the Python core computed - no hidden axis swaps.
import * as THREE from 'three';
import { OrbitControls } from '../vendor/OrbitControls.js';
import { Line2 } from '../vendor/lines/Line2.js';
import { LineGeometry } from '../vendor/lines/LineGeometry.js';
import { LineMaterial } from '../vendor/lines/LineMaterial.js';

export const COL = {
  water: 0xe7f0f1, deep: 0xc6dce4, ink: 0x14304a, soft: 0x4d6679, rule: 0xa9c1c9, land: 0xead9ac,
  landEdge: 0xb89b52, magenta: 0xb0127a, stable: 0x1c7557, unstable: 0xc0432b, paper: 0xf8fbfb,
};
export const css = (n) => '#' + n.toString(16).padStart(6, '0');

export class Stage {
  constructor(canvas) {
    this.canvas = canvas;
    this.renderer = new THREE.WebGLRenderer({ canvas, antialias: true, preserveDrawingBuffer: true });
    this.renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
    this.renderer.setClearColor(COL.water);
    this.scene = new THREE.Scene();
    this.camera = new THREE.PerspectiveCamera(38, 1, 0.001, 500);
    this.camera.up.set(0, 0, 1);
    this.camera.position.set(0.3, -2.6, 2.2);
    this.controls = new OrbitControls(this.camera, canvas);
    this.controls.enableDamping = true;
    this.controls.dampingFactor = 0.12;
    this.controls.target.set(0.3, 0, 0);
    this.scene.add(new THREE.HemisphereLight(0xffffff, 0x9fb6c0, 2.2));
    const sun = new THREE.DirectionalLight(0xffffff, 1.4);
    sun.position.set(-3, -2, 4);
    this.scene.add(sun);
    this.root = new THREE.Group();          // everything belonging to the current sheet
    this.scene.add(this.root);
    this.labels = [];                       // {el, pos:Vector3 | fn}
    this.labelLayer = document.getElementById('labels');
    this.fat = [];                          // wide-line materials need the viewport size
    this.onFrame = () => {};
    this.onPick = () => {};
    this.pickables = [];
    addEventListener('resize', () => this.resize());
    this.resize();
    canvas.addEventListener('pointerdown', (e) => { this._down = [e.clientX, e.clientY]; });
    canvas.addEventListener('pointerup', (e) => {
      if (!this._down || Math.hypot(e.clientX - this._down[0], e.clientY - this._down[1]) > 4) return;
      const m = new THREE.Vector2((e.clientX / innerWidth) * 2 - 1, -(e.clientY / innerHeight) * 2 + 1);
      const rc = new THREE.Raycaster();
      rc.setFromCamera(m, this.camera);
      const hit = rc.intersectObjects(this.pickables, false)[0];
      this.onPick(hit ? hit.object : null, hit ? hit.point : null);
    });
    this.clock = new THREE.Clock();
    const loop = () => {
      requestAnimationFrame(loop);
      const dt = Math.min(this.clock.getDelta(), 0.1);
      this.controls.update();
      this.onFrame(dt);
      this.renderer.render(this.scene, this.camera);
      this.placeLabels();
    };
    loop();
  }

  resize() {
    const w = innerWidth, h = innerHeight;
    this.renderer.setSize(w, h, false);
    this.camera.aspect = w / h;
    this.camera.updateProjectionMatrix();
    for (const m of this.fat) m.resolution.set(w, h);
  }

  clear() {
    this.root.traverse((o) => { o.geometry?.dispose?.(); });
    this.root.clear();
    for (const l of this.labels) l.el.remove();
    this.labels = [];
    this.pickables = [];
    this.fat = [];
  }

  group(name) {
    const g = new THREE.Group();
    g.name = name;
    this.root.add(g);
    return g;
  }

  // thin 1-pixel line: used for the many curves of a tube or a family
  line(points, color, parent, opts = {}) {
    const geo = new THREE.BufferGeometry();
    geo.setAttribute('position', new THREE.Float32BufferAttribute(points.flat(), 3));
    const mat = opts.dash
      ? new THREE.LineDashedMaterial({ color, dashSize: opts.dash, gapSize: opts.dash * 0.7, transparent: true, opacity: opts.opacity ?? 1 })
      : new THREE.LineBasicMaterial({ color, transparent: true, opacity: opts.opacity ?? 1 });
    const l = opts.loop ? new THREE.LineLoop(geo, mat) : new THREE.Line(geo, mat);
    if (opts.dash) l.computeLineDistances();
    (parent || this.root).add(l);
    return l;
  }

  // wide line: used for the few curves that matter most
  wide(points, color, width, parent, opts = {}) {
    const geo = new LineGeometry();
    const flat = points.flat();
    if (opts.loop) flat.push(flat[0], flat[1], flat[2]);
    geo.setPositions(flat);
    const mat = new LineMaterial({ color, linewidth: width, transparent: true, opacity: opts.opacity ?? 1,
      dashed: !!opts.dash, dashSize: opts.dash || 1, gapSize: (opts.dash || 1) * 0.6 });
    mat.resolution.set(innerWidth, innerHeight);
    this.fat.push(mat);
    const l = new Line2(geo, mat);
    l.computeLineDistances();
    (parent || this.root).add(l);
    return l;
  }

  sphere(radius, color, parent, opts = {}) {
    const m = new THREE.Mesh(new THREE.SphereGeometry(radius, opts.wire ? 18 : 40, opts.wire ? 10 : 24),
      opts.flat ? new THREE.MeshBasicMaterial({ color, transparent: true, opacity: opts.opacity ?? 1, wireframe: !!opts.wire })
        : new THREE.MeshStandardMaterial({ color, roughness: 0.9, metalness: 0, transparent: (opts.opacity ?? 1) < 1, opacity: opts.opacity ?? 1 }));
    (parent || this.root).add(m);
    return m;
  }

  // a small chart symbol (cross in a circle) that always faces the viewer
  mark(pos, color, parent, size = 0.012) {
    const c = document.createElement('canvas');
    c.width = c.height = 64;
    const g = c.getContext('2d');
    g.strokeStyle = css(color); g.lineWidth = 5;
    g.beginPath(); g.arc(32, 32, 20, 0, 7); g.stroke();
    g.beginPath(); g.moveTo(32, 4); g.lineTo(32, 60); g.moveTo(4, 32); g.lineTo(60, 32); g.stroke();
    const s = new THREE.Sprite(new THREE.SpriteMaterial({ map: new THREE.CanvasTexture(c), depthTest: false, sizeAttenuation: false }));
    s.scale.set(size * 2.6, size * 2.6, 1);
    s.position.copy(pos);
    (parent || this.root).add(s);
    return s;
  }

  label(text, pos, cls = '', onClick = null) {
    const el = document.createElement(onClick ? 'button' : 'div');
    el.className = 'lab ' + cls;
    el.textContent = text;
    if (onClick) el.addEventListener('click', onClick);
    this.labelLayer.appendChild(el);
    const rec = { el, pos, visible: true };
    this.labels.push(rec);
    return rec;
  }

  dropLabel(rec) {
    rec.el.remove();
    this.labels = this.labels.filter((l) => l !== rec);
  }

  project(v) {
    const p = v.clone().project(this.camera);
    return { x: (p.x * 0.5 + 0.5) * innerWidth, y: (-p.y * 0.5 + 0.5) * innerHeight, front: p.z < 1 && p.z > -1 };
  }

  placeLabels() {
    for (const l of this.labels) {
      const v = typeof l.pos === 'function' ? l.pos() : l.pos;
      if (!v || !l.visible) { l.el.style.display = 'none'; continue; }
      const p = this.project(v);
      if (!p.front) { l.el.style.display = 'none'; continue; }
      l.el.style.display = '';
      l.el.style.left = p.x + 'px';
      l.el.style.top = p.y + 'px';
    }
  }

  frameOn(center, radius) {
    this.controls.target.copy(center);
    const dir = new THREE.Vector3(0.12, -0.72, 0.68).normalize();
    this.camera.position.copy(center).addScaledVector(dir, radius * 2.9);
    this.camera.near = radius / 2000;
    this.camera.far = radius * 400;
    this.camera.updateProjectionMatrix();
  }

  png() {
    this.renderer.render(this.scene, this.camera);
    return this.canvas.toDataURL('image/png');
  }
}

export const V = (a) => new THREE.Vector3(a[0], a[1], a[2] ?? 0);
export { THREE };
