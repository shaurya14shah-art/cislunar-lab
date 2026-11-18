"""Smoke tests of every server endpoint (the front end's contract)."""
import numpy as np
import pytest
from fastapi.testclient import TestClient
from cislab.server.app import app
from cislab.constants import PRESETS, GM_EARTH, GM_MOON

client = TestClient(app)


def bodies():
    from cislab.nbody.system import two_body
    s = two_body()
    return [dict(name=b.name, gm=b.gm, radius=b.radius, pos=list(b.pos), vel=list(b.vel)) for b in s.bodies]


def test_meta_and_lagrange():
    m = client.get("/api/meta").json()
    assert "Earth" in m["presets"] and m["c"] == 299792458.0
    r = client.get("/api/cr3bp/lagrange").json()
    assert abs(r["points"]["L1"]["position"][0] - 0.836915) < 1e-6
    assert abs(r["points"]["L2"]["lambda"] - 2.1587) < 1e-4
    assert client.get("/api/cr3bp/lagrange?mu=0.9").status_code == 400


def test_hill_family_halo_manifold_poincare():
    h = client.get("/api/cr3bp/hill?C=3.16&n=60&zlevels=3").json()
    assert len(h["grid"]) == 3600 and len(h["contours"]) > 0
    f = client.get("/api/cr3bp/family").json()
    assert len(f["table"]) > 100
    o = client.get("/api/cr3bp/halo?z=0.08").json()
    assert o["stability_index"] > 100 and len(o["points"]) == 240 and o["richardson"]
    m = client.get("/api/cr3bp/manifold?z=0.08&n=12").json()
    assert len(m["tubes"]) == 12 and m["convergence"]["ok"]
    p = client.get("/api/cr3bp/poincare?n=6&t_max=10").json()
    assert len(p["orbits"]) >= 1


def test_propagate_and_live_monte_carlo():
    o = client.get("/api/cr3bp/halo?z=0.08").json()
    r = client.post("/api/cr3bp/propagate", json={"state": o["state0"], "t": o["period"], "compare_pn": True}).json()
    assert r["jacobi_drift"] < 1e-11
    assert r["pn_gap_m"][-1] < 100 and 1e-12 < r["pn_ratio"] < 1e-9
    r = client.post("/api/cr3bp/montecarlo", json={"z": 0.10, "sigma_r_km": 0.001, "n": 40, "n_draw": 4}).json()
    assert r["A"]["captured"] > 0.8 and len(r["tracks"]) == 4
    r = client.post("/api/cr3bp/montecarlo", json={"z": 0.10, "sigma_r_km": 300.0, "n": 40, "n_draw": 2}).json()
    assert r["A"]["captured"] < 0.2


def test_nbody_endpoints():
    b = bodies()
    a = client.post("/api/nbody/analyze", json={"bodies": b, "pair": [0, 1]}).json()
    assert a["regime"][1]["band"] == "newtonian" and abs(a["pair"]["mu"] - 0.01215) < 1e-4
    r = client.post("/api/nbody/run", json={"bodies": b, "t_end": 2.4e6, "compare": True, "frame": "corotating",
                                             "pair": [0, 1], "n_frames": 50}).json()
    assert r["max_energy_drift"] < 1e-9 and len(r["pos"]) == 51 and len(r["divergence_m"]) == 51
    r = client.post("/api/nbody/run", json={"bodies": b, "t_end": 1e5, "method": "yoshida4", "pn": True, "dt": 100})
    assert r.status_code == 400
    f = client.post("/api/nbody/field", json={"bodies": b, "half_width": 5e8, "n": 24, "probe": [1e8, 0, 0],
                                               "n_lines": 6}).json()
    assert len(f["phi"]) == 576 and abs(sum(f["tidal"]["values"])) < 1e-20 and len(f["lines"]) == 6
    sun = PRESETS["Sun"]
    r = client.post("/api/nbody/rays", json={"bodies": [dict(name="Sun", gm=sun.gm, radius=sun.radius, pos=[0, 0, 0], vel=[0, 0, 0])],
                                              "origin": [-5e11, 0, 0], "direction": [1, 0, 0], "spread": 3e9,
                                              "n_rays": 5, "length": 1e12}).json()
    assert len(r["rays"]) == 5 and r["rays"][2]["blocked"] and 0.3 < r["rays"][0]["deflection_arcsec"] < 0.5


def test_schwarzschild_endpoints_and_checks():
    o = client.get("/api/schwarzschild/orbit?rp=8&ra=20").json()
    assert abs(o["isco"] - 6) < 1e-6 and o["precession_deg_per_orbit"] > o["weak_field_deg_per_orbit"]
    r = client.get("/api/schwarzschild/rays?n=5").json()
    assert r["rays"][0]["captured"] and not r["rays"][-1]["captured"]
    c = client.get("/api/checks").json()
    assert abs(c["gps_us_per_day"] - 38.6) < 0.2 and abs(c["mercury_arcsec_per_century"] - 42.98) < 0.1


def test_websocket_live_session():
    with client.websocket_connect("/ws/nbody") as ws:
        import json
        ws.send_text(json.dumps({"cmd": "set", "bodies": bodies()}))
        assert json.loads(ws.receive_text())["set"] == 2
        ws.send_text(json.dumps({"cmd": "advance", "t": 1e5, "n_frames": 10}))
        a = json.loads(ws.receive_text())
        assert a["ok"] and len(a["pos"]) == 11 and abs(a["energy_drift"]) < 1e-9
        ws.send_text(json.dumps({"cmd": "advance", "t": -1e5, "n_frames": 10}))     # reverse
        b = json.loads(ws.receive_text())
        assert np.allclose(b["pos"][-1], a["pos"][0], rtol=1e-8, atol=1.0)
