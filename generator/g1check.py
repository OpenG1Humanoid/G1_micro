"""Pose interference check for the G1 micro (joint angles in degrees, G1 names and signs)."""
import os, json, numpy as np
from scipy.ndimage import map_coordinates
from g1assemble import *

THRESH = -0.4
B = lambda j, l, r=None: {f"left_{j}": l, f"right_{j}": (l if r is None else r)}
POSES = {
    "zero (G1 zero, elbows 90)": {},
    "arms down (shoulder roll 8)": {**B("elbow", 90), **B("shoulder_roll", 8, -8)},
    "squat": {**B("hip_pitch", -90), **B("knee", 115), **B("ankle_pitch", -25), **B("elbow", 60), **B("shoulder_roll", 15, -15)},
    "T-pose": {**B("shoulder_roll", 90, -90), **B("elbow", 90)},
    "arms up (roll 50)": {**B("shoulder_pitch", -160), **B("shoulder_roll", 50, -50), **B("elbow", 60)},
    "arms forward": {**B("shoulder_pitch", -80), **B("elbow", 90), **B("shoulder_roll", 8, -8)},
    "walk sway": {**B("hip_roll", 10), **B("ankle_roll", -10), **B("elbow", 90), **B("shoulder_roll", 12, -12)},
    "step": {"left_hip_pitch": -40, "left_knee": 60, "left_ankle_pitch": -20, "right_hip_pitch": 15, "right_ankle_pitch": -15,
             **B("elbow", 70), **B("shoulder_roll", 12, -12), "left_shoulder_pitch": 25, "right_shoulder_pitch": -25},
    "turn": {**B("hip_yaw", 30), "waist_yaw": 45, **B("elbow", 60)},
    "push-up arms": {**B("shoulder_pitch", -90), **B("elbow", 90), **B("shoulder_roll", 8, -8)},
    # G1 limits kept for reference: these also collide on the untouched G1 parts
    "REF arms straight down, no roll": B("elbow", 90),
}
_sdf = {}


def sdf_sample(link, P):
    src = link.replace("right_", "left_")
    if src not in _sdf:
        z = np.load(os.path.join(CACHE, src + "_sdf.npz")); _sdf[src] = (z["vol"].astype(np.float32), z["lo"], float(z["h"]))
    vol, lo, h = _sdf[src]
    if link.startswith("right_"): P = P * MIR
    return map_coordinates(vol, ((P - lo) / h).T, order=1, mode="constant", cval=8.0)


def check(q):
    T = fk(q); hits = []
    W = {l: apply(T[l], np.concatenate([v[0] for v in pieces(l).values()])) for l in LINKS}
    for a in LINKS:
        Ti = np.linalg.inv(T[a])
        vol_lo = None
        for b in LINKS:
            if a == b: continue
            P = apply(Ti, W[b])
            dv = sdf_sample(a if a != "head" else "head", P)
            n = int((dv < THRESH).sum())
            if n > 3: hits.append((b, a, n, round(float(dv.min()), 2)))
    return hits


if __name__ == "__main__":
    out = {}
    for k, q in POSES.items():
        h = check(q); out[k] = h
        print(f"{k:28s}", "OK" if not h else h, flush=True)
    os.makedirs(os.path.join(ROOT, "report"), exist_ok=True)
    json.dump(out, open(os.path.join(ROOT, "report", "pose_check.json"), "w"), indent=1)
