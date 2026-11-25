import os, numpy as np
from PIL import Image, ImageDraw
from g1assemble import *
import g1build as gb

R = os.path.join(ROOT, "renders"); os.makedirs(R, exist_ok=True)


def shot(ms, eye, lab, W=520, H=980, sc=2.05, target=None):
    allV = np.concatenate([m[0] for m in ms]); c = (allV.min(0) + allV.max(0)) / 2 if target is None else np.asarray(target)
    render(ms, "/tmp/_s.png", eye=c + np.array(eye), target=c, scale=sc, W=W, H=H)
    im = Image.open("/tmp/_s.png").convert("RGB"); ImageDraw.Draw(im).text((10, 10), lab, fill=(0, 0, 0)); return im


def row(ims, path):
    out = Image.new("RGB", (sum(i.width for i in ims), max(i.height for i in ims)), "white"); x = 0
    for i in ims: out.paste(i, (x, 0)); x += i.width
    out.save(path)


def rods(q):
    """Ankle push rods (thin boxes) between crank tips and heel studs, for the renders."""
    T = fk(q); ms = []
    for side, mir in (("left", 1), ("right", -1)):
        for a in gb.ANK:
            t = gb.crank_tip(a); h = gb.heel_stud(a)
            if mir == -1: t = t * MIR; h = h * MIR
            t = apply(T[f"{side}_knee_link"], t[None])[0]; h = apply(T[f"{side}_ankle_roll_link"], h[None])[0]
            v = h - t; L = np.linalg.norm(v); z = v / L; x = np.cross(z, [0.3, 0.7, 0.1]); x /= np.linalg.norm(x); y = np.cross(z, x)
            cs = np.array([[sx, sy, sz] for sx in (-0.9, 0.9) for sy in (-0.9, 0.9) for sz in (0, L)])
            V = t + cs[:, 0:1] * x + cs[:, 1:2] * y + cs[:, 2:3] * z
            F = np.array([[0, 1, 3], [0, 3, 2], [4, 6, 7], [4, 7, 5], [0, 4, 5], [0, 5, 1], [2, 3, 7], [2, 7, 6], [0, 2, 6], [0, 6, 4], [1, 5, 7], [1, 7, 3]])
            ms.append((V, F, (60, 60, 60)))
    return ms


if __name__ == "__main__":
    ref, mic = g1_reference(), scene({}) + rods({})
    row([shot(ref, (600, 0, 0), "G1 parts x0.35, untouched"), shot(mic, (600, 0, 0), "G1 micro"),
         shot(ref, (380, 420, 150), "G1 parts x0.35, untouched"), shot(mic, (380, 420, 150), "G1 micro")],
        os.path.join(R, "01_g1_micro_vs_g1.png"))
    B = lambda j, l, r=None: {f"left_{j}": l, f"right_{j}": (l if r is None else r)}
    P = [("wave", {"left_shoulder_roll": 100, "left_elbow": -40, "left_shoulder_pitch": -30, **B("elbow", 60)} | {"left_elbow": -40}),
         ("T-pose", {**B("shoulder_roll", 90, -90), **B("elbow", 90)}),
         ("squat", {**B("hip_pitch", -90), **B("knee", 115), **B("ankle_pitch", -25), **B("shoulder_pitch", -60), **B("elbow", 0)}),
         ("step", {"left_hip_pitch": -40, "left_knee": 60, "left_ankle_pitch": -20, "right_hip_pitch": 15, "right_ankle_pitch": -15,
                   "left_shoulder_pitch": 25, "right_shoulder_pitch": -25, **B("elbow", 50)})]
    row([shot(scene(q) + rods(q), (450, 330, 120), n, W=520, H=900, sc=1.75) for n, q in P], os.path.join(R, "02_poses.png"))
    # exploded lids (left leg + arm)
    T = fk({}); ms = []
    for l in LINKS:
        for k, (V, F) in pieces(l).items():
            off = np.zeros(3)
            if k.startswith("lid") and not l.startswith("right_"):
                srcl = l
                Sfs = [s["S"] for s in gb.SERVOS.values() if s["host"] == l] + ([a["S"] for a in gb.ANK] if l == "left_knee_link" else [])
                i = int(k[3:]) - 1
                if i < len(Sfs): off = Sfs[i].R[:, 1] * 14
            ms.append((apply(T[l], V) + off, F, (240, 150, 60) if k.startswith("lid") else colour(l, k)))
    row([shot(ms, (380, 420, 150), "lids pulled off (orange): every servo drops in under its lid", W=760, H=980)], os.path.join(R, "03_lids.png"))
    # ankle close-up
    cl = [m for m in scene({}, exclude=[l for l in LINKS if not any(t in l for t in ("knee", "ankle"))]) + rods({})]
    row([shot(cl, (260, 300, 60), "G1-style ankle: two calf servos, cranks, push rods, passive cross", W=700, H=760, sc=3.2)],
        os.path.join(R, "04_ankle.png"))
