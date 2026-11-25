"""Servo placement inside the scaled G1 links.
For every joint the servo axis is the G1 joint axis. The case (+1.2 mm wall) should sit in the
host link and the horn (spline end) inside the driven link, so the horn couples straight into
the G1 part on the other side of the joint. Writes placements.json."""
import json, os, numpy as np
import g1, design as d
from fit import vox_link, frame, WALL, H

S = 0.35
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "placements.json")


def pts_case():
    xs = np.arange(d.X0 - WALL, d.X1 + WALL + 1e-6, 1.0)
    ys = np.arange(-d.HW - WALL, d.HW + WALL + 1e-6, 1.0)
    zs = np.arange(d.FLOOR, d.SV["TOP"] + WALL + 1e-6, 1.0)
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij"); return np.stack([X.ravel(), Y.ravel(), Z.ravel()], 1)


def pts_horn():
    r = np.arange(0, 8.01, 1.0); a = np.radians(np.arange(0, 360, 20)); z = np.arange(-0.8, 2.61, 0.85)
    R, A, Z = np.meshgrid(r, a, z, indexing="ij")
    return np.stack([(R * np.cos(A)).ravel(), (R * np.sin(A)).ravel(), Z.ravel()], 1)


def inside(occ, lo, W):
    idx = np.round((W - lo) / H).astype(int)
    ok = np.all((idx >= 0) & (idx < np.array(occ.shape)), axis=1)
    r = np.zeros(len(W), bool); r[ok] = occ[idx[ok, 0], idx[ok, 1], idx[ok, 2]]
    return r


def place(host, other, p, a, vox):
    C, Hn = pts_case(), pts_horn()
    oh, lh = vox[host]; oo, lo_ = vox[other]
    best = None
    for sgn in (1, -1):
        for th in np.radians(np.arange(0, 360, 10)):
            x, y, z = frame(a * sgn, th); R = np.stack([x, y, z], 1)
            Cl, Hl = C @ R.T, Hn @ R.T
            for t in np.arange(-35, 35.1, 1.0):
                o = p + z * t
                c_out = 1 - inside(oh, lh, Cl + o).mean()
                h_out = 1 - inside(oo, lo_, Hl + o).mean()
                h_inhost = inside(oh, lh, Hl + o).mean()
                cost = c_out + 1.5 * h_out + 0.5 * h_inhost
                if best is None or cost < best[0]:
                    best = (cost, dict(o=o.tolist(), x=x.tolist(), z=z.tolist(), case_out=round(100 * c_out, 1),
                                       horn_out=round(100 * h_out, 1), horn_in_host=round(100 * h_inhost, 1)))
    return best[1]


if __name__ == "__main__":
    J = {k: v for k, v in g1.micro_joints().items() if not k.startswith("right") and "ankle" not in k}
    HOST = {"left_hip_pitch": "pelvis", "waist_yaw": "pelvis", "left_hip_roll": "left_hip_pitch_link",
            "left_hip_yaw": "left_hip_yaw_link", "left_knee": "left_hip_yaw_link",
            "left_shoulder_pitch": "left_shoulder_pitch_link", "left_shoulder_roll": "left_shoulder_pitch_link",
            "left_shoulder_yaw": "left_shoulder_yaw_link", "left_elbow": "left_shoulder_yaw_link"}
    vox = {}
    res = {}
    for k, (par, chi, p, a, lim) in J.items():
        host = HOST[k]; other = chi if host == par else par
        for n in (host, other):
            if n not in vox: vox[n] = vox_link(n, S)
        r = place(host, other, p * S, a, vox)
        r.update(host=host, other=other, axis_pt=(p * S).tolist(), axis=a.tolist(), lim=lim)
        res[k] = r
        print(f"{k:22s} host {host:26s} case out {r['case_out']:5.1f}%  horn out of {other}: {r['horn_out']:5.1f}%  horn in host {r['horn_in_host']:5.1f}%", flush=True)
    json.dump(res, open(OUT, "w"), indent=1)
