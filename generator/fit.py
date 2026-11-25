"""For each G1 joint, search where a micro servo (case + 1.2 mm wall) fits best inside the
scaled G1 links: host link, spline direction, axial position, rotation about the axis."""
import sys, json, numpy as np
import g1, design as d

WALL = 1.2
H = 0.5


def vox_link(name, s, pad=6):
    V, F = g1.link_mesh(name, s)
    lo = V.min(0) - pad; hi = V.max(0) + pad
    n = np.ceil((hi - lo) / H).astype(int) + 1
    return g1.voxelize(V, F, lo, H, n), lo


def envelope_pts():
    xs = np.arange(d.X0 - WALL, d.X1 + WALL + 1e-6, 1.0)
    ys = np.arange(-d.HW - WALL, d.HW + WALL + 1e-6, 1.0)
    zs = np.arange(d.FLOOR, d.SV["TOP"] + WALL + 1e-6, 1.0)
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    return np.stack([X.ravel(), Y.ravel(), Z.ravel()], 1)


def frame(a, theta, ref=np.array([0.3, 0.2, 0.93])):
    z = a / np.linalg.norm(a)
    x0 = np.cross(ref, z); x0 /= np.linalg.norm(x0); y0 = np.cross(z, x0)
    x = np.cos(theta) * x0 + np.sin(theta) * y0
    return x, np.cross(z, x), z


def best_fit(occ, lo, p, a, rng):
    E = envelope_pts()
    best = (1e9, None)
    for sgn in (1, -1):
        for th in np.radians(np.arange(0, 360, 15)):
            x, y, z = frame(a * sgn, th)
            R = np.stack([x, y, z], 1)
            local = E @ R.T
            for t in rng:
                W = local + p + z * t
                idx = np.round((W - lo) / H).astype(int)
                ok = np.all((idx >= 0) & (idx < np.array(occ.shape)), axis=1)
                ins = np.zeros(len(W), bool)
                ins[ok] = occ[idx[ok, 0], idx[ok, 1], idx[ok, 2]]
                out = 1 - ins.mean()
                if out < best[0]:
                    best = (out, dict(sign=sgn, theta=float(np.degrees(th)), t=float(t)))
    return best


if __name__ == "__main__":
    scales = [float(x) for x in sys.argv[1:]] or [0.3, 0.35, 0.4, 0.45]
    J = {k: v for k, v in g1.micro_joints().items() if not k.startswith("right")}
    res = {}
    for s in scales:
        vox = {}
        res[s] = {}
        for k, (par, chi, p, a, lim) in J.items():
            row = {}
            for host in (par, chi):
                if host not in vox: vox[host] = vox_link(host, s)
                occ, lo = vox[host]
                out, cfg = best_fit(occ, lo, p * s, a, np.arange(-30, 30.1, 1.5))
                row[host] = (round(100 * out, 1), cfg)
            res[s][k] = row
            bh = min(row, key=lambda h: row[h][0])
            print(f"s={s:.2f} {k:20s} best host {bh:26s} outside {row[bh][0]:5.1f}%   (other {min(v[0] for h, v in row.items() if h != bh):5.1f}%)", flush=True)
    json.dump({str(k): v for k, v in res.items()}, open("/tmp/fit.json", "w"), indent=1)
