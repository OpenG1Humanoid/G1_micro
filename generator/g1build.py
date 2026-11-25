"""Build the G1 micro from the real G1 parts, shrunk to S = 0.35 (45 cm).

Every link starts as the scaled G1 mesh (solid). Edits only:
  * host link: + housing around the servo (a bulge only where the servo pokes out of the G1 shape),
               - servo pocket, - clearance for the driven link's horn boss
  * driven link: + horn boss (only where the horn region is outside the G1 shape),
                 - horn pocket (cross horn, arms snipped to 13 mm), - centre-screw access hole,
                 - clearance for the host's servo housing over the joint range
  * a screw-on lid over every servo (the G1 shell cut along the servo's mid plane) so the servo can go in
  * ankles like the real G1: two servos in the calf with cranks, push rods to the foot heel,
    the ankle cross is passive on two M2 pins.
"""
import os, sys, json, time
import numpy as np
from scipy import ndimage
import g1, design as d
from sdf import *
from mesh import Grid, mc, decimate, check_mesh, write_stl, undercut_report
from fit import vox_link

S = 0.35
H = 0.25
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
CACHE = os.path.join(ROOT, "g1_cache")
os.makedirs(CACHE, exist_ok=True)
PL = json.load(open(os.path.join(os.path.dirname(__file__), "placements.json")))
WALL = 1.2

# micro joint ranges (deg): G1 limits, clipped to what the servos and the shapes allow
RANGE = {"hip_pitch": (-130, 90), "hip_roll": (-25, 90), "hip_yaw": (-90, 90), "knee": (-5, 140),
         "waist_yaw": (-90, 90), "shoulder_pitch": (-175, 120), "shoulder_roll": (-20, 120),
         "shoulder_yaw": (-120, 120), "elbow": (-60, 120), "ankle_pitch": (-45, 30), "ankle_roll": (-15, 15)}


def jbase(k):
    return k.replace("left_", "").replace("right_", "")


def servo_frames():
    out = {}
    for k, p in PL.items():
        fr = Frame(np.array(p["o"]), x=np.array(p["x"]), z=np.array(p["z"]))
        out[k] = dict(S=fr, host=p["host"], other=p["other"], axis_pt=np.array(p["axis_pt"]),
                      axis=np.array(p["axis"]), rng=RANGE[jbase(k)])
    return out


SERVOS = servo_frames()
MIRV = np.array([1, -1, 1.0])
for k in [k for k in SERVOS if k.startswith("left_")]:
    s = SERVOS[k]; Sf = s["S"]
    SERVOS[k.replace("left_", "right_")] = dict(
        S=Frame(Sf.o * MIRV, x=Sf.R[:, 0] * MIRV, z=Sf.R[:, 2] * MIRV),
        host=s["host"].replace("left_", "right_"), other=s["other"].replace("left_", "right_"),
        axis_pt=s["axis_pt"] * MIRV, axis=s["axis"] * MIRV, rng=(-s["rng"][1], -s["rng"][0]))


# ---------------------------------------------------------------- servo-side features
def housing(Sf):
    return U(rboxr(Sf, (d.X0 - WALL, -d.HW - WALL, d.FLOOR), (d.X1 + WALL, d.HW + WALL, d.SV["TOP"] + WALL), 1.2),
             *[cyl(Frame(Sf.world(p), x=Sf.R[:, 0], z=Sf.R[:, 1]), 2.9, -d.HW - WALL, d.HW + WALL) for p in lid_screws()])


def lid_screws():
    zm = (d.ZB + d.SV["TOP"]) / 2
    return [(d.X0 - 2.2, 0, zm), (d.X1 + 2.2, 0, zm)]


def pocket(Sf):
    e = (d.SV["EARLEN"] - d.SV["L"]) / 2
    return U(boxr(Sf, (d.X0, -d.HW, d.ZB), (d.X1, d.HW, d.SV["TOP"] + 0.01)),
             boxr(Sf, (d.X0 - e - 0.3, -d.HW, d.SV["EAR"][0] - 0.2), (d.X1 + e + 0.3, d.HW, d.SV["EAR"][1] + 0.2)),
             cyl(Sf, 3.6, d.ZB - 0.8, d.ZB + 0.1),
             boxr(Sf, (d.X1 - 0.6, -2.2, d.ZB), (d.X1 + WALL + 3.0, 2.2, d.SV["EAR"][0])),     # cable exit
             cyl(Sf, 6.5, d.SV["TOP"] - 0.1, 1.0),
             boxr(Sf, (0, -4.0, d.SV["TOP"] - 0.1), (8.5, 4.0, d.SV["TOP"] + WALL + 0.6)))


def lid_region(Sf):
    return boxr(Sf, (d.X0 - 5.2, 0.0, d.FLOOR - 4), (d.X1 + 5.2, 80, -3.4))


def lid_screw_cuts(Sf):
    c = []
    for p in lid_screws():
        Fr = Frame(Sf.world(p), x=Sf.R[:, 0], z=Sf.R[:, 1])    # +z = lid side
        c += [cyl(Fr, d.M2_CLR, -0.01, 3.0), cyl(Fr, d.M2_HEAD, 3.0, 60), cyl(Fr, d.M2_PILOT, -7, 0.01)]
    return U(*c)


# ---------------------------------------------------------------- driven-side features
def horn_boss(Sf):
    """Material for the horn socket, with a stem running into the driven G1 part."""
    return U(rcyl(Sf, 9.3, -1.0, 3.6, 0.8), rcyl(Sf, 8.0, 2.0, 26.0, 1.5))


def horn_cuts(Sf):
    return U(d.horn_cut(Sf, -0.8, 1.4), cyl(Sf, 4.1, -3.6, 1.4), cyl(Sf, 2.45, 1.4, 60),
             *[cyl(Sf, 0.7, -1, 4.6, c=(5.0 * np.cos(a), 5.0 * np.sin(a))) for a in np.radians([45, 225])])


def other_clear(Sf):
    return cyl(Sf, 9.9, d.FLOOR - 2, -1.1)


def host_clear(Sf):
    return cyl(Sf, 10.0, -3.4, 20)


# ---------------------------------------------------------------- ankle (G1 style linkage)
def ankle_layout():
    """Two crank servos in the calf (axis along x, horn behind the calf), rods to the heel."""
    J = g1.micro_joints()
    _, _, pp, _, _ = J["left_ankle_pitch"]; _, _, pr, _, _ = J["left_ankle_roll"]
    pp = pp * S; pr = pr * S
    occ, lo = vox_link("left_knee_link", S)
    best = []
    E = []
    xs = np.arange(d.X0 - WALL, d.X1 + WALL + 1e-6, 1.0); ys = np.arange(-d.HW - WALL, d.HW + WALL + 1e-6, 1.0)
    zs = np.arange(d.FLOOR, d.SV["TOP"] + WALL + 1e-6, 1.0)
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij"); E = np.stack([X.ravel(), Y.ravel(), Z.ravel()], 1)
    def ins(W):
        idx = np.round((W - lo) / 0.5).astype(int); ok = np.all((idx >= 0) & (idx < occ.shape), 1)
        r = np.zeros(len(W), bool); r[ok] = occ[idx[ok, 0], idx[ok, 1], idx[ok, 2]]; return r.mean()
    cands = []
    kz = SERVOS["left_knee"]["axis_pt"][2]
    for up in (1, -1):                 # case long side up or down
        Fr0 = (np.array([0, 0, up], float), np.array([-1, 0, 0], float))   # x_s, z_s (spline backwards)
        for zc in np.arange(pp[2] + 30, kz - 15, 1.0):
            for xo in np.arange(-6, 14, 1.0):
                for yo in np.arange(-6, 6.1, 1.0):
                    o = np.array([xo, pp[1] + yo, zc])
                    Fr = Frame(o, x=Fr0[0], z=Fr0[1])
                    cands.append((ins(E @ Fr.R.T + o), up, zc, xo, yo))
    cands.sort(key=lambda c: -c[0])
    a = cands[0]
    # second servo: same x, no overlap (case 23.2 + walls along z)
    b = max((c for c in cands if abs(c[2] - a[2]) >= 26.5 and c[1] == a[1]), key=lambda c: c[0])
    out = []
    for c, side in ((a, +1), (b, -1)):
        o = np.array([c[3], pp[1] + c[4], c[2]])
        out.append(dict(S=Frame(o, x=np.array([0, 0, c[1]], float), z=np.array([-1, 0, 0], float)), fill=c[0], side=side))
    print("ankle servos in calf, fraction inside G1 calf:", [round(float(c[0]), 2) for c in (a, b)])
    return out, pp, pr


ANK, PP, PR = ankle_layout()
CRANK = 9.0     # horn arm hole radius used as crank


def crank_tip(a):
    return a["S"].o + a["S"].R[:, 2] * 2.0 + np.array([0, a["side"] * CRANK, 0])


def heel_stud(a):
    t = crank_tip(a)
    return np.array([t[0], t[1], PR[2] + 1.0])


def g1assemble_parent(link):
    """Parent micro link of a micro link (left or right)."""
    for j, (p, c, *_ ) in g1.micro_joints().items():
        if c == link: return p
    return None


# ---------------------------------------------------------------- cached-shape obstacles
_obs = {}


def cached_shape(name, clearance=0.5):
    """SDF of an already built link (from g1_cache), right side by mirroring."""
    from scipy.ndimage import map_coordinates
    src = name.replace("right_", "left_")
    if src not in _obs:
        z = np.load(os.path.join(CACHE, src + "_sdf.npz")); _obs[src] = (z["vol"].astype(np.float32), z["lo"], float(z["h"]))
    vol, lo, h = _obs[src]
    mir = name.startswith("right_")
    def f(P):
        Q = P.T.astype(np.float64)
        if mir: Q = Q * np.array([1, -1, 1.0])
        return (map_coordinates(vol, ((Q - lo) / h).T, order=1, mode="constant", cval=8.0) - clearance).astype(np.float32)
    return f


def near_axis(fn, pt, ax, r):
    a = np.asarray(ax, float); a /= np.linalg.norm(a)
    ref = np.array([0.3, 0.5, 0.8]); x = np.cross(a, ref); x /= np.linalg.norm(x)
    return I(fn, cyl(Frame(pt, x=x, z=a), r, -60, 60))


# ---------------------------------------------------------------- per-link edit list
def link_edits(name):
    adds, cuts, lids, sweeps = [], [], [], []
    for k, s in SERVOS.items():
        Sf = s["S"]
        if s["host"] == name:
            adds.append(housing(Sf)); cuts += [pocket(Sf), host_clear(Sf), lid_screw_cuts(Sf)]; lids.append(Sf)
        if s["other"] == name:
            adds.append(horn_boss(Sf)); cuts += [horn_cuts(Sf), other_clear(Sf)]
            # rotate the host relative to this link: host is the child -> +angle, host is the parent -> -angle
            host_is_child = g1assemble_parent(s["host"]) == name
            rng = s["rng"] if not host_is_child else (-s["rng"][1], -s["rng"][0])
            sweeps.append((near_axis(cached_shape(s["host"]), s["axis_pt"], s["axis"], 24.0), s["axis_pt"], s["axis"], rng))
    if name == "left_knee_link":
        for a in ANK:
            Sf = a["S"]
            adds.append(housing(Sf)); cuts += [pocket(Sf), lid_screw_cuts(Sf), cyl(Sf, 10.5, -3.4, 30)]; lids.append(Sf)
        cuts.append(cyl(Frame(PP, x=(1, 0, 0), z=(0, 1, 0)), d.M2_PILOT, -30, 30))          # ankle pitch pin
    if name == "head":
        cuts.append(cached_shape("torso", 0.3))
    if name == "torso":
        # the upper arm (shoulder roll + yaw links) swings past the torso side: clear its path
        J = g1.micro_joints()
        for side in ("left", "right"):
            pt, ax = J[f"{side}_shoulder_pitch"][2] * S, J[f"{side}_shoulder_pitch"][3]
            for l in (f"{side}_shoulder_roll_link", f"{side}_shoulder_yaw_link"):
                sweeps.append((cached_shape(l, 0.6), pt, ax, (-60, 120)))   # arm pitch -120 (forward) .. +60 (back)
    if name == "left_ankle_pitch_link":
        J = g1.micro_joints()
        for k, host in (("left_ankle_pitch", "left_knee_link"), ("left_ankle_roll", "left_ankle_roll_link")):
            pt, ax = J[k][2] * S, J[k][3]
            sweeps.append((near_axis(cached_shape(host), pt, ax, 15.0), pt, ax, RANGE[jbase(k)]))
        cuts.append(cyl(Frame(PP, x=(1, 0, 0), z=(0, 1, 0)), 1.15, -30, 30))
        cuts.append(cyl(Frame(PR, x=(0, 1, 0), z=(1, 0, 0)), 1.15, -30, 30))
    if name == "left_ankle_roll_link":
        cuts.append(cyl(Frame(PR, x=(0, 1, 0), z=(1, 0, 0)), d.M2_PILOT, -30, 30))          # ankle roll pin
        for a in ANK:
            hs = heel_stud(a)
            adds.append(rcyl(Frame(hs, x=(1, 0, 0), z=(0, 0, 1)), 3.2, -6, 0.5, 0.8))
            cuts.append(cyl(Frame(hs, x=(1, 0, 0), z=(0, 0, 1)), d.M2_PILOT, -6, 2))
    return adds, cuts, lids, sweeps


def g1_sdf(name, grid):
    """Signed distance of the scaled G1 link: sign from a solid voxelisation, exact distance to the
    G1 surface in a band near it (no voxel steps), loose fragments dropped."""
    from mesh import read_stl
    meshes = []
    if name == "head":
        for f in ["frame_mm", "panel_0_mm", "panel_1_mm", "panel_2_mm", "panel_3_mm", "panel_4_mm", "brow_left_mm", "brow_right_mm"]:
            Vp, Fp = read_stl(os.path.join(g1.REPO, "head", f + ".stl")); meshes.append(((Vp + np.array([0, 0, g1.PELVIS_Z])) * S, Fp))
    else:
        meshes.append(g1.link_mesh(name, S))
    occ = np.zeros(tuple(grid.n), bool)
    for V, F in meshes:
        occ |= g1.voxelize(V, F, grid.lo, grid.h, grid.n)
    if name == "head":
        occ = ndimage.binary_closing(occ, iterations=2)
    occ = ndimage.binary_fill_holes(occ)
    lab, n = ndimage.label(occ)
    if n > 1:
        sizes = ndimage.sum(occ, lab, range(1, n + 1)); keep = np.flatnonzero(sizes >= 0.03 * sizes.max()) + 1
        occ = np.isin(lab, keep)
    sd = (ndimage.distance_transform_edt(~occ) - ndimage.distance_transform_edt(occ)) * grid.h
    band = np.abs(sd) < 1.6
    idx = np.argwhere(band)
    P = grid.lo + idx * grid.h
    V = np.concatenate([m[0] for m in meshes]); off = 0; Fs = []
    for m in meshes: Fs.append(m[1] + off); off += len(m[0])
    dist = g1.mesh_distance(V, np.concatenate(Fs), P, 2.0)
    sgn = np.where(occ[band], -1.0, 1.0)
    sd = sd.astype(np.float32)
    sd[band] = sgn * dist
    return sd


def bbox_of(name, adds):
    if name == "head":
        V, _ = g1.head_mesh(S)
    else:
        V, _ = g1.link_mesh(name, S)
    lo, hi = V.min(0), V.max(0)
    return lo - 4, hi + 4


def build(name):
    t0 = time.time()
    adds, cuts, lids, sweeps = link_edits(name)
    lo, hi = bbox_of(name, adds)
    for a in adds:   # grow box to include bulges (coarse probe)
        pass
    # enlarge box by servo envelopes hosted here
    for Sf in lids:
        c = np.array([Sf.world(p) for p in [(x, y, z) for x in (d.X0 - 5, d.X1 + 5) for y in (-9, 9) for z in (d.FLOOR - 2, 4)]])
        lo = np.minimum(lo, c.min(0)); hi = np.maximum(hi, c.max(0))
    for k, s in SERVOS.items():
        if s["other"] == name:
            c = np.array([s["S"].world((x, y, z)) for x in (-10, 10) for y in (-10, 10) for z in (-1.5, 4)])
            lo = np.minimum(lo, c.min(0)); hi = np.maximum(hi, c.max(0))
    if name == "left_ankle_roll_link":
        for a in ANK:
            hs = heel_stud(a); lo = np.minimum(lo, hs - 5); hi = np.maximum(hi, hs + 5)
    h = H if np.prod(hi - lo) < 2.5e5 else 0.3
    grid = Grid(lo, hi, h)
    base = g1_sdf(name, grid)
    vol = base
    if adds:
        A = grid.eval(U(*adds))
        k = 1.5
        hh = np.clip(0.5 + 0.5 * (A - vol) / k, 0, 1)
        vol = A * (1 - hh) + vol * hh - k * hh * (1 - hh)
    if cuts:
        vol = np.maximum(vol, -grid.eval(U(*cuts)))
    for (hous, pt, ax, rng) in sweeps:
        lo_, hi_ = rng
        sw = None
        for ang in np.unique(np.concatenate([np.arange(lo_, hi_, 6.0), [hi_]])):
            v = grid.eval(rotated(hous, pt, ax, -ang))
            sw = v if sw is None else np.minimum(sw, v)
        vol = np.maximum(vol, -sw)
    # pieces: sequential lids, then body
    pieces = {}
    taken = None
    for i, Sf in enumerate(lids):
        R = grid.eval(lid_region(Sf))
        if taken is not None: R = np.maximum(R, -taken)
        pieces[f"lid{i + 1}"] = np.maximum(vol, R)
        taken = R if taken is None else np.minimum(taken, R)
    pieces["body"] = vol if taken is None else np.maximum(vol, -taken)
    res = {}
    for tag, pv in pieces.items():
        occp = pv < 0
        lab, n = ndimage.label(occp)
        if n > 1:
            sizes = ndimage.sum(occp, lab, range(1, n + 1))
            drop = np.isin(lab, np.flatnonzero(sizes < 0.03 * sizes.max()) + 1)
            frac = drop.sum() / max(occp.sum(), 1)
            if frac > 0.002: print(f"   WARNING {name}/{tag}: dropped loose pieces = {100*frac:.2f}% of volume", flush=True)
            pv = np.where(drop, np.float32(grid.h), pv)
        V, F = mc(pv, grid)
        if len(F) == 0: continue
        # drop tiny islands
        V0, F0 = V, F
        V, F = decimate(V0, F0, 0.002, 0)
        chk = check_mesh(V, F)
        if not chk["watertight"]:
            V, F = decimate(V0, F0, 0.0004, 0); chk = check_mesh(V, F)
        if not chk["watertight"]:
            for off in (0.03, -0.03, 0.08, -0.08):
                V1, F1 = mc(pv + np.float32(off), grid)
                V, F = decimate(V1, F1, 0.002, 0); chk = check_mesh(V, F)
                if chk["watertight"]: break
        res[tag] = (V, F, chk)
    np.savez_compressed(os.path.join(CACHE, name + ".npz"), **{f"{k}_V": v[0] for k, v in res.items()},
                        **{f"{k}_F": v[1] for k, v in res.items()})
    np.savez_compressed(os.path.join(CACHE, name + "_sdf.npz"), vol=np.clip(vol, -8, 8).astype(np.float16), lo=grid.lo, h=grid.h)
    print(f"{name:26s} {time.time() - t0:5.1f}s  " + "  ".join(f"{k}: {v[2]['triangles']} tri wt={v[2]['watertight']} shells={v[2]['shells']} vol={v[2]['volume_mm3']}" for k, v in res.items()), flush=True)
    return {k: v[2] for k, v in res.items()}


LEFT = ["pelvis", "torso", "head", "left_hip_pitch_link", "left_hip_roll_link", "left_hip_yaw_link", "left_knee_link",
        "left_ankle_pitch_link", "left_ankle_roll_link", "left_shoulder_pitch_link", "left_shoulder_roll_link",
        "left_shoulder_yaw_link", "left_elbow"]

if __name__ == "__main__":
    names = sys.argv[1:] or LEFT
    rep_p = os.path.join(CACHE, "report.json")
    rep = json.load(open(rep_p)) if os.path.exists(rep_p) else {}
    for n in names:
        rep[n] = build(n)
        json.dump(rep, open(rep_p, "w"), indent=1)
