"""Export the G1 micro: print STLs, assembly STLs, OpenSCAD assembly, masses, parts list, sim."""
import os, json, shutil
import numpy as np
import g1, design as d
from mesh import write_stl, check_mesh
from g1assemble import JOINTS, LINKS, pieces, S, ROOT, MIR

OUT = os.path.join(ROOT, "STL")
PLA = 1.24e-3
STS, SCS = 20.0, 13.2


def flat_down(V, F):
    """Rotate so the biggest flat face lies on z = 0."""
    T = V[F]; n = np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0]); a = np.linalg.norm(n, axis=1) / 2
    n = n / np.maximum(2 * a[:, None], 1e-12)
    key = np.round(n, 2); off = np.round(np.einsum("ij,ij->i", n, T[:, 0]), 1)
    groups = {}
    for i in np.flatnonzero(a > 0):
        k = (tuple(key[i]), off[i]); groups[k] = groups.get(k, 0) + a[i]
    (nk, _), area = max(groups.items(), key=lambda kv: kv[1])
    nn = np.array(nk, float); nn /= np.linalg.norm(nn)
    z = np.array([0, 0, -1.0])
    v = np.cross(nn, z); s = np.linalg.norm(v); c = nn @ z
    if s < 1e-9:
        R = np.eye(3) if c > 0 else np.diag([1, -1, -1.0])
    else:
        K = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]]); R = np.eye(3) + K + K @ K * ((1 - c) / s ** 2)
    W = V @ R.T; W[:, 2] -= W[:, 2].min()
    W[:, :2] -= (W[:, :2].min(0) + W[:, :2].max(0)) / 2
    return W, area


def components(V, F):
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    e = np.concatenate([F[:, [0, 1]], F[:, [1, 2]]])
    A = coo_matrix((np.ones(len(e)), (e[:, 0], e[:, 1])), shape=(len(V), len(V)))
    _, lab = connected_components(A, directed=False)
    fl = lab[F[:, 0]]; T = V[F]
    vol = np.einsum("ij,ij->i", T[:, 0], np.cross(T[:, 1], T[:, 2])) / 6
    v = np.bincount(fl, weights=vol)
    v = v[np.unique(fl)]
    return int((v > 0).sum() - 1), int((v < 0).sum())


def printed_mass(V, F):
    T = V[F]; area = float(np.linalg.norm(np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0]), axis=1).sum() / 2)
    vol = abs(check_mesh(V, F)["volume_mm3"]); shell = min(vol, area * 0.8)
    return PLA * (shell + 0.15 * max(vol - shell, 0)), vol


def servos_in(link):
    pl = json.load(open(os.path.join(os.path.dirname(__file__), "placements.json")))
    n = sum(1 for k, p in pl.items() if p["host"] == link.replace("right_", "left_"))
    if link in ("pelvis",):
        n = 3
    if link.endswith("knee_link"):
        n += 2
    return n


def main():
    for sub in ("print", "assembly"):
        shutil.rmtree(os.path.join(OUT, sub), ignore_errors=True)
    rows, link_mass, rep = [], {}, {}
    for l in LINKS:
        m_l = 0
        for k, (V, F) in pieces(l).items():
            write_stl(os.path.join(OUT, "assembly", f"{l}_{k}.stl"), V, F, l)
            m, vol = printed_mass(V, F); m_l += m
            W, area = flat_down(V, F)
            write_stl(os.path.join(OUT, "print", f"{l}_{k}.stl"), W, F, l)
            c = check_mesh(V, F)
            c["loose_pieces"], c["internal_voids"] = components(V, F)
            rep[f"{l}/{k}"] = dict(c, flat_face_mm2=round(float(area)), est_mass_g=round(m, 1))
            if not l.startswith("right_"):
                rows.append((l, k, round(vol / 1000, 1), round(m, 1), c["watertight"], c["loose_pieces"]))
        arm = any(t in l for t in ("shoulder_roll", "shoulder_yaw", "elbow"))
        sv = servos_in(l)
        servo_g = sv * (SCS if (arm and "shoulder_pitch" not in l) else STS)
        if l == "left_shoulder_pitch_link" or l == "right_shoulder_pitch_link":
            servo_g = STS + SCS
        if l.endswith("shoulder_yaw_link"):
            servo_g = 2 * SCS
        link_mass[l] = round(m_l + servo_g + (65 if l == "torso" else 0), 1)
    # the dome head's brows are separate parts in the OpenG1 design: export them as they are, scaled
    from mesh import read_stl
    for b in ("brow_left_mm", "brow_right_mm"):
        Vb, Fb = read_stl(os.path.join(g1.REPO, "head", b + ".stl")); Vb = (Vb + np.array([0, 0, g1.PELVIS_Z])) * S
        nm = "head_" + b.replace("_mm", "")
        write_stl(os.path.join(OUT, "assembly", nm + ".stl"), Vb, Fb, nm)
        W, _ = flat_down(Vb, Fb); write_stl(os.path.join(OUT, "print", nm + ".stl"), W, Fb, nm)
    os.makedirs(os.path.join(ROOT, "report"), exist_ok=True)
    json.dump(rep, open(os.path.join(ROOT, "report", "mesh_check.json"), "w"), indent=1)
    json.dump(link_mass, open(os.path.join(ROOT, "report", "link_mass_g.json"), "w"), indent=1)
    plastic = sum(r[3] for r in rows) + sum(r[3] for r in rows if r[0].startswith("left_"))
    lines = ["| Part (left shown; right is the mirror) | Piece | Solid volume (cm3) | Est. printed mass (g) |", "|---|---|---|---|"]
    for r in rows:
        lines.append(f"| {r[0]} | {r[1]} | {r[2]} | {r[3]} |")
    lines.append(f"\nEstimated plastic for one robot: **{plastic:.0f} g** (PLA, 4 walls, 15 % infill). "
                 f"Estimated total with servos and electronics: **{sum(link_mass.values()):.0f} g**.")
    open(os.path.join(ROOT, "parts.md"), "w").write("\n".join(lines) + "\n")
    print("plastic", round(plastic), "total", round(sum(link_mass.values())))
    write_scad()
    import g1sim
    print(g1sim.main({k: v for k, v in link_mass.items()} | {"left_wrist_roll_rubber_hand": 0.0, "right_wrist_roll_rubber_hand": 0.0}
                     | {"torso_link": link_mass["torso"], "pelvis": link_mass["pelvis"],
                        "left_elbow_link": link_mass["left_elbow"], "right_elbow_link": link_mass["right_elbow"]}))


def write_scad():
    f = lambda v: "[" + ",".join(f"{x:.4f}" for x in v) + "]"
    parent = {c: (p, j) for j, (p, c, *_ ) in JOINTS.items()}
    children = {}
    for j, (p, c, *_ ) in JOINTS.items(): children.setdefault(p, []).append((c, j))
    sliders = []
    for side in ("left", "right"):
        sliders.append(f"/* [{side.capitalize()}] */")
        for j, (p, c, pt, ax, lim) in JOINTS.items():
            if j.startswith(side):
                lo, hi = np.degrees(lim)
                sliders.append(f"{j} = 0; // [{lo:.0f}:1:{hi:.0f}]")
    lo, hi = np.degrees(JOINTS["waist_yaw"][4])
    sliders = ["/* [Waist] */", f"waist_yaw = 0; // [{lo:.0f}:1:{hi:.0f}]"] + sliders
    def body(link, ind):
        sp = "  " * ind
        out = []
        for k in pieces(link):
            col = "c_dark" if link in ("head",) or "ankle_pitch" in link else ("c_lid" if k.startswith("lid") else "c_white")
            show = " if (show_lids)" if k.startswith("lid") else ""
            out.append(f'{sp}{show} color({col}) import("STL/assembly/{link}_{k}.stl", convexity=10);')
        for c, j in children.get(link, []):
            _, _, pt, ax, _ = JOINTS[j]
            out.append(f"{sp}translate({f(pt * S)}) rotate(a={j}, v={f(ax)}) translate({f(-pt * S)}) {{")
            out += body(c, ind + 1)
            out.append(f"{sp}}}")
        if link == "torso":
            out += body("head", ind)
            out += [f'{sp}color([0.95,0.75,0.2]) import("STL/assembly/head_brow_{s}.stl");' for s in ("left", "right")]
        return out
    txt = ["// G1 micro: the Unitree G1 parts shrunk to 0.35 and edited for micro servos. Keep next to STL/. F5.",
           "// Joint angles in degrees, same names, zero pose and signs as the G1 (elbow 0 = forearm forward)."] + sliders + [
           "/* [Display] */", "show_lids = true;", "/* [Hidden] */",
           "c_white = [0.92,0.92,0.90]; c_lid = [0.85,0.86,0.88]; c_dark = [0.30,0.31,0.34];"] + body("pelvis", 0)
    open(os.path.join(ROOT, "g1_micro_assembly.scad"), "w").write("\n".join(txt) + "\n")


def fit_test():
    """Servo pocket (body + lid, exactly as in the robot) and a horn socket disc."""
    from sdf import Frame, U, D, I, rboxr, rcyl, halfspace
    from mesh import Grid, mc, decimate
    import g1build as gb
    Sf = Frame((0, 0, 0))
    solid = gb.housing(Sf)
    cut = U(gb.pocket(Sf), gb.lid_screw_cuts(Sf))
    g = Grid((-14, -12, -36), (24, 12, 2), 0.2)
    for tag, keep in (("body", halfspace(Sf, 1, 0.0, keep_above=False)), ("lid", halfspace(Sf, 1, 0.0))):
        V, F = mc(g.eval(I(D(solid, cut), keep)), g); V, F = decimate(V, F, 0.002, 0)
        W, _ = flat_down(V, F); write_stl(os.path.join(OUT, "print", "fit_test", f"servo_pocket_{tag}.stl"), W, F)
    disc = D(rcyl(Sf, 11.0, -1.0, 4.0, 0.6), gb.horn_cuts(Sf))
    g = Grid((-12, -12, -2), (12, 12, 5), 0.15)
    V, F = mc(g.eval(disc), g); V, F = decimate(V, F, 0.002, 0)
    W, _ = flat_down(V, F); write_stl(os.path.join(OUT, "print", "fit_test", "horn_socket.stl"), W, F)


if __name__ == "__main__":
    import sys
    if "fit_test" in sys.argv:
        fit_test()
    else:
        main(); fit_test()
