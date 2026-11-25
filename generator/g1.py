"""The real G1 geometry: URDF kinematics + link meshes (Unitree g1_23dof_rev_1_0), with the
user's OpenG1 dome head. Everything returned in mm, zero pose, robot frame = G1 pelvis frame
shifted so the feet stand on z = 0."""
import os, subprocess, struct, tempfile
import numpy as np
import xml.etree.ElementTree as ET
from mesh import read_stl

G1 = os.environ.get("G1_DESCRIPTION", "/home/claude/src/unitree_ros/robots/g1_description")   # unitree_ros/robots/g1_description
REPO = os.environ.get("OPENG1_STL", "/home/claude/src/humanoid/STL")                      # OpenG1Humanoid/humanoid/STL
URDF = os.path.join(G1, "g1_23dof_rev_1_0.urdf")
HERE = os.path.dirname(os.path.abspath(__file__))


def rpy_mat(r, p, y):
    cr, sr, cp, sp, cy, sy = np.cos(r), np.sin(r), np.cos(p), np.sin(p), np.cos(y), np.sin(y)
    return np.array([[cy * cp, cy * sp * sr - sy * cr, cy * sp * cr + sy * sr],
                     [sy * cp, sy * sp * sr + cy * cr, sy * sp * cr - cy * sr],
                     [-sp, cp * sr, cp * cr]])


def tf(xyz, rpy):
    T = np.eye(4); T[:3, :3] = rpy_mat(*rpy); T[:3, 3] = xyz; return T


def parse():
    r = ET.parse(URDF).getroot()
    links, joints = {}, []
    for l in r.findall("link"):
        vis = []
        for v in l.findall("visual"):
            m = v.find("geometry/mesh")
            if m is None: continue
            o = v.find("origin")
            xyz = [float(x) for x in (o.get("xyz", "0 0 0") if o is not None else "0 0 0").split()]
            rpy = [float(x) for x in (o.get("rpy", "0 0 0") if o is not None else "0 0 0").split()]
            vis.append((m.get("filename").split("/")[-1], tf(xyz, rpy)))
        links[l.get("name")] = vis
    for j in r.findall("joint"):
        o = j.find("origin"); a = j.find("axis"); lim = j.find("limit")
        if o is None: o = ET.Element("origin")
        joints.append(dict(name=j.get("name"), type=j.get("type"), parent=j.find("parent").get("link"),
                           child=j.find("child").get("link"),
                           T=tf([float(x) for x in o.get("xyz", "0 0 0").split()], [float(x) for x in o.get("rpy", "0 0 0").split()]),
                           axis=np.array([float(x) for x in a.get("xyz").split()]) if a is not None else None,
                           lim=(float(lim.get("lower")), float(lim.get("upper"))) if lim is not None and lim.get("lower") else None))
    return links, joints


LINKS, JOINTS = parse()
PELVIS_Z = 793.0     # G1 pelvis frame height at zero pose (feet on the floor)

# micro robot: which G1 links are fused into one rigid micro link
FUSE = {"pelvis": ["pelvis", "pelvis_contour_link"], "torso": ["torso_link", "logo_link"]}
for s in ("left", "right"):
    FUSE[f"{s}_elbow"] = [f"{s}_elbow_link", f"{s}_wrist_roll_rubber_hand"]   # no wrist servo
SKIP = {"head_link", "imu_in_torso", "imu_in_pelvis", "d435_link", "mid360_link"}


def world_T():
    """Zero-pose world transform of every URDF link (mm)."""
    T = {"pelvis": np.eye(4)}
    T["pelvis"][2, 3] = PELVIS_Z / 1000
    changed = True
    while changed:
        changed = False
        for j in JOINTS:
            if j["parent"] in T and j["child"] not in T:
                T[j["child"]] = T[j["parent"]] @ j["T"]; changed = True
    out = {}
    for k, v in T.items():
        M = v.copy(); M[:3, 3] *= 1000; out[k] = M
    return out


WT = world_T()


def micro_links():
    """name -> list of URDF links fused into it."""
    used = set(x for v in FUSE.values() for x in v)
    out = dict(FUSE)
    for l in LINKS:
        if l in used or l in SKIP or not LINKS[l]: continue
        out[l] = [l]
    return out


def micro_joints():
    """Revolute joints kept in the micro robot: name -> (parent micro link, child micro link, point, axis, lim)."""
    owner = {u: m for m, us in micro_links().items() for u in us}
    out = {}
    for j in JOINTS:
        if j["type"] != "revolute" or "wrist" in j["name"]: continue
        Tc = WT[j["child"]]
        p = Tc[:3, 3]; a = Tc[:3, :3] @ j["axis"]
        out[j["name"].replace("_joint", "")] = (owner[j["parent"]], owner[j["child"]], p, a / np.linalg.norm(a), j["lim"])
    return out


_mesh_cache = {}


def link_mesh(name, s=1.0):
    """Zero-pose world mesh (mm, scaled by s about the origin) of a micro link."""
    key = (name, s)
    if key in _mesh_cache: return _mesh_cache[key]
    Vs, Fs, n = [], [], 0
    for u in micro_links()[name]:
        for fn, Tv in LINKS[u]:
            V, F = read_stl(os.path.join(G1, "meshes", fn))
            M = WT[u] @ Tv
            W = (V * 1000) @ M[:3, :3].T + M[:3, 3]
            Vs.append(W); Fs.append(F + n); n += len(W)
    V = np.concatenate(Vs) * s; F = np.concatenate(Fs)
    _mesh_cache[key] = (V, F)
    return V, F


def head_mesh(s=1.0):
    """The user's OpenG1 dome head (frame, panels, brows), placed as in robot_assembly.scad."""
    Vs, Fs, n = [], [], 0
    for f in ["frame_mm", "panel_0_mm", "panel_1_mm", "panel_2_mm", "panel_3_mm", "panel_4_mm", "brow_left_mm", "brow_right_mm"]:
        V, F = read_stl(os.path.join(REPO, "head", f + ".stl"))
        Vs.append(V + np.array([0, 0, PELVIS_Z])); Fs.append(F + n); n += len(V)
    return np.concatenate(Vs) * s, np.concatenate(Fs)


def voxelize(V, F, lo, h, n):
    with tempfile.TemporaryDirectory() as td:
        a = os.path.join(td, "a"); b = os.path.join(td, "b")
        with open(a, "wb") as fp:
            fp.write(struct.pack("i", len(F))); fp.write(np.ascontiguousarray(V[F].reshape(-1, 9), np.float64).tobytes())
            fp.write(np.asarray(lo, np.float64).tobytes()); fp.write(struct.pack("d", h)); fp.write(np.asarray(n, np.int32).tobytes())
        subprocess.run([os.path.join(HERE, "voxelize"), a, b], check=True)
        return np.fromfile(b, np.uint8).reshape(n).astype(bool)


def mesh_distance(V, F, P, rmax=1.5):
    with tempfile.TemporaryDirectory() as td:
        a = os.path.join(td, "a"); b = os.path.join(td, "b")
        with open(a, "wb") as fp:
            fp.write(struct.pack("i", len(F))); fp.write(np.ascontiguousarray(V[F].reshape(-1, 9), np.float64).tobytes())
            fp.write(struct.pack("i", len(P))); fp.write(np.ascontiguousarray(P, np.float64).tobytes()); fp.write(struct.pack("d", rmax))
        subprocess.run([os.path.join(HERE, "dist"), a, b], check=True)
        return np.fromfile(b, np.float32)
