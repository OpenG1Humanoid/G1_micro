"""Assembly of the G1 micro: G1 kinematics scaled by S, pieces from g1_cache, right side = mirror
of the left pieces. Joint angles are G1 joint angles (radians in G1, degrees here)."""
import os, json
import numpy as np
import g1
from mesh import render

S = 0.35
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
CACHE = os.path.join(ROOT, "g1_cache")
MIR = np.array([1, -1, 1.0])

JOINTS = g1.micro_joints()                      # name -> (parent, child, pt, axis, lim)
PARENT = {c: (p, j) for j, (p, c, *_ ) in JOINTS.items()}


def rot(axis, deg):
    a = np.radians(deg); k = np.asarray(axis, float); k /= np.linalg.norm(k)
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) + np.sin(a) * K + (1 - np.cos(a)) * K @ K


def T_about(pt, axis, deg):
    T = np.eye(4); R = rot(axis, deg); T[:3, :3] = R; T[:3, 3] = pt - R @ pt; return T


def fk(q):
    T = {}
    def get(link):
        if link in T: return T[link]
        if link not in PARENT:
            T[link] = np.eye(4); return T[link]
        p, j = PARENT[link]
        _, _, pt, ax, _ = JOINTS[j]
        T[link] = get(p) @ T_about(pt * S, ax, q.get(j, 0.0)); return T[link]
    for l in list(g1.micro_links()) + ["head"]:
        if l == "head":
            T["head"] = get("torso")
        else:
            get(l)
    return T


_cache = {}


def pieces(link):
    """dict piece -> (V, F) in zero-pose robot coordinates for any link (right = mirrored left)."""
    if link in _cache: return _cache[link]
    src = link.replace("right_", "left_")
    z = np.load(os.path.join(CACHE, src + ".npz"))
    out = {}
    for k in sorted(set(f[:-2] for f in z.files)):
        V, F = z[k + "_V"], z[k + "_F"]
        if link.startswith("right_"):
            V = V * MIR; F = F[:, ::-1]
        out[k] = (V, F)
    _cache[link] = out
    return out


def apply(T, V):
    return V @ T[:3, :3].T + T[:3, 3]


LINKS = [l for l in g1.micro_links()] + ["head"]
WHITE, GREY, DARK, ORANGE = (236, 236, 232), (175, 178, 186), (70, 72, 78), (240, 150, 50)


def colour(link, piece):
    if link == "head": return DARK
    if "ankle_pitch" in link: return DARK
    if piece.startswith("lid"): return (226, 226, 222)
    return WHITE


def scene(q, exclude=()):
    T = fk(q); ms = []
    for l in LINKS:
        if l in exclude: continue
        for k, (V, F) in pieces(l).items():
            ms.append((apply(T[l], V), F, colour(l, k)))
    return ms


def g1_reference(q=None):
    """Unmodified G1 parts at the same scale, for comparison renders."""
    T = fk(q or {}); ms = []
    for l in g1.micro_links():
        V, F = g1.link_mesh(l, S); ms.append((apply(T[l], V), F, (200, 205, 215)))
    V, F = g1.head_mesh(S); ms.append((apply(T["torso"], V), F, DARK))
    return ms
