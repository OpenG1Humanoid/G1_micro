"""Grid evaluation, meshing, decimation, STL I/O, checks and rendering helpers."""
import os, struct, subprocess, tempfile
import numpy as np
from skimage.measure import marching_cubes

HERE = os.path.dirname(os.path.abspath(__file__))
F32 = np.float32


class Grid:
    def __init__(self, lo, hi, h):
        self.h = h
        lo = np.floor(np.asarray(lo, float) / h) * h - 2 * h + np.array([0.0371, 0.0529, 0.0613]) * h
        hi = np.ceil(np.asarray(hi, float) / h) * h + 2 * h
        self.lo = lo
        self.n = np.round((hi - lo) / h).astype(int) + 1
        self.axes = [lo[i] + h * np.arange(self.n[i]) for i in range(3)]

    def points(self, sl=None):
        X, Y, Z = np.meshgrid(*[a.astype(F32) for a in self.axes], indexing="ij")
        return np.stack([X.ravel(), Y.ravel(), Z.ravel()])

    def eval(self, f, chunk=4_000_000):
        """Evaluate in z-slabs to bound memory."""
        nx, ny, nz = self.n
        out = np.empty((nx, ny, nz), F32)
        per = max(1, chunk // (nx * ny))
        ax = [a.astype(F32) for a in self.axes]
        for k0 in range(0, nz, per):
            k1 = min(nz, k0 + per)
            X, Y, Z = np.meshgrid(ax[0], ax[1], ax[2][k0:k1], indexing="ij")
            P = np.stack([X.ravel(), Y.ravel(), Z.ravel()])
            out[:, :, k0:k1] = f(P).reshape(nx, ny, k1 - k0)
        return out


def mc(vol, grid):
    v = vol.astype(np.float32, copy=True)
    v[v == 0] = np.float32(1e-5)
    v[0, :, :] = v[-1, :, :] = 1; v[:, 0, :] = v[:, -1, :] = 1; v[:, :, 0] = v[:, :, -1] = 1
    if (v < 0).sum() == 0:
        return np.zeros((0, 3)), np.zeros((0, 3), int)
    verts, faces, _, _ = marching_cubes(v, level=0.0, spacing=(grid.h,) * 3, allow_degenerate=False)
    verts = verts + grid.lo
    pass
    return verts.astype(np.float64), faces.astype(np.int32)


def weld(V, F, tol=1e-6):
    key = np.round(V / tol).astype(np.int64)
    _, idx, inv = np.unique(key, axis=0, return_index=True, return_inverse=True)
    V2 = V[idx]; F2 = inv.reshape(-1)[F]
    good = (F2[:, 0] != F2[:, 1]) & (F2[:, 1] != F2[:, 2]) & (F2[:, 0] != F2[:, 2])
    return V2, F2[good]


def decimate(V, F, maxerr=0.002, target=0):
    with tempfile.TemporaryDirectory() as td:
        a = os.path.join(td, "a.bin"); b = os.path.join(td, "b.bin")
        with open(a, "wb") as fp:
            fp.write(struct.pack("ii", len(V), len(F)))
            fp.write(np.ascontiguousarray(V, np.float64).tobytes())
            fp.write(np.ascontiguousarray(F, np.int32).tobytes())
        subprocess.run([os.path.join(HERE, "decimate"), a, b, str(maxerr), str(target)], check=True)
        with open(b, "rb") as fp:
            nv, nf = struct.unpack("ii", fp.read(8))
            V2 = np.frombuffer(fp.read(nv * 24), np.float64).reshape(nv, 3).copy()
            F2 = np.frombuffer(fp.read(nf * 12), np.int32).reshape(nf, 3).copy()
    return V2, F2


def write_stl(path, V, F, name="part"):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    T = V[F].astype(np.float32)
    n = np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0])
    n /= np.maximum(np.linalg.norm(n, axis=1, keepdims=True), 1e-12)
    data = np.zeros(len(F), dtype=[("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")])
    data["n"] = n; data["v"] = T
    with open(path, "wb") as fp:
        hdr = (name.encode() + b" " * 80)[:80]
        fp.write(hdr); fp.write(struct.pack("<I", len(F))); fp.write(data.tobytes())


def read_stl(path):
    with open(path, "rb") as fp:
        fp.read(80); n = struct.unpack("<I", fp.read(4))[0]
        data = np.frombuffer(fp.read(n * 50), dtype=[("n", "<f4", 3), ("v", "<f4", (3, 3)), ("a", "<u2")])
    T = data["v"].astype(np.float64).reshape(-1, 3)
    V, inv = np.unique(np.round(T, 5), axis=0, return_inverse=True)
    return V, inv.reshape(-1, 3).astype(np.int32)


def check_mesh(V, F):
    """Watertight / manifold / orientation / volume."""
    e = np.concatenate([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]])
    und = np.sort(e, axis=1)
    _, cnt = np.unique(und, axis=0, return_counts=True)
    manifold = bool((cnt == 2).all())
    # directed edges should each appear once if consistently oriented
    _, dcnt = np.unique(e, axis=0, return_counts=True)
    oriented = bool((dcnt == 1).all())
    T = V[F]
    vol = float(np.einsum("ij,ij->i", T[:, 0], np.cross(T[:, 1], T[:, 2])).sum() / 6.0)
    # connected components (shells)
    parent = np.arange(len(V))
    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]; a = parent[a]
        return a
    for a, b in und[::1]:
        ra, rb = find(a), find(b)
        if ra != rb: parent[ra] = rb
    roots = len({find(i) for i in np.unique(F)})
    return dict(watertight=manifold, oriented=oriented, volume_mm3=round(vol, 1), shells=roots,
                triangles=int(len(F)))


def undercut_report(occ, axis):
    """Count columns (along axis) that hit the solid in more than one interval."""
    o = np.moveaxis(occ, axis, -1).astype(np.int8)
    starts = np.diff(np.concatenate([np.zeros(o.shape[:-1] + (1,), np.int8), o], axis=-1), axis=-1) == 1
    n = starts.sum(-1)
    bad = int((n > 1).sum()); tot = int((n > 0).sum())
    return bad, tot


def monotone_cut(occ_cut, axis, toward_low):
    """Extend a cut region along `axis` toward the low (or high) end so the remaining
    solid stays a single interval per column."""
    c = np.moveaxis(occ_cut, axis, -1)
    if toward_low:
        ext = np.flip(np.logical_or.accumulate(np.flip(c, -1), axis=-1), -1)
    else:
        ext = np.logical_or.accumulate(c, axis=-1)
    return np.moveaxis(ext, -1, axis)


# ------------------------------------------------------------------ render
def look(eye, target, up=(0, 0, 1)):
    eye = np.asarray(eye, float); target = np.asarray(target, float)
    f = target - eye; f /= np.linalg.norm(f)
    r = np.cross(f, up); r /= np.linalg.norm(r)
    u = np.cross(r, f)
    return r, u, f, eye


def render(meshes, path, eye, target, scale, W=1000, H=1200, bg=(1, 1, 1), up=(0, 0, 1)):
    """meshes: list of (V, F, rgba 0-255). Orthographic, `scale` px per mm."""
    r, u, f, e = look(eye, target, up)
    tris, cols = [], []
    for V, F, c in meshes:
        tris.append(V[F].reshape(-1, 9).astype(np.float32))
        cc = np.zeros((len(F), 4), np.uint8); cc[:] = c if len(c) == 4 else tuple(c) + (255,)
        cols.append(cc)
    T = np.concatenate(tris); C = np.concatenate(cols)
    # draw opaque first, transparent last
    order = np.argsort(C[:, 3] < 250, kind="stable")
    T, C = T[order], C[order]
    with tempfile.TemporaryDirectory() as td:
        a = os.path.join(td, "in.bin"); b = os.path.join(td, "o.ppm")
        with open(a, "wb") as fp:
            fp.write(struct.pack("i", len(T))); fp.write(T.tobytes()); fp.write(C.tobytes())
            cam = np.array(list(r) + list(u) + list(f) + list(e) + [scale] + list(bg), np.float32)
            fp.write(cam.tobytes()); fp.write(struct.pack("ii", W, H))
        subprocess.run([os.path.join(HERE, "raster"), a, b], check=True)
        from PIL import Image
        im = Image.open(b)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        im.save(path)
    return path


def undercut_where(occ, axis, grid, max_clusters=8):
    """Bounding boxes (world) of clusters of undercut columns, for debugging."""
    from scipy import ndimage as ndi
    o = np.moveaxis(occ, axis, -1).astype(np.int8)
    starts = np.diff(np.concatenate([np.zeros(o.shape[:-1] + (1,), np.int8), o], axis=-1), axis=-1) == 1
    bad = starts.sum(-1) > 1
    lab, n = ndi.label(bad)
    others = [i for i in range(3) if i != axis]
    out = []
    for k in range(1, n + 1):
        idx = np.argwhere(lab == k)
        if len(idx) < 3: continue
        lo = [grid.axes[others[j]][idx[:, j].min()] for j in range(2)]
        hi = [grid.axes[others[j]][idx[:, j].max()] for j in range(2)]
        # depth intervals of the first column
        c = tuple(idx[len(idx) // 2])
        col = o[c]; ch = np.flatnonzero(np.diff(np.concatenate([[0], col, [0]])))
        zs = [round(float(grid.axes[axis][min(i, len(grid.axes[axis]) - 1)]), 2) for i in ch]
        out.append((len(idx), [round(float(v), 1) for v in lo], [round(float(v), 1) for v in hi], zs))
    out.sort(key=lambda t: -t[0])
    return out[:max_clusters]


def keep_largest_interval(occ, axis):
    """Per column along axis keep only the longest solid run (removes sliver undercuts)."""
    o = np.moveaxis(occ, axis, -1)
    n = o.shape[-1]
    flat = o.reshape(-1, n)
    starts = np.diff(np.concatenate([np.zeros((len(flat), 1), bool), flat], 1).astype(np.int8), axis=1) == 1
    multi = np.flatnonzero(starts.sum(1) > 1)
    out = flat.copy()
    for i in multi:
        col = flat[i]
        e = np.flatnonzero(np.diff(np.concatenate([[0], col.astype(np.int8), [0]])))
        runs = list(zip(e[0::2], e[1::2]))
        a, b = max(runs, key=lambda r: r[1] - r[0])
        out[i] = False; out[i, a:b] = True
    removed = flat & ~out
    return np.moveaxis(removed.reshape(o.shape), -1, axis), len(multi)
