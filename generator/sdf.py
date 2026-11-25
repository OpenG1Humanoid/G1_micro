"""Tiny signed-distance-field toolkit (numpy).

Every shape is a function f(P) -> distance, where P is a (3, N) float32 array of
world points in millimetres. Negative = inside. Shapes are built in a local
Frame and combined with min / max. Good enough for printable parts at 0.2-0.3 mm.
"""
import numpy as np

F32 = np.float32


class Frame:
    """Origin o and rotation R whose columns are the local x, y, z axes in world."""

    def __init__(self, o=(0, 0, 0), x=(1, 0, 0), y=None, z=None):
        o = np.asarray(o, dtype=np.float64)
        x = np.asarray(x, dtype=np.float64)
        if z is not None and y is None:
            z = np.asarray(z, dtype=np.float64)
            y = np.cross(z, x)
        elif y is not None and z is None:
            y = np.asarray(y, dtype=np.float64)
            z = np.cross(x, y)
        elif y is None and z is None:
            y = np.array([0, 1, 0.0]); z = np.array([0, 0, 1.0])
        else:
            y = np.asarray(y, float); z = np.asarray(z, float)
        R = np.stack([x, y, z], axis=1)
        assert np.allclose(R.T @ R, np.eye(3), atol=1e-9), "frame not orthonormal"
        assert np.linalg.det(R) > 0, "frame not right handed"
        self.o, self.R = o, R

    def local(self, P):
        d = P - self.o.astype(F32)[:, None]
        return (self.R.T.astype(F32) @ d)

    def world(self, p):
        return self.o + self.R @ np.asarray(p, float)

    def wdir(self, v):
        return self.R @ np.asarray(v, float)

    def shifted(self, d):
        """Frame moved by d given in LOCAL coords."""
        f = Frame.__new__(Frame)
        f.o = self.o + self.R @ np.asarray(d, float)
        f.R = self.R.copy()
        return f


WORLD = Frame()


# ---------------------------------------------------------------- primitives
def box(fr, c, h):
    """Box centre c (local), half sizes h."""
    def f(P):
        p = fr.local(P)
        q = np.abs(p - np.asarray(c, F32)[:, None]) - np.asarray(h, F32)[:, None]
        out = np.linalg.norm(np.maximum(q, 0), axis=0)
        ins = np.minimum(np.max(q, axis=0), 0)
        return out + ins
    return f


def boxr(fr, lo, hi):
    lo = np.asarray(lo, float); hi = np.asarray(hi, float)
    return box(fr, (lo + hi) / 2, (hi - lo) / 2)


def rbox(fr, c, h, r):
    """Rounded box: outer half size h, corner radius r."""
    h = np.asarray(h, float)
    r = min(r, float(h.min()) - 1e-3)
    def f(P):
        p = fr.local(P)
        q = np.abs(p - np.asarray(c, F32)[:, None]) - (h - r).astype(F32)[:, None]
        out = np.linalg.norm(np.maximum(q, 0), axis=0)
        ins = np.minimum(np.max(q, axis=0), 0)
        return out + ins - F32(r)
    return f


def rboxr(fr, lo, hi, r):
    lo = np.asarray(lo, float); hi = np.asarray(hi, float)
    return rbox(fr, (lo + hi) / 2, (hi - lo) / 2, r)


def cyl(fr, r, z0, z1, c=(0, 0)):
    """Capped cylinder along local z, radius r, from z0 to z1, centred at local (cx, cy)."""
    def f(P):
        p = fr.local(P)
        d = np.sqrt((p[0] - F32(c[0])) ** 2 + (p[1] - F32(c[1])) ** 2) - F32(r)
        zc = F32((z0 + z1) / 2); zh = F32(abs(z1 - z0) / 2)
        e = np.abs(p[2] - zc) - zh
        out = np.sqrt(np.maximum(d, 0) ** 2 + np.maximum(e, 0) ** 2)
        return out + np.minimum(np.maximum(d, e), 0)
    return f


def rcyl(fr, r, z0, z1, rr, c=(0, 0)):
    """Cylinder with rounded rims (radius rr)."""
    rr = min(rr, r - 1e-3, abs(z1 - z0) / 2 - 1e-3)
    def f(P):
        p = fr.local(P)
        d = np.sqrt((p[0] - F32(c[0])) ** 2 + (p[1] - F32(c[1])) ** 2) - F32(r - rr)
        zc = F32((z0 + z1) / 2); zh = F32(abs(z1 - z0) / 2 - rr)
        e = np.abs(p[2] - zc) - zh
        out = np.sqrt(np.maximum(d, 0) ** 2 + np.maximum(e, 0) ** 2)
        return out + np.minimum(np.maximum(d, e), 0) - F32(rr)
    return f


def cone(fr, r0, r1, z0, z1, c=(0, 0)):
    """Truncated cone along local z (approximate SDF, fine for small slopes)."""
    def f(P):
        p = fr.local(P)
        t = np.clip((p[2] - F32(z0)) / F32(z1 - z0), 0, 1)
        rad = F32(r0) + (F32(r1) - F32(r0)) * t
        k = np.sqrt(1 + ((r1 - r0) / (z1 - z0)) ** 2)
        d = (np.sqrt((p[0] - F32(c[0])) ** 2 + (p[1] - F32(c[1])) ** 2) - rad) / F32(k)
        zc = F32((z0 + z1) / 2); zh = F32(abs(z1 - z0) / 2)
        e = np.abs(p[2] - zc) - zh
        out = np.sqrt(np.maximum(d, 0) ** 2 + np.maximum(e, 0) ** 2)
        return out + np.minimum(np.maximum(d, e), 0)
    return f


def sphere(fr, r, c=(0, 0, 0)):
    def f(P):
        p = fr.local(P) - np.asarray(c, F32)[:, None]
        return np.linalg.norm(p, axis=0) - F32(r)
    return f


def ellipsoid(fr, rads, c=(0, 0, 0)):
    rads = np.asarray(rads, F32)
    def f(P):
        p = (fr.local(P) - np.asarray(c, F32)[:, None])
        k0 = np.linalg.norm(p / rads[:, None], axis=0)
        k1 = np.linalg.norm(p / (rads[:, None] ** 2), axis=0)
        return k0 * (k0 - 1) / np.maximum(k1, 1e-6)
    return f


def capsule(fr, a, b, r):
    a = np.asarray(a, F32); b = np.asarray(b, F32)
    def f(P):
        p = fr.local(P)
        pa = p - a[:, None]; ba = (b - a)[:, None]
        h = np.clip((pa * ba).sum(0) / float((ba * ba).sum()), 0, 1)
        return np.linalg.norm(pa - ba * h, axis=0) - F32(r)
    return f


def halfspace(fr, axis, c, keep_above=True):
    """Solid where local[axis] >= c (or <= c)."""
    def f(P):
        p = fr.local(P)
        return (F32(c) - p[axis]) if keep_above else (p[axis] - F32(c))
    return f


def prism(fr, poly2d, z0, z1):
    """Extrude a convex 2D polygon (CCW, local xy) between z0 and z1 (approx SDF)."""
    poly = np.asarray(poly2d, float)
    n = len(poly)
    normals, offs = [], []
    for i in range(n):
        a, b = poly[i], poly[(i + 1) % n]
        e = b - a
        nrm = np.array([e[1], -e[0]]); nrm /= np.linalg.norm(nrm)
        normals.append(nrm); offs.append(nrm @ a)
    normals = np.asarray(normals, F32); offs = np.asarray(offs, F32)
    def f(P):
        p = fr.local(P)
        d2 = normals @ p[:2] - offs[:, None]
        d = d2.max(0)
        zc = F32((z0 + z1) / 2); zh = F32(abs(z1 - z0) / 2)
        e = np.abs(p[2] - zc) - zh
        out = np.sqrt(np.maximum(d, 0) ** 2 + np.maximum(e, 0) ** 2)
        return out + np.minimum(np.maximum(d, e), 0)
    return f


# ---------------------------------------------------------------- operators
def U(*fs):
    fs = [f for f in fs if f is not None]
    def f(P):
        d = fs[0](P)
        for g in fs[1:]:
            d = np.minimum(d, g(P))
        return d
    return f


def I(*fs):
    def f(P):
        d = fs[0](P)
        for g in fs[1:]:
            d = np.maximum(d, g(P))
        return d
    return f


def D(a, *bs):
    bs = [b for b in bs if b is not None]
    if not bs:
        return a
    def f(P):
        d = a(P)
        for b in bs:
            d = np.maximum(d, -b(P))
        return d
    return f


def SU(a, b, k):
    """Smooth union with blend radius k."""
    def f(P):
        da, db = a(P), b(P)
        h = np.clip(0.5 + 0.5 * (db - da) / k, 0, 1)
        return db * (1 - h) + da * h - F32(k) * h * (1 - h)
    return f


def offset(a, r):
    return lambda P: a(P) - F32(r)


def shell(a, t):
    """Hollow version of a: keeps a band of thickness t inside the surface."""
    return lambda P: np.maximum(a(P), -(a(P) + F32(t)))


def rotated(fn, axis_point, axis_dir, ang_deg):
    """Shape fn rotated by ang about the given world axis."""
    a = np.radians(ang_deg)
    k = np.asarray(axis_dir, float); k = k / np.linalg.norm(k)
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    Rm = np.eye(3) + np.sin(a) * K + (1 - np.cos(a)) * K @ K
    Rinv = Rm.T.astype(F32)
    c = np.asarray(axis_point, F32)[:, None]
    return lambda P: fn(Rinv @ (P - c) + c)


def mirrored_y(fn):
    def f(P):
        Q = P.copy(); Q[1] = -Q[1]
        return fn(Q)
    return f
