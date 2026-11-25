// Exact unsigned distance from query points to a triangle mesh (bucketed), capped at rmax.
// in: int32 nf, double[nf*9], int32 np, double[np*3], double rmax ; out: float[np]
#include <cstdio>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cstdint>
struct V { double x, y, z; };
static inline V sub(V a, V b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
static inline double dot(V a, V b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static double tri_d2(V p, V a, V b, V c) {   // Ericson, closest point on triangle
  V ab = sub(b, a), ac = sub(c, a), ap = sub(p, a);
  double d1 = dot(ab, ap), d2 = dot(ac, ap); if (d1 <= 0 && d2 <= 0) return dot(ap, ap);
  V bp = sub(p, b); double d3 = dot(ab, bp), d4 = dot(ac, bp); if (d3 >= 0 && d4 <= d3) return dot(bp, bp);
  double vc = d1 * d4 - d3 * d2;
  if (vc <= 0 && d1 >= 0 && d3 <= 0) { double v = d1 / (d1 - d3); V q = {a.x + v * ab.x - p.x, a.y + v * ab.y - p.y, a.z + v * ab.z - p.z}; return dot(q, q); }
  V cp = sub(p, c); double d5 = dot(ab, cp), d6 = dot(ac, cp); if (d6 >= 0 && d5 <= d6) return dot(cp, cp);
  double vb = d5 * d2 - d1 * d6;
  if (vb <= 0 && d2 >= 0 && d6 <= 0) { double w = d2 / (d2 - d6); V q = {a.x + w * ac.x - p.x, a.y + w * ac.y - p.y, a.z + w * ac.z - p.z}; return dot(q, q); }
  double va = d3 * d6 - d5 * d4;
  if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0) { double w = (d4 - d3) / ((d4 - d3) + (d5 - d6)); V bc = sub(c, b);
    V q = {b.x + w * bc.x - p.x, b.y + w * bc.y - p.y, b.z + w * bc.z - p.z}; return dot(q, q); }
  double den = 1.0 / (va + vb + vc), v = vb * den, w = vc * den;
  V q = {a.x + ab.x * v + ac.x * w - p.x, a.y + ab.y * v + ac.y * w - p.y, a.z + ab.z * v + ac.z * w - p.z}; return dot(q, q);
}
int main(int argc, char** argv) {
  FILE* f = fopen(argv[1], "rb"); int nf, np; double rmax;
  if (fread(&nf, 4, 1, f) != 1) return 1; std::vector<double> T(nf * 9); if (fread(T.data(), 8, nf * 9, f) != (size_t)nf * 9) return 1;
  if (fread(&np, 4, 1, f) != 1) return 1; std::vector<double> P(np * 3); if (fread(P.data(), 8, np * 3, f) != (size_t)np * 3) return 1;
  if (fread(&rmax, 8, 1, f) != 1) return 1; fclose(f);
  double cs = std::max(rmax, 0.5);
  auto key = [&](long i, long j, long k) { return (i * 73856093L) ^ (j * 19349663L) ^ (k * 83492791L); };
  std::unordered_map<long, std::vector<int>> grid;
  for (int t = 0; t < nf; t++) {
    double* p = &T[9 * t]; double lo[3], hi[3];
    for (int a = 0; a < 3; a++) { lo[a] = std::min({p[a], p[3 + a], p[6 + a]}); hi[a] = std::max({p[a], p[3 + a], p[6 + a]}); }
    for (long i = (long)std::floor(lo[0] / cs); i <= (long)std::floor(hi[0] / cs); i++)
      for (long j = (long)std::floor(lo[1] / cs); j <= (long)std::floor(hi[1] / cs); j++)
        for (long k = (long)std::floor(lo[2] / cs); k <= (long)std::floor(hi[2] / cs); k++) grid[key(i, j, k)].push_back(t);
  }
  std::vector<float> out(np);
  for (int q = 0; q < np; q++) {
    V p = {P[3 * q], P[3 * q + 1], P[3 * q + 2]}; double best = rmax * rmax;
    long ci = (long)std::floor(p.x / cs), cj = (long)std::floor(p.y / cs), ck = (long)std::floor(p.z / cs);
    for (long i = ci - 1; i <= ci + 1; i++) for (long j = cj - 1; j <= cj + 1; j++) for (long k = ck - 1; k <= ck + 1; k++) {
      auto it = grid.find(key(i, j, k)); if (it == grid.end()) continue;
      for (int t : it->second) { double* a = &T[9 * t]; double d2 = tri_d2(p, {a[0], a[1], a[2]}, {a[3], a[4], a[5]}, {a[6], a[7], a[8]}); if (d2 < best) best = d2; }
    }
    out[q] = (float)std::sqrt(best);
  }
  FILE* o = fopen(argv[2], "wb"); fwrite(out.data(), 4, np, o); fclose(o); return 0;
}
