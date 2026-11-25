// Solid voxelisation of a triangle mesh by ray parity along x, y and z (2-of-3 vote,
// tolerant of small holes). in: int32 nf, double[nf*9], double lo[3], double h, int32 n[3]
// out: uint8[n0*n1*n2] (index i*n1*n2 + j*n2 + k)
#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cstdint>
int main(int argc, char** argv) {
  FILE* f = fopen(argv[1], "rb"); int nf; if (fread(&nf, 4, 1, f) != 1) return 1;
  std::vector<double> T(nf * 9); if (fread(T.data(), 8, nf * 9, f) != (size_t)nf * 9) return 1;
  double lo[3], h; int n[3];
  if (fread(lo, 8, 3, f) != 3 || fread(&h, 8, 1, f) != 1 || fread(n, 4, 3, f) != 3) return 1; fclose(f);
  size_t N = (size_t)n[0] * n[1] * n[2];
  std::vector<uint8_t> votes(N, 0);
  for (int ax = 0; ax < 3; ax++) {
    int a = (ax + 1) % 3, b = (ax + 2) % 3;
    std::vector<std::vector<float>> hits((size_t)n[a] * n[b]);
    for (int t = 0; t < nf; t++) {
      double* p = &T[9 * t];
      double A[3] = {p[a], p[3 + a], p[6 + a]}, B[3] = {p[b], p[3 + b], p[6 + b]}, C[3] = {p[ax], p[3 + ax], p[6 + ax]};
      double amin = std::min({A[0], A[1], A[2]}), amax = std::max({A[0], A[1], A[2]});
      double bmin = std::min({B[0], B[1], B[2]}), bmax = std::max({B[0], B[1], B[2]});
      int i0 = std::max(0, (int)std::ceil((amin - lo[a]) / h)), i1 = std::min(n[a] - 1, (int)std::floor((amax - lo[a]) / h));
      int j0 = std::max(0, (int)std::ceil((bmin - lo[b]) / h)), j1 = std::min(n[b] - 1, (int)std::floor((bmax - lo[b]) / h));
      double den = (A[1] - A[0]) * (B[2] - B[0]) - (A[2] - A[0]) * (B[1] - B[0]);
      if (std::fabs(den) < 1e-14) continue;
      for (int i = i0; i <= i1; i++) for (int j = j0; j <= j1; j++) {
        double x = lo[a] + i * h + 1.3e-7, y = lo[b] + j * h + 2.9e-7;
        double w1 = ((x - A[0]) * (B[2] - B[0]) - (A[2] - A[0]) * (y - B[0])) / den;
        double w2 = ((A[1] - A[0]) * (y - B[0]) - (x - A[0]) * (B[1] - B[0])) / den;
        double w0 = 1 - w1 - w2;
        if (w0 < 0 || w1 < 0 || w2 < 0) continue;
        hits[(size_t)i * n[b] + j].push_back((float)(w0 * C[0] + w1 * C[1] + w2 * C[2]));
      }
    }
    for (int i = 0; i < n[a]; i++) for (int j = 0; j < n[b]; j++) {
      auto& v = hits[(size_t)i * n[b] + j]; if (v.size() < 2) continue;
      std::sort(v.begin(), v.end());
      for (size_t q = 0; q + 1 < v.size(); q += 2) {
        int k0 = std::max(0, (int)std::ceil((v[q] - lo[ax]) / h)), k1 = std::min(n[ax] - 1, (int)std::floor((v[q + 1] - lo[ax]) / h));
        for (int k = k0; k <= k1; k++) {
          int idx[3]; idx[ax] = k; idx[a] = i; idx[b] = j;
          votes[(size_t)idx[0] * n[1] * n[2] + (size_t)idx[1] * n[2] + idx[2]]++;
        }
      }
    }
  }
  std::vector<uint8_t> out(N); for (size_t q = 0; q < N; q++) out[q] = votes[q] >= 2;
  FILE* o = fopen(argv[2], "wb"); fwrite(out.data(), 1, N, o); fclose(o); return 0;
}
