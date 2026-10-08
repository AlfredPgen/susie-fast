// Kernels for the susie-fast post-processing paths (see R/susie_fast_post.R).

#include <cpp11.hpp>
#include <cmath>

using namespace cpp11;

#if defined(__clang__)
#pragma clang fp contract(off)
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif

// Largest |M[pos[a], pos[b]]| over the block M[pos, pos] (pos 1-based and
// in range, checked in R), or Inf as soon as an entry is NaN, NA or
// infinite. Reads the block in place instead of copying it.
[[cpp11::register]]
double block_maxabs_cpp(const doubles_matrix<>& M, const integers& pos) {
  const R_xlen_t m = pos.size();
  const R_xlen_t nr = M.nrow();
  const double* x = REAL(M.data());
  double mx = 0.0;
  for (R_xlen_t b = 0; b < m; b++) {
    const double* col = x + static_cast<R_xlen_t>(pos[b] - 1) * nr;
    for (R_xlen_t a = 0; a < m; a++) {
      const double v = col[pos[a] - 1];
      if (!std::isfinite(v)) return R_PosInf;
      const double av = std::fabs(v);
      if (av > mx) mx = av;
    }
  }
  return mx;
}
