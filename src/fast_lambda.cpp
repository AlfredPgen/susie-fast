// Kernels for the susie_rss_lambda() fast paths (see R/susie_fast_lambda.R).
//
// Same IEEE operations, in the same order, as the R expression replaced;
// R/susie_fast_lambda.R checks this at run time (fast_lambda_self_test) and
// falls back to the R code if a platform ever disagrees.

#include <cpp11.hpp>
#include <vector>

using namespace cpp11;

#if defined(__clang__)
#pragma clang fp contract(off)
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif

// Diagonal of R S^{-1} R in get_ER2.rss_lambda:
//
//   rowSums((V^2) * rep(w, each = nrow(V)))
//
// i.e. out[i] = sum_j (V[i, j] * V[i, j]) * w[j], each term rounded to
// double, accumulated column by column in long double as R's rowSums()
// does, without the two n x p temporaries.
[[cpp11::register]]
doubles rss_diag_rsinvr_cpp(const doubles_matrix<>& V, const doubles& w) {
  const int n = V.nrow();
  const int p = V.ncol();
  std::vector<long double> acc(n, 0.0L);
  const double* v = REAL(V.data());
  const double* wp = REAL(w.data());
  for (int j = 0; j < p; j++) {
    const double wj = wp[j];
    const double* vc = v + static_cast<R_xlen_t>(j) * n;
    for (int i = 0; i < n; i++) {
      const double t = vc[i] * vc[i];
      const double u = t * wj;
      acc[i] += u;
    }
  }
  writable::doubles out(n);
  double* o = REAL(out.data());
  for (int i = 0; i < n; i++) o[i] = static_cast<double>(acc[i]);
  return out;
}
