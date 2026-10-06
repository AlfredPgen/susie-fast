// Kernels for the susie-fast code paths (see R/susie_fast.R).
//
// Each kernel performs the same IEEE operations, in the same order, as
// the R expressions it replaces, so its results are bit-identical to
// them. R/susie_fast.R checks this at run time (fast_self_test) and
// falls back to the R code if a platform ever disagrees.

#include <cpp11.hpp>
#include <cfloat>
#include <cmath>
#include <vector>

using namespace cpp11;

#if defined(__clang__)
#pragma clang fp contract(off)
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif

// Log of the model-level Bayes factor of one single-effect regression,
// as computed by gaussian_ser_lbf() + lbf_stabilization() +
// compute_posterior_weights():
//
//   lbf  <- -0.5 * log(1 + V / s) + h * V / (s * (V + s))
//   lbf[zero] <- 0
//   lpo  <- lbf + logpw
//   m    <- max(lpo)
//   log(sum(exp(lpo - m))) + m
//
// with h = 0.5 * betahat^2 and s = pmax(shat2, eps) precomputed in R.
// The sum is accumulated in long double, as R's sum() does. Returns NA
// if any lpo is NaN, so the caller can use the R code instead.
[[cpp11::register]]
double ser_lbf_model_cpp(const doubles& h, const doubles& s,
                         const logicals& zero, const doubles& logpw,
                         double V) {
  const R_xlen_t p = h.size();
  std::vector<double> lpo(p);
  double m = R_NegInf;
  for (R_xlen_t j = 0; j < p; j++) {
    double lbf;
    if (zero[j] == TRUE) {
      lbf = 0.0;
    } else {
      const double sj = s[j];
      const double a = -0.5 * std::log(1.0 + V / sj);
      const double b = h[j] * V / (sj * (V + sj));
      lbf = a + b;
    }
    const double x = lbf + logpw[j];
    if (std::isnan(x)) return NA_REAL;
    lpo[j] = x;
    if (x > m) m = x;
  }
  long double acc = 0.0;
  for (R_xlen_t j = 0; j < p; j++) acc += std::exp(lpo[j] - m);
  double total;
  if (acc > DBL_MAX) total = R_PosInf;
  else if (acc < -DBL_MAX) total = R_NegInf;
  else total = (double) acc;
  return std::log(total) + m;
}

// Standardisation of a p x p cross-product matrix, as computed by
//   t((1 / csd) * XtX) / csd
// i.e. out[i, j] = ((1 / csd[j]) * XtX[j, i]) / csd[i], without the
// three full-size temporaries. Works in tiles for the transposed reads.
[[cpp11::register]]
doubles scale_xtx_cpp(const doubles_matrix<>& XtX, const doubles& csd) {
  const int p = XtX.nrow();
  writable::doubles out(static_cast<R_xlen_t>(p) * p);
  std::vector<double> inv(p);
  for (int k = 0; k < p; k++) inv[k] = 1 / csd[k];
  const double* x = REAL(XtX.data());
  double* o = REAL(out.data());
  const double* c = REAL(csd.data());
  const int B = 64;
  for (int jj = 0; jj < p; jj += B) {
    const int jend = jj + B < p ? jj + B : p;
    for (int ii = 0; ii < p; ii += B) {
      const int iend = ii + B < p ? ii + B : p;
      for (int j = jj; j < jend; j++) {
        const double ij = inv[j];
        double* oc = o + static_cast<R_xlen_t>(j) * p;
        for (int i = ii; i < iend; i++)
          oc[i] = (ij * x[j + static_cast<R_xlen_t>(i) * p]) / c[i];
      }
    }
  }
  return out;
}
