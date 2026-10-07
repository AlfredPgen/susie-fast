// Single-effect regression kernels (see R/susie_fast_ser.R).
//
// As in susie_fast.cpp, each kernel performs the same IEEE operations, in
// the same order, as the R expressions it replaces. R/susie_fast_ser.R
// checks this at run time (fast_ser_self_test) and falls back to the R
// code if a platform ever disagrees.

#include <cpp11.hpp>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>

using namespace cpp11;
namespace writable = cpp11::writable;
using namespace cpp11::literals;

#if defined(__clang__)
#pragma clang fp contract(off)
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif

// R's sum() of doubles: long double accumulator, clamped to +-Inf.
static double ld_to_double(long double acc) {
  if (acc > DBL_MAX) return R_PosInf;
  if (acc < -DBL_MAX) return R_NegInf;
  return (double) acc;
}

// Groups the entries of s by bit pattern. Returns list(u, sidx): u holds
// the distinct values among the entries with zero[j] FALSE, and sidx[j]
// is 0 where zero[j] is TRUE and otherwise the 1-based index of s[j] in
// u. Returns NULL when there are more than `cap` distinct values.
[[cpp11::register]]
SEXP ser_s_table_cpp(const doubles& s, const logicals& zero, int cap) {
  const R_xlen_t p = s.size();
  const double* sp = REAL(s.data());
  const int* zp = LOGICAL(zero.data());
  writable::integers sidx(p);
  int* ip = INTEGER(sidx.data());
  std::vector<std::uint64_t> bits;
  std::vector<double> u;
  int last = 0;
  for (R_xlen_t j = 0; j < p; j++) {
    if (zp[j] == TRUE) {
      ip[j] = 0;
      continue;
    }
    std::uint64_t b;
    std::memcpy(&b, sp + j, sizeof b);
    if (last == 0 || bits[last - 1] != b) {
      int k = 0;
      const int K = (int) bits.size();
      while (k < K && bits[k] != b) k++;
      if (k == K) {
        if (K >= cap) return R_NilValue;
        bits.push_back(b);
        u.push_back(sp[j]);
      }
      last = k + 1;
    }
    ip[j] = last;
  }
  writable::doubles uu(u.begin(), u.end());
  return writable::list({"u"_nm = uu, "sidx"_nm = sidx});
}

// ser_lbf_model_cpp (susie_fast.cpp) for s taking few distinct values:
// log(1 + V / s) and s * (V + s) are computed once per distinct s, from
// the same operands, so every lbf is bit-identical to the per-variant
// kernel. u and sidx come from ser_s_table_cpp.
[[cpp11::register]]
double ser_lbf_model_tab_cpp(const doubles& h, const doubles& u,
                             const integers& sidx, const doubles& logpw,
                             double V) {
  const R_xlen_t p = h.size();
  const R_xlen_t K = u.size();
  const double* hp = REAL(h.data());
  const double* lp = REAL(logpw.data());
  const int* ip = INTEGER(sidx.data());
  std::vector<double> A(K), D(K);
  for (R_xlen_t k = 0; k < K; k++) {
    const double sk = u[k];
    A[k] = -0.5 * std::log(1.0 + V / sk);
    D[k] = sk * (V + sk);
  }
  std::vector<double> lpo(p);
  double m = R_NegInf;
  for (R_xlen_t j = 0; j < p; j++) {
    double lbf;
    const int k = ip[j];
    if (k == 0) {
      lbf = 0.0;
    } else {
      const double a = A[k - 1];
      const double b = hp[j] * V / D[k - 1];
      lbf = a + b;
    }
    const double x = lbf + lp[j];
    if (std::isnan(x)) return NA_REAL;
    lpo[j] = x;
    if (x > m) m = x;
  }
  long double acc = 0.0;
  for (R_xlen_t j = 0; j < p; j++) acc += std::exp(lpo[j] - m);
  return std::log(ld_to_double(acc)) + m;
}

// Log Bayes factors, posterior weights and model log Bayes factor of one
// SER at its final V, as computed by gaussian_ser_lbf() +
// lbf_stabilization() + compute_posterior_weights():
//
//   s     <- pmax(shat2, eps)
//   lbf   <- -0.5 * log(1 + V / s) + 0.5 * betahat^2 * V / (s * (V + s))
//   lbf[!is.finite(betahat) | !is.finite(shat2)] <- 0
//   lpo   <- lbf + log(pi + sqrt(eps))
//   lbf[is.infinite(shat2)] <- 0; lpo[is.infinite(shat2)] <- log(pi + sqrt(eps))
//   w     <- exp(lpo - max(lpo)); alpha <- w / sum(w)
//   lbf_model <- log(sum(w)) + max(lpo)
//
// Returns list(lbf, alpha, lbf_model), or NULL when an lpo is not finite
// (the caller then uses the R code). The caller checks pi >= 0, so the
// R code would not warn.
[[cpp11::register]]
SEXP ser_lbf_l_cpp(const doubles& betahat, const doubles& shat2,
                   const doubles& pi, double V) {
  const R_xlen_t p = betahat.size();
  const double* bp = REAL(betahat.data());
  const double* sp = REAL(shat2.data());
  const double* pp = REAL(pi.data());
  const double eps = DBL_EPSILON;
  const double seps = std::sqrt(DBL_EPSILON);
  writable::doubles lbf(p), alpha(p);
  double* lb = REAL(lbf.data());
  double* al = REAL(alpha.data());
  std::vector<double> lpo(p);
  double m = R_NegInf;
  for (R_xlen_t j = 0; j < p; j++) {
    const double bj = bp[j], sj = sp[j];
    double l;
    if (!std::isfinite(bj) || !std::isfinite(sj)) {
      l = 0.0;
    } else {
      const double s = sj < eps ? eps : sj;
      const double a = -0.5 * std::log(1.0 + V / s);
      const double b = 0.5 * (bj * bj) * V / (s * (V + s));
      l = a + b;
    }
    const double logpw = std::log(pp[j] + seps);
    double x = l + logpw;
    if (std::isinf(sj)) {
      l = 0.0;
      x = logpw;
    }
    if (!std::isfinite(x)) return R_NilValue;
    lb[j] = l;
    lpo[j] = x;
    if (x > m) m = x;
  }
  long double acc = 0.0;
  for (R_xlen_t j = 0; j < p; j++) {
    al[j] = std::exp(lpo[j] - m);
    acc += al[j];
  }
  const double S = ld_to_double(acc);
  for (R_xlen_t j = 0; j < p; j++) al[j] = al[j] / S;
  const double lbf_model = std::log(S) + m;
  return writable::list({"lbf"_nm = lbf, "alpha"_nm = alpha,
                         "lbf_model"_nm = lbf_model});
}

// Posterior moments of one SER for ss data, as computed by
// calculate_posterior_moments.ss() + gaussian_ser_moments():
//
//   shat2 <- rv / pw; betahat <- r / pw
//   post_var  <- V * shat2 / (V + shat2)
//   post_mean <- post_var / shat2 * betahat
//   post_var[no_info] <- 0; post_mean[no_info] <- 0
//   post_mean2 <- post_var + post_mean^2
//
// Returns list(post_mean, post_mean2), or NULL when a value is NaN.
[[cpp11::register]]
SEXP ser_moments_cpp(const doubles& r, const doubles& pw, double rv,
                     double V) {
  const R_xlen_t p = r.size();
  const double* rp = REAL(r.data());
  const double* wp = REAL(pw.data());
  writable::doubles mu(p), mu2(p);
  double* m1 = REAL(mu.data());
  double* m2 = REAL(mu2.data());
  for (R_xlen_t j = 0; j < p; j++) {
    const double s2 = rv / wp[j];
    const double b = rp[j] / wp[j];
    double pv = V * s2 / (V + s2);
    double pm = pv / s2 * b;
    if (!std::isfinite(b) || !std::isfinite(s2)) {
      pv = 0.0;
      pm = 0.0;
    }
    const double pm2 = pv + pm * pm;
    if (std::isnan(pm) || std::isnan(pm2)) return R_NilValue;
    m1[j] = pm;
    m2[j] = pm2;
  }
  return writable::list({"post_mean"_nm = mu, "post_mean2"_nm = mu2});
}

// Posterior expected log-likelihood of one SER for ss data, as computed
// by SER_posterior_e_loglik.ss() + gaussian_ser_posterior_e_loglik():
//
//   shat2 <- rv / pw; betahat <- r / pw; fi <- is.finite(shat2)
//   Eb <- alpha * mu; Eb2 <- alpha * mu2
//   -0.5 * sum((-2 * Eb[fi] * betahat[fi] + Eb2[fi]) / shat2[fi])
//
// Returns NA when the result is not finite (the caller then uses the R
// code).
[[cpp11::register]]
double ser_e_loglik_cpp(const doubles& alpha, const doubles& mu,
                        const doubles& mu2, const doubles& r,
                        const doubles& pw, double rv) {
  const R_xlen_t p = alpha.size();
  const double* ap = REAL(alpha.data());
  const double* m1 = REAL(mu.data());
  const double* m2 = REAL(mu2.data());
  const double* rp = REAL(r.data());
  const double* wp = REAL(pw.data());
  long double acc = 0.0;
  for (R_xlen_t j = 0; j < p; j++) {
    const double s2 = rv / wp[j];
    if (!std::isfinite(s2)) continue;
    const double b = rp[j] / wp[j];
    const double Eb = ap[j] * m1[j];
    const double Eb2 = ap[j] * m2[j];
    const double t = (-2 * Eb * b + Eb2) / s2;
    acc += t;
  }
  const double e = -0.5 * ld_to_double(acc);
  return std::isfinite(e) ? e : NA_REAL;
}
