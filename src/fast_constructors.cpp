// Kernels for the data constructors (see R/susie_fast_constructors.R).
//
// The checks compare bits and values only; the scaling kernels perform
// the same IEEE operations, in the same order, as the R expressions they
// replace. R/susie_fast_constructors.R checks them against the R code at
// run time and falls back to the R code if a platform ever disagrees.

#include <cpp11.hpp>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>

using namespace cpp11;

#if defined(__clang__)
#pragma clang fp contract(off)
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif

static const int TILE = 64;
static const uint64_t EXP_MASK = 0x7ff0000000000000ULL;

static inline uint64_t double_bits(double a) {
  uint64_t u;
  std::memcpy(&u, &a, sizeof u);
  return u;
}

// Non-finite (Inf, NaN or NA) from the exponent bits, so that compiler
// flags such as -ffinite-math-only cannot fold the test away.
static inline bool bits_nonfinite(uint64_t u) {
  return (u & EXP_MASK) == EXP_MASK;
}

// One tiled pass over a square double matrix x. Returns
//   0  some entry is not finite (the caller uses is_symmetric_matrix()
//      and anyNA() instead);
//   1  all entries finite, and some pair x[i, j], x[j, i] differs in
//      value, i.e. Rfast::is.symmetric(x) is FALSE;
//   2  all entries finite and symmetric in value, but some pair differs
//      in bits (+0 against -0);
//   3  all entries finite and x is bitwise symmetric, so t(x) has the
//      same data as x.
// For finite doubles, different bits mean different values except for a
// pair of zeros.
[[cpp11::register]]
int xtx_check_cpp(const doubles_matrix<>& XtX) {
  const R_xlen_t p = XtX.nrow();
  const double* x = REAL_RO(XtX.data());
  bool asym = false, bitsym = true;
  for (R_xlen_t k = 0; k < p; k++)
    if (bits_nonfinite(double_bits(x[k + k * p]))) return 0;
  for (R_xlen_t jj = 0; jj < p; jj += TILE) {
    const R_xlen_t jend = jj + TILE < p ? jj + TILE : p;
    for (R_xlen_t ii = 0; ii <= jj; ii += TILE) {
      const R_xlen_t iend = ii + TILE < p ? ii + TILE : p;
      for (R_xlen_t j = jj; j < jend; j++) {
        const double* col = x + j * p;
        const R_xlen_t ie = iend < j ? iend : j;
        for (R_xlen_t i = ii; i < ie; i++) {
          const uint64_t a = double_bits(col[i]);
          const uint64_t b = double_bits(x[j + i * p]);
          if (bits_nonfinite(a) || bits_nonfinite(b)) return 0;
          if (a != b) {
            if (((a | b) << 1) != 0) asym = true;
            else bitsym = false;
          }
        }
      }
    }
  }
  if (asym) return 1;
  return bitsym ? 3 : 2;
}

// t((1 / csd) * XtX) / csd for a bitwise symmetric XtX: the transposed
// operand XtX[j, i] has the same bits as XtX[i, j], so the matrix is read
// column by column. out[i, j] = ((1 / csd[j]) * XtX[i, j]) / csd[i].
[[cpp11::register]]
doubles scale_xtx_sym_cpp(const doubles_matrix<>& XtX, const doubles& csd) {
  const R_xlen_t p = XtX.nrow();
  writable::doubles out(p * p);
  const double* x = REAL_RO(XtX.data());
  const double* c = REAL_RO(csd.data());
  double* o = REAL(out.data());
  for (R_xlen_t j = 0; j < p; j++) {
    const double ij = 1 / c[j];
    const double* xc = x + j * p;
    double* oc = o + j * p;
    for (R_xlen_t i = 0; i < p; i++) oc[i] = (ij * xc[i]) / c[i];
  }
  return out;
}

// The original-scale XtX of summary_stats_constructor():
//   XtX <- t(R * s) * s
//   XtX <- (XtX + t(XtX)) / 2
// i.e. A[i, j] = (R[j, i] * s[j]) * s[i] and
// out[i, j] = (A[i, j] + A[j, i]) / 2, every element (the diagonal
// included) computed in that order. Works in tile pairs for the
// transposed reads.
[[cpp11::register]]
doubles orig_scale_xtx_cpp(const doubles_matrix<>& R, const doubles& s) {
  const R_xlen_t p = R.nrow();
  writable::doubles out(p * p);
  const double* r = REAL_RO(R.data());
  const double* sv = REAL_RO(s.data());
  double* o = REAL(out.data());
  for (R_xlen_t jj = 0; jj < p; jj += TILE) {
    const R_xlen_t jend = jj + TILE < p ? jj + TILE : p;
    for (R_xlen_t ii = 0; ii <= jj; ii += TILE) {
      const R_xlen_t iend = ii + TILE < p ? ii + TILE : p;
      for (R_xlen_t j = jj; j < jend; j++) {
        const double sj = sv[j];
        const R_xlen_t ie = iend < j + 1 ? iend : j + 1;
        for (R_xlen_t i = ii; i < ie; i++) {
          const double si = sv[i];
          const double aij = (r[j + i * p] * sj) * si;
          const double aji = (r[i + j * p] * si) * sj;
          o[i + j * p] = (aij + aji) / 2;
          o[j + i * p] = (aji + aij) / 2;
        }
      }
    }
  }
  return out;
}

// Columns of X (n >= 1 rows) whose variance is not decided by their range. A column
// with every |x| <= 1e150 and max - min >= 1e-100 has var() > 0 and not
// NA (some entry lies at least (max - min) / 2 from the mean, and its
// squared deviation neither underflows to 0 nor makes the sum NaN). The
// other columns (constant, tiny range, very large or non-finite entries)
// are flagged TRUE so that the caller evaluates var() on them.
[[cpp11::register]]
logicals const_col_screen_cpp(const doubles_matrix<>& X) {
  const R_xlen_t n = X.nrow();
  const int p = X.ncol();
  writable::logicals out(p);
  const double* x = REAL_RO(X.data());
  for (int j = 0; j < p; j++) {
    const double* col = x + static_cast<R_xlen_t>(j) * n;
    bool undecided = false;
    double mn = col[0], mx = col[0];
    for (R_xlen_t k = 0; k < n; k++) {
      const double v = col[k];
      if (!(std::fabs(v) <= 1e150)) { undecided = true; break; }
      if (v < mn) mn = v;
      if (v > mx) mx = v;
    }
    if (!undecided) undecided = !(mx - mn >= 1e-100);
    out[j] = undecided ? TRUE : FALSE;
  }
  return out;
}
