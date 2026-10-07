// Matrix-vector products for the susie-fast code paths (see
// R/susie_fast_products.R).
//
// R's %*%, crossprod() and tcrossprod() first scan both operands for
// NaN/Inf (mayHaveNaNOrInf in R's array.c) and then call BLAS dgemv.
// The data matrix does not change during a fit, so the R code scans it
// once per fit with fast_mayhave_cpp() and then calls fast_gemv_cpp(),
// which scans only the vector and calls the same dgemv with the same
// arguments. With R's reference BLAS it can instead use the kernels
// below, which perform the reference dgemv's operations in the same
// order for every output element (checked at run time against R's own
// products before use).

#define USE_FC_LEN_T
#include <Rconfig.h>
#include <cpp11.hpp>
#include <R_ext/BLAS.h>
#ifndef FCONE
#define FCONE
#endif
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif

#if defined(__clang__)
#pragma clang fp contract(off)
#pragma STDC FP_CONTRACT OFF
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif

// R's mayHaveNaNOrInf(): if the length is odd, test x[0]; then test
// x[i] + x[i+1] for each following pair. A pair sum is non-finite
// exactly when s - s is NaN, so the pairs are reduced as sums of s - s
// (0 when all pair sums are finite, NaN otherwise), in blocks so that a
// non-finite value stops the scan early.
static bool mayhave_nan_inf(const double* x, R_xlen_t n) {
  if ((n & 1) != 0 && !std::isfinite(x[0])) return true;
  R_xlen_t i = n & 1;
  const R_xlen_t B = 4096;
  while (i < n) {
    const R_xlen_t end = (n - i > B) ? i + B : n;
    double a0 = 0.0, a1 = 0.0, a2 = 0.0, a3 = 0.0;
    for (; i + 8 <= end; i += 8) {
      const double s0 = x[i] + x[i + 1], s1 = x[i + 2] + x[i + 3];
      const double s2 = x[i + 4] + x[i + 5], s3 = x[i + 6] + x[i + 7];
      a0 += s0 - s0; a1 += s1 - s1; a2 += s2 - s2; a3 += s3 - s3;
    }
    for (; i < end; i += 2) {
      const double s = x[i] + x[i + 1];
      a0 += s - s;
    }
    if (std::isnan(a0 + a1 + a2 + a3)) return true;
  }
  return false;
}

// TRUE when R's product code would not call BLAS for this double vector
// or matrix (it may contain NaN or Inf).
[[cpp11::register]]
bool fast_mayhave_cpp(SEXP x) {
  if (TYPEOF(x) != REALSXP) cpp11::stop("fast_mayhave_cpp: double required");
  return mayhave_nan_inf(REAL_RO(x), XLENGTH(x));
}

// TRUE when a and b are the same R object.
[[cpp11::register]]
bool fast_same_cpp(SEXP a, SEXP b) {
  return a == b;
}

// TRUE when the square matrix A equals its transpose exactly.
[[cpp11::register]]
bool fast_symmetric_cpp(SEXP A) {
  const int p = Rf_nrows(A);
  if (Rf_ncols(A) != p) return false;
  const double* a = REAL_RO(A);
  const int T = 64;
  for (int jj = 0; jj < p; jj += T) {
    const int jend = jj + T < p ? jj + T : p;
    for (int ii = 0; ii <= jj; ii += T) {
      const int iend = ii + T < p ? ii + T : p;
      for (int j = jj; j < jend; j++)
        for (int i = ii; i < iend && i < j; i++)
          if (!(a[i + static_cast<R_xlen_t>(j) * p] ==
                a[j + static_cast<R_xlen_t>(i) * p]))
            return false;
    }
  }
  return true;
}

[[cpp11::register]]
int fast_max_threads_cpp() {
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
}

// Reference dgemv 'N' with alpha = 1, beta = 0, for output rows
// [i0, i1): y(i) = 0, then y(i) = y(i) + x(j) * A(i, j) for j ascending.
// Four columns per sweep; each y(i) still receives its terms one at a
// time in ascending j.
static void ref_gemv_n_rows(const double* a, R_xlen_t m, int n,
                            const double* x, double* y,
                            R_xlen_t i0, R_xlen_t i1) {
  for (R_xlen_t i = i0; i < i1; i++) y[i] = 0.0;
  int j = 0;
  for (; j + 4 <= n; j += 4) {
    const double* a0 = a + static_cast<R_xlen_t>(j) * m;
    const double* a1 = a0 + m;
    const double* a2 = a1 + m;
    const double* a3 = a2 + m;
    const double t0 = 1.0 * x[j], t1 = 1.0 * x[j + 1];
    const double t2 = 1.0 * x[j + 2], t3 = 1.0 * x[j + 3];
    for (R_xlen_t i = i0; i < i1; i++) {
      double s = y[i];
      s = s + t0 * a0[i];
      s = s + t1 * a1[i];
      s = s + t2 * a2[i];
      s = s + t3 * a3[i];
      y[i] = s;
    }
  }
  for (; j < n; j++) {
    const double* a0 = a + static_cast<R_xlen_t>(j) * m;
    const double t0 = 1.0 * x[j];
    for (R_xlen_t i = i0; i < i1; i++) y[i] = y[i] + t0 * a0[i];
  }
}

// Reference dgemv 'T' with alpha = 1, beta = 0, for output columns
// [j0, j1): temp = 0, temp = temp + A(i, j) * x(i) for i ascending, then
// y(j) = 0 + 1 * temp. Four columns at a time as independent chains.
static void ref_gemv_t_cols(const double* a, R_xlen_t m, const double* x,
                            double* y, int j0, int j1) {
  int j = j0;
  for (; j + 4 <= j1; j += 4) {
    const double* a0 = a + static_cast<R_xlen_t>(j) * m;
    const double* a1 = a0 + m;
    const double* a2 = a1 + m;
    const double* a3 = a2 + m;
    double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
    for (R_xlen_t i = 0; i < m; i++) {
      const double xi = x[i];
      s0 = s0 + a0[i] * xi;
      s1 = s1 + a1[i] * xi;
      s2 = s2 + a2[i] * xi;
      s3 = s3 + a3[i] * xi;
    }
    y[j] = 0.0 + 1.0 * s0;
    y[j + 1] = 0.0 + 1.0 * s1;
    y[j + 2] = 0.0 + 1.0 * s2;
    y[j + 3] = 0.0 + 1.0 * s3;
  }
  for (; j < j1; j++) {
    const double* a0 = a + static_cast<R_xlen_t>(j) * m;
    double s0 = 0.0;
    for (R_xlen_t i = 0; i < m; i++) s0 = s0 + a0[i] * x[i];
    y[j] = 0.0 + 1.0 * s0;
  }
}

static void ref_gemv(bool trans, const double* a, int m, int n,
                     const double* x, double* y, int nthreads) {
  if (!trans) {
    // Row tiles; the tile size does not change any element's operations.
    R_xlen_t tile = 1024;
    if (nthreads > 1) {
      const R_xlen_t per = (m + nthreads - 1) / nthreads;
      if (per < tile) tile = per < 64 ? 64 : per;
    }
    const R_xlen_t ntile = (m + tile - 1) / tile;
    if (nthreads > 1) {
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) schedule(static)
#endif
      for (R_xlen_t t = 0; t < ntile; t++) {
        const R_xlen_t i0 = t * tile, i1 = i0 + tile < m ? i0 + tile : m;
        ref_gemv_n_rows(a, m, n, x, y, i0, i1);
      }
    } else {
      for (R_xlen_t t = 0; t < ntile; t++) {
        const R_xlen_t i0 = t * tile, i1 = i0 + tile < m ? i0 + tile : m;
        ref_gemv_n_rows(a, m, n, x, y, i0, i1);
      }
    }
  } else {
    if (nthreads > 1) {
      const int nblock = (n + 3) / 4;
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads) schedule(static)
#endif
      for (int b = 0; b < nblock; b++) {
        const int j0 = 4 * b, j1 = j0 + 4 < n ? j0 + 4 : n;
        ref_gemv_t_cols(a, m, x, y, j0, j1);
      }
    } else {
      ref_gemv_t_cols(a, m, x, y, 0, n);
    }
  }
}

// A %*% x (trans = FALSE) or t(A) %*% x (trans = TRUE) for a plain double
// matrix A, as R's product code computes it when neither operand may hold
// NaN/Inf: dgemv(trans, nrow(A), ncol(A), 1, A, nrow(A), x, 1, 0, z, 1).
// The caller has checked A (plain double matrix, both dimensions >= 2,
// no NaN/Inf). Returns NULL, so the caller uses R's own product, when x
// is not a plain double vector of the right length or may hold NaN/Inf.
// ref = TRUE uses the reference-order kernels instead of BLAS.
[[cpp11::register]]
SEXP fast_gemv_cpp(SEXP A, SEXP x, bool trans, bool ref, int nthreads) {
  const int m = Rf_nrows(A), n = Rf_ncols(A);
  const R_xlen_t len_x = trans ? m : n, len_z = trans ? n : m;
  if (TYPEOF(x) != REALSXP || OBJECT(x) ||
      Rf_getAttrib(x, R_DimSymbol) != R_NilValue || XLENGTH(x) != len_x)
    return R_NilValue;
  const double* xv = REAL_RO(x);
  if (mayhave_nan_inf(xv, len_x)) return R_NilValue;
  const double* a = REAL_RO(A);
  SEXP z = PROTECT(cpp11::safe[Rf_allocVector](REALSXP, len_z));
  double* zv = REAL(z);
  if (ref) {
    ref_gemv(trans, a, m, n, xv, zv, nthreads);
  } else {
    const char* tr = trans ? "T" : "N";
    const double one = 1.0, zero = 0.0;
    const int ione = 1;
    F77_CALL(dgemv)(tr, &m, &n, &one, a, &m, xv, &ione, &zero, zv, &ione
                    FCONE);
  }
  UNPROTECT(1);
  return z;
}
