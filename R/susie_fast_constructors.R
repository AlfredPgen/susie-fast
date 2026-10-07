# =============================================================================
# SUSIE-FAST: data constructors
#
# Faster input checks and scaling in the constructors. The data objects,
# warnings and errors are the same as upstream, in the same order.
#
#   1. sufficient_stats_constructor: one tiled pass checks XtX for symmetry
#      and for non-finite entries, instead of Rfast::is.symmetric (which
#      reads half the matrix row by row) and anyNA. When XtX is bitwise
#      symmetric, t((1 / csd) * XtX) / csd reads it column by column, and
#      when every csd is exactly 1 (diag(R) == 1 in susie_rss) it returns
#      XtX itself: (1 * x) / 1 == x for every non-NaN double.
#   2. summary_stats_constructor, original scale (bhat, shat, var_y):
#      t(R * s) * s and (XtX + t(XtX)) / 2 in one pass instead of four
#      full-size temporaries.
#   3. individual_data_constructor: the constant-column warning calls var()
#      only on the columns whose range does not already prove var() > 0.
#
# Each path is used only for a plain double matrix, outside
# options(susieR.fast = "off"), and after its kernel has reproduced the R
# code on this platform (fast_constructors_ok); otherwise the upstream
# expression runs.
# =============================================================================

# Symmetry and finiteness of XtX in one pass. Returns NULL (use
# is_symmetric_matrix() and anyNA() as upstream) or list(sym, bitsym) for a
# matrix whose entries are all finite: sym is is_symmetric_matrix(XtX),
# bitsym says t(XtX) has the same data as XtX.
#' @keywords internal
fast_xtx_check <- function(XtX) {
  if (fast_mode() == "off" || !is.matrix(XtX) || !is.double(XtX) ||
      is.object(XtX) || nrow(XtX) != ncol(XtX) || nrow(XtX) < 1 ||
      !fast_kernel_ok() || !fast_constructors_ok())
    return(NULL)
  code <- xtx_check_cpp(XtX)
  if (code == 0L) return(NULL)
  list(sym = code >= 2L, bitsym = code == 3L)
}

# t((1 / csd) * XtX) / csd, given the flags of fast_xtx_check() for this
# XtX (NULL when unknown).
#' @keywords internal
fast_scale_xtx_sym <- function(XtX, csd, fl) {
  if (is.null(fl) || !fl$bitsym || fast_mode() == "off" ||
      !is.matrix(XtX) || !is.double(XtX) || is.object(XtX) ||
      nrow(XtX) != ncol(XtX) || length(csd) != nrow(XtX) ||
      !is.double(csd) || !fast_plain_attributes(XtX) ||
      !fast_constructors_ok())
    return(fast_scale_xtx(XtX, csd))
  fast_scale_xtx_sym_kernel(XtX, csd)
}

#' @keywords internal
fast_scale_xtx_sym_kernel <- function(XtX, csd) {
  dn <- dimnames(XtX)
  if (!anyNA(csd) && all(csd == 1) && (is.null(dn) || identical(dn, dn[2:1])))
    return(XtX)
  at <- attributes(XtX)
  if (!is.null(at$dimnames)) at$dimnames <- at$dimnames[2:1]
  out <- scale_xtx_sym_cpp(XtX, csd)
  attributes(out) <- at
  out
}

# Original-scale XtX: same values and attributes as
# XtX <- t(R * sqrt(XtXdiag)) * sqrt(XtXdiag); (XtX + t(XtX)) / 2.
#' @keywords internal
fast_orig_scale_xtx <- function(R, XtXdiag) {
  if (fast_mode() == "off" || !is.matrix(R) || !is.double(R) ||
      is.object(R) || nrow(R) != ncol(R) || length(XtXdiag) != nrow(R) ||
      !is.double(XtXdiag) || !fast_plain_attributes(R) ||
      !all(names(attributes(XtXdiag)) == "names") ||
      !fast_constructors_ok()) {
    XtX <- t(R * sqrt(XtXdiag)) * sqrt(XtXdiag)
    return((XtX + t(XtX)) / 2)
  }
  fast_orig_scale_kernel(R, sqrt(XtXdiag))
}

#' @keywords internal
fast_orig_scale_kernel <- function(R, s) {
  at <- attributes(R)
  if (!is.null(at$dimnames)) at$dimnames <- at$dimnames[2:1]
  out <- orig_scale_xtx_cpp(R, s)
  attributes(out) <- at
  out
}

# Indices of the columns of X with var() equal to 0 or NA, as
# which(col_vars == 0 | is.na(col_vars)) with col_vars <- apply(X, 2, var).
#' @keywords internal
fast_const_cols <- function(X) {
  if (fast_mode() == "off" || !is.matrix(X) || !is.double(X) ||
      is.object(X) || !is.null(attr(X, "matrix.type")) ||
      nrow(X) < 2 || ncol(X) < 1 || !fast_constructors_ok()) {
    col_vars <- apply(X, 2, var)
    return(which(col_vars == 0 | is.na(col_vars)))
  }
  fast_const_cols_screened(X)
}

#' @keywords internal
fast_const_cols_screened <- function(X) {
  flags <- logical(ncol(X))
  for (j in which(const_col_screen_cpp(X))) {
    v <- var(X[, j])
    flags[j] <- v == 0 | is.na(v)
  }
  which(flags)
}

# Arithmetic drops a matrix's names attribute and treats tsp specially, so
# the kernels, which copy the attributes, are used only without them.
#' @keywords internal
fast_plain_attributes <- function(M)
  is.null(attr(M, "names")) && is.null(attr(M, "tsp"))

# The constructor kernels are used only after they have reproduced the R
# code bit for bit on this platform (checked once per session).
#' @keywords internal
fast_constructors_ok <- local({
  ok <- NA
  function() {
    if (is.na(ok))
      ok <<- isTRUE(tryCatch(fast_constructors_self_test(),
                             error = function(e) FALSE))
    ok
  }
})

#' @keywords internal
fast_constructors_self_test <- function() {
  # Deterministic test matrices (no use of the RNG).
  wave <- function(n, p, k) matrix(sin(seq_len(n * p) * k + k), n, p)
  ref_code <- function(M) {
    if (!all(is.finite(M))) return(0L)
    if (!is_symmetric_matrix(M)) return(1L)
    if (identical(unname(M), unname(t(M)), num.eq = FALSE)) 3L else 2L
  }
  bit_same <- function(a, b) identical(a, b, num.eq = FALSE)

  # 1. Symmetry / finiteness check.
  checks <- list()
  for (p in c(1, 2, 5, 70, 150)) {
    M <- crossprod(wave(p + 3, p, 0.7))
    checks <- c(checks, list(M))
    if (p >= 2) {
      A <- M; A[1, p] <- A[1, p] + 1; checks <- c(checks, list(A))
      Z <- M; Z[p, 1] <- 0; Z[1, p] <- -0; checks <- c(checks, list(Z))
      for (bad in list(NaN, NA_real_, Inf, -Inf)) {
        B <- M; B[p, 1] <- bad; checks <- c(checks, list(B))
        B <- M; B[p, 1] <- bad; B[1, p] <- bad; checks <- c(checks, list(B))
        B <- A; B[p, p] <- bad; checks <- c(checks, list(B))
      }
    }
    D <- M; D[p, p] <- Inf; checks <- c(checks, list(D))
  }
  for (M in checks)
    if (!identical(xtx_check_cpp(M), ref_code(M))) return(FALSE)

  # 2. Column-by-column scaling of a bitwise symmetric matrix, and the
  #    identity case.
  for (p in c(1, 3, 130)) {
    M <- crossprod(wave(p + 5, p, 1.3))
    if (!identical(xtx_check_cpp(M), 3L)) return(FALSE)
    for (csd in list(sqrt(diag(M) / 9), rep(1, p),
                     c(NaN, sqrt(diag(M) / 9))[seq_len(p)])) {
      if (!bit_same(scale_xtx_sym_cpp(M, csd), as.vector(t((1 / csd) * M) / csd)))
        return(FALSE)
    }
    dimnames(M) <- list(paste0("v", seq_len(p)), paste0("v", seq_len(p)))
    attr(M, "extra") <- "x"
    if (!identical(fast_scale_xtx_sym_kernel(M, rep(1, p)),
                   t((1 / rep(1, p)) * M) / rep(1, p)))
      return(FALSE)
    dimnames(M) <- list(NULL, paste0("v", seq_len(p)))
    if (!identical(fast_scale_xtx_sym_kernel(M, rep(1, p)),
                   t((1 / rep(1, p)) * M) / rep(1, p)))
      return(FALSE)
  }

  # 3. Original-scale XtX, including a non-symmetric R, signed zeros,
  #    infinities, NaN and a near-overflow diagonal.
  for (p in c(1, 2, 67, 140)) {
    R <- wave(p, p, 0.9)
    s <- sqrt(seq_len(p) * 0.37 + 0.5)
    names(s) <- paste0("s", seq_len(p))
    if (p >= 2) {
      R[1, 2] <- -0; R[2, 1] <- 0
      R[p, 1] <- Inf; R[2, 2] <- NaN
      R[1, 1] <- 1e308
    }
    dimnames(R) <- list(paste0("r", seq_len(p)), NULL)
    attr(R, "eigen") <- list(values = 1)
    up <- t(R * s) * s
    up <- (up + t(up)) / 2
    if (!bit_same(fast_orig_scale_kernel(R, s), up)) return(FALSE)
  }

  # 4. Constant-column screen.
  n <- 40
  cols <- list(sin(1:n), rep(2, n), rep(0, n), c(1e-170, rep(0, n - 1)),
               rep(1 / 3, n), c(Inf, rep(1, n - 1)), rep(-0, n),
               c(0, 1e-100, rep(0, n - 2)), c(-1e150, 1e150, rep(0, n - 2)),
               c(1, 1 + 2^-52, rep(1, n - 2)), c(1e160, rep(1e160, n - 1)),
               c(rep(5, n - 1), 5 + 1e-14), cos(1:n) * 1e-120)
  X <- do.call(cbind, cols)
  col_vars <- apply(X, 2, var)
  if (!identical(fast_const_cols_screened(X),
                 unname(which(col_vars == 0 | is.na(col_vars)))))
    return(FALSE)
  TRUE
}
