# =============================================================================
# SUSIE-FAST: MATRIX-VECTOR PRODUCTS
#
# compute_Rv, compute_Xb and compute_Xty multiply the data matrix (XtX, or
# X) by a vector with %*%, crossprod() or tcrossprod(). For each such
# product R first scans both operands for NaN/Inf and, when neither may
# hold one, calls BLAS dgemv. The scan of the matrix costs about as much as
# the product, and the matrix does not change during a fit.
#
#   1. The matrix is scanned once per fit (the same pairwise test as R's
#      mayHaveNaNOrInf). Each product then scans only its vector and calls
#      the same dgemv with the same arguments. Whenever R would not call
#      BLAS (options(matprod) other than "default", a vector that may hold
#      NaN/Inf, a matrix that is not a plain double matrix), the upstream
#      expression is used.
#   2. With R's reference BLAS, C++ kernels that perform the reference
#      dgemv's operations in the same order for every output element are
#      used instead (several columns per pass, independent sums). For an
#      exactly symmetric XtX, XtX %*% v is computed as t(XtX) %*% v, whose
#      column sums read memory in order. The kernels are used only after a
#      session self-test against R's own products, and only with Rblas.
#   3. The first product of each kind in each fit is also computed with the
#      upstream expression and compared bit for bit; on any difference the
#      fit goes back to the upstream expressions.
#   4. initialize_fitted: without model_init the initial fitted values are
#      a product with an all-zero vector, which is +0 in every entry when
#      the matrix is finite; it is not computed.
#
# options(susieR.threads = n): number of threads for the kernels in 2
# (default 1). Results are bit-identical for every n, because each output
# element is computed by one thread in the same order. Keep 1 when fits
# already run in parallel (e.g. mclapply over regions). Other BLAS
# libraries (OpenBLAS, MKL) use their own threading and are not affected.
# =============================================================================

# Per-fit state of the matrix that compute_Rv / compute_Xb / compute_Xty
# multiply, or NULL when the fast products do not apply. Created by
# susie_workhorse before ibss_initialize and kept in the fast cache; it
# holds a reference to the matrix, which therefore cannot be modified in
# place while the state exists.
#' @keywords internal
fast_products_new <- function(data) {
  if (fast_mode() == "off" || !class(data)[1] %in% c("ss", "individual"))
    return(NULL)
  ind <- class(data)[1] == "individual"
  M <- if (ind || !is.null(data$X)) data$X else data$XtX
  if (!is.matrix(M) || !is.double(M) || is.object(M) ||
      !is.null(attr(M, "matrix.type")) || any(dim(M) < 2))
    return(NULL)
  zero_ok <- TRUE
  if (ind) {
    cm  <- attr(M, "scaled:center")
    csd <- attr(M, "scaled:scale")
    if (!is.double(cm) || !is.double(csd) || length(cm) != ncol(M) ||
        length(csd) != ncol(M) || !is.null(dim(cm)) || !is.null(dim(csd)))
      return(NULL)
    zero_ok <- isTRUE(data$p == ncol(M)) &&
      all(is.finite(csd) & csd != 0) && all(is.finite(cm))
  }
  if (fast_mayhave_cpp(M) || !fast_products_ok())
    return(NULL)
  st <- new.env(parent = emptyenv())
  st$M       <- M
  st$zero_ok <- zero_ok
  st$threads <- fast_threads()
  st$ref     <- fast_refblas_ok(st$threads)
  st$sym     <- st$ref && !ind && is.null(data$X) && fast_symmetric_cpp(M)
  st$ok      <- TRUE
  st$checked <- character(0)
  st
}

# Adds the product state to the fast cache. A matrix without NaN/Inf is
# finite, which is what the cache's zero_ok test establishes, so the test
# is not repeated.
#' @keywords internal
fast_products_attach <- function(cache, st) {
  if (!is.null(cache) && !is.null(st)) {
    cache$products <- st
    if (st$zero_ok) cache$zero_ok <- TRUE
  }
  cache
}

#' @keywords internal
fast_threads <- function() {
  nt <- suppressWarnings(as.integer(getOption("susieR.threads", 1L))[1])
  if (is.na(nt) || nt <= 1L) return(1L)
  min(nt, fast_max_threads_cpp())
}

#' @keywords internal
fast_products_active <- function(st)
  !is.null(st) && st$ok && identical(getOption("matprod"), "default") &&
    fast_mode() != "off"

# Returns the fast result z after the first product of each kind in each
# fit has been checked against the upstream result (computed by upstream()
# and returned in its place while unchecked). A check counts only for a
# nonzero input vector.
#' @keywords internal
fast_checked <- function(st, kind, input, z, upstream) {
  if (kind %in% st$checked) return(z)
  ref <- upstream()
  if (!identical(z, ref, num.eq = FALSE)) st$ok <- FALSE
  else if (any(input != 0)) st$checked <- c(st$checked, kind)
  ref
}

# compute_Rv(data, v, Rv_matrix) with the fit's product state.
#' @keywords internal
fast_Rv <- function(model, data, v, Rv_matrix = NULL)
  fast_Rv_with(model$runtime$fast_cache$products, data, v, Rv_matrix)

#' @keywords internal
fast_Rv_with <- function(st, data, v, Rv_matrix = NULL) {
  z <- NULL
  if (is.null(Rv_matrix) && fast_products_active(st)) {
    if (!is.null(data$X)) {
      if (fast_same_cpp(st$M, data$X)) {
        Xv <- fast_gemv_cpp(st$M, v, FALSE, st$ref, st$threads)
        if (!is.null(Xv))
          z <- fast_gemv_cpp(st$M, Xv, TRUE, st$ref, st$threads)
      }
    } else if (fast_same_cpp(st$M, data$XtX)) {
      z <- fast_gemv_cpp(st$M, v, st$sym, st$ref, st$threads)
    }
  }
  if (is.null(z)) return(compute_Rv(data, v, Rv_matrix))
  fast_checked(st, "Rv", v, z, function() compute_Rv(data, v))
}

# compute_Xb(X, b) with the fit's product state.
#' @keywords internal
fast_Xb <- function(model, X, b)
  fast_Xb_with(model$runtime$fast_cache$products, X, b)

#' @keywords internal
fast_Xb_with <- function(st, X, b) {
  if (fast_products_active(st) && fast_same_cpp(st$M, X) && is.double(b) &&
      is.null(dim(b)) && length(b) == ncol(X)) {
    cm  <- attr(X, "scaled:center")
    csd <- attr(X, "scaled:scale")
    z <- fast_gemv_cpp(X, b / csd, FALSE, st$ref, st$threads)
    if (!is.null(z))
      return(fast_checked(st, "Xb", b, as.numeric(z - sum(cm * b / csd)),
                          function() compute_Xb(X, b)))
  }
  compute_Xb(X, b)
}

# compute_Xty(X, y) with the fit's product state.
#' @keywords internal
fast_Xty <- function(model, X, y)
  fast_Xty_with(model$runtime$fast_cache$products, X, y)

#' @keywords internal
fast_Xty_with <- function(st, X, y) {
  if (fast_products_active(st) && fast_same_cpp(st$M, X)) {
    z <- fast_gemv_cpp(X, y, TRUE, st$ref, st$threads)
    if (!is.null(z)) {
      cm  <- attr(X, "scaled:center")
      csd <- attr(X, "scaled:scale")
      return(fast_checked(st, "Xty", y, as.numeric(z / csd - cm / csd * sum(y)),
                          function() compute_Xty(X, y)))
    }
  }
  compute_Xty(X, y)
}

# initialize_fitted.ss / .individual. The product state is found in the
# calling susie_workhorse frame (initialize_fitted has no model argument).
# Without model_init, b = colSums(alpha * mu) is all zero, and the product
# of a finite matrix with a zero vector is +0 in every entry (BLAS and R's
# own loops both sum from +0), so it is not computed.
#' @keywords internal
fast_init_Rv <- function(data, b) {
  st <- fast_products_find()
  if (!is.null(st) && class(data)[1] == "ss" && is.null(data$X) &&
      fast_mode() != "off" && fast_same_cpp(st$M, data$XtX) &&
      fast_zero_vector(b, ncol(st$M)))
    return(numeric(length(b)))
  fast_Rv_with(st, data, b)
}

#' @keywords internal
fast_init_Xb <- function(data, b) {
  st <- fast_products_find()
  if (!is.null(st) && class(data)[1] == "individual" && st$zero_ok &&
      fast_mode() != "off" && fast_same_cpp(st$M, data$X) &&
      fast_zero_vector(b, ncol(st$M)))
    return(numeric(nrow(st$M)))
  fast_Xb_with(st, data$X, b)
}

#' @keywords internal
fast_products_find <- function() {
  st <- dynGet("fast_products", ifnotfound = NULL)
  if (is.environment(st) && !is.null(st$M)) st else NULL
}

#' @keywords internal
fast_zero_vector <- function(b, p)
  is.double(b) && is.null(dim(b)) && length(b) == p && !anyNA(b) &&
    !any(b != 0)

# -----------------------------------------------------------------------------
# Session self-tests
# -----------------------------------------------------------------------------

# The direct dgemv calls and the NaN/Inf scan are used only after they
# have reproduced R's products and R's scan on this platform (once per
# session).
#' @keywords internal
fast_products_ok <- local({
  ok <- NA
  function() {
    if (is.na(ok))
      ok <<- isTRUE(tryCatch(fast_products_self_test(FALSE, 1L),
                             error = function(e) FALSE))
    ok
  }
})

# The reference-order kernels are used only with R's reference BLAS, and
# only after reproducing R's products at this thread count (once per
# session, BLAS library and thread count).
#' @keywords internal
fast_refblas_ok <- local({
  ok <- list()
  function(nthreads) {
    blas <- extSoftVersion()[["BLAS"]]
    if (!(blas == "" || grepl("^(lib)?Rblas", basename(blas)))) return(FALSE)
    key <- paste(blas, nthreads)
    if (is.null(ok[[key]]))
      ok[[key]] <<- isTRUE(tryCatch(fast_products_self_test(TRUE, nthreads),
                                    error = function(e) FALSE))
    ok[[key]]
  }
})

# Compares the scan with an R version of R's mayHaveNaNOrInf on edge cases,
# and the products (direct dgemv when ref = FALSE, the kernels when
# ref = TRUE) with R's %*%, crossprod(y, X), tcrossprod(X, t(b)) and
# crossprod(X, X %*% v), called as the package calls them.
#' @keywords internal
fast_products_self_test <- function(ref, nthreads) {
  old <- options(matprod = "default")
  st <- if (exists(".Random.seed", globalenv(), inherits = FALSE))
          get(".Random.seed", globalenv()) else NULL
  on.exit({
    options(old)
    if (is.null(st)) suppressWarnings(rm(".Random.seed", envir = globalenv()))
    else assign(".Random.seed", st, envir = globalenv())
  })
  set.seed(20261007)

  mayhave <- function(x) {
    n <- length(x)
    if (n %% 2 == 1 && !is.finite(x[1])) return(TRUE)
    if (n < 2) return(FALSE)
    i <- seq(n %% 2 + 1, n - 1, by = 2)
    any(!is.finite(x[i] + x[i + 1]))
  }
  big <- rnorm(10001)
  vs <- list(numeric(0), 5, NaN, Inf, c(1, 2), c(NaN, 1), c(1, NA), c(-Inf, 1),
             c(Inf, -Inf), c(1e308, 1e308), c(-1e308, -1e308), c(1e308, -1e308),
             c(Inf, 1, 2), c(1, Inf, 2), c(1e308, 1e308, 1), c(1, 1e308, 1e308),
             c(-0, 0), c(-0, -0, -0), big, replace(big, 10001, Inf),
             replace(big, 1, NaN), replace(big, 8193, NA),
             replace(big[-1], 4097, -Inf), replace(big[-1], 9999:10000, 1e308),
             matrix(rnorm(12), 3, 4), matrix(c(rnorm(11), Inf), 3, 4))
  for (x in vs)
    if (!identical(fast_mayhave_cpp(x), mayhave(x))) return(FALSE)

  # Odd shapes, and shapes that are not multiples of the kernels' 4 columns
  # or 1024-row tiles. Products at the size of the data are checked in each
  # fit (fast_checked).
  shapes <- list(c(2, 2), c(3, 2), c(2, 5), c(7, 3), c(5, 7), c(50, 40),
                 c(101, 37), c(131, 131), c(1025, 9), c(9, 1031), c(2050, 6))
  for (s in shapes) {
    n <- s[1]; p <- s[2]
    X <- matrix(runif(n * p, -1, 1), n, p) * 10^runif(n, -3, 3)
    X <- X * rep(10^runif(p, -3, 3), each = n)
    X[sample(n * p, n)] <- 0
    v <- rnorm(p) * 10^runif(p, -8, 8); v[sample(p, 1)] <- 0; v[sample(p, 1)] <- -0
    y <- rnorm(n) * 10^runif(n, -8, 8); y[sample(n, 1)] <- -0
    if (!identical(fast_gemv_cpp(X, v, FALSE, ref, nthreads),
                   as.vector(tcrossprod(X, t(v))), num.eq = FALSE) ||
        !identical(fast_gemv_cpp(X, y, TRUE, ref, nthreads),
                   as.vector(crossprod(y, X)), num.eq = FALSE) ||
        !identical(fast_gemv_cpp(X, fast_gemv_cpp(X, v, FALSE, ref, nthreads),
                                 TRUE, ref, nthreads),
                   as.vector(crossprod(X, X %*% v)), num.eq = FALSE))
      return(FALSE)
    if (n == p) {
      S <- X + t(X)
      if (!identical(fast_gemv_cpp(X, v, FALSE, ref, nthreads),
                     as.vector(X %*% v), num.eq = FALSE) ||
          !identical(fast_gemv_cpp(S, v, ref, ref, nthreads),
                     as.vector(S %*% v), num.eq = FALSE) ||
          !fast_symmetric_cpp(S) || fast_symmetric_cpp(X))
        return(FALSE)
    }
  }
  # Vectors the products must refuse (R would not call BLAS for them).
  X <- matrix(rnorm(20), 5, 4)
  if (!is.null(fast_gemv_cpp(X, c(1, NaN, 2, 3), FALSE, ref, nthreads)) ||
      !is.null(fast_gemv_cpp(X, c(1, 2, 1e308, 1e308), FALSE, ref, nthreads)) ||
      !is.null(fast_gemv_cpp(X, 1:4 + 0.5, TRUE, ref, nthreads)) ||
      !is.null(fast_gemv_cpp(X, matrix(1, 4, 1), FALSE, ref, nthreads)))
    return(FALSE)
  TRUE
}
