# =============================================================================
# SUSIE-FAST: post-processing (credible sets, z-scores, refinement)
#
#   1. Purity early exit. susie_get_cs keeps a credible set only if its
#      min |corr| reaches min_abs_corr. A diffuse effect (V > 0, flat alpha)
#      gives a CS of most variables, whose full purity costs O(m^2) (O(n m^2)
#      from X) and is then thrown away. A few rows of the correlation matrix,
#      O(m) (O(n m)), usually show a pair below the threshold; the CS is then
#      certified as dropped and its purity row, which never reaches the
#      output, is not computed. In every other case get_purity runs as before.
#   2. z-score memo. With compute_univariate_zscore = TRUE every workhorse
#      fit (main fit, refine candidates, greedy rounds) recomputes the same
#      univariate z-scores from the same X and y. They are computed once per
#      top-level fit.
#   3. Parallel refine candidates (opt-in, options(susieR.refine_cores = k),
#      Linux/macOS). The candidate fits of one refine step are independent;
#      they run in forked processes and are collected in the serial order.
#
# All three apply only when fast_mode() != "off".
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Purity early exit
# -----------------------------------------------------------------------------

# Returns c(v, v, v) with v < min_abs_corr when the credible set pos is
# certain to be dropped by susie_get_cs (its purity row is then discarded),
# or NULL when get_purity must be called. The certificate holds only when
# median_abs_corr is NULL and squared is FALSE: the keep rule is then
# purity[, 1] >= min_abs_corr, and purity[, 1] is the min over all pairs.
# If it draws the X-path subsample and returns NULL, it restores the RNG
# state first, so get_purity draws the same subsample.
#' @keywords internal
fast_purity_drop <- function(pos, X, Xcorr, squared, n, min_abs_corr,
                             median_abs_corr) {
  if (length(pos) < 2 || !isFALSE(squared) || !is.null(median_abs_corr) ||
      !is.numeric(min_abs_corr) || length(min_abs_corr) != 1 ||
      !isTRUE(min_abs_corr > 0) || fast_mode() == "off")
    return(NULL)
  if (!is.numeric(pos) || anyNA(pos) || any(pos < 1) || any(pos != trunc(pos)))
    return(NULL)
  if (!is.null(Xcorr))
    return(fast_purity_drop_corr(as.integer(pos), Xcorr, min_abs_corr))

  # X path: draw the subsample exactly as get_purity does.
  if (!(identical(n, "auto") ||
        (is.numeric(n) && length(n) == 1 && !is.na(n))))
    return(NULL)
  if (!((is.matrix(X) && is.double(X) && !is.object(X)) ||
        identical(as.character(class(X)), "dgCMatrix")) ||
      max(pos) > ncol(X) || !fast_cora_ok())
    return(NULL)
  n <- resolve_n_purity(n, nrow(X), length(pos))
  seed <- NULL
  if (length(pos) > n) {
    seed <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)
    if (is.null(seed) || fast_rng_user(seed))
      return(NULL)
    pos <- sample(pos, n)
  }
  v <- fast_purity_drop_X(pos, X, min_abs_corr)
  if (is.null(v) && !is.null(seed))
    assign(".Random.seed", seed, envir = globalenv())
  v
}

# Xcorr path (a correlation matrix, or a scaled_XtX object). Each probed
# value is the element of the correlation block that get_purity computes,
# by the same IEEE expression, so upstream's min |corr| is at most v.
#' @keywords internal
fast_purity_drop_corr <- function(pos, Xcorr, thr) {
  scaled <- identical(class(Xcorr), "scaled_XtX")
  M <- if (scaled) Xcorr$XtX else if (is.null(oldClass(Xcorr))) Xcorr
  if (!is.matrix(M) || !is.double(M) || isS4(M) ||
      max(pos) > min(dim(M)))
    return(NULL)
  m <- length(pos)
  if (scaled) {
    # As safe_cov2cor(M[pos, pos]): R[a, b] = (M[pos[b], pos[a]] * dinv[b]) * dinv[a]
    dg <- M[cbind(pos, pos)]
    if (anyNA(dg) || any(dg < 0)) return(NULL)
    d <- sqrt(dg)
    dinv <- 1 / d
    dinv[d == 0] <- 0
  }
  for (r in unique(c(1L, m %/% 2L, m))) {
    # pairs (lo, hi), lo < hi, of upper_tri() that involve row r
    o  <- seq_len(m)[-r]
    lo <- pmin(r, o)
    hi <- pmax(r, o)
    v <- if (scaled) abs((M[cbind(pos[hi], pos[lo])] * dinv[hi]) * dinv[lo])
         else abs(M[cbind(pos[lo], pos[hi])])
    k <- which(v < thr)[1]
    if (!is.na(k)) {
      # get_purity stops on a NaN anywhere in the block. With every entry
      # finite (and, for a zero diagonal, no overflow before the 0 factor)
      # no element can be NaN, so upstream returns, with min <= v[k].
      mx <- block_maxabs_cpp(M, pos)
      if (!is.finite(mx)) return(NULL)
      if (scaled && any(d == 0) && !is.finite(mx * max(dinv))) return(NULL)
      return(rep(v[k], 3))
    }
  }
  NULL
}

# X path. get_purity computes safe_cor(X_sub), which is Rfast::cora(X_sub)
# when no column is constant: mat <- t(x) - colmeans(x); mat <- mat /
# sqrt(rowsums(mat^2)); tcrossprod(mat). mat is rebuilt here by the same
# lines (fast_cora_ok checks cora's body), and the entries (r, j) of
# tcrossprod(mat) are estimated with an error bound that holds for any
# summation order: |fl(x.y) - x.y| <= gam_n sum|x_k y_k| (Higham, Thm 3.4),
# for upstream's value and ours, so |b - dd| <= 2 gam_n A / (1 - gam_n).
# The bound is doubled again to cover the rounding of the test itself.
#' @keywords internal
fast_purity_drop_X <- function(pos, X, thr) {
  X_sub <- X[, pos]
  X_sub <- as.matrix(X_sub)
  if (!is.double(X_sub)) return(NULL)
  # safe_cor's branch: cora is used only when no column is constant (NA:
  # safe_cor stops, so get_purity must run)
  n <- nrow(X_sub)
  cm <- colMeans(X_sub)
  css <- colSums(X_sub^2) - n * cm^2
  if (!identical(any(css == 0), FALSE)) return(NULL)
  mat <- base::t(X_sub) - Rfast::colmeans(X_sub)
  mat <- mat / sqrt(Rfast::rowsums(mat^2))
  if (!all(is.finite(mat))) return(NULL)
  u <- .Machine$double.eps / 2
  gam <- n * u / (1 - n * u)
  amat <- abs(mat)
  m <- nrow(mat)
  for (r in unique(c(1L, m %/% 2L, m))) {
    dd <- drop(mat %*% mat[r, ])
    A  <- drop(amat %*% amat[r, ])
    bound <- 2 * gam * A / (1 - gam) * 2 + 1e-300
    below <- abs(dd) + bound < thr
    below[r] <- FALSE
    k <- which(below)[1]
    if (!is.na(k)) return(rep(abs(dd[k]), 3))
  }
  NULL
}

# TRUE when .Random.seed says the RNG is "user-supplied": its state is then
# not (all) in .Random.seed. Read from the seed's kind code (RNGkind()
# would call into the RNG to answer).
#' @keywords internal
fast_rng_user <- function(seed) {
  !is.null(seed) && (!is.integer(seed) || length(seed) < 1 || is.na(seed[1]) ||
                     seed[1] %% 100L == 5L)
}

# TRUE when Rfast::cora is the two-line R function the X-path certificate
# rebuilds, with base t, sqrt, arithmetic and tcrossprod. Checked once per
# session; any other Rfast version keeps the upstream path.
#' @keywords internal
fast_cora_ok <- local({
  ok <- NA
  function() {
    if (is.na(ok)) ok <<- isTRUE(tryCatch({
      f <- Rfast::cora
      ns <- environment(f)
      same <- function(nm, fun) identical(get(nm, envir = ns), fun)
      identical(deparse(body(f)), c(
        "{",
        "    mat <- t(x) - Rfast::colmeans(x)",
        "    mat <- mat/sqrt(Rfast::rowsums(mat^2))",
        "    if (large) {",
        "        Rfast::Tcrossprod(mat, mat)",
        "    }",
        "    else tcrossprod(mat)",
        "}")) &&
        identical(formals(f), as.pairlist(alist(x = , large = FALSE))) &&
        same("t", base::t) && same("sqrt", base::sqrt) &&
        same("-", base::`-`) && same("/", base::`/`) && same("^", base::`^`) &&
        same("tcrossprod", base::tcrossprod)
    }, error = function(e) FALSE))
    ok
  }
})

# -----------------------------------------------------------------------------
# 2. z-score memo
# -----------------------------------------------------------------------------

# Called at the top of susie_workhorse. The outermost call attaches a memo
# environment to its (local) data object; nested workhorse calls (refine
# candidates, greedy rounds) receive that data object and share the memo.
#' @keywords internal
fast_zscore_memo_attach <- function(data, params) {
  if (!identical(class(data), "individual") ||
      !isTRUE(params$compute_univariate_zscore) ||
      !is.null(data$fast_zscore_memo) || fast_mode() == "off")
    return(data)
  data$fast_zscore_memo <- new.env(parent = emptyenv())
  data
}

# Returns the z-scores z (a promise for the calc_z call in
# get_zscore.individual). A stored vector is returned only for the same X,
# y, number of columns and centring/scaling flags. A result is stored only
# when calc_z signalled no condition and set no error message (its try()
# fallback), so a later call that would emit something still recomputes.
#' @keywords internal
fast_zscore_memo <- function(data, params, X, z) {
  memo <- data$fast_zscore_memo
  if (is.null(memo) || !is.environment(memo))
    return(z)
  key <- list(ncol(X), params$intercept, params$standardize)
  if (!is.null(memo$z) && identical(memo$key, key) &&
      identical(memo$X, data$X) &&
      identical(memo$y, data$y, num.eq = FALSE)) {
    memo$hits <- memo$hits + 1L
    return(memo$z)
  }
  signalled <- FALSE
  err <- geterrmessage()
  z <- withCallingHandlers(z, condition = function(cond) signalled <<- TRUE)
  if (!signalled && identical(geterrmessage(), err)) {
    memo$key  <- key
    memo$X    <- data$X
    memo$y    <- data$y
    memo$z    <- z
    memo$hits <- 0L
  }
  z
}

# -----------------------------------------------------------------------------
# 3. Parallel refine candidates
# -----------------------------------------------------------------------------

# One refine candidate: the body of the loop in run_refine's cs_refine_step
# (step 1 with the CS's prior weights zeroed, step 2 warm-started from it).
#' @keywords internal
fast_refine_candidate <- function(model, data, params, pw_s, cs_idx) {
  pw_cs <- pw_s
  pw_cs[model$sets$cs[[cs_idx]]] <- 0

  p1 <- params
  p1$prior_weights <- reconstruct_full_weights(pw_cs, model$null_weight)
  p1$null_weight   <- model$null_weight
  p1$model_init    <- NULL
  p1$verbose       <- FALSE
  p1$track_fit     <- FALSE
  p1$refine        <- FALSE
  m1 <- susie_workhorse(data, p1)

  init <- list(alpha = m1$alpha, mu = m1$mu, mu2 = m1$mu2)
  class(init) <- "susie"
  p2 <- params
  p2$prior_weights <- reconstruct_full_weights(pw_s, model$null_weight)
  p2$null_weight   <- model$null_weight
  p2$model_init    <- init
  p2$verbose       <- FALSE
  p2$track_fit     <- FALSE
  p2$refine        <- FALSE
  susie_workhorse(data, p2)
}

# The candidate list of one refine step, computed in forked processes, or
# NULL when the serial loop must run. Used only when
# options(susieR.refine_cores) >= 2, on a platform with fork(), for the
# package's own data classes, outside a forked worker, and with an RNG
# whose whole state is in .Random.seed.
#
# Each candidate is a deterministic function of (data, params, CS, prior
# weights). A child inherits the parent's .Random.seed (mc.set.seed =
# FALSE) and records whether its fit used the RNG (only get_purity's
# subsample does); if any did, or any child failed, the results are
# discarded and the serial loop runs, from the parent's untouched RNG
# state. Messages and warnings are recorded in the children and replayed
# here in the serial order. The BLAS must be fork-safe (single-threaded,
# or OpenBLAS with pthreads); Apple's Accelerate is excluded.
#' @keywords internal
fast_refine_parallel <- function(model, data, params, pw_s) {
  cores <- getOption("susieR.refine_cores", 1)
  if (!is.numeric(cores) || length(cores) != 1 || is.na(cores) ||
      cores < 2 || .Platform$OS.type != "unix" ||
      !class(data)[1] %in% c("ss", "individual", "rss_lambda") ||
      length(class(data)) != 1 ||
      fast_rng_user(get0(".Random.seed", envir = globalenv(), inherits = FALSE)) ||
      fast_mode() == "off" ||
      grepl("Accelerate|vecLib", extSoftVersion()[["BLAS"]]) ||
      !requireNamespace("parallel", quietly = TRUE) ||
      isTRUE(get("isChild", envir = asNamespace("parallel"))()))
    return(NULL)

  # The CS indices the serial loop visits (it stops at the first CS whose
  # removal leaves no prior weight).
  idx <- integer(0)
  for (cs_idx in seq_along(model$sets$cs)) {
    pw_cs <- pw_s
    pw_cs[model$sets$cs[[cs_idx]]] <- 0
    all_zero <- all(pw_cs == 0)
    if (is.na(all_zero)) return(NULL)
    if (all_zero) break
    idx <- c(idx, cs_idx)
  }
  if (length(idx) < 2) return(NULL)

  worker <- function(cs_idx) {
    seed  <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)
    conds <- list()
    fit <- tryCatch(
      withCallingHandlers(
        fast_refine_candidate(model, data, params, pw_s, cs_idx),
        message = function(cond) {
          conds[[length(conds) + 1L]] <<- cond
          invokeRestart("muffleMessage")
        },
        warning = function(cond) {
          conds[[length(conds) + 1L]] <<- cond
          invokeRestart("muffleWarning")
        }),
      error = function(e) NULL)
    list(fit = fit, conds = conds,
         rng = !identical(seed, get0(".Random.seed", envir = globalenv(),
                                     inherits = FALSE)))
  }
  res <- tryCatch(
    suppressWarnings(parallel::mclapply(idx, worker,
                                        mc.cores = as.integer(cores),
                                        mc.preschedule = FALSE,
                                        mc.set.seed = FALSE)),
    error = function(e) NULL)
  if (!is.list(res) || length(res) != length(idx)) return(NULL)
  for (r in res)
    if (!is.list(r) || is.null(r$fit) || !identical(r$rng, FALSE))
      return(NULL)

  for (r in res)
    for (cond in r$conds)
      if (inherits(cond, "warning")) warning(cond) else message(cond)
  lapply(res, `[[`, "fit")
}
