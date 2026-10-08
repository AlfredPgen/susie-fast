# =============================================================================
# SUSIE-FAST: SuSiE-ash (unmappable_effects = "ash", "ash_filter_archived")
#
#   1. Correlation matrix. get_xcorr() caches cor(X) (or cov2cor(XtX)) on
#      its local copy of the data object, which is dropped when the ash
#      update returns, so upstream recomputes it in every iteration: an
#      O(n p^2) cor(X) for individual data. It is kept in the per-workhorse
#      cache instead (only while the workhorse runs), together with the
#      safe_cov2cor(XtX) that compute_ash_from_summary_stats() passes to
#      mr.ash.rss and the LD adjacency matrix of the neighbourhood PIP.
#      These are deterministic functions of the data, which does not change
#      within a workhorse call, so the stored copies are the values upstream
#      recomputes. A computation that signals a condition (e.g. the 'NaNs
#      produced' warning of sqrt() on a negative diagonal of XtX) is not
#      stored: every later call then runs the upstream code again, so the
#      warnings repeat as often as upstream's do. Memory: besides the data,
#      the cache holds two p x p double matrices (Xcorr or cov2cor(XtX),
#      and the LD adjacency; p = 6000: 576 MB) for the whole workhorse
#      call, where upstream allocates them anew in each ash update.
#   2. mr.ash.rss. fast_mr_ash_rss_cpp (src/fast_ash.cpp) is mr_ash_rss_cpp
#      without the per-coordinate containers. It is used only after it has
#      reproduced mr_ash_rss_cpp bit for bit on this platform (checked once
#      per session).
#
# Both are off under options(susieR.fast = "off"); 1 also needs the
# per-workhorse cache, i.e. data of class "ss" or "individual".
# =============================================================================

# get_xcorr(data) for the ash updates, computed once per workhorse call.
#' @keywords internal
fast_get_xcorr <- function(data, model) {
  cache <- model$runtime$fast_cache
  if (is.null(cache) || isTRUE(cache$ash_nocache)) return(get_xcorr(data))
  Xcorr <- cache$xcorr
  if (is.null(Xcorr)) {
    # The safe_cov2cor() branch of get_xcorr(), shared with fast_ash_R().
    if (is.null(data$Xcorr_cache) && !is.null(data$XtX) &&
        any(!(diag(data$XtX) %in% c(0, 1))))
      Xcorr <- fast_ash_cov2cor(data, cache)
    else
      Xcorr <- fast_ash_once(get_xcorr(data)$Xcorr, cache)
    if (!isTRUE(cache$ash_nocache)) cache$xcorr <- Xcorr
  }
  data$Xcorr_cache <- Xcorr
  list(Xcorr = Xcorr, data = data)
}

# Evaluates expr; if that signalled a condition (passed on, not muffled),
# the ash entries of the cache are not used for the rest of the workhorse
# call, so that later calls recompute and signal it again as upstream does.
#' @keywords internal
fast_ash_once <- function(expr, cache) {
  signalled <- FALSE
  value <- withCallingHandlers(expr, condition = function(c) signalled <<- TRUE)
  if (signalled) cache$ash_nocache <- TRUE
  value
}

#' @keywords internal
fast_ash_cov2cor <- function(data, cache) {
  if (is.null(cache$ash_R)) {
    R <- fast_ash_once(safe_cov2cor(data$XtX), cache)
    if (isTRUE(cache$ash_nocache)) return(R)
    cache$ash_R <- R
  }
  cache$ash_R
}

# safe_cov2cor(data$XtX) in compute_ash_from_summary_stats().
#' @keywords internal
fast_ash_R <- function(data, model) {
  cache <- model$runtime$fast_cache
  if (is.null(cache) || isTRUE(cache$ash_nocache)) return(safe_cov2cor(data$XtX))
  fast_ash_cov2cor(data, cache)
}

# abs(Xcorr) > threshold, for LD_adj %*% pip. The cached copy is stored as
# 0/1 doubles (NA stays NA), which is what %*% coerces the logical matrix
# to, so the product is the same. Only for a base double matrix: anything
# else (e.g. the S4 Matrix that safe_cor() returns for sparse X with a
# constant column) takes the upstream expression, and its errors.
#' @keywords internal
fast_ld_adj <- function(Xcorr, threshold, model) {
  cache <- model$runtime$fast_cache
  if (is.null(cache) || isS4(Xcorr) || !is.matrix(Xcorr) || !is.double(Xcorr) ||
      !identical(Xcorr, cache$xcorr))
    return(abs(Xcorr) > threshold)
  if (is.null(cache$ld_adj) || !identical(cache$ld_adj_threshold, threshold)) {
    A <- abs(Xcorr) > threshold
    storage.mode(A) <- "double"
    cache$ld_adj <- A
    cache$ld_adj_threshold <- threshold
  }
  cache$ld_adj
}

# The mr.ash.rss kernel to call: fast_mr_ash_rss_cpp or mr_ash_rss_cpp.
# Inputs with NA/NaN go to the upstream kernel: the two propagate missing
# values differently (NA vs NaN), and identical() tells them apart.
#' @keywords internal
fast_mr_ash_rss_fn <- function(...) {
  if (fast_mode() != "off" &&
      !any(vapply(list(...), function(x) is.numeric(x) && anyNA(x), TRUE)) &&
      fast_mr_ash_rss_ok()) fast_mr_ash_rss_cpp
  else mr_ash_rss_cpp
}

#' @keywords internal
fast_mr_ash_rss_ok <- local({
  ok <- NA
  function() {
    if (is.na(ok))
      ok <<- isTRUE(tryCatch(fast_mr_ash_rss_self_test(), error = function(e) FALSE))
    ok
  }
})

# Runs mr_ash_rss_cpp and fast_mr_ash_rss_cpp on the same input and returns
# TRUE when the results and the in-place update of w0 are identical. Each
# call gets its own copy of w0, because both overwrite it.
#' @keywords internal
fast_mr_ash_rss_same <- function(bhat, shat, z, R, var_y, n, sigma2_e, s0, w0,
                                 mu1_init, tol, max_iter, update_w0,
                                 update_sigma, compute_ELBO, standardize) {
  run <- function(f) {
    w <- w0 + 0
    fit <- f(bhat, shat, z, R, var_y, n, sigma2_e, s0, w, mu1_init, tol,
             max_iter, update_w0, update_sigma, compute_ELBO, standardize)
    list(fit, w)
  }
  identical(run(mr_ash_rss_cpp), run(fast_mr_ash_rss_cpp))
}

#' @keywords internal
fast_mr_ash_rss_self_test <- function() {
  # Inputs come from a local Park-Miller generator (exact in doubles), so
  # the user's random number stream is left alone.
  state <- 20261008
  draw <- function(k, sd = 1) {
    u <- numeric(k)
    for (i in seq_len(k)) {
      state <<- (16807 * state) %% 2147483647
      u[i] <- state / 2147483647
    }
    sd * qnorm(u)
  }
  for (trial in 1:12) {
    p <- c(1, 7, 30)[trial %% 3 + 1]
    K <- c(1, 5, 25)[(trial %/% 3) %% 3 + 1]
    n <- 200L
    X <- matrix(draw(n * p), n, p)
    if (p > 1) X[, 2] <- X[, 1] + 0.3 * X[, 2]
    R <- cor(X)
    b <- draw(p, sd = 0.05); b[1] <- 0.3
    y <- as.vector(X %*% b + draw(n))
    r <- as.vector(cor(X, y))
    z <- r * sqrt((n - 2) / (1 - r^2))
    shat <- rep(1 / sqrt(n), p) * exp(draw(p, sd = 0.1))
    bhat <- z * shat
    s0 <- (2^((0:(K - 1)) / K) - 1)^2 / (n - 1) * n
    if (trial %% 4 == 0) s0[1] <- 0
    w0 <- rep(1 / K, K)
    if (K > 1 && trial %% 5 == 0) w0 <- c(0, rep(1 / (K - 1), K - 1))
    mu1 <- if (trial %% 2 == 0) numeric(0) else draw(p, sd = 0.01)
    fl <- as.logical(bitwAnd(trial, c(1, 2, 4, 8)))
    var_y <- if (trial %% 6 == 0) Inf else var(y)
    if (!fast_mr_ash_rss_same(bhat, shat, z, R, var_y, n, 1.2, s0, w0, mu1,
                              1e-4, 200L, fl[1], fl[2], fl[3], fl[4]))
      return(FALSE)
  }
  TRUE
}
