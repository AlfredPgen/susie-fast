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
#      recomputes.
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
  if (is.null(cache)) return(get_xcorr(data))
  if (is.null(cache$xcorr)) {
    # The safe_cov2cor() branch of get_xcorr(), shared with fast_ash_R().
    if (is.null(data$Xcorr_cache) && !is.null(data$XtX) &&
        any(!(diag(data$XtX) %in% c(0, 1))))
      cache$xcorr <- fast_ash_cov2cor(data, cache)
    else
      cache$xcorr <- get_xcorr(data)$Xcorr
  }
  data$Xcorr_cache <- cache$xcorr
  list(Xcorr = cache$xcorr, data = data)
}

#' @keywords internal
fast_ash_cov2cor <- function(data, cache) {
  if (is.null(cache$ash_R)) cache$ash_R <- safe_cov2cor(data$XtX)
  cache$ash_R
}

# safe_cov2cor(data$XtX) in compute_ash_from_summary_stats().
#' @keywords internal
fast_ash_R <- function(data, model) {
  cache <- model$runtime$fast_cache
  if (is.null(cache)) return(safe_cov2cor(data$XtX))
  fast_ash_cov2cor(data, cache)
}

# abs(Xcorr) > threshold, for LD_adj %*% pip. The cached copy is stored as
# 0/1 doubles (NA stays NA), which is what %*% coerces the logical matrix
# to, so the product is the same.
#' @keywords internal
fast_ld_adj <- function(Xcorr, threshold, model) {
  cache <- model$runtime$fast_cache
  if (is.null(cache) || !identical(Xcorr, cache$xcorr))
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
#' @keywords internal
fast_mr_ash_rss_fn <- function() {
  if (fast_mode() != "off" && fast_mr_ash_rss_ok()) fast_mr_ash_rss_cpp
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
  st <- if (exists(".Random.seed", globalenv(), inherits = FALSE))
          get(".Random.seed", globalenv()) else NULL
  on.exit(if (is.null(st)) suppressWarnings(rm(".Random.seed", envir = globalenv()))
          else assign(".Random.seed", st, envir = globalenv()))
  set.seed(20261008)
  for (trial in 1:24) {
    p <- c(1, 7, 60)[trial %% 3 + 1]
    K <- c(1, 5, 25)[(trial %/% 3) %% 3 + 1]
    n <- 200L
    X <- matrix(rnorm(n * p), n, p)
    if (p > 1) X[, 2] <- X[, 1] + 0.3 * X[, 2]
    R <- cor(X)
    b <- rnorm(p, sd = 0.05); b[1] <- 0.3
    y <- as.vector(X %*% b + rnorm(n))
    r <- as.vector(cor(X, y))
    z <- r * sqrt((n - 2) / (1 - r^2))
    shat <- rep(1 / sqrt(n), p) * exp(rnorm(p, sd = 0.1))
    bhat <- z * shat
    s0 <- (2^((0:(K - 1)) / K) - 1)^2 / (n - 1) * n
    if (trial %% 4 == 0) s0[1] <- 0
    w0 <- rep(1 / K, K)
    if (K > 1 && trial %% 5 == 0) w0 <- c(0, rep(1 / (K - 1), K - 1))
    mu1 <- if (trial %% 2 == 0) numeric(0) else rnorm(p, sd = 0.01)
    fl <- as.logical(bitwAnd(trial, c(1, 2, 4, 8)))
    var_y <- if (trial %% 6 == 0) Inf else var(y)
    if (!fast_mr_ash_rss_same(bhat, shat, z, R, var_y, n, 1.2, s0, w0, mu1,
                              1e-4, 1000L, fl[1], fl[2], fl[3], fl[4]))
      return(FALSE)
  }
  TRUE
}
