# =============================================================================
# SUSIE-FAST: susie_rss_lambda()
#
# The round-1 techniques (R/susie_fast.R) applied to the rss_lambda class,
# plus three repeats specific to it. As in round 1, the fit is bit-identical
# to upstream susieR: every value is either computed by the same R
# expression or reused from an earlier evaluation on bit-identical inputs.
#
#   1. Per-effect product cache and SER Bayes factor kernel (round 1), with
#      R supplied directly (data$R) as well as X.
#   2. signal = crossprod(SinvRj, residuals) is computed by
#      compute_ser_statistics and again by calculate_posterior_moments on the
#      same residuals; the second is looked up.
#   3. SER_posterior_e_loglik for a null effect (b_l = 0): S^{-1} b_l is +0,
#      so the products R r and S^{-1} b_l are not needed.
#   4. SinvRj = V %*% (Dinv * D * t(V)) is a p^3 product. refine and L_greedy
#      run the workhorse several times on the same data and residual
#      variance; within one susie_rss_lambda() call it is computed once.
#   5. Residual-variance update: the parts of the expected log-likelihood
#      that do not depend on sigma2 are computed once per optimize() call,
#      diag(R S^{-1} R) in one C++ pass, and repeated sigma2 looked up.
#
# Only data of class exactly "rss_lambda" takes these paths, and only when
# options(susieR.fast) is not "off".
# =============================================================================

#' @keywords internal
fast_lambda_class <- function(data)
  identical(class(data)[1], "rss_lambda")

# -----------------------------------------------------------------------------
# 2. signal = crossprod(SinvRj, residuals)
# -----------------------------------------------------------------------------

#' @keywords internal
fast_lambda_signal <- function(data, model) {
  cache <- model$runtime$fast_cache
  if (is.null(cache) || !fast_lambda_class(data))
    return(as.vector(crossprod(model$SinvRj, model$residuals)))
  if (!is.null(cache$sig_value) &&
      identical(cache$sig_residuals, model$residuals, num.eq = FALSE) &&
      identical(cache$sig_SinvRj, model$SinvRj, num.eq = FALSE))
    return(cache$sig_value)
  signal <- as.vector(crossprod(model$SinvRj, model$residuals))
  cache$sig_residuals <- model$residuals
  cache$sig_SinvRj    <- model$SinvRj
  cache$sig_value     <- signal
  signal
}

# -----------------------------------------------------------------------------
# 3. Null effect in SER_posterior_e_loglik
# -----------------------------------------------------------------------------

# TRUE when Eb = 0 and the upstream expression
#   -2 * sum(rR * SinvEb),  rR = R %*% r,  SinvEb = V %*% (Dinv * crossprod(V, Eb))
# is -2 * (+0): with V finite and Dinv finite, crossprod(V, 0) and every
# product after it are +/-0; with rR finite (|rR| <= max|R| * sum|r|),
# rR * SinvEb is +/-0, and sum() (long double accumulator starting at +0)
# returns +0. Only when R itself is used for rR (no X, a dense double R).
#' @keywords internal
fast_lambda_null_effect <- function(data, model, Eb, V, Dinv) {
  cache <- model$runtime$fast_cache
  if (is.null(cache) || !fast_lambda_class(data) ||
      anyNA(Eb) || any(Eb != 0) ||
      !is.null(data$X) || !is.null(data$XtX) ||
      !is.matrix(data$R) || !is.double(data$R))
    return(FALSE)
  if (!identical(cache$null_V, V, num.eq = FALSE)) {
    cache$null_ok   <- is.finite(sum(V)) && is.finite(sum(data$R))
    cache$null_Rmax <- max(abs(data$R))
    cache$null_V    <- V
  }
  bnd <- cache$null_Rmax * sum(abs(model$residuals))
  cache$null_ok && all(is.finite(Dinv)) && isTRUE(is.finite(bnd) && bnd < 1e300)
}

# -----------------------------------------------------------------------------
# 4. SinvRj and RjSinvRj within one susie_rss_lambda() call
# -----------------------------------------------------------------------------

# The memo lives only while fast_lambda_scope() evaluates its argument; it
# is never stored in the data or model objects.
.fast_lambda <- new.env(parent = emptyenv())

#' @keywords internal
fast_lambda_scope <- function(expr) {
  old <- .fast_lambda$memo
  .fast_lambda$memo <- new.env(parent = emptyenv())
  on.exit(.fast_lambda$memo <- old)
  expr
}

# list(SinvRj, RjSinvRj): from the memo when V, D and Dinv are bitwise the
# ones it was computed from, otherwise by forcing `compute` (the upstream
# expressions, evaluated in the caller). The first result in the scope is
# kept: it is the one for the initial residual variance, which every
# workhorse call of refine and L_greedy starts from.
#' @keywords internal
fast_lambda_sinv <- function(data, V, D, Dinv, compute) {
  memo <- .fast_lambda$memo
  if (is.null(memo) || fast_mode() == "off" || !fast_lambda_class(data))
    return(compute)
  e <- memo$sinv
  if (!is.null(e) &&
      identical(e$V, V, num.eq = FALSE) &&
      identical(e$D, D, num.eq = FALSE) &&
      identical(e$Dinv, Dinv, num.eq = FALSE))
    return(e$value)
  value <- compute
  if (is.null(e))
    memo$sinv <- list(V = V, D = D, Dinv = Dinv, value = value)
  value
}

# -----------------------------------------------------------------------------
# 5. Residual-variance objective
# -----------------------------------------------------------------------------

# Returns a function equal to `objective` (Eloglik.rss_lambda at sigma2 with
# the rest of the model fixed), or `objective` itself when the fast path
# does not apply. Z %*% V, crossprod(V, zbar) and the other sigma2-free
# terms are computed once; the rest uses the expressions of compute_Dinv,
# get_ER2.rss_lambda and Eloglik.rss_lambda in the same order.
#' @keywords internal
fast_lambda_objective <- function(data, model, objective) {
  if (fast_mode() == "off" || !fast_lambda_class(data) ||
      is.null(model$Z) || is.null(model$zbar) || is.null(model$diag_postb2))
    return(objective)

  eigen_R <- get_eigen_R(data, model)
  D       <- eigen_R$values
  V       <- eigen_R$vectors
  lambda  <- data$lambda
  Vtz     <- get_Vtz(data, model)
  zbar    <- model$zbar
  postb2  <- model$diag_postb2
  z_null_norm2 <- if (!is.null(model$z_null_norm2)) model$z_null_norm2 else data$z_null_norm2

  D2      <- D^2
  DVtz    <- D * Vtz
  Vtzbar2 <- crossprod(V, zbar)^2
  VtZ2    <- (model$Z %*% V)^2
  use_cpp <- is.matrix(V) && is.double(V) && length(D) == ncol(V) &&
             is.finite(sum(V)) && fast_lambda_kernel_ok()
  V2      <- if (use_cpp) NULL else V^2

  f <- function(sigma2) {
    d    <- sigma2 * D + lambda
    Dinv <- 1 / d
    Dinv[is.infinite(Dinv)] <- 0

    zSinvz <- sum((Dinv * Vtz) * Vtz)
    if (lambda > 0) zSinvz <- zSinvz + z_null_norm2 / lambda
    tmp   <- V %*% (Dinv * DVtz)
    term2 <- -2 * sum(tmp * zbar)
    w     <- Dinv * D2
    term3 <- sum(Vtzbar2 * w)
    term4 <- sum(VtZ2 %*% w)
    diag_RSinvR <- if (use_cpp && all(is.finite(w))) rss_diag_rsinvr_cpp(V, w)
                   else rowSums((if (is.null(V2)) V^2 else V2) * rep(w, each = nrow(V)))
    term5 <- sum(diag_RSinvR * postb2)
    ER2   <- zSinvz + term2 + term3 - term4 + term5

    d_pos <- d[d > 0]
    r_eff <- length(d_pos)
    -(r_eff / 2) * log(2 * pi) - 0.5 * sum(log(d_pos)) - 0.5 * ER2
  }

  # optimize() evaluates the objective again at its result, and the caller
  # evaluates it once more there; look those up.
  memo_s <- numeric(0)
  memo_f <- numeric(0)
  function(sigma2) {
    k <- match(sigma2, memo_s)
    if (!is.na(k) && identical(memo_s[k], sigma2, num.eq = FALSE))
      return(memo_f[k])
    out <- f(sigma2)
    memo_s <<- c(memo_s, sigma2)
    memo_f <<- c(memo_f, out)
    out
  }
}

# rowSums((V^2) * rep(w, each = nrow(V))) in one pass, used only after it
# has reproduced rowSums bit for bit on this platform (checked once per
# session).
#' @keywords internal
fast_lambda_kernel_ok <- local({
  ok <- NA
  function() {
    if (is.na(ok))
      ok <<- isTRUE(tryCatch(fast_lambda_self_test(), error = function(e) FALSE))
    ok
  }
})

#' @keywords internal
fast_lambda_self_test <- function() {
  st <- if (exists(".Random.seed", globalenv(), inherits = FALSE))
          get(".Random.seed", globalenv()) else NULL
  on.exit(if (is.null(st)) suppressWarnings(rm(".Random.seed", envir = globalenv()))
          else assign(".Random.seed", st, envir = globalenv()))
  set.seed(20261008)
  for (trial in 1:12) {
    n <- sample(c(1, 2, 7, 64, 300), 1)
    p <- if (trial %% 3 == 0) n else sample(c(5, 64, 301), 1)
    V <- matrix(rnorm(n * p, sd = 10^runif(1, -3, 3)), n, p)
    if (n * p > 4) V[sample(n * p, 2)] <- 0
    # w on one scale, so that the long double accumulation matters
    w <- rexp(p) * 10^runif(1, -8, 8)
    if (trial %% 4 == 0) w[1:2] <- c(0, 1e300)
    if (!identical(rss_diag_rsinvr_cpp(V, w),
                   rowSums((V^2) * rep(w, each = nrow(V)))))
      return(FALSE)
  }
  TRUE
}
