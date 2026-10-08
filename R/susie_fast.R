# =============================================================================
# SUSIE-FAST
#
# Faster evaluation of the standard IBSS code paths. The mathematics is
# unchanged; the speed comes from not repeating work.
#
#   1. Per-effect product cache. Each effect update needs R %*% b_l (or
#      X %*% b_l) for the old b_l (to remove effect l from the fit) and for
#      the new b_l (to add it back). The new product is stored and reused
#      as the old product at the next visit, and within one update the
#      individual-data KL term reuses it too. Hits are exact: a product is
#      reused only when b_l is bit-for-bit the vector it was computed from.
#      Null effects (b_l = 0, e.g. prior variance estimated as 0) need no
#      product at all.
#   2. Single-effect Bayes factor kernel. The prior-variance optimiser
#      evaluates the SER log Bayes factor 20-40 times per effect. The terms
#      that do not depend on V are computed once, the rest in one C++ pass,
#      and repeated evaluations at the same V are looked up, not redone.
#   3. Expected residual sum of squares (ELBO, residual variance). It is
#      memoised within an iteration. When it feeds only the ELBO (residual
#      variance fixed, no NIG prior), full mode assembles it from the cached
#      per-effect products instead of L more matrix products.
#   4. Constructor: XtX is standardised in one pass instead of three
#      full-size copies.
#
# Controlled by options(susieR.fast = ...):
#   "full"  (default) everything above. The fit is bit-identical to upstream
#           susieR; the ELBO trace may differ in the last digit or two,
#           because its terms in 3 are summed in a different order.
#   "exact" 1, 2, 4 and the memo in 3. Bit-identical, ELBO included.
#   "off"   the upstream code paths.
# =============================================================================

#' @keywords internal
fast_mode <- function() {
  m <- getOption("susieR.fast", "full")
  if (isFALSE(m)) return("off")
  if (isTRUE(m)) return("full")
  m <- as.character(m)[1]
  if (!m %in% c("full", "exact", "off"))
    stop("options(susieR.fast) must be \"full\", \"exact\" or \"off\".")
  m
}

# Per-workhorse runtime cache. The data object is constant for the duration
# of one susie_workhorse call, so products of it can be cached there.
#' @keywords internal
fast_cache_new <- function(data, params, L) {
  mode <- fast_mode()
  if (mode == "off" || !class(data)[1] %in% c("ss", "individual", "rss_lambda"))
    return(NULL)
  cache <- new.env(parent = emptyenv())
  # get_ER2 may be re-summed from the cache only when it feeds just the ELBO.
  cache$er2_from_products <- mode == "full" && class(data)[1] != "rss_lambda" &&
    !isTRUE(params$estimate_residual_variance) && !isTRUE(params$use_NIG)
  cache$b    <- vector("list", L)
  cache$prod <- vector("list", L)
  cache$zero_ok <- NA
  cache
}

# R %*% b (ss data) or standardized X %*% b (individual data) for effect l,
# reusing the cached product when b is unchanged. A null effect (b = 0) has
# product 0: BLAS returns +0 in every entry when the matrix is finite, so
# the product is skipped.
#' @keywords internal
fast_effect_product <- function(data, model, l, b, compute, out_len, finite) {
  cache <- model$runtime$fast_cache
  if (is.null(cache) || l > length(cache$b))
    return(compute(b))
  if (!is.null(cache$prod[[l]]) &&
      identical(cache$b[[l]], b, num.eq = FALSE))
    return(cache$prod[[l]])
  if (!anyNA(b) && !any(b != 0)) {
    if (is.na(cache$zero_ok)) cache$zero_ok <- isTRUE(finite())
    v <- if (cache$zero_ok) numeric(out_len) else compute(b)
  } else {
    v <- compute(b)
  }
  cache$b[[l]]    <- b
  cache$prod[[l]] <- v
  v
}

#' @keywords internal
fast_Rv_l <- function(data, model, l, b)
  fast_effect_product(data, model, l, b,
    compute = function(v) fast_Rv(model, data, v),
    out_len = length(b),
    finite  = function() {
      M <- if (!is.null(data$X)) data$X else if (!is.null(data$XtX)) data$XtX
           else if (is.matrix(data$R) && is.double(data$R)) data$R
      !is.null(M) && is.finite(sum(M))
    })

#' @keywords internal
fast_Xb_l <- function(data, model, l, b)
  fast_effect_product(data, model, l, b,
    compute = function(v) fast_Xb(model, data$X, v),
    out_len = nrow(data$X),
    finite  = function() {
      csd <- attr(data$X, "scaled:scale")
      cm  <- attr(data$X, "scaled:center")
      is.null(attr(data$X, "matrix.type")) &&
        is.finite(sum(data$X)) && length(csd) == length(b) &&
        all(is.finite(csd) & csd != 0) && all(is.finite(cm))
    })

# All L cached products as an L x m matrix, or NULL when any is missing or
# stale. Used by full mode to assemble the expected residual sum of squares.
#' @keywords internal
fast_all_products <- function(model, B) {
  cache <- model$runtime$fast_cache
  if (is.null(cache) || !cache$er2_from_products || nrow(B) != length(cache$b))
    return(NULL)
  for (l in seq_len(nrow(B)))
    if (is.null(cache$prod[[l]]) ||
        !identical(cache$b[[l]], B[l, ], num.eq = FALSE))
      return(NULL)
  do.call(rbind, cache$prod)
}

# Memo for get_ER2: the ELBO and the residual-variance update evaluate it on
# the same model state in each iteration.
#' @keywords internal
fast_er2_key <- function(model)
  list(model$alpha, model$mu, model$mu2, model$slot_weights,
       model$predictor_weights, model$Xr, model$XtXr, model$X_theta,
       model$XtX_theta)

#' @keywords internal
fast_er2_lookup <- function(model) {
  cache <- model$runtime$fast_cache
  if (is.null(cache) || is.null(cache$er2_key)) return(NULL)
  if (identical(cache$er2_key, fast_er2_key(model), num.eq = FALSE))
    return(cache$er2_value)
  NULL
}

#' @keywords internal
fast_er2_store <- function(model, value) {
  cache <- model$runtime$fast_cache
  if (!is.null(cache)) {
    cache$er2_key   <- fast_er2_key(model)
    cache$er2_value <- value
  }
  value
}

# -----------------------------------------------------------------------------
# Single-effect log Bayes factor evaluator
# -----------------------------------------------------------------------------

# Returns list(neg = function(V_param), ll = function(V)) reproducing
# neg_loglik() / loglik() (with l = NULL) for the standard Gaussian SER on
# ss and individual data, or NULL when another method must be used.
#' @keywords internal
fast_ser_evaluator <- function(data, params, model, ser_stats) {
  if (fast_mode() == "off" ||
      !class(data)[1] %in% c("ss", "individual", "rss_lambda") ||
      isTRUE(params$use_NIG) ||
      identical(params$unmappable_effects, "inf") ||
      !identical(ser_stats$optim_scale, "log"))
    return(NULL)

  betahat <- ser_stats$betahat
  shat2   <- ser_stats$shat2
  # Same expressions as gaussian_ser_lbf() and lbf_stabilization().
  h     <- 0.5 * betahat^2
  s     <- pmax(shat2, .Machine$double.eps)
  zero  <- !is.finite(betahat) | !is.finite(shat2)
  logpw <- log(model$pi + sqrt(.Machine$double.eps))
  if (length(h) != length(logpw)) return(NULL)

  use_cpp <- fast_kernel_ok()
  kernel  <- if (use_cpp) fast_ser_kernel(h, s, zero, logpw)
  lbf_model <- function(V) {
    if (use_cpp) {
      out <- kernel(V)
      if (!is.na(out)) return(out)
    }
    lbf <- -0.5 * log(1 + V / s) + h * V / (s * (V + s))
    lbf[zero] <- 0
    lpo <- lbf + logpw
    compute_posterior_weights(lpo)$lbf_model
  }

  # Brent re-evaluates the objective at its minimiser, and the caller
  # evaluates it again at the chosen V; look those up instead.
  memo_V <- numeric(0)
  memo_f <- numeric(0)
  ll <- function(V) {
    k <- match(V, memo_V)
    if (!is.na(k) && identical(memo_V[k], V, num.eq = FALSE))
      return(memo_f[k])
    f <- lbf_model(V)
    memo_V <<- c(memo_V, V)
    memo_f <<- c(memo_f, f)
    f
  }
  list(neg = function(V_param) -ll(exp(V_param)), ll = ll)
}

# The C++ kernel is used only after it has reproduced the R arithmetic bit
# for bit on this platform (checked once per session).
#' @keywords internal
fast_kernel_ok <- local({
  ok <- NA
  function() {
    if (is.na(ok)) ok <<- isTRUE(tryCatch(fast_self_test(), error = function(e) FALSE))
    ok
  }
})

#' @keywords internal
fast_self_test <- function() {
  ref <- function(betahat, shat2, pw, V) {
    lbf <- gaussian_ser_lbf(betahat, shat2, V)
    stable <- lbf_stabilization(lbf, pw, shat2)
    compute_posterior_weights(stable$lpo)$lbf_model
  }
  st <- if (exists(".Random.seed", globalenv(), inherits = FALSE))
          get(".Random.seed", globalenv()) else NULL
  on.exit(if (is.null(st)) suppressWarnings(rm(".Random.seed", envir = globalenv()))
          else assign(".Random.seed", st, envir = globalenv()))
  set.seed(20261006)
  for (trial in 1:20) {
    p <- sample(c(1, 2, 7, 100, 1000, 5000), 1)
    betahat <- rnorm(p, sd = sample(c(0.01, 1, 30), 1))
    shat2 <- rexp(p) * 10^runif(1, -6, 2)
    if (p > 2 && trial %% 3 == 0) {
      betahat[2] <- NaN; shat2[3] <- Inf; shat2[1] <- 0
    }
    pw <- if (trial %% 2 == 0) rep(1 / p, p) else { w <- rexp(p); w / sum(w) }
    h <- 0.5 * betahat^2
    s <- pmax(shat2, .Machine$double.eps)
    zero <- !is.finite(betahat) | !is.finite(shat2)
    logpw <- log(pw + sqrt(.Machine$double.eps))
    for (V in c(0, exp(-30), 1e-3, 0.37, 1, 50, exp(15))) {
      a <- ser_lbf_model_cpp(h, s, zero, logpw, V)
      b <- ref(betahat, shat2, pw, V)
      if (!identical(a, b)) return(FALSE)
    }
  }
  for (p in c(1, 3, 130)) {
    M <- crossprod(matrix(rnorm(p * (p + 5)), p + 5, p))
    csd <- sqrt(diag(M) / 9)
    if (!identical(scale_xtx_cpp(M, csd), as.vector(t((1 / csd) * M) / csd)))
      return(FALSE)
  }
  TRUE
}

# Standardise XtX by csd: same values and attributes as
# t((1 / csd) * XtX) / csd.
#' @keywords internal
fast_scale_xtx <- function(XtX, csd) {
  if (fast_mode() == "off" || !is.matrix(XtX) || !is.double(XtX) ||
      isS4(XtX) || nrow(XtX) != ncol(XtX) || !is.null(attr(XtX, "names")) ||
      !is.null(attr(XtX, "tsp")) || !fast_kernel_ok())
    return(t((1 / csd) * XtX) / csd)
  at <- attributes(XtX)
  if (!is.null(at$dimnames)) at$dimnames <- at$dimnames[2:1]
  out <- scale_xtx_cpp(XtX, csd)
  attributes(out) <- at
  out
}
