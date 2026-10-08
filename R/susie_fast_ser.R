# =============================================================================
# SUSIE-FAST: SINGLE-EFFECT REGRESSION
#
# Further savings in the single-effect regression (SER), on top of
# R/susie_fast.R. As there, the mathematics is unchanged.
#
#   1. SER memo. A null effect leaves the fitted values bit-for-bit
#      unchanged, so the next effect often sees exactly the same residuals
#      and prior variance as the previous one. Its SER is then the same
#      computation and the stored rows are reused. For individual data the
#      X'r product of the residuals is reused in the same way.
#   2. Table kernel. In the standard susie_rss and standardised susie and
#      susie_ss fits, shat2 takes one or a few distinct values. The terms
#      log(1 + V / s) and s * (V + s) are computed once per distinct value
#      instead of once per variable.
#   3. Final SER step for ss data. The log Bayes factors, posterior
#      weights, posterior moments and expected log-likelihood at the chosen
#      V are computed in C++ instead of about 40 R vector operations.
#
# All three apply in modes "full" and "exact"; mode "off" runs the
# upstream code.
# =============================================================================

# -----------------------------------------------------------------------------
# SER memo
# -----------------------------------------------------------------------------

# single_effect_regression(), reusing the previous SER of this IBSS run when
# every input it reads is bit-for-bit the same.
#' @keywords internal
fast_single_effect_regression <- function(data, params, model, l) {
  cache <- model$runtime$fast_cache
  if (is.null(cache) || !fast_ser_memo_ok(cache, data, params, model))
    return(single_effect_regression(data, params, model, l))

  # Everything the SER reads, under the default methods for ss and
  # individual data. l enters only through V[l] and the rows written.
  key <- list(data, params, model$residuals, model$residual_variance,
              model$shat2_inflation, model$predictor_weights, model$pi,
              model$V[l], model$sigma2, model$raw_residuals)
  out <- cache$ser_out
  if (!is.null(out) && identical(cache$ser_key, key, num.eq = FALSE)) {
    # The writes of apply_ser_lbf, record_R_bf_attenuation,
    # store_ser_moments, compute_kl and set_prior_variance_l, in order.
    model$alpha[l, ]        <- out$alpha
    model$lbf[l]            <- out$lbf
    model$lbf_variable[l, ] <- out$lbf_variable
    if (!is.null(out$atten)) {
      if (is.null(model$R_bf_attenuation))
        model$R_bf_attenuation <- matrix(NA_real_, nrow(model$alpha),
                                         length(out$atten))
      model$R_bf_attenuation[l, ] <- out$atten
    }
    model$mu[l, ]  <- out$mu
    model$mu2[l, ] <- out$mu2
    model$KL[l]    <- out$KL
    model$V[l]     <- out$V
    return(model)
  }

  # An SER that signals a warning or message is not stored, so a repeat
  # signals it again.
  cache$ser_out <- NULL
  quiet <- TRUE
  model <- withCallingHandlers(single_effect_regression(data, params, model, l),
                               condition = function(c) quiet <<- FALSE)
  if (quiet) {
    cache$ser_key <- key
    cache$ser_out <- list(
      alpha = model$alpha[l, ], lbf = model$lbf[l],
      lbf_variable = model$lbf_variable[l, ],
      atten = if (fast_ser_atten_recorded(model)) model$R_bf_attenuation[l, ],
      mu = model$mu[l, ], mu2 = model$mu2[l, ], KL = model$KL[l],
      V = model$V[l])
  }
  model
}

# TRUE when the SER wrote row l of R_bf_attenuation: the condition in
# record_R_bf_attenuation(), with length(shat2) from the recycling rules of
# compute_ser_statistics(). The inputs are unchanged by the SER.
#' @keywords internal
fast_ser_atten_recorded <- function(model) {
  infl <- model$shat2_inflation
  if (is.null(infl)) return(FALSE)
  n <- c(length(model$residual_variance), length(model$predictor_weights),
         length(infl))
  n_shat2 <- if (any(n == 0)) 0 else max(n)
  length(infl) == n_shat2 && any(is.finite(infl) & infl > 1)
}

# The memo is used only for the standard Gaussian SER on ss and individual
# data, with susieR's own methods for every generic the SER dispatches.
#' @keywords internal
fast_ser_memo_ok <- function(cache, data, params, model) {
  if (is.null(cache$ser_memo_ok)) {
    cls <- class(data)
    cache$ser_memo_ok <-
      (identical(cls, "ss") || identical(cls, "individual")) &&
      isTRUE(params$estimate_prior_method %in%
               c("optim", "uniroot", "simple", "none", "EM")) &&
      !isTRUE(params$use_NIG) &&
      identical(params$unmappable_effects, "none") &&
      all(vapply(c("compute_ser_statistics", "loglik", "neg_loglik",
                   "calculate_posterior_moments", "compute_kl",
                   "SER_posterior_e_loglik", "pre_loglik_prior_hook",
                   "post_loglik_prior_hook", "optimize_prior_variance",
                   "em_update_prior_variance"),
                 fast_own_method, TRUE, cls = cls)) &&
      all(vapply(c("get_prior_variance_l", "set_prior_variance_l",
                   "get_alpha_l", "get_posterior_moments_l"),
                 fast_own_method, TRUE, cls = "susie"))
  }
  cache$ser_memo_ok && identical(class(model), "susie") &&
    is.null(model$slot_weights) && is.null(model$c_hat_state)
}

# TRUE when S3 dispatch of `generic` on class `cls`, called from within
# susieR, finds susieR's own method (or its default method): methods
# defined in the namespace come first, then registered ones, then those on
# the search path.
#' @keywords internal
fast_own_method <- function(generic, cls) {
  ns <- asNamespace("susieR")
  table <- get0(".__S3MethodsTable__.", envir = ns, inherits = FALSE)
  for (k in c(cls, "default")) {
    name <- paste(generic, k, sep = ".")
    if (exists(name, envir = ns, inherits = FALSE)) return(TRUE)
    if ((!is.null(table) && exists(name, envir = table, inherits = FALSE)) ||
        exists(name, envir = parent.env(ns), mode = "function"))
      return(FALSE)
  }
  FALSE
}

# fast_Xty(model, data$X, y), reusing the previous product when X and y are
# unchanged (consecutive null effects give the same residuals).
#' @keywords internal
fast_Xty_memo <- function(data, model, y) {
  cache <- model$runtime$fast_cache
  if (is.null(cache))
    return(fast_Xty(model, data$X, y))
  if (!is.null(cache$xty_v) &&
      identical(cache$xty_X, data$X, num.eq = FALSE) &&
      identical(cache$xty_y, y, num.eq = FALSE))
    return(cache$xty_v)
  v <- fast_Xty(model, data$X, y)
  cache$xty_X <- data$X
  cache$xty_y <- y
  cache$xty_v <- v
  v
}

# -----------------------------------------------------------------------------
# Table kernel
# -----------------------------------------------------------------------------

# The model log Bayes factor as a function of V for fast_ser_evaluator():
# the table kernel when s has at most 32 distinct values, otherwise the
# per-variable kernel.
#' @keywords internal
fast_ser_kernel <- function(h, s, zero, logpw) {
  tab <- if (fast_ser_ok()) ser_s_table_cpp(s, zero, 32L)
  if (is.null(tab))
    return(function(V) ser_lbf_model_cpp(h, s, zero, logpw, V))
  u <- tab$u
  sidx <- tab$sidx
  function(V) ser_lbf_model_tab_cpp(h, u, sidx, logpw, V)
}

# -----------------------------------------------------------------------------
# Final SER step for ss data
# -----------------------------------------------------------------------------

#' @keywords internal
fast_ser_final_ok <- function(data, params, model) {
  fast_mode() != "off" && identical(class(data), "ss") &&
    !isTRUE(params$use_NIG) && identical(params$unmappable_effects, "none") &&
    is.null(model$shat2_inflation) && fast_ser_ok()
}

# loglik.ss() with l given, or NULL when the R code must be used.
#' @keywords internal
fast_ser_loglik_ss <- function(data, params, model, V, ser_stats, l) {
  if (is.null(l) || !fast_ser_final_ok(data, params, model)) return(NULL)
  betahat <- ser_stats$betahat
  shat2   <- ser_stats$shat2
  pi      <- model$pi
  p       <- length(betahat)
  if (p == 0 || !is.double(betahat) || !is.double(shat2) || !is.double(pi) ||
      length(shat2) != p || length(pi) != p ||
      !is.double(V) || length(V) != 1 || !is.finite(V) || V < 0 ||
      anyNA(pi) || any(pi < 0))
    return(NULL)
  out <- ser_lbf_l_cpp(betahat, shat2, pi, V)
  if (is.null(out)) return(NULL)
  model$alpha[l, ]        <- out$alpha
  model$lbf[l]            <- out$lbf_model
  model$lbf_variable[l, ] <- out$lbf
  model
}

# calculate_posterior_moments.ss(), or NULL when the R code must be used.
#' @keywords internal
fast_ser_moments_ss <- function(data, params, model, V, l) {
  if (!fast_ser_final_ok(data, params, model)) return(NULL)
  r  <- model$residuals
  pw <- model$predictor_weights
  rv <- model$residual_variance
  if (!is.double(r) || !is.double(pw) || length(r) != length(pw) ||
      !is.double(rv) || length(rv) != 1 || !is.double(V) || length(V) != 1)
    return(NULL)
  moments <- ser_moments_cpp(r, pw, rv, V)
  if (is.null(moments)) return(NULL)
  store_ser_moments(model, l, moments)
}

# SER_posterior_e_loglik.ss(), or NULL when the R code must be used.
#' @keywords internal
fast_ser_e_loglik_ss <- function(data, params, model, l) {
  if (!fast_ser_final_ok(data, params, model)) return(NULL)
  r  <- model$residuals
  pw <- model$predictor_weights
  rv <- model$residual_variance
  alpha <- model$alpha[l, ]
  mu    <- model$mu[l, ]
  mu2   <- model$mu2[l, ]
  p <- length(pw)
  if (!is.double(r) || !is.double(pw) || length(r) != p ||
      !is.double(rv) || length(rv) != 1 || !is.double(alpha) ||
      !is.double(mu) || !is.double(mu2) ||
      length(alpha) != p || length(mu) != p || length(mu2) != p)
    return(NULL)
  e <- ser_e_loglik_cpp(alpha, mu, mu2, r, pw, rv)
  if (is.na(e)) NULL else e
}

# -----------------------------------------------------------------------------
# Self-test
# -----------------------------------------------------------------------------

# The kernels in src/fast_ser.cpp are used only after they have reproduced
# the R code bit for bit on this platform (checked once per session), and
# only together with the round-1 kernel.
#' @keywords internal
fast_ser_ok <- local({
  ok <- NA
  function() {
    if (is.na(ok))
      ok <<- fast_kernel_ok() &&
        isTRUE(tryCatch(fast_ser_self_test(), error = function(e) FALSE))
    ok
  }
})

#' @keywords internal
fast_ser_self_test <- function() {
  eps <- .Machine$double.eps
  ref_lbf_model <- function(betahat, shat2, pw, V) {
    lbf <- gaussian_ser_lbf(betahat, shat2, V)
    stable <- lbf_stabilization(lbf, pw, shat2)
    compute_posterior_weights(stable$lpo)$lbf_model
  }
  st <- if (exists(".Random.seed", globalenv(), inherits = FALSE))
          get(".Random.seed", globalenv()) else NULL
  on.exit(if (is.null(st)) suppressWarnings(rm(".Random.seed", envir = globalenv()))
          else assign(".Random.seed", st, envir = globalenv()))
  set.seed(20261007)
  Vs <- c(0, exp(-30), 1e-3, 0.37, 1, 50, exp(15))
  for (trial in 1:20) {
    p <- c(1, 2, 7, 60, 300)[trial %% 5 + 1]
    K <- min(p, c(1, 3, 7)[trial %% 3 + 1])
    # shat2 with K distinct values, some of them below eps.
    vals <- rexp(K) * 10^runif(K, -6, 2)
    if (trial %% 4 == 0) vals[1] <- eps / 4
    shat2 <- vals[sample.int(K, p, replace = TRUE)]
    betahat <- rnorm(p, sd = c(0.01, 1, 30, 200)[trial %% 4 + 1]) * sqrt(shat2)
    if (p > 6 && trial %% 2 == 0) {
      betahat[2] <- NaN; shat2[3] <- Inf; shat2[4] <- 0; shat2[5] <- -0
      shat2[6] <- eps / 2; betahat[7] <- NA
    }
    pw <- if (trial %% 3 == 0) rep(1 / p, p) else {
      w <- rexp(p); if (p > 2) w[1] <- 0; w / sum(w) }
    if (trial == 20) betahat[1] <- 1e160
    h <- 0.5 * betahat^2
    s <- pmax(shat2, eps)
    zero <- !is.finite(betahat) | !is.finite(shat2)
    logpw <- log(pw + sqrt(eps))
    tab <- ser_s_table_cpp(s, zero, 32L)
    if (is.null(tab) || !identical(is.null(ser_s_table_cpp(s, zero, 0L)), !all(zero)))
      return(FALSE)
    rv <- c(1, 0.37, 2999)[trial %% 3 + 1]
    pwt <- rexp(p) * 1000
    if (p > 2) pwt[2] <- 0
    r <- rnorm(p, sd = 30)
    if (p > 6 && trial %% 2 == 1) { r[3] <- NaN; r[4] <- NA }
    for (V in Vs) {
      a <- ser_lbf_model_tab_cpp(h, tab$u, tab$sidx, logpw, V)
      if (!identical(a, ser_lbf_model_cpp(h, s, zero, logpw, V))) return(FALSE)
      if (!is.na(a) && !identical(a, ref_lbf_model(betahat, shat2, pw, V)))
        return(FALSE)
      # Final step: log Bayes factors and posterior weights.
      out <- ser_lbf_l_cpp(betahat, shat2, pw, V)
      lbf <- gaussian_ser_lbf(betahat, shat2, V)
      stable <- lbf_stabilization(lbf, pw, shat2)
      if (!is.null(out)) {
        w <- compute_posterior_weights(stable$lpo)
        if (!identical(out, list(lbf = stable$lbf, alpha = w$alpha,
                                 lbf_model = w$lbf_model)))
          return(FALSE)
      } else if (all(is.finite(stable$lpo))) {
        return(FALSE)
      }
      # Posterior moments.
      mom <- ser_moments_cpp(r, pwt, rv, V)
      ref <- gaussian_ser_moments(r / pwt, rv / pwt, V)
      if (!is.null(mom) && !identical(mom, ref)) return(FALSE)
      if (is.null(mom) && !anyNA(unlist(ref))) return(FALSE)
      # Expected log-likelihood.
      if (!is.null(mom)) {
        alpha <- rexp(p); alpha <- alpha / sum(alpha)
        e <- ser_e_loglik_cpp(alpha, mom$post_mean, mom$post_mean2, r, pwt, rv)
        ref <- gaussian_ser_posterior_e_loglik(alpha, mom$post_mean,
                                               mom$post_mean2, r / pwt, rv / pwt)
        if (!identical(e, ref) && !(is.na(e) && !is.finite(ref))) return(FALSE)
      }
    }
  }
  # More distinct values than the cap: no table.
  s <- seq(1, 2, length.out = 40)
  if (!is.null(ser_s_table_cpp(s, logical(40), 32L))) return(FALSE)
  TRUE
}
