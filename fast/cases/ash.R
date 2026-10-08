# SuSiE-ash (work package "ash"): the cached correlation, cov2cor and LD
# adjacency matrices of the ash updates, and the mr.ash.rss kernel.
ash_ss01 <- local({
  # Unscaled correlation matrix passed straight to susie_ss: diag(XtX) is
  # all 1, so get_xcorr() returns XtX itself while mr.ash.rss gets
  # safe_cov2cor(XtX). Variant 5 has zero variance (row, column and
  # diagonal 0), which safe_cov2cor turns into 0 / 1.
  i <- 1:120
  R <- ds$R[i, i]; R[5, ] <- 0; R[, 5] <- 0; R[5, 5] <- 0
  list(XtX = R, Xty = ds$z[i] / sqrt(ds$n), yty = 1, n = 200)
})
ash_mr <- function(...) {
  # Direct mr.ash.rss call on the p = 300 region; returns the fit together
  # with the in-place updated w0.
  K <- 25; w0 <- rep(1 / K, K) + 0
  s0 <- (2^((0:24) / 25) - 1)^2 / (ds$n - 1) * ds$n
  args <- modifyList(list(bhat = ds$z / sqrt(ds$n), shat = rep(1 / sqrt(ds$n), length(ds$z)),
                          R = ds$R, var_y = 1, n = ds$n, s0 = s0, w0 = w0, tol = 1e-4,
                          max_iter = 1000), list(...))
  w0 <- args$w0
  fit <- do.call(mr.ash.rss, args)
  list(fit = fit, w0 = w0)
}
ash_conds <- function(f) {
  # Result (or error message) of f() together with its warnings, in order.
  w <- character(0)
  r <- tryCatch(withCallingHandlers(f(), warning = function(c) {
         w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") }),
       error = function(e) paste("ERROR:", conditionMessage(e)))
  if (is.list(r)) r$.diag_env <- NULL
  list(result = r, warnings = w)
}
ash_negdiag <- local({
  # Centred, unscaled XtX with a negative diagonal entry far from the
  # signal: safe_cov2cor() warns 'NaNs produced' in every call upstream.
  Xc <- scale(ds$X, scale = FALSE); yc <- ds$y - mean(ds$y)
  XtX <- crossprod(Xc); j <- order(abs(ds$z))[1]
  XtX[j, j] <- -XtX[j, j]
  list(XtX = XtX, Xty = as.vector(crossprod(Xc, yc)), yty = sum(yc^2), n = ds$n)
})
cases <- c(cases, list(
  ash_rss_p1000      = function() susie_rss(d1$z, d1$R, n = d1$n, unmappable_effects = "ash"),
  ash_rss_p1000_erv  = function() susie_rss(d1$z, d1$R, n = d1$n, unmappable_effects = "ash",
                                            estimate_residual_variance = TRUE),
  ash_rss_p3000      = function() susie_rss(d2$z, d2$R, n = d2$n, unmappable_effects = "ash"),
  ash_rss_archived   = function() susie_rss(ds$z, ds$R, n = ds$n, unmappable_effects = "ash_filter_archived",
                                            estimate_residual_variance = TRUE),
  ash_ss_p1000       = function() susie_ss(ss$XtX, ss$Xty, ss$yty, ss$n, unmappable_effects = "ash"),
  ash_ss_diag01      = function() susie_ss(ash_ss01$XtX, ash_ss01$Xty, ash_ss01$yty, ash_ss01$n,
                                           standardize = FALSE, unmappable_effects = "ash",
                                           estimate_residual_variance = TRUE),
  ash_ss_diag01_arch = function() susie_ss(ash_ss01$XtX, ash_ss01$Xty, ash_ss01$yty, ash_ss01$n,
                                           standardize = FALSE, unmappable_effects = "ash_filter_archived",
                                           estimate_residual_variance = TRUE),
  ash_ind_p1000      = function() susie(d1$X, d1$y, unmappable_effects = "ash"),
  ash_ind_archived   = function() susie(ds$X, ds$y, unmappable_effects = "ash_filter_archived"),
  ash_ind_const      = function() { X <- ds$X; X[, 7] <- 1; susie(X, ds$y, unmappable_effects = "ash") },
  ash_ss_negdiag     = function() ash_conds(function()
                         susie_ss(ash_negdiag$XtX, ash_negdiag$Xty, ash_negdiag$yty, ash_negdiag$n,
                                  standardize = FALSE, unmappable_effects = "ash",
                                  estimate_residual_variance = TRUE)),
  ash_ss_negdiag_arch = function() ash_conds(function()
                         susie_ss(ash_negdiag$XtX, ash_negdiag$Xty, ash_negdiag$yty, ash_negdiag$n,
                                  standardize = FALSE, unmappable_effects = "ash_filter_archived")),
  ash_mr_default     = function() ash_mr(),
  ash_mr_std         = function() ash_mr(standardize = TRUE, sigma2_e = 0.9, mu1_init = ds$z / 100),
  ash_mr_flags_off   = function() ash_mr(update_w0 = FALSE, update_sigma = FALSE, compute_ELBO = FALSE),
  ash_mr_varyInf     = function() ash_mr(var_y = Inf, z = ds$z),
  ash_mr_w0zero_s00  = function() { w <- c(0, rep(1 / 24, 24)); s <- (2^((0:24) / 25) - 1)^2
                                    ash_mr(w0 = w, s0 = s) },
  ash_mr_nan_maxit   = function() { z <- ds$z; z[3] <- NaN; ash_mr(z = z, bhat = z / sqrt(ds$n), max_iter = 5) },
  ash_mr_zero_diag   = function() { R <- ds$R; R[4, 4] <- 0; ash_mr(R = R, max_iter = 5) }
))
