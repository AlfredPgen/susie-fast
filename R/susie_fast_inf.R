# =============================================================================
# SUSIE-FAST: SuSiE-inf (unmappable_effects = "inf")
#
# The inf path recomputes several eigenspace products that are either never
# read or identical to one computed moments earlier. Each helper below does
# the same IEEE operations as the upstream function it replaces, minus the
# repeated or dead ones, and is used only inside susie_workhorse on the
# standard ss and individual classes (a non-NULL model$runtime$fast_cache)
# with options(susieR.fast) not "off". Otherwise the upstream code runs.
#
#   a. Residuals: compute_shat2_inflation returns NULL before reading
#      XtXr_without_l and r unless an R_finite reference size is set, so
#      those products are only formed when it is.
#   b. Fitted values: on this path nothing reads model$XtXr / model$Xr
#      between effect updates (see fast_inf_skip_fitted), so only the last
#      effect of a sweep forms the product, on the same state as upstream.
#   c. Consecutive null effects see the same b_minus_l, so the two
#      eigenspace products of the Omega-weighted residual are reused.
#   d. ELBO and method-of-moments update: t(V) is formed once per fit, not
#      L + 1 times per ELBO, and a null effect (b_l = 0) contributes
#      t(V) %*% b_l = +0 when V is finite, so its product is skipped.
#      rowSums(sweep(t(V^2), 2, tmpD, `*`)) is computed as
#      colSums(V^2 * tmpD): the same products (same operand order) summed in
#      the same order in long double (R's do_colsum, both branches).
#   e. compute_omega_quantities is looked up when tau2 and sigma2 are those
#      of the previous call (the ELBO follows update_derived_quantities).
#   f. compute_theta_blup reuses the ELBO's t(V) %*% b for the same b.
#
# Every lookup is keyed on bit-for-bit identical inputs.
# =============================================================================

# The per-workhorse cache, or NULL when the upstream code must run.
#' @keywords internal
fast_inf_cache <- function(model) {
  if (fast_mode() == "off") return(NULL)
  model$runtime$fast_cache
}

#' @keywords internal
fast_inf_lookup <- function(cache, name, key) {
  hit <- cache[[name]]
  if (!is.null(hit) && identical(hit$key, key, num.eq = FALSE)) hit$value
  else NULL
}

#' @keywords internal
fast_inf_store <- function(cache, name, key, value) {
  assign(name, list(key = key, value = value), envir = cache)
  value
}

# t(V), formed once per fit.
#' @keywords internal
fast_inf_tV <- function(cache, eigen_vectors) {
  if (is.null(cache$inf_tV)) cache$inf_tV <- t(eigen_vectors)
  cache$inf_tV
}

# TRUE when t(V) %*% b_l and crossprod(V, b_l) are +0 in every entry, so
# Vtbl^2 = +0 and diagVtMV - Vtbl^2 = diagVtMV: b_l is all +-0 and V is
# finite (a NaN or Inf in V would turn 0 * V into NaN).
#' @keywords internal
fast_inf_null <- function(cache, eigen_vectors, bl) {
  if (anyNA(bl) || any(bl != 0)) return(FALSE)
  if (is.null(cache$inf_V_finite))
    cache$inf_V_finite <- is.finite(sum(eigen_vectors))
  cache$inf_V_finite
}

# rowSums(sweep(t(V_sq), 2, tmpD, `*`)) without the transpose, or NULL when
# the shapes would make sweep() warn (the caller then runs it verbatim).
#' @keywords internal
fast_inf_diag_term <- function(V_sq, tmpD) {
  if (!is.matrix(V_sq) || !is.double(V_sq) || !is.double(tmpD) ||
      nrow(V_sq) != length(tmpD))
    return(NULL)
  colSums(V_sq * tmpD)
}

# -----------------------------------------------------------------------------
# a, c: compute_residuals (inf branch)
# -----------------------------------------------------------------------------

# Returns the updated model, or NULL to run the upstream branch.
#' @keywords internal
fast_residuals_inf <- function(data, model, b_minus_l) {
  cache <- fast_inf_cache(model)
  if (is.null(cache)) return(NULL)
  individual <- class(data)[1] == "individual"

  key <- list(b_minus_l, model$omega_var)
  hit <- fast_inf_lookup(cache, "inf_residual", key)
  if (is.null(hit)) {
    if (individual) {
      Vtb_minus_l <- as.vector(crossprod(data$eigen_vectors, b_minus_l))
      XtOmegaXb   <- as.vector(data$eigen_vectors %*%
                                 (Vtb_minus_l * data$eigen_values / model$omega_var))
    } else {
      Vtb_minus_l <- NULL
      XtOmegaXb   <- as.vector(data$eigen_vectors %*%
                                 ((crossprod(data$eigen_vectors, b_minus_l)) *
                                    data$eigen_values / model$omega_var))
    }
    hit <- fast_inf_store(cache, "inf_residual", key,
                          list(Vtb = Vtb_minus_l, XtOmegaXb = XtOmegaXb))
  }

  model$residuals         <- model$XtOmegay - hit$XtOmegaXb
  model$residual_variance <- 1

  # compute_shat2_inflation would return NULL without reading its inputs.
  if (is.null(get_current_R_finite_B(data, model)) ||
      isTRUE(model$sigma2 <= .Machine$double.eps))
    return(apply_inflation_state(model, NULL))

  XtXr_without_l <- if (individual)
                      as.vector(data$eigen_vectors %*%
                                  (data$eigen_values * hit$Vtb))
                    else compute_Rv(data, b_minus_l)
  r <- data$Xty - XtXr_without_l
  infl_state <- compute_shat2_inflation(data, model, XtXr_without_l,
                                        b_minus_l, r)
  apply_inflation_state(model, infl_state)
}

# -----------------------------------------------------------------------------
# b: update_fitted_values (inf branch)
# -----------------------------------------------------------------------------

# TRUE when the fitted-value product of effect l < L can be skipped. On the
# inf path model$XtXr / model$Xr are read only by
#   - compute_R_mismatch_state, initialize_R_mismatch (R_mismatch != "none"),
#   - adjust_fitted_for_c_hat, recompute_fitted_weighted (c_hat, slot prior),
#   - compute_shat2_inflation callers, which use their own b_minus_l product,
#   - get_ER2 (not called: the inf ELBO is compute_elbo_inf, the variance
#     update is MoM),
#   - the NIG residual sum of squares (excluded below as well),
#   - get_fitted and the final fit, after effect L has refreshed them.
# Effect L always runs last (no c_hat skipping) and recomputes the value
# from the whole state, so the final XtXr / Xr are upstream's.
#' @keywords internal
fast_inf_skip_fitted <- function(data, params, model, l) {
  !is.null(fast_inf_cache(model)) &&
    l < nrow(model$alpha) &&
    is.null(model$c_hat_state) && is.null(model$slot_weights) &&
    (is.null(params$R_mismatch) || identical(params$R_mismatch, "none")) &&
    is.null(get_current_R_finite_B(data, model)) &&
    !isTRUE(params$use_NIG)
}

# -----------------------------------------------------------------------------
# e: compute_omega_quantities
# -----------------------------------------------------------------------------

#' @keywords internal
fast_omega_quantities <- function(data, model) {
  cache <- fast_inf_cache(model)
  if (is.null(cache))
    return(compute_omega_quantities(data, model$tau2, model$sigma2))
  key <- list(model$tau2, model$sigma2)
  hit <- fast_inf_lookup(cache, "inf_omega", key)
  if (!is.null(hit)) return(hit)
  fast_inf_store(cache, "inf_omega", key,
                 compute_omega_quantities(data, model$tau2, model$sigma2))
}

# -----------------------------------------------------------------------------
# d, f: compute_elbo_inf, mom_unmappable, compute_theta_blup
# -----------------------------------------------------------------------------

# compute_elbo_inf with the per-workhorse cache as first argument.
#' @keywords internal
fast_elbo_inf <- function(cache, alpha, mu, omega, lbf, sigma2, tau2, n, p,
                          eigen_vectors, eigen_values, VtXty, yty,
                          eigen_vectors_sq = NULL) {
  if (is.null(cache) || fast_mode() == "off")
    return(compute_elbo_inf(alpha, mu, omega, lbf, sigma2, tau2, n, p,
                            eigen_vectors, eigen_values, VtXty, yty,
                            eigen_vectors_sq = eigen_vectors_sq))
  L <- nrow(mu)
  r <- length(eigen_values)
  if (is.null(eigen_vectors_sq))
    eigen_vectors_sq <- eigen_vectors^2

  tV  <- fast_inf_tV(cache, eigen_vectors)
  b   <- colSums(mu * alpha)
  Vtb <- tV %*% b
  fast_inf_store(cache, "inf_Vtb", b, Vtb)
  diagVtMV <- Vtb^2
  tmpD <- rep(0, p)

  for (l in seq_len(L)) {
    bl <- mu[l, ] * alpha[l, ]
    if (!fast_inf_null(cache, eigen_vectors, bl)) {
      Vtbl <- tV %*% bl
      diagVtMV <- diagVtMV - Vtbl^2
    }
    tmpD <- tmpD + alpha[l, ] * (mu[l, ]^2 + 1 / omega[l, ])
  }

  dterm <- fast_inf_diag_term(eigen_vectors_sq, tmpD)
  if (is.null(dterm))
    dterm <- rowSums(sweep(t(eigen_vectors_sq), 2, tmpD, `*`))
  diagVtMV <- diagVtMV + dterm

  var <- tau2 * eigen_values + sigma2
  neg_elbo <- 0.5 * (n - r) * log(sigma2) + 0.5 / sigma2 * yty +
    sum(0.5 * log(var) -
          0.5 * tau2 / sigma2 * VtXty^2 / var -
          Vtb * VtXty / var +
          0.5 * eigen_values / var * diagVtMV)
  -neg_elbo
}

# compute_theta_blup, reusing t(V) and the ELBO's t(V) %*% b.
#' @keywords internal
fast_theta_blup <- function(data, model) {
  cache <- fast_inf_cache(model)
  if (is.null(cache) || is.null(model$omega_var) || is.null(model$XtOmegay))
    return(compute_theta_blup(data, model))
  b   <- colSums(model$mu * model$alpha)
  Vtb <- fast_inf_lookup(cache, "inf_Vtb", b)
  if (is.null(Vtb)) Vtb <- fast_inf_tV(cache, data$eigen_vectors) %*% b
  XtOmegaXb <- as.vector(data$eigen_vectors %*%
                           (Vtb * data$eigen_values / model$omega_var))
  XtOmegar  <- model$XtOmegay - XtOmegaXb
  model$tau2 * XtOmegar
}

# mom_unmappable with null effects skipped and the transpose-free diagonal.
#' @keywords internal
fast_mom_unmappable <- function(data, params, model, omega, tau2,
                                est_tau2 = TRUE, est_sigma2 = TRUE) {
  cache <- fast_inf_cache(model)
  dterm_ok <- !is.null(cache) && is.matrix(data$eigen_vectors_sq)
  if (!dterm_ok)
    return(mom_unmappable(data, params, model, omega, tau2,
                          est_tau2 = est_tau2, est_sigma2 = est_sigma2))
  L <- nrow(model$mu)

  A <- matrix(0, nrow = 2, ncol = 2)
  A[1, 1] <- data$n
  A[1, 2] <- sum(data$eigen_values)
  A[2, 1] <- A[1, 2]
  A[2, 2] <- sum(data$eigen_values^2)

  b <- colSums(model$mu * model$alpha)
  Vtb <- crossprod(data$eigen_vectors, b)
  diagVtMV <- Vtb^2
  tmpD <- rep(0, data$p)

  for (l in seq_len(L)) {
    bl <- model$mu[l, ] * model$alpha[l, ]
    if (!fast_inf_null(cache, data$eigen_vectors, bl)) {
      Vtbl <- crossprod(data$eigen_vectors, bl)
      diagVtMV <- diagVtMV - Vtbl^2
    }
    tmpD <- tmpD + model$alpha[l, ] * (model$mu[l, ]^2 + 1 / omega[l, ])
  }

  dterm <- fast_inf_diag_term(data$eigen_vectors_sq, tmpD)
  if (is.null(dterm))
    dterm <- rowSums(sweep(t(data$eigen_vectors_sq), 2, tmpD, `*`))
  diagVtMV <- diagVtMV + dterm

  x <- rep(0, 2)
  x[1] <- data$yty - 2 * sum(b * data$Xty) + sum(data$eigen_values * diagVtMV)
  x[2] <- sum(data$Xty^2) - 2 * sum(Vtb * data$VtXty * data$eigen_values) +
    sum(data$eigen_values^2 * diagVtMV)

  if (est_tau2) {
    sol <- solve(A, x)
    if (sol[1] > 0 && sol[2] > 0) {
      sigma2 <- sol[1]
      tau2   <- sol[2]
    } else {
      sigma2 <- x[1] / data$n
      tau2   <- 0
    }
    if (params$verbose) {
      message(sprintf("Update (sigma^2,tau^2) to (%f,%e)\n", sigma2, tau2))
    }
  } else if (est_sigma2) {
    sigma2 <- (x[1] - A[1, 2] * tau2) / data$n
    if (params$verbose) {
      message(sprintf("Update sigma^2 to %f\n", sigma2))
    }
  }
  return(list(sigma2 = sigma2, tau2 = tau2))
}
