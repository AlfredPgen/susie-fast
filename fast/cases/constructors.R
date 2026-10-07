# Data constructors: symmetry/NA check and scaling of XtX, original-scale
# XtX (bhat, shat, var_y), constant-column screen. Each case returns the
# fit or the constructed data object together with every message, warning
# and error raised, in order.
con_capture <- function(expr) {
  msgs <- character(0)
  val <- withCallingHandlers(
    tryCatch(expr, error = function(e) structure(list(msg = conditionMessage(e)), class = "err")),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") },
    warning = function(w) { msgs <<- c(msgs, paste("W:", conditionMessage(w))); invokeRestart("muffleWarning") })
  if (inherits(val, "err")) val <- list(error = val$msg)
  if (inherits(val, "susie")) { val$captured_msgs <- msgs; return(val) }
  c(unclass(val), list(captured_msgs = msgs))
}
con_ss  <- function(...) con_capture(susieR:::summary_stats_constructor(...))
# As con_capture, but each warning also records its call.
con_capture_calls <- function(expr) {
  msgs <- character(0)
  val <- withCallingHandlers(
    tryCatch(expr, error = function(e) structure(list(msg = conditionMessage(e)), class = "err")),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") },
    warning = function(w) {
      msgs <<- c(msgs, paste("W:", conditionMessage(w), "| call:", paste(deparse(conditionCall(w)), collapse = " ")))
      invokeRestart("muffleWarning")
    })
  if (inherits(val, "err")) val <- list(error = val$msg)
  if (inherits(val, "susie")) { val$captured_msgs <- msgs; return(val) }
  c(unclass(val), list(captured_msgs = msgs))
}
con_suf <- function(...) con_capture(susieR:::sufficient_stats_constructor(...))
con_ind <- function(...) con_capture(susieR:::individual_data_constructor(...))

con_R <- d1$R
con_R_dn <- con_R; dimnames(con_R_dn) <- list(paste0("rs", 1:ncol(con_R)), paste0("rs", 1:ncol(con_R)))
attr(con_R_dn, "foo") <- "bar"
con_R_cn <- con_R; colnames(con_R_cn) <- paste0("rs", 1:ncol(con_R))
con_R_asym <- con_R; con_R_asym[3, 7] <- con_R_asym[3, 7] + 1e-3
con_R_z <- con_R; con_R_z[5, 9] <- 0; con_R_z[9, 5] <- -0
con_R_na <- con_R; con_R_na[4, 4] <- NA
con_R_nan <- con_R; con_R_nan[2, 6] <- NaN; con_R_nan[6, 2] <- NaN
con_R_inf <- con_R; con_R_inf[2, 6] <- Inf; con_R_inf[6, 2] <- Inf
con_R_diag <- con_R; diag(con_R_diag) <- 1 + (seq_len(ncol(con_R)) %% 7) * 1e-3
con_bh <- setNames(uni$betahat, paste0("rs", seq_along(uni$betahat)))
con_sh <- setNames(uni$sebetahat, paste0("rs", seq_along(uni$sebetahat)))
con_XtX_neg <- ss$XtX; con_XtX_neg[10, 10] <- -con_XtX_neg[10, 10]
con_XtX_names <- ss$XtX; names(con_XtX_names) <- paste0("e", seq_along(con_XtX_names))
con_maf <- { set.seed(11); runif(ncol(ss$XtX)) }
con_X_const <- d1$X[, 1:200]
con_X_const[, 5] <- 2; con_X_const[, 9] <- 0; con_X_const[, 11] <- c(1e-170, rep(0, nrow(con_X_const) - 1))
con_X_const[, 12] <- 1 / 3; con_X_const[, 13] <- -0; con_X_const[, 20] <- con_X_const[, 20] * 1e-120
con_X_inf <- con_X_const; con_X_inf[1, 30] <- Inf
con_R200 <- d1$R[1:200, 1:200]
con_bh200 <- uni$betahat[1:200]; con_sh200 <- uni$sebetahat[1:200]
con_XtX_tsp <- ss$XtX; attr(con_XtX_tsp, "tsp") <- c(1, nrow(con_XtX_tsp), 1)
con_XtX_dnn <- ss$XtX
dimnames(con_XtX_dnn) <- list(a = paste0("v", 1:ncol(ss$XtX)), b = paste0("v", 1:ncol(ss$XtX)))
con_X_dn <- con_X_const; dimnames(con_X_dn) <- list(paste0("i", 1:nrow(con_X_dn)), paste0("v", 1:ncol(con_X_dn)))

cases <- c(cases, list(
  # sufficient_stats_constructor via summary_stats_constructor (z path)
  con_rss_obj      = function() con_ss(z = d1$z, R = d1$R, n = d1$n),
  con_rss_obj3000  = function() con_ss(z = d2$z, R = d2$R, n = d2$n),
  con_rss_dn       = function() con_ss(z = d1$z, R = con_R_dn, n = d1$n),
  con_rss_cn       = function() con_ss(z = d1$z, R = con_R_cn, n = d1$n),
  con_rss_asym     = function() con_ss(z = d1$z, R = con_R_asym, n = d1$n),
  con_rss_negzero  = function() con_ss(z = d1$z, R = con_R_z, n = d1$n),
  con_rss_na       = function() con_ss(z = d1$z, R = con_R_na, n = d1$n),
  con_rss_nan      = function() con_ss(z = d1$z, R = con_R_nan, n = d1$n),
  con_rss_inf      = function() con_ss(z = d1$z, R = con_R_inf, n = d1$n),
  con_rss_diag     = function() con_ss(z = d1$z, R = con_R_diag, n = d1$n),
  con_rss_nullw_dn = function() con_ss(z = d1$z, R = con_R_dn, n = d1$n, null_weight = 0.1),
  con_rss_maf      = function() con_ss(z = d1$z, R = con_R_dn, n = d1$n, maf = con_maf, maf_thresh = 0.3),
  con_rss_no_n     = function() con_ss(z = d1$z, R = con_R_dn),
  con_rss_nostd    = function() con_ss(z = d1$z, R = d1$R, n = d1$n, standardize = FALSE),
  con_rss_chk      = function() con_ss(z = d1$z, R = d1$R, n = d1$n, check_input = TRUE),
  con_fit_dn       = function() con_capture(susie_rss(d1$z, con_R_dn, n = d1$n)),
  con_fit_asym     = function() con_capture(susie_rss(d1$z, con_R_asym, n = d1$n)),
  con_fit_diag     = function() con_capture(susie_rss(d1$z, con_R_diag, n = d1$n)),
  con_fit_negzero  = function() con_capture(susie_rss(d1$z, con_R_z, n = d1$n)),
  con_fit_nullw_cn = function() con_capture(susie_rss(d1$z, con_R_cn, n = d1$n, null_weight = 0.1)),
  # original scale (bhat, shat, var_y)
  con_os_obj       = function() con_ss(bhat = uni$betahat, shat = uni$sebetahat, R = d1$R, n = d1$n, var_y = var(d1$y)),
  con_os_named     = function() con_ss(bhat = con_bh, shat = con_sh, R = con_R_dn, n = d1$n, var_y = var(d1$y)),
  con_os_asym      = function() con_ss(bhat = uni$betahat, shat = uni$sebetahat, R = con_R_asym, n = d1$n, var_y = var(d1$y)),
  con_os_nan       = function() con_ss(bhat = uni$betahat, shat = uni$sebetahat, R = con_R_nan, n = d1$n, var_y = var(d1$y)),
  con_os_maf       = function() con_ss(bhat = con_bh, shat = con_sh, R = con_R_cn, n = d1$n, var_y = var(d1$y), maf = con_maf, maf_thresh = 0.5),
  con_os_one       = function() con_ss(bhat = uni$betahat[1], shat = uni$sebetahat[1], R = matrix(1), n = d1$n, var_y = var(d1$y)),
  con_fit_os_named = function() con_capture(susie_rss(bhat = con_bh, shat = con_sh, R = con_R_dn, n = d1$n, var_y = var(d1$y))),
  # original scale with a negative PVE-adjusted diagonal: sqrt() warns
  # twice, from the call sqrt(working$XtXdiag)
  con_os_negvar    = function() con_capture_calls(susieR:::summary_stats_constructor(bhat = con_bh200, shat = con_sh200, R = con_R200, n = d1$n, var_y = -1)),
  con_os_n15       = function() con_capture_calls(susieR:::summary_stats_constructor(bhat = con_bh200, shat = con_sh200, R = con_R200, n = 1.5, var_y = 1)),
  con_os_neginf    = function() con_capture_calls(susieR:::summary_stats_constructor(bhat = con_bh200, shat = con_sh200, R = con_R200, n = d1$n, var_y = -Inf)),
  con_os_varNA     = function() con_capture_calls(susieR:::summary_stats_constructor(bhat = con_bh200, shat = con_sh200, R = con_R200, n = d1$n, var_y = NA_real_)),
  con_os_varInf    = function() con_capture_calls(susieR:::summary_stats_constructor(bhat = con_bh200, shat = con_sh200, R = con_R200, n = d1$n, var_y = Inf)),
  con_fit_os_negvar = function() con_capture_calls(susie_rss(bhat = con_bh200, shat = con_sh200, R = con_R200, n = d1$n, var_y = -1)),
  con_fit_os_n15   = function() con_capture_calls(susie_rss(bhat = con_bh200, shat = con_sh200, R = con_R200, n = 1.5, var_y = 1)),
  # tile-boundary sizes
  con_os_p65       = function() con_ss(bhat = uni$betahat[1:65], shat = uni$sebetahat[1:65], R = d1$R[1:65, 1:65], n = d1$n, var_y = var(d1$y)),
  con_rss_p64      = function() con_ss(z = d1$z[1:64], R = d1$R[1:64, 1:64], n = d1$n),
  con_rss_p129     = function() con_ss(z = d1$z[1:129], R = d1$R[1:129, 1:129], n = d1$n),
  con_fit_os_nullw = function() con_capture(susie_rss(bhat = uni$betahat, shat = uni$sebetahat, R = d1$R, n = d1$n, var_y = var(d1$y), null_weight = 0.2)),
  # sufficient_stats_constructor directly (susie_ss)
  con_suf_obj      = function() con_suf(XtX = ss$XtX, Xty = ss$Xty, yty = ss$yty, n = ss$n),
  con_suf_maf      = function() con_suf(XtX = ss$XtX, Xty = ss$Xty, yty = ss$yty, n = ss$n, maf = con_maf, maf_thresh = 0.4),
  con_suf_maf1     = function() con_suf(XtX = ss$XtX, Xty = ss$Xty, yty = ss$yty, n = ss$n, maf = c(1, rep(0, ncol(ss$XtX) - 1)), maf_thresh = 0.5),
  con_suf_negdiag  = function() con_suf(XtX = con_XtX_neg, Xty = ss$Xty, yty = ss$yty, n = ss$n),
  con_suf_names    = function() con_suf(XtX = con_XtX_names, Xty = ss$Xty, yty = ss$yty, n = ss$n),
  con_suf_tsp      = function() con_capture_calls(susieR:::sufficient_stats_constructor(XtX = con_XtX_tsp, Xty = ss$Xty, yty = ss$yty, n = ss$n)),
  con_suf_dnn      = function() con_suf(XtX = con_XtX_dnn, Xty = ss$Xty, yty = ss$yty, n = ss$n),
  con_suf_p63      = function() con_suf(XtX = ss$XtX[1:63, 1:63], Xty = ss$Xty[1:63], yty = ss$yty, n = ss$n),
  con_suf_nullw    = function() con_suf(XtX = ss$XtX, Xty = ss$Xty, yty = ss$yty, n = ss$n, null_weight = 0.3),
  con_suf_nostd    = function() con_suf(XtX = ss$XtX, Xty = ss$Xty, yty = ss$yty, n = ss$n, standardize = FALSE),
  con_suf_chk      = function() con_suf(XtX = ss$XtX, Xty = ss$Xty, yty = ss$yty, n = ss$n, check_input = TRUE),
  con_suf_sparse   = function() con_suf(XtX = Matrix::Matrix(ss$XtX, sparse = TRUE), Xty = ss$Xty, yty = ss$yty, n = ss$n),
  con_fit_ss_maf   = function() con_capture(susie_ss(ss$XtX, ss$Xty, ss$yty, ss$n, maf = con_maf, maf_thresh = 0.4)),
  con_fit_ss_neg   = function() con_capture(susie_ss(con_XtX_neg, ss$Xty, ss$yty, ss$n)),
  # individual_data_constructor: constant-column warning
  con_ind_obj      = function() con_ind(d1$X, d1$y),
  con_ind_const    = function() con_ind(con_X_const, d1$y),
  con_ind_inf      = function() con_ind(con_X_inf, d1$y),
  con_ind_dn       = function() con_ind(con_X_dn, d1$y),
  con_ind_scaled   = function() con_ind(scale(con_X_const[, 1:50]), d1$y),
  con_ind_na       = function() { X <- con_X_const; X[3, 3] <- NA; con_ind(X, d1$y) },
  con_ind_row2     = function() con_ind(con_X_const[1:2, ], d1$y[1:2]),
  con_ind_wide     = function() con_ind(cbind(d1$X[1:100, ], 7), d1$y[1:100]),
  con_ind_sparse   = function() con_ind(Matrix::Matrix(con_X_const, sparse = TRUE), d1$y),
  con_fit_ind_const = function() con_capture(susie(con_X_const, d1$y)),
  con_fit_ind_dn   = function() con_capture(susie(con_X_dn, d1$y, null_weight = 0.1))
))
