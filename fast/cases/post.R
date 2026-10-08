# Post-processing configurations: purity early exit (weak signals give
# diffuse credible sets that are dropped), z-score memo
# (compute_univariate_zscore with refine / greedy / null weight) and the
# opt-in parallel refine (options(susieR.refine_cores); serial on Windows).
post_y <- function(X, seed, k, sd) {
  set.seed(seed)
  Xs <- scale(X)
  j <- sample(ncol(X), k)
  drop(Xs[, j, drop = FALSE] %*% rnorm(k, sd = sd)) + rnorm(nrow(X))
}
post_z <- function(X, y) {
  r <- as.vector(cor(X, y))
  r * sqrt((nrow(X) - 2) / (1 - r^2))
}
post_w1 <- post_y(d1$X, 2, 2, 0.05)          # weak: diffuse CSs, all dropped
post_w2 <- post_y(ds$X, 1, 2, 0.05)
post_s3 <- post_y(ds$X, 5, 3, 0.15)          # three signals, for refine
post_zw1 <- post_z(d1$X, post_w1)
post_zw2 <- post_z(ds$X, post_w2)
post_ssw <- list(XtX = crossprod(d1$X), Xty = as.vector(crossprod(d1$X, post_w1)),
                 yty = sum(post_w1^2), n = nrow(d1$X))
post_Xc <- ds$X; post_Xc[, 150] <- 1          # constant column
post_cores <- function(expr) {
  op <- options(susieR.refine_cores = as.integer(Sys.getenv("SUSIE_REFINE_CORES", "3")))
  on.exit(options(op))
  expr
}

cases <- c(cases, list(
  # purity early exit: scaled_XtX, dense Xcorr, low-rank X, ss, X, sparse X
  post_rss_weak     = function() susie_rss(post_zw1, d1$R, n = d1$n),
  post_rss_weak_L15 = function() susie_rss(post_zw1, d1$R, n = d1$n, L = 15),
  post_rss_weak_non = function() susie_rss(post_zw1, d1$R),
  post_rss_weak_med = function() susie_rss(post_zw1, d1$R, n = d1$n, median_abs_corr = 0.5),
  post_rss_weak_ref = function() susie_rss(post_zw1, d1$R, n = d1$n, refine = TRUE),
  post_rss_weak_nw  = function() susie_rss(post_zw1, d1$R, n = d1$n, null_weight = 0.1),
  post_rss_weak_ext = function() susie_rss(post_zw1, d1$R, n = d1$n, cs_extension_corr = 0.9),
  post_rss_lowrank  = function() susie_rss(post_zw1, X = d1$X[1:600, ], n = d1$n),
  post_rss_small    = function() susie_rss(post_zw2, ds$R, n = ds$n),
  post_ss_weak      = function() susie_ss(post_ssw$XtX, post_ssw$Xty, post_ssw$yty, post_ssw$n),
  post_ind_weak     = function() susie(ds$X, post_w2),
  post_ind_weak_np  = function() susie(ds$X, post_w2, n_purity = 100),
  post_ind_weak_spr = function() susie(Matrix::Matrix(ds$X, sparse = TRUE), post_w2),
  post_ind_weak_nw  = function() susie(ds$X, post_w2, null_weight = 0.1),
  post_ind_weak_mac = function() susie(ds$X, post_w2, min_abs_corr = 0.05),
  post_ind_constcol = function() susie(post_Xc, post_w2),
  post_ind_ref_np   = function() susie(ds$X, ds$y, refine = TRUE, n_purity = 50),
  post_ind_p1000_wk = function() susie(d1$X, post_w1),
  # z-score memo
  post_z_refine     = function() susie(ds$X, ds$y, refine = TRUE, compute_univariate_zscore = TRUE),
  post_z_refine_nw  = function() susie(ds$X, ds$y, refine = TRUE, null_weight = 0.1, compute_univariate_zscore = TRUE),
  post_z_greedy     = function() susie(ds$X, post_s3, L = 12, L_greedy = 4, compute_univariate_zscore = TRUE),
  post_z_nostd      = function() susie(ds$X, ds$y, standardize = FALSE, intercept = FALSE, refine = TRUE, compute_univariate_zscore = TRUE),
  post_z_pw         = function() susie(ds$X, ds$y, prior_weights = pw[idx] / sum(pw[idx]), refine = TRUE, compute_univariate_zscore = TRUE),
  post_z_sparse     = function() susie(Matrix::Matrix(ds$X, sparse = TRUE), ds$y, refine = TRUE, compute_univariate_zscore = TRUE),
  post_z_constcol   = function() susie(post_Xc, ds$y, refine = TRUE, compute_univariate_zscore = TRUE),
  # parallel refine (opt-in; candidates in forked processes on Linux/macOS)
  post_par_rss      = function() post_cores(susie_rss(d1$z, d1$R, n = d1$n, refine = TRUE)),
  post_par_rss_weak = function() post_cores(susie_rss(post_zw1, d1$R, n = d1$n, L = 15, refine = TRUE)),
  post_par_ss       = function() post_cores(susie_ss(ss$XtX, ss$Xty, ss$yty, ss$n, refine = TRUE)),
  post_par_ind      = function() post_cores(susie(ds$X, post_s3, refine = TRUE)),
  post_par_ind_z    = function() post_cores(susie(ds$X, post_s3, refine = TRUE, compute_univariate_zscore = TRUE)),
  post_par_ind_np   = function() post_cores(susie(ds$X, post_s3, refine = TRUE, n_purity = 20)),
  post_par_ind_nw   = function() post_cores(susie(ds$X, post_s3, refine = TRUE, null_weight = 0.1)),
  post_par_ind_mi   = function() post_cores(susie(ds$X, post_s3, refine = TRUE, max_iter = 3))
))
