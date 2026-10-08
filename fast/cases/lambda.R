# susie_rss_lambda() configurations (work package "lambda").
lam_rd <- local({
  # rank-deficient R (200 samples, 300 variants) for lambda = 0 / "estimate"
  X <- ds$X[1:200, ]
  X <- X[, apply(X, 2, sd) > 0]
  r <- as.vector(cor(X, d1$y[1:200]))
  list(R = cor(X), z = r * sqrt(198 / (1 - r^2)), n = 200)
})
lam_pw <- { set.seed(5); w <- runif(ncol(ds$R)); w / sum(w) }
lam_zna <- { z <- ds$z; z[c(4, 150)] <- NA; z }
cases <- c(cases, list(
  rlam_p1000      = function() susie_rss_lambda(d1$z, d1$R, n = d1$n, lambda = 1e-3),
  rlam_erv_p1000  = function() susie_rss_lambda(d1$z, d1$R, n = d1$n, lambda = 1e-3,
                                                estimate_residual_variance = TRUE),
  rlam_refine     = function() susie_rss_lambda(ds$z, ds$R, n = ds$n, lambda = 1e-3, refine = TRUE),
  rlam_erv_refine = function() susie_rss_lambda(ds$z, ds$R, n = ds$n, lambda = 0.01,
                                                estimate_residual_variance = TRUE, refine = TRUE),
  rlam_greedy     = function() susie_rss_lambda(ds$z, ds$R, n = ds$n, lambda = 1e-3,
                                                L = 15, L_greedy = 5),
  rlam_greedy_erv = function() susie_rss_lambda(ds$z, ds$R, n = ds$n, lambda = 1e-3, L = 9,
                                                L_greedy = 3, estimate_residual_variance = TRUE),
  rlam_init_s2    = function() {
    f0 <- susie_rss_lambda(ds$z, ds$R, n = ds$n, lambda = 0.02, L = 3,
                           estimate_residual_variance = TRUE)
    susie_rss_lambda(ds$z, ds$R, n = ds$n, lambda = 0.02, L = 5, model_init = f0,
                     estimate_residual_variance = TRUE)
  },
  rlam_nullw_pw   = function() susie_rss_lambda(ds$z, ds$R, n = ds$n, lambda = 1e-3,
                                                null_weight = 0.1, prior_weights = lam_pw,
                                                estimate_residual_variance = TRUE),
  rlam_slot       = function() susie_rss_lambda(ds$z, ds$R, n = ds$n, lambda = 1e-3,
                                                slot_prior = slot_prior_betabinom()),
  rlam_X          = function() susie_rss_lambda(ds$z, X = ds$X, n = ds$n, lambda = 1e-3,
                                                estimate_residual_variance = TRUE),
  rlam_X_sparse   = function() susie_rss_lambda(ds$z, X = Matrix::Matrix(ds$X, sparse = TRUE),
                                                n = ds$n, lambda = 1e-3),
  rlam_R_sparse   = function() susie_rss_lambda(ds$z, Matrix::Matrix(ds$R, sparse = TRUE),
                                                n = ds$n, lambda = 1e-3),
  rlam_lam0_rd    = function() susie_rss_lambda(lam_rd$z, lam_rd$R, n = lam_rd$n, lambda = 0,
                                                estimate_residual_variance = TRUE),
  rlam_estlam     = function() susie_rss_lambda(lam_rd$z, lam_rd$R, n = lam_rd$n,
                                                lambda = "estimate", refine = TRUE),
  rlam_NAz        = function() susie_rss_lambda(lam_zna, ds$R, n = ds$n, lambda = "estimate",
                                                refine = TRUE),
  rlam_EM         = function() susie_rss_lambda(-ds$z, ds$R, n = ds$n, lambda = 1e-3,
                                                estimate_prior_method = "EM"),
  rlam_simple     = function() susie_rss_lambda(ds$z, ds$R, n = ds$n, lambda = 1e-3, L = 15,
                                                estimate_prior_method = "simple"),
  rlam_noest      = function() susie_rss_lambda(ds$z, ds$R, n = ds$n, lambda = 1e-3,
                                                estimate_prior_variance = FALSE),
  rlam_mixture    = function() susie_rss_lambda(ds$z, ds$R, n = ds$n, lambda = 1e-3,
                                                prior_variance_grid = c(0.1, 1, 10),
                                                mixture_weights = c(0.3, 0.4, 0.3)),
  rlam_chknull    = function() susie_rss_lambda(ds$z, ds$R, n = ds$n, lambda = 1e-3,
                                                check_null_threshold = 0.1, track_fit = TRUE),
  rlam_no_n_rv    = function() susie_rss_lambda(ds$z, ds$R, lambda = 0.05, residual_variance = 0.9,
                                                estimate_residual_variance = TRUE, max_iter = 20)
))
