# Matrix-vector products (work package "products"): direct dgemv, the
# reference-order kernels and threads, and the zero initial product.
with_opts <- function(opts, f) { old <- options(opts); on.exit(options(old)); f() }
Xu <- d1$X[, 1:400]                          # unstandardised crossprod: XtX scaled by csd != 1
ssu <- list(XtX = crossprod(Xu), Xty = as.vector(crossprod(Xu, d1$y)), yty = sum(d1$y^2), n = d1$n)
Xn <- d1$X[, 1:300]; colnames(Xn) <- paste0("v", 1:300)
Rn <- cor(Xn); zn <- d1$z[1:300]; names(zn) <- colnames(Xn)
Xc <- d1$X[, 1:200]; Xc[, 7] <- 1           # a constant column
init_ind <- { set.seed(5); susie(d1$X[, 1:300], d1$y, L = 5) }
init_rss <- { set.seed(5); susie_rss(d1$z, d1$R, n = d1$n, L = 5) }
zero_init <- local({ m <- init_rss; m$mu[] <- 0; m$mu2[] <- 0; m })
cases <- c(cases, list(
  pr_rss_threads4   = function() with_opts(list(susieR.threads = 4), function() susie_rss(d1$z, d1$R, n = d1$n)),
  pr_rss_p3000_t3   = function() with_opts(list(susieR.threads = 3), function() susie_rss(d2$z, d2$R, n = d2$n)),
  pr_ind_threads4   = function() with_opts(list(susieR.threads = 4), function() susie(d1$X, d1$y)),
  pr_ind_refine_t2  = function() with_opts(list(susieR.threads = 2), function() susie(d1$X[, 1:300], d1$y, refine = TRUE)),
  pr_rss_internal   = function() with_opts(list(matprod = "internal"), function() susie_rss(d1$z, d1$R, n = d1$n)),
  pr_rss_blas       = function() with_opts(list(matprod = "blas"), function() susie_rss(d1$z, d1$R, n = d1$n)),
  pr_ind_internal   = function() with_opts(list(matprod = "internal"), function() susie(d1$X[, 1:300], d1$y)),
  pr_ss_unstd       = function() susie_ss(ssu$XtX, ssu$Xty, ssu$yty, ssu$n),
  pr_ss_unstd_erv   = function() susie_ss(ssu$XtX, ssu$Xty, ssu$yty, ssu$n, estimate_residual_method = "MLE"),
  pr_ss_unstd_t4    = function() with_opts(list(susieR.threads = 4), function() susie_ss(ssu$XtX, ssu$Xty, ssu$yty, ssu$n)),
  pr_rss_lowrank    = function() susie_rss(d1$z, X = d1$X[1:600, ], n = d1$n),
  pr_rss_lowrank_t2 = function() with_opts(list(susieR.threads = 2), function() susie_rss(d1$z, X = d1$X[1:600, ], n = d1$n, estimate_residual_variance = TRUE)),
  pr_rss_names      = function() susie_rss(zn, Rn, n = d1$n),
  pr_rss_nullw      = function() susie_rss(d1$z, d1$R, n = d1$n, null_weight = 0.1, estimate_residual_variance = TRUE),
  pr_rss_na_z       = function() susie_rss(replace(d1$z, c(3, 50), NA), d1$R, n = d1$n),
  pr_rss_init       = function() susie_rss(d1$z, d1$R, n = d1$n, model_init = init_rss),
  pr_rss_init_zero  = function() susie_rss(d1$z, d1$R, n = d1$n, model_init = zero_init),
  pr_rss_p2         = function() susie_rss(d1$z[1:2], d1$R[1:2, 1:2], n = d1$n, L = 1),
  pr_ind_names      = function() susie(Xn, d1$y, L = 5),
  pr_ind_const      = function() susie(Xc, d1$y, L = 5),
  pr_ind_nullw      = function() susie(d1$X[, 1:300], d1$y, null_weight = 0.2),
  pr_ind_init       = function() susie(d1$X[, 1:300], d1$y, model_init = init_ind),
  pr_ind_nostd_int  = function() susie(d1$X[, 1:300], d1$y, standardize = FALSE),
  pr_ind_p2         = function() susie(d1$X[, 1:2], d1$y, L = 1),
  pr_ind_y_na       = function() susie(d1$X[, 1:300], replace(d1$y, 5, NA), na.rm = TRUE),
  pr_ind_inf_t2     = function() with_opts(list(susieR.threads = 2), function() susie(ds$X, ds$y, unmappable_effects = "inf")),
  pr_rss_inf_t2     = function() with_opts(list(susieR.threads = 2), function() susie_rss(ds$z, ds$R, n = ds$n, unmappable_effects = "inf")),
  pr_rss_ash_t2     = function() with_opts(list(susieR.threads = 2), function() susie_rss(ds$z, ds$R, n = ds$n, unmappable_effects = "ash")),
  pr_ind_ash_t2     = function() with_opts(list(susieR.threads = 2), function() susie(ds$X, ds$y, unmappable_effects = "ash")),
  pr_ss_NIG         = function() susie_ss(ss$XtX, ss$Xty, ss$yty, ss$n, L = 1, estimate_residual_method = "NIG"),
  pr_rss_inf_matrix = function() susie_rss(d1$z, replace(d1$R, 5, Inf), n = d1$n),
  pr_rss_huge       = function() susie_rss(d1$z[1:50], replace(d1$R[1:50, 1:50], c(2, 3), 1e308), n = d1$n)
))
