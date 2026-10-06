# Usage: Rscript run_suite.R <lib> <fast-mode|-> <out.rds> [filter-regex]
# Fits 35 configurations on data/ made by make_data.R (seeds 1 and 2) and
# saves the fits; compare two runs with compare.R.
args <- commandArgs(TRUE)
lib <- args[1]; mode <- args[2]; out <- args[3]
filt <- if (length(args) > 3) args[4] else "."
if (mode != "-") options(susieR.fast = mode)
suppressPackageStartupMessages(library(susieR, lib.loc = lib))
d1 <- readRDS("data/region_n3000_p1000_s1.rds")
d2 <- readRDS("data/region_n3000_p3000_s2.rds")
# small region for slow configurations
set.seed(7); idx <- 301:600
ds <- list(X = d1$X[, idx], y = d1$y, n = d1$n)
ds$R <- cor(ds$X); r <- as.vector(cor(ds$X, ds$y)); ds$z <- r * sqrt((ds$n - 2) / (1 - r^2))
data(N3finemapping, package = "susieR", lib.loc = lib)
N3 <- N3finemapping
Xs <- scale(d1$X); ss <- list(XtX = crossprod(Xs), Xty = as.vector(crossprod(Xs, d1$y - mean(d1$y))),
                             yty = sum((d1$y - mean(d1$y))^2), n = d1$n)
uni <- susieR::univariate_regression(d1$X, d1$y)
pw <- { set.seed(3); w <- runif(ncol(d1$X)); w / sum(w) }

cases <- list(
  rss_p1000        = function() susie_rss(d1$z, d1$R, n = d1$n),
  rss_p3000        = function() susie_rss(d2$z, d2$R, n = d2$n),
  rss_erv          = function() susie_rss(d1$z, d1$R, n = d1$n, estimate_residual_variance = TRUE),
  rss_L20_pw       = function() susie_rss(d1$z, d1$R, n = d1$n, L = 20, prior_weights = pw),
  rss_nullw        = function() susie_rss(d1$z, d1$R, n = d1$n, null_weight = 0.2),
  rss_EM           = function() susie_rss(d1$z, d1$R, n = d1$n, estimate_prior_method = "EM"),
  rss_simple       = function() susie_rss(d1$z, d1$R, n = d1$n, estimate_prior_method = "simple"),
  rss_noest        = function() susie_rss(d1$z, d1$R, n = d1$n, estimate_prior_variance = FALSE),
  rss_refine       = function() susie_rss(d1$z, d1$R, n = d1$n, refine = TRUE),
  rss_pipconv      = function() susie_rss(d1$z, d1$R, n = d1$n, convergence_method = "pip"),
  rss_no_n         = function() susie_rss(d1$z, d1$R),
  rss_bhat_shat    = function() susie_rss(bhat = uni$betahat, shat = uni$sebetahat, R = d1$R, n = d1$n, var_y = var(d1$y)),
  rss_greedy       = function() susie_rss(d1$z, d1$R, n = d1$n, L = 15, L_greedy = 5),
  rss_chknull      = function() susie_rss(d1$z, d1$R, n = d1$n, check_null_threshold = 0.1),
  rss_track        = function() susie_rss(d1$z, d1$R, n = d1$n, track_fit = TRUE),
  rss_small_inf    = function() susie_rss(ds$z, ds$R, n = ds$n, unmappable_effects = "inf"),
  rss_small_ash    = function() susie_rss(ds$z, ds$R, n = ds$n, unmappable_effects = "ash"),
  rss_small_mism   = function() susie_rss(ds$z, ds$R, n = ds$n, R_mismatch = "eb", R_finite = 500),
  rss_small_slot   = function() susie_rss(ds$z, ds$R, n = ds$n, slot_prior = slot_prior_betabinom()),
  ss_p1000         = function() susie_ss(ss$XtX, ss$Xty, ss$yty, ss$n),
  ss_erv_MLE       = function() susie_ss(ss$XtX, ss$Xty, ss$yty, ss$n, estimate_residual_method = "MLE"),
  ind_p1000        = function() susie(d1$X, d1$y),
  ind_p3000        = function() susie(d2$X, d2$y),
  ind_nostd        = function() susie(d1$X, d1$y, standardize = FALSE, intercept = FALSE),
  ind_EM           = function() susie(d1$X, d1$y, estimate_prior_method = "EM"),
  ind_refine       = function() susie(d1$X, d1$y, refine = TRUE),
  ind_noerv        = function() susie(d1$X, d1$y, estimate_residual_variance = FALSE),
  ind_L1           = function() susie(d1$X, d1$y, L = 1),
  ind_greedy       = function() susie(d1$X, d1$y, L = 12, L_greedy = 4),
  ind_N3           = function() susie(N3$X, N3$Y[, 1], L = 10),
  ind_N3_NIG       = function() susie(N3$X, N3$Y[, 1], L = 1, estimate_residual_method = "NIG"),
  ind_small_inf    = function() susie(ds$X, ds$y, unmappable_effects = "inf"),
  ind_small_ash    = function() susie(ds$X, ds$y, unmappable_effects = "ash"),
  ind_small_slot   = function() susie(ds$X, ds$y, slot_prior = slot_prior_poisson(C = 3)),
  ind_sparse       = function() susie(Matrix::Matrix(d1$X, sparse = TRUE), d1$y)
)
cases <- cases[grepl(filt, names(cases))]
res <- list()
for (nm in names(cases)) {
  set.seed(1)
  t0 <- proc.time()[[3]]
  fit <- tryCatch(suppressWarnings(suppressMessages(cases[[nm]]())),
                  error = function(e) structure(list(msg = conditionMessage(e)), class = "err"))
  el <- proc.time()[[3]] - t0
  res[[nm]] <- list(fit = fit, time = el)
  cat(sprintf("%-16s %7.2fs %s\n", nm, el, if (inherits(fit, "err")) paste("ERROR:", fit$msg) else paste("niter", fit$niter)))
}
saveRDS(res, out)
