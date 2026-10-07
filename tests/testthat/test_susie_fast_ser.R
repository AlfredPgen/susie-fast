context("susie-fast SER code paths")

# Fits in modes exact and full are compared with the upstream code path
# (options(susieR.fast = "off")) on the same input.
ser_fit_modes <- function(fun) {
  old <- options(susieR.fast = "off")
  on.exit(options(old))
  set.seed(1); ref <- suppressWarnings(suppressMessages(fun()))
  options(susieR.fast = "exact")
  set.seed(1); exact <- suppressWarnings(suppressMessages(fun()))
  options(susieR.fast = "full")
  set.seed(1); full <- suppressWarnings(suppressMessages(fun()))
  list(ref = ref, exact = exact, full = full)
}

ser_test_data <- function() {
  set.seed(42)
  n <- 400; p <- 120
  X <- matrix(as.double(rbinom(n * p, 2, 0.3)), n, p)
  for (j in 2:p) X[, j] <- ifelse(runif(n) < 0.7, X[, j - 1], X[, j])
  X <- X[, apply(X, 2, sd) > 0]
  b <- rep(0, ncol(X)); b[c(10, 60)] <- c(0.5, -0.4)
  y <- as.vector(X %*% b + rnorm(n))
  R <- cor(X)
  r <- as.vector(cor(X, y))
  z <- r * sqrt((n - 2) / (1 - r^2))
  list(X = X, y = y, R = R, z = z, n = n)
}

# Counts calls of the namespace functions in `fns` while `expr` runs.
ser_count_calls <- function(fns, expr) {
  ns <- asNamespace("susieR")
  counter <- new.env()
  for (f in fns) {
    counter[[f]] <- 0
    suppressMessages(trace(f, bquote(assign(.(f), get(.(f), .(counter)) + 1,
                                            envir = .(counter))),
                           print = FALSE, where = ns))
  }
  on.exit(for (f in fns) suppressMessages(untrace(f, where = ns)))
  force(expr)
  mget(fns, counter)
}

test_that("the SER kernels reproduce the R arithmetic bit for bit", {
  expect_true(susieR:::fast_ser_self_test())
  expect_true(susieR:::fast_ser_ok())
})

test_that("SER fast paths give fits identical to the upstream code path", {
  d <- ser_test_data()
  Rd <- d$R; diag(Rd) <- 1 + (seq_len(nrow(Rd)) %% 3) * 1e-3
  Xc <- scale(d$X, scale = FALSE)
  XtX <- crossprod(Xc); Xty <- as.vector(crossprod(Xc, d$y - mean(d$y)))
  yty <- sum((d$y - mean(d$y))^2)
  pw <- rep(1, ncol(d$X)); pw[1:10] <- 0; pw <- pw / sum(pw)
  cases <- list(
    function() susie_rss(d$z, d$R, n = d$n, L = 15),
    function() susie_rss(d$z * 6, d$R, n = d$n),
    function() susie_rss(d$z, Rd, n = d$n),
    function() susie_rss(d$z, d$R, n = d$n, L = 12, null_weight = 0.4),
    function() susie_rss(d$z, d$R, n = d$n, prior_weights = pw),
    function() susie_rss(d$z, d$R, n = d$n, estimate_prior_method = "simple"),
    function() susie_rss(rep(0, length(d$z)), d$R, n = d$n),
    function() susie_rss(d$z, d$R, n = d$n, L = 6, refine = TRUE),
    function() susie_rss(d$z, d$R, n = d$n, R_mismatch = "eb", R_finite = 200),
    function() susie_rss(d$z * 0.3, d$R, n = d$n, L = 4, R_mismatch = "eb_mix",
                         R_finite = 150),
    function() susie_ss(XtX, Xty, yty, d$n, L = 12),
    function() susie_ss(XtX, Xty, yty, d$n, standardize = FALSE),
    function() susie(d$X, d$y, L = 15),
    function() susie(d$X, d$y, standardize = FALSE, intercept = FALSE),
    function() susie(d$X, d$y, estimate_prior_method = "EM"),
    function() susie(d$X, d$y, L = 6, L_greedy = 2))
  for (fun in cases) {
    m <- ser_fit_modes(fun)
    expect_identical(m$exact, m$ref)
    m$full$elbo <- m$ref$elbo <- NULL
    expect_identical(m$full, m$ref)
  }
})

test_that("repeated SERs and X'r products are reused", {
  d <- ser_test_data()
  old <- options(susieR.fast = "exact")
  on.exit(options(old))
  n <- ser_count_calls(c("single_effect_regression", "compute_Xty"),
                       fit <- suppressWarnings(susie(d$X, d$y, L = 15)))
  expect_lt(n$single_effect_regression, 15 * fit$niter)
  expect_lt(n$compute_Xty, 15 * fit$niter)
  n <- ser_count_calls("single_effect_regression",
                       fit <- suppressWarnings(susie_rss(d$z, d$R, n = d$n, L = 15)))
  expect_lt(n$single_effect_regression, 15 * fit$niter)
})

test_that("the SER memo reproduces R_bf_attenuation rows exactly", {
  d <- ser_test_data()
  old <- options(susieR.fast = "exact")
  on.exit(options(old))
  ns <- asNamespace("susieR")
  # Capture a mid-fit state of an R_finite fit.
  cap <- new.env()
  suppressMessages(trace("fast_single_effect_regression",
    bquote(if (is.null(.(cap)$model) && l == 2 && !is.null(model$runtime)) {
      assign("data", data, envir = .(cap)); assign("params", params, envir = .(cap))
      assign("model", model, envir = .(cap)) }),
    print = FALSE, where = ns))
  fit <- suppressWarnings(suppressMessages(
    susie_rss(d$z, cor(d$X[1:100, ]), n = d$n, R_mismatch = "eb", R_finite = 100)))
  suppressMessages(untrace("fast_single_effect_regression", where = ns))
  expect_false(is.null(cap$model))
  for (infl in c(1, 1.5)) {
    m <- cap$model
    L <- nrow(m$alpha); p <- ncol(m$alpha)
    # With inflation exactly 1 (all other effects null) no row is recorded,
    # so rows left from earlier sweeps must stay as they are.
    m$shat2_inflation  <- rep(infl, p)
    m$R_bf_attenuation <- matrix(seq_len(L * p) / 7, L, p)
    m$V[2] <- m$V[3] <- 0
    m$runtime$fast_cache <- ns$fast_cache_new(cap$data, cap$params, L)
    n <- ser_count_calls("single_effect_regression", {
      a2 <- ns$fast_single_effect_regression(cap$data, cap$params, m, 2)
      a3 <- ns$fast_single_effect_regression(cap$data, cap$params, a2, 3)
    })
    expect_equal(n$single_effect_regression, 1)
    b2 <- ns$single_effect_regression(cap$data, cap$params, m, 2)
    b3 <- ns$single_effect_regression(cap$data, cap$params, b2, 3)
    a3$runtime <- b3$runtime <- NULL
    expect_identical(a3, b3)
  }
})

test_that("the SER memo is off when a generic has a non-susieR method", {
  d <- ser_test_data()
  expect_true(susieR:::fast_own_method("pre_loglik_prior_hook", "ss"))
  assign("pre_loglik_prior_hook.ss", function(data, params, model, ser_stats, l, V_init)
    susieR:::pre_loglik_prior_hook.default(data, params, model, ser_stats, l, V_init),
    envir = globalenv())
  on.exit(rm("pre_loglik_prior_hook.ss", envir = globalenv()))
  expect_false(susieR:::fast_own_method("pre_loglik_prior_hook", "ss"))
  old <- options(susieR.fast = "exact")
  on.exit(options(old), add = TRUE)
  n <- ser_count_calls("single_effect_regression",
                       fit <- suppressWarnings(susie_rss(d$z, d$R, n = d$n, L = 15)))
  expect_equal(n$single_effect_regression, 15 * fit$niter)
})
