context("susie-fast: SuSiE-ash")

ash_test_data <- function() {
  set.seed(7)
  n <- 300; p <- 80
  X <- matrix(as.double(rbinom(n * p, 2, 0.3)), n, p)
  for (j in 2:p) X[, j] <- ifelse(runif(n) < 0.7, X[, j - 1], X[, j])
  X <- X[, apply(X, 2, sd) > 0]
  b <- rep(0, ncol(X)); b[c(10, 50)] <- c(0.6, -0.5)
  b <- b + rnorm(ncol(X), sd = 0.03)
  y <- as.vector(X %*% b + rnorm(n))
  R <- cor(X)
  r <- as.vector(cor(X, y))
  z <- r * sqrt((n - 2) / (1 - r^2))
  list(X = X, y = y, R = R, z = z, n = n)
}

# Fits fun() under "off", "exact" and "full". The ash diagnostics live in an
# environment (.diag_env), which identical() compares by reference, so it is
# split off and its history compared separately.
ash_fit_modes <- function(fun) {
  old <- options(susieR.fast = "off")
  on.exit(options(old))
  run <- function(mode) {
    options(susieR.fast = mode)
    set.seed(1)
    fit <- suppressWarnings(suppressMessages(fun()))
    hist <- if (is.environment(fit$.diag_env)) fit$.diag_env$history
    fit$.diag_env <- NULL
    list(fit = fit, history = hist)
  }
  list(ref = run("off"), exact = run("exact"), full = run("full"))
}

test_that("fast mr.ash.rss kernel reproduces mr_ash_rss_cpp bit for bit", {
  expect_true(susieR:::fast_mr_ash_rss_self_test())
})

test_that("fast mr.ash.rss kernel matches on edge-case inputs", {
  d <- ash_test_data()
  p <- length(d$z); K <- 25
  s0 <- (2^((0:24) / 25) - 1)^2 / (d$n - 1) * d$n
  same <- function(..., z = d$z, R = d$R, shat = rep(1 / sqrt(d$n), p),
                   var_y = 1, s0. = s0, w0 = rep(1 / K, K), mu1 = numeric(0),
                   max_iter = 200L, flags = c(TRUE, TRUE, TRUE, FALSE))
    susieR:::fast_mr_ash_rss_same(z * shat, shat, z, R, var_y, d$n, 1, s0.,
                                  w0, mu1, 1e-4, max_iter, flags[1], flags[2],
                                  flags[3], flags[4])
  expect_true(same())
  expect_true(same(flags = c(TRUE, TRUE, TRUE, TRUE), mu1 = d$z / 50))
  expect_true(same(flags = c(FALSE, FALSE, FALSE, FALSE)))
  expect_true(same(var_y = Inf))
  expect_true(same(s0. = replace(s0, 3, 0), w0 = c(0, rep(1 / 24, 24))))
  expect_true(same(s0. = 0.1, w0 = 1))
  expect_true(same(z = d$z[1], R = d$R[1, 1, drop = FALSE],
                   shat = 1 / sqrt(d$n)))
  # Non-finite input runs to max_iter (and prints the upstream notice).
  expect_true(same(z = replace(d$z, 3, NaN), max_iter = 4L))
  R0 <- d$R; R0[4, 4] <- 0
  expect_true(same(R = R0, max_iter = 4L))
  expect_true(same(shat = replace(rep(1 / sqrt(d$n), p), 6, Inf), max_iter = 4L))
})

test_that("mr.ash.rss updates w0 in place as upstream does", {
  d <- ash_test_data()
  p <- length(d$z); K <- 10
  call <- function(mode) {
    old <- options(susieR.fast = mode); on.exit(options(old))
    w0 <- rep(1 / K, K) + 0
    fit <- mr.ash.rss(d$z / sqrt(d$n), rep(1 / sqrt(d$n), p), d$R, 1, d$n,
                      s0 = (1:K) / 1000, w0 = w0, tol = 1e-4)
    list(fit = fit, w0 = w0)
  }
  ref <- call("off")
  expect_false(identical(ref$w0, rep(1 / K, K)))
  expect_identical(call("exact"), ref)
  expect_identical(call("full"), ref)
})

test_that("SuSiE-ash fits are bit-identical to the upstream code path", {
  skip_on_cran()
  d <- ash_test_data()
  Xs <- scale(d$X)
  yc <- d$y - mean(d$y)
  R01 <- d$R; R01[5, ] <- 0; R01[, 5] <- 0; R01[5, 5] <- 0
  cases <- list(
    function() susie_rss(d$z, d$R, n = d$n, unmappable_effects = "ash"),
    function() susie_rss(d$z, d$R, n = d$n, unmappable_effects = "ash",
                         estimate_residual_variance = TRUE),
    function() susie_rss(d$z, d$R, n = d$n,
                         unmappable_effects = "ash_filter_archived",
                         estimate_residual_variance = TRUE),
    function() susie_ss(crossprod(Xs), as.vector(crossprod(Xs, yc)),
                        sum(yc^2), d$n, unmappable_effects = "ash"),
    # diag(XtX) all 0 or 1: get_xcorr() returns XtX, mr.ash.rss gets
    # safe_cov2cor(XtX), which differs in the zero-variance variant.
    function() susie_ss(R01, d$z / sqrt(d$n), 1, d$n, standardize = FALSE,
                        unmappable_effects = "ash", max_iter = 4,
                        estimate_residual_variance = TRUE),
    function() susie_ss(R01, d$z / sqrt(d$n), 1, d$n, standardize = FALSE,
                        unmappable_effects = "ash_filter_archived",
                        max_iter = 4, estimate_residual_variance = TRUE),
    function() susie(d$X, d$y, unmappable_effects = "ash"),
    function() susie(d$X, d$y, unmappable_effects = "ash_filter_archived"),
    function() { X <- d$X; X[, 3] <- 1; susie(X, d$y, unmappable_effects = "ash") })
  for (fun in cases) {
    m <- ash_fit_modes(fun)
    expect_identical(m$exact, m$ref)
    expect_identical(m$full$fit[names(m$full$fit) != "elbo"],
                     m$ref$fit[names(m$ref$fit) != "elbo"])
    expect_identical(m$full$history, m$ref$history)
    expect_null(m$exact$fit$runtime)
  }
})

test_that("the cached LD adjacency is keyed on the threshold", {
  d <- ash_test_data()
  old <- options(susieR.fast = "full"); on.exit(options(old))
  cache <- new.env(parent = emptyenv())
  model <- list(runtime = list(fast_cache = cache))
  data <- structure(list(XtX = crossprod(scale(d$X))), class = "ss")
  Xc <- susieR:::fast_get_xcorr(data, model)$Xcorr
  expect_identical(Xc, susieR:::get_xcorr(data)$Xcorr)
  expect_identical(susieR:::fast_ash_R(data, model), susieR:::safe_cov2cor(data$XtX))
  pip <- runif(ncol(Xc)); pip[2] <- NA
  for (t in c(0.5, 0.3, 0.5, 0.95)) {
    A <- susieR:::fast_ld_adj(Xc, t, model)
    expect_identical(as.vector(A %*% pip), as.vector((abs(Xc) > t) %*% pip))
    expect_identical(A == 1, abs(Xc) > t)
  }
  # A matrix other than the cached one is not served from the cache.
  expect_identical(susieR:::fast_ld_adj(Xc * 0.5, 0.3, model), abs(Xc * 0.5) > 0.3)
})

# Result (or error message) and warnings of fun() under "off", "exact" and
# "full".
ash_conds_modes <- function(fun) {
  old <- options(susieR.fast = "off")
  on.exit(options(old))
  run <- function(mode) {
    options(susieR.fast = mode)
    w <- character(0)
    r <- tryCatch(withCallingHandlers(fun(), warning = function(c) {
           w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") }),
         error = function(e) c(conditionMessage(e), deparse(conditionCall(e))))
    if (is.list(r)) r$.diag_env <- NULL
    list(result = r, warnings = w)
  }
  list(ref = run("off"), exact = run("exact"), full = run("full"))
}

# Near-constant column whose centred sum of squares rounds to a negative
# number, next to a constant column: safe_cor() takes its fallback and
# sqrt() warns.
ash_negcss_X <- function(X) {
  X <- X[1:50, ]
  X[, 2] <- 0
  X[, 4] <- 5.7712482971837744
  X[1, 4] <- 5.7712482971837593
  X
}

test_that("a correlation matrix that warned is recomputed, warning again", {
  d <- ash_test_data()
  Xc <- scale(d$X, scale = FALSE)
  XtX <- crossprod(Xc); XtX[7, 7] <- -XtX[7, 7]
  X <- ash_negcss_X(d$X)
  expect_lt(sum(X[, 4]^2) - 50 * mean(X[, 4])^2, 0)
  conds <- function(f) {
    w <- character(0)
    v <- withCallingHandlers(f(), warning = function(c) {
      w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") })
    list(v, w)
  }
  ss <- structure(list(XtX = XtX), class = "ss")
  ind <- structure(list(X = X), class = "individual")
  model <- list(runtime = list(fast_cache = new.env(parent = emptyenv())))
  # Two ash iterations: get_xcorr() then (ss only) safe_cov2cor().
  ss_fast <- conds(function() for (i in 1:2)
    list(susieR:::fast_get_xcorr(ss, model), susieR:::fast_ash_R(ss, model)))
  ss_ref <- conds(function() for (i in 1:2)
    list(susieR:::get_xcorr(ss), susieR:::safe_cov2cor(ss$XtX)))
  expect_length(ss_ref[[2]], 4)
  expect_identical(ss_fast, ss_ref)
  model <- list(runtime = list(fast_cache = new.env(parent = emptyenv())))
  ind_fast <- conds(function() lapply(1:2, function(i) susieR:::fast_get_xcorr(ind, model)))
  ind_ref <- conds(function() lapply(1:2, function(i) susieR:::get_xcorr(ind)))
  expect_length(ind_ref[[2]], 2)
  expect_identical(ind_fast, ind_ref)
})

test_that("SuSiE-ash warnings and errors match the upstream code path", {
  d <- ash_test_data()
  Xc <- scale(d$X, scale = FALSE); yc <- d$y - mean(d$y)
  XtX <- crossprod(Xc); j <- order(abs(d$z))[1]
  XtX[j, j] <- -XtX[j, j]
  Xty <- as.vector(crossprod(Xc, yc))
  X <- ash_negcss_X(d$X)
  cases <- list(
    function() susie_ss(XtX, Xty, sum(yc^2), d$n, standardize = FALSE,
                        unmappable_effects = "ash",
                        estimate_residual_variance = TRUE),
    function() susie_ss(XtX, Xty, sum(yc^2), d$n, standardize = FALSE,
                        unmappable_effects = "ash_filter_archived"),
    function() susie(X, d$y[1:50], unmappable_effects = "ash"))
  for (fun in cases) {
    m <- ash_conds_modes(fun)
    expect_gt(length(m$ref$warnings), 0)
    expect_identical(m$exact, m$ref)
    expect_identical(m$full$warnings, m$ref$warnings)
    if (is.character(m$ref$result))
      expect_identical(m$full$result, m$ref$result)
  }
  # Sparse X with a zero column: the correlation matrix is an S4 Matrix,
  # which upstream passes on until caisa_cpp rejects it.
  Xs <- d$X; Xs[, 9] <- 0
  Xs <- Matrix::Matrix(Xs, sparse = TRUE)
  for (ue in c("ash", "ash_filter_archived")) {
    m <- ash_conds_modes(function() susie(Xs, d$y, unmappable_effects = ue))
    expect_identical(m$exact, m$ref)
    expect_identical(m$full, m$ref)
  }
})

test_that("the mr.ash.rss self-test leaves the random number stream alone", {
  old <- RNGkind(normal.kind = "Box-Muller")
  on.exit(RNGkind(old[1], old[2], old[3]))
  set.seed(7); rnorm(1)
  ref <- rnorm(4)
  set.seed(7); rnorm(1)
  expect_true(susieR:::fast_mr_ash_rss_self_test())
  expect_identical(rnorm(4), ref)
})
