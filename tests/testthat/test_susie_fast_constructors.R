context("susie-fast data constructors")

# Each constructor call is run with options(susieR.fast = "off") (upstream
# code) and in the fast modes; the result and every message, warning and
# error, in order, must be identical.
con_run <- function(fun) {
  msgs <- character(0)
  val <- withCallingHandlers(
    tryCatch(fun(), error = function(e) paste("E:", conditionMessage(e))),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage")
    },
    warning = function(w) {
      msgs <<- c(msgs, paste("W:", conditionMessage(w)))
      invokeRestart("muffleWarning")
    })
  list(val = val, msgs = msgs)
}

con_modes <- function(fun) {
  old <- options(susieR.fast = "off")
  on.exit(options(old))
  ref <- con_run(fun)
  options(susieR.fast = "exact")
  exact <- con_run(fun)
  options(susieR.fast = "full")
  full <- con_run(fun)
  list(ref = ref, exact = exact, full = full)
}

# In full mode the ELBO trace of a fit may differ in its last digits.
expect_same_modes <- function(fun) {
  r <- con_modes(fun)
  expect_identical(r$exact, r$ref)
  no_elbo <- function(x) {
    if (inherits(x$val, "susie")) x$val$elbo <- NULL
    x
  }
  expect_identical(no_elbo(r$full), no_elbo(r$ref))
  invisible(r)
}

con_test_data <- function() {
  set.seed(7)
  n <- 300; p <- 90
  X <- matrix(as.double(rbinom(n * p, 2, 0.3)), n, p)
  for (j in 2:p) X[, j] <- ifelse(runif(n) < 0.7, X[, j - 1], X[, j])
  X <- X[, apply(X, 2, sd) > 0]
  b <- rep(0, ncol(X)); b[c(10, 40)] <- c(0.5, -0.4)
  y <- as.vector(X %*% b + rnorm(n))
  R <- cor(X)
  r <- as.vector(cor(X, y))
  z <- r * sqrt((n - 2) / (1 - r^2))
  uni <- univariate_regression(X, y)
  Xs <- scale(X)
  XtX <- crossprod(Xs)
  Xty <- as.vector(crossprod(Xs, y - mean(y)))
  list(X = X, y = y, R = R, z = z, n = n, uni = uni, XtX = XtX, Xty = Xty,
       yty = sum((y - mean(y))^2))
}

test_that("the constructor kernels pass their self-test", {
  expect_true(susieR:::fast_constructors_self_test())
  expect_true(susieR:::fast_constructors_ok())
})

test_that("xtx_check_cpp classifies symmetry and non-finite entries", {
  M <- crossprod(matrix(rnorm(60), 10, 6))
  expect_identical(susieR:::xtx_check_cpp(M), 3L)
  Z <- M; Z[1, 2] <- 0; Z[2, 1] <- -0
  expect_identical(susieR:::xtx_check_cpp(Z), 2L)
  A <- M; A[1, 6] <- A[1, 6] + 1
  expect_identical(susieR:::xtx_check_cpp(A), 1L)
  for (bad in list(NA_real_, NaN, Inf, -Inf)) {
    B <- M; B[3, 2] <- bad
    expect_identical(susieR:::xtx_check_cpp(B), 0L)
    B <- A; B[6, 6] <- bad
    expect_identical(susieR:::xtx_check_cpp(B), 0L)
  }
})

test_that("summary_stats_constructor (z, R) matches upstream", {
  d <- con_test_data()
  p <- ncol(d$R)
  R_dn <- d$R
  dimnames(R_dn) <- list(paste0("v", 1:p), paste0("v", 1:p))
  attr(R_dn, "foo") <- "bar"
  R_cn <- d$R; colnames(R_cn) <- paste0("v", 1:p)
  R_asym <- d$R; R_asym[2, 5] <- R_asym[2, 5] + 0.01
  R_z <- d$R; R_z[3, 4] <- 0; R_z[4, 3] <- -0
  R_na <- d$R; R_na[2, 2] <- NA
  R_nan <- d$R; R_nan[1, 3] <- NaN; R_nan[3, 1] <- NaN
  R_inf <- d$R; R_inf[1, 3] <- Inf; R_inf[3, 1] <- Inf
  R_diag <- d$R; diag(R_diag) <- 1 + (1:p %% 5) * 1e-3
  maf <- seq(0.01, 0.5, length.out = p)
  ssc <- function(...) susieR:::summary_stats_constructor(...)
  r <- expect_same_modes(function() ssc(z = d$z, R = d$R, n = d$n))
  # diag(R) == 1: the standardised XtX is (n - 1) * R itself.
  expect_identical(unclass(r$full$val$data$XtX)[1:5],
                   ((d$n - 1) * d$R)[1:5])
  for (R in list(R_dn, R_cn, R_asym, R_z, R_na, R_nan, R_inf, R_diag)) {
    expect_same_modes(function() ssc(z = d$z, R = R, n = d$n))
    expect_same_modes(function() ssc(z = d$z, R = R, n = d$n, null_weight = 0.2))
  }
  expect_same_modes(function() ssc(z = d$z, R = R_dn, n = d$n, maf = maf, maf_thresh = 0.2))
  expect_same_modes(function() ssc(z = d$z, R = R_dn, n = d$n, maf = c(0.4, rep(0.01, p - 1)), maf_thresh = 0.2))
  expect_same_modes(function() ssc(z = d$z, R = R_dn))
  expect_same_modes(function() ssc(z = d$z, R = d$R, n = d$n, standardize = FALSE))
  expect_same_modes(function() ssc(z = d$z, R = d$R, n = d$n, check_input = TRUE))
  r <- expect_same_modes(function() ssc(z = d$z, R = R_asym, n = d$n))
  expect_true(any(grepl("not symmetric", r$full$msgs)))
  r <- expect_same_modes(function() ssc(z = d$z, R = R_na, n = d$n))
  expect_match(r$full$val, "contains NAs")
})

test_that("original-scale XtX (bhat, shat, var_y) matches upstream", {
  d <- con_test_data()
  p <- ncol(d$R)
  bh <- setNames(d$uni$betahat, paste0("v", 1:p))
  sh <- setNames(d$uni$sebetahat, paste0("v", 1:p))
  R_dn <- d$R
  dimnames(R_dn) <- list(paste0("v", 1:p), NULL)
  attr(R_dn, "eigen") <- list(values = 1)
  R_asym <- d$R; R_asym[2, 5] <- R_asym[2, 5] + 0.01
  R_nan <- d$R; R_nan[1, 3] <- NaN
  ssc <- function(...) susieR:::summary_stats_constructor(...)
  vy <- var(d$y)
  for (R in list(d$R, R_dn, R_asym, R_nan)) {
    expect_same_modes(function() ssc(bhat = d$uni$betahat, shat = d$uni$sebetahat, R = R, n = d$n, var_y = vy))
    expect_same_modes(function() ssc(bhat = bh, shat = sh, R = R, n = d$n, var_y = vy, null_weight = 0.1))
  }
  expect_same_modes(function() ssc(bhat = bh, shat = sh, R = R_dn, n = d$n, var_y = vy,
                                   maf = seq(0.01, 0.5, length.out = p), maf_thresh = 0.3))
  expect_same_modes(function() ssc(bhat = bh[1], shat = sh[1], R = matrix(1), n = d$n, var_y = vy))
  # The kernel against the upstream expression on its own.
  s <- sqrt(seq_len(p) + 0.5)
  up <- t(R_dn * s) * s
  up <- (up + t(up)) / 2
  expect_identical(susieR:::fast_orig_scale_kernel(R_dn, s), up)
})

test_that("sufficient_stats_constructor (XtX) matches upstream", {
  d <- con_test_data()
  p <- ncol(d$XtX)
  suf <- function(...) susieR:::sufficient_stats_constructor(...)
  XtX_neg <- d$XtX; XtX_neg[4, 4] <- -XtX_neg[4, 4]
  XtX_names <- d$XtX; names(XtX_names) <- paste0("e", seq_along(XtX_names))
  XtX_dn <- d$XtX; dimnames(XtX_dn) <- list(NULL, paste0("v", 1:p))
  maf <- seq(0.01, 0.5, length.out = p)
  for (XtX in list(d$XtX, XtX_neg, XtX_names, XtX_dn)) {
    expect_same_modes(function() suf(XtX = XtX, Xty = d$Xty, yty = d$yty, n = d$n))
    expect_same_modes(function() suf(XtX = XtX, Xty = d$Xty, yty = d$yty, n = d$n, null_weight = 0.3))
  }
  expect_same_modes(function() suf(XtX = d$XtX, Xty = d$Xty, yty = d$yty, n = d$n, maf = maf, maf_thresh = 0.2))
  expect_same_modes(function() suf(XtX = d$XtX, Xty = d$Xty, yty = d$yty, n = d$n,
                                   maf = c(0.4, rep(0.01, p - 1)), maf_thresh = 0.2))
  expect_same_modes(function() suf(XtX = d$XtX, Xty = d$Xty, yty = d$yty, n = d$n, standardize = FALSE))
  expect_same_modes(function() suf(XtX = d$XtX, Xty = d$Xty, yty = d$yty, n = d$n, check_input = TRUE))
  expect_same_modes(function() suf(XtX = Matrix::Matrix(d$XtX, sparse = TRUE), Xty = d$Xty, yty = d$yty, n = d$n))
  # The user's matrix is not modified.
  XtX <- d$XtX + 0
  keep <- XtX + 0
  options(susieR.fast = "full"); on.exit(options(susieR.fast = NULL))
  res <- suf(XtX = XtX, Xty = d$Xty, yty = d$yty, n = d$n)
  expect_identical(XtX, keep)
})

test_that("individual_data_constructor constant-column warning matches upstream", {
  d <- con_test_data()
  X <- d$X[, 1:60]
  n <- nrow(X)
  X[, 5] <- 2; X[, 9] <- 0; X[, 11] <- c(1e-170, rep(0, n - 1))
  X[, 12] <- 1 / 3; X[, 13] <- -0; X[, 20] <- X[, 20] * 1e-120
  X[, 21] <- c(-1e150, 1e150, rep(0, n - 2)); X[, 22] <- 1e160
  X_inf <- X; X_inf[1, 30] <- Inf
  X_dn <- X; dimnames(X_dn) <- list(paste0("i", 1:n), paste0("v", 1:ncol(X)))
  X_na <- X; X_na[3, 3] <- NA
  idc <- function(...) susieR:::individual_data_constructor(...)
  r <- expect_same_modes(function() idc(X, d$y))
  expect_true(any(grepl("constant columns", r$full$msgs)))
  for (Xi in list(X_inf, X_dn, X_na, scale(X[, 30:60]), X[1:2, ], X[, 1, drop = FALSE]))
    expect_same_modes(function() idc(Xi, d$y[seq_len(nrow(Xi))]))
  expect_same_modes(function() idc(Matrix::Matrix(X, sparse = TRUE), d$y))
  expect_same_modes(function() idc(d$X, d$y))
})

test_that("fits through the constructors match upstream", {
  d <- con_test_data()
  R_dn <- d$R
  dimnames(R_dn) <- list(colnames(R_dn), colnames(R_dn))
  R_diag <- d$R; diag(R_diag) <- 1 + (seq_len(ncol(d$R)) %% 5) * 1e-3
  expect_same_modes(function() susie_rss(d$z, d$R, n = d$n))
  expect_same_modes(function() susie_rss(d$z, R_diag, n = d$n, null_weight = 0.1))
  expect_same_modes(function() susie_rss(bhat = d$uni$betahat, shat = d$uni$sebetahat,
                                         R = d$R, n = d$n, var_y = var(d$y)))
  expect_same_modes(function() susie_ss(d$XtX, d$Xty, d$yty, d$n))
  X <- d$X; X[, 3] <- 1
  expect_same_modes(function() susie(X, d$y))
})
