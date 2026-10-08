context("susie-fast code paths: susie_rss_lambda")

# Every fast-mode fit is compared with the upstream code path
# (options(susieR.fast = "off")) on the same input.
lambda_fit_modes <- function(fun) {
  old <- options(susieR.fast = "off")
  on.exit(options(old))
  set.seed(1); ref <- suppressWarnings(suppressMessages(fun()))
  options(susieR.fast = "exact")
  set.seed(1); exact <- suppressWarnings(suppressMessages(fun()))
  options(susieR.fast = "full")
  set.seed(1); full <- suppressWarnings(suppressMessages(fun()))
  list(ref = ref, exact = exact, full = full)
}

lambda_test_data <- function() {
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
  # rank-deficient R from 80 samples
  X80 <- X[1:80, ]
  X80 <- X80[, apply(X80, 2, sd) > 0]
  r80 <- as.vector(cor(X80, y[1:80]))
  list(X = X, y = y, R = R, z = z, n = n,
       R80 = cor(X80), z80 = r80 * sqrt(78 / (1 - r80^2)))
}

test_that("the diag(R S^-1 R) kernel reproduces rowSums bit for bit", {
  expect_true(susieR:::fast_lambda_self_test())
})

test_that("susie_rss_lambda fits are bit-identical to the upstream code path", {
  d <- lambda_test_data()
  pw <- { set.seed(3); w <- runif(ncol(d$R)); w / sum(w) }
  cases <- list(
    function() susie_rss_lambda(d$z, d$R, n = d$n, lambda = 1e-3),
    function() susie_rss_lambda(d$z, d$R, n = d$n, lambda = 0.01,
                                estimate_residual_variance = TRUE),
    function() susie_rss_lambda(d$z, d$R, n = d$n, lambda = 1e-3, L = 5, refine = TRUE),
    function() susie_rss_lambda(d$z, d$R, n = d$n, lambda = 1e-3, L = 6, L_greedy = 2,
                                estimate_residual_variance = TRUE),
    function() susie_rss_lambda(d$z, d$R, n = d$n, lambda = 1e-3, null_weight = 0.1,
                                prior_weights = pw, estimate_residual_variance = TRUE),
    function() susie_rss_lambda(d$z, X = d$X, n = d$n, lambda = 1e-3),
    function() susie_rss_lambda(d$z, d$R, n = d$n, lambda = 1e-3,
                                slot_prior = slot_prior_betabinom()),
    function() susie_rss_lambda(d$z80, d$R80, n = 80, lambda = 0,
                                estimate_residual_variance = TRUE),
    function() susie_rss_lambda(d$z80, d$R80, n = 80, lambda = "estimate", refine = TRUE),
    function() susie_rss_lambda(d$z, d$R, n = d$n, lambda = 1e-3, L = 12,
                                estimate_prior_method = "simple"),
    function() susie_rss_lambda(d$z, d$R, n = d$n, lambda = 1e-3,
                                estimate_prior_method = "EM"))
  for (fun in cases) {
    m <- lambda_fit_modes(fun)
    expect_identical(m$exact, m$ref)
    expect_identical(m$full, m$ref)
    expect_true(identical(m$exact, m$ref, num.eq = FALSE))
  }
})

test_that("SinvRj is computed once per susie_rss_lambda call", {
  data <- structure(list(), class = "rss_lambda")
  V <- diag(2); D <- c(1, 0.5); Dinv <- 1 / (D + 0.1)
  old <- options(susieR.fast = "exact")
  on.exit(options(old))
  n_eval <- 0
  compute <- function() { n_eval <<- n_eval + 1; list(SinvRj = n_eval) }
  susieR:::fast_lambda_scope({
    a <- susieR:::fast_lambda_sinv(data, V, D, Dinv, compute())
    b <- susieR:::fast_lambda_sinv(data, V, D, Dinv, compute())
    x <- susieR:::fast_lambda_sinv(data, V, D, Dinv * 2, compute())
    e <- susieR:::fast_lambda_sinv(data, V, D, Dinv, compute())
  })
  expect_equal(n_eval, 2)
  expect_identical(b, a)
  expect_identical(e, a)
  # No memo outside a susie_rss_lambda() call, in "off" mode, or for
  # subclasses.
  susieR:::fast_lambda_sinv(data, V, D, Dinv, compute())
  expect_equal(n_eval, 3)
  options(susieR.fast = "off")
  susieR:::fast_lambda_scope({
    susieR:::fast_lambda_sinv(data, V, D, Dinv, compute())
    susieR:::fast_lambda_sinv(data, V, D, Dinv, compute())
  })
  expect_equal(n_eval, 5)
  options(susieR.fast = "exact")
  sub <- structure(list(), class = c("my_lambda", "rss_lambda"))
  susieR:::fast_lambda_scope({
    susieR:::fast_lambda_sinv(sub, V, D, Dinv, compute())
    susieR:::fast_lambda_sinv(sub, V, D, Dinv, compute())
  })
  expect_equal(n_eval, 7)
})

test_that("the product cache and null-effect shortcut are used for rss_lambda", {
  d <- lambda_test_data()
  old <- options(susieR.fast = "exact")
  on.exit(options(old))
  counter <- new.env()
  counter$calls <- 0
  suppressMessages(trace("compute_Rv",
                         bquote(assign("calls", get("calls", .(counter)) + 1,
                                       envir = .(counter))),
                         print = FALSE, where = asNamespace("susieR")))
  on.exit(suppressMessages(untrace("compute_Rv", where = asNamespace("susieR"))),
          add = TRUE)
  fit <- suppressWarnings(susie_rss_lambda(d$z, d$R, n = d$n, lambda = 1e-3, L = 10))
  # Upstream: three products with R per effect update. Here: at most one
  # per non-null effect update, none for a null effect.
  expect_lt(counter$calls, 10 * fit$niter)
})
