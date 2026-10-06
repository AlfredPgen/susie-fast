context("susie-fast code paths")

# Every fast-mode fit is compared with the upstream code path
# (options(susieR.fast = "off")) on the same input.
fit_modes <- function(fun) {
  old <- options(susieR.fast = "off")
  on.exit(options(old))
  set.seed(1); ref <- suppressWarnings(suppressMessages(fun()))
  options(susieR.fast = "exact")
  set.seed(1); exact <- suppressWarnings(suppressMessages(fun()))
  options(susieR.fast = "full")
  set.seed(1); full <- suppressWarnings(suppressMessages(fun()))
  list(ref = ref, exact = exact, full = full)
}

without_elbo <- function(fit) { fit$elbo <- NULL; fit }

fast_test_data <- function() {
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

test_that("the C++ kernels reproduce the R arithmetic bit for bit", {
  expect_true(susieR:::fast_self_test())
})

test_that("exact mode is bit-identical to the upstream code path", {
  d <- fast_test_data()
  cases <- list(
    function() susie_rss(d$z, d$R, n = d$n),
    function() susie_rss(d$z, d$R, n = d$n, estimate_residual_variance = TRUE),
    function() susie_rss(d$z, d$R, n = d$n, L = 5, refine = TRUE),
    function() susie_rss(d$z, d$R, n = d$n, estimate_prior_method = "EM"),
    function() susie(d$X, d$y),
    function() susie(d$X, d$y, standardize = FALSE, intercept = FALSE),
    function() susie(d$X, d$y, L = 6, L_greedy = 2))
  for (fun in cases) {
    m <- fit_modes(fun)
    expect_identical(m$exact, m$ref)
  }
})

test_that("full mode changes at most the last digits of the ELBO", {
  d <- fast_test_data()
  cases <- list(
    function() susie_rss(d$z, d$R, n = d$n),
    function() susie_rss(d$z, d$R, n = d$n, estimate_residual_variance = TRUE),
    function() susie(d$X, d$y),
    function() susie(d$X, d$y, estimate_residual_variance = FALSE))
  for (fun in cases) {
    m <- fit_modes(fun)
    expect_identical(without_elbo(m$full), without_elbo(m$ref))
    expect_equal(m$full$elbo, m$ref$elbo, tolerance = 1e-12)
  }
})

test_that("the per-effect product cache is used and stays exact", {
  d <- fast_test_data()
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
  fit <- suppressWarnings(susie_rss(d$z, d$R, n = d$n, L = 10))
  # Upstream: two products per effect update plus L + 1 per ELBO. With the
  # cache: at most one per non-null effect update plus L + 1 per ELBO.
  expect_lt(counter$calls, 2 * 10 * fit$niter)
})
