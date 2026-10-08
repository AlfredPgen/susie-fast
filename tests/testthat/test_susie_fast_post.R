context("susie-fast post-processing paths")

# Runs fun() in the given fast mode from a fixed RNG state and returns the
# value (or the error message), the messages and warnings in order, and
# the RNG state afterwards.
post_run <- function(fun, mode, seed = 1) {
  old <- options(susieR.fast = mode)
  on.exit(options(old))
  set.seed(seed)
  conds <- list()
  value <- tryCatch(
    withCallingHandlers(fun(),
      message = function(c) {
        conds[[length(conds) + 1]] <<- list("message", conditionMessage(c))
        invokeRestart("muffleMessage")
      },
      warning = function(c) {
        conds[[length(conds) + 1]] <<- list("warning", conditionMessage(c),
                                            conditionCall(c))
        invokeRestart("muffleWarning")
      }),
    error = function(e) list(error = conditionMessage(e),
                             call = deparse(conditionCall(e))))
  list(value = value, conds = conds, seed = .Random.seed)
}

expect_same_run <- function(fun, modes = c("exact", "full"), seed = 1) {
  ref <- post_run(fun, "off", seed)
  for (m in modes) {
    out <- post_run(fun, m, seed)
    if (m == "full" && is.list(out$value) && !is.null(out$value$elbo)) {
      out$value$elbo <- NULL
      ref$value$elbo <- NULL
    }
    expect_identical(out, ref)
  }
  invisible(ref)
}

post_test_data <- function(n = 400, p = 150, seed = 42, b = c(0.5, -0.4)) {
  set.seed(seed)
  X <- matrix(as.double(rbinom(n * p, 2, 0.3)), n, p)
  for (j in 2:p) X[, j] <- ifelse(runif(n) < 0.7, X[, j - 1], X[, j])
  X <- X[, apply(X, 2, sd) > 0]
  beta <- rep(0, ncol(X))
  beta[round(seq(10, ncol(X) - 10, length.out = length(b)))] <- b
  y <- as.vector(X %*% beta + rnorm(n))
  R <- cor(X)
  r <- as.vector(cor(X, y))
  z <- r * sqrt((n - 2) / (1 - r^2))
  list(X = X, y = y, R = R, z = z, n = n)
}

# A susie-like object whose effects have the given credible sets.
cs_res <- function(sets, p) {
  alpha <- t(vapply(sets, function(s) {
    a <- rep(1e-12, p); a[s] <- 1; a / sum(a)
  }, numeric(p)))
  list(alpha = alpha, V = rep(1, length(sets)))
}

# Number of non-NULL results of the internal function fname while expr runs.
count_hits <- function(fname, expr) {
  ns <- asNamespace("susieR")
  env <- new.env()
  env$n <- 0L
  f <- get(fname, ns)
  g <- function(...) {
    v <- f(...)
    if (!is.null(v)) env$n <- env$n + 1L
    v
  }
  unlockBinding(fname, ns)
  assign(fname, g, ns)
  on.exit({
    assign(fname, f, ns)
    lockBinding(fname, ns)
  })
  force(expr)
  env$n
}

test_that("the purity early exit leaves susie_get_cs and fits unchanged", {
  d <- post_test_data(b = c(0.12, 0.1))
  fit <- suppressMessages(susie_rss(d$z, d$R, n = d$n))
  res <- list(alpha = fit$alpha, V = fit$V)
  XtX <- crossprod(d$X)
  for (args in list(list(Xcorr = d$R), list(X = d$X),
                    list(Xcorr = structure(list(XtX = XtX), class = "scaled_XtX"),
                         check_symmetric = FALSE),
                    list(X = d$X, n_purity = 30),
                    list(X = Matrix::Matrix(d$X, sparse = TRUE)),
                    list(Xcorr = d$R, min_abs_corr = 0.2),
                    list(Xcorr = d$R, median_abs_corr = 0.5),
                    list(Xcorr = d$R, squared = TRUE))) {
    expect_same_run(function() do.call(susie_get_cs, c(list(res), args)))
  }
  # weak signals give diffuse sets, which are dropped early
  ff <- function() susie_rss(d$z, d$R, n = d$n)
  expect_same_run(ff)
  expect_gt(count_hits("fast_purity_drop", post_run(ff, "exact")), 0)
  expect_same_run(function() susie(d$X, d$y))
  expect_same_run(function() susie(d$X, d$y, n_purity = 20))
  expect_same_run(function() susie(d$X, d$y, null_weight = 0.2))
  expect_same_run(function() susie_ss(XtX, as.vector(crossprod(d$X, d$y - mean(d$y))),
                                      sum((d$y - mean(d$y))^2), d$n))
  expect_same_run(function() susie_rss(d$z, X = d$X[1:60, ], n = d$n))
})

test_that("the purity early exit preserves upstream errors", {
  p <- 6
  res <- cs_res(list(1:p), p)
  # zero diagonal next to an overflowing entry: NaN in safe_cov2cor
  V <- diag(p); V[1, 1] <- 0; V[2, 2] <- 1e-300; V[1, 2] <- V[2, 1] <- 1e300
  sx <- structure(list(XtX = V), class = "scaled_XtX")
  r <- expect_same_run(function() susie_get_cs(res, Xcorr = sx, check_symmetric = FALSE))
  expect_match(r$value$error, "NaN")
  # zero diagonal without overflow: dropped, as upstream
  V <- diag(p) * 2; V[1, 1] <- 0; V[4, 5] <- V[5, 4] <- 1.9
  sx <- structure(list(XtX = V), class = "scaled_XtX")
  get_sx <- function() susie_get_cs(res, Xcorr = sx, check_symmetric = FALSE)
  expect_same_run(get_sx)
  expect_gt(count_hits("fast_purity_drop", post_run(get_sx, "exact")), 0)
  # negative diagonal
  V <- diag(p); V[3, 3] <- -1
  expect_same_run(function() susie_get_cs(res, Xcorr = structure(list(XtX = V),
                                                                 class = "scaled_XtX"),
                                          check_symmetric = FALSE))
  # NaN in an entry the probe does not read
  V <- diag(p) * 0.01 + 0; V[2, 3] <- V[3, 2] <- NaN
  expect_same_run(function() susie_get_cs(res, Xcorr = V, check_symmetric = FALSE))
  V[2, 3] <- V[3, 2] <- Inf
  expect_same_run(function() susie_get_cs(res, Xcorr = V, check_symmetric = FALSE))
  # sparse Xcorr
  expect_same_run(function() susie_get_cs(res, Xcorr = Matrix::Matrix(diag(p) * 0.1 + 0,
                                                                      sparse = TRUE)))
  # NA in X, constant column in X
  set.seed(3)
  X <- matrix(rnorm(40 * p), 40, p)
  Xna <- X; Xna[5, 2] <- NA
  expect_same_run(function() susie_get_cs(res, X = Xna))
  Xc <- X; Xc[, 4] <- 2
  expect_same_run(function() susie_get_cs(res, X = Xc))
  # n_purity sampling with and without an RNG state
  expect_same_run(function() susie_get_cs(res, X = X, n_purity = 3))
  # without an RNG state the subsample is drawn by get_purity itself
  ns <- asNamespace("susieR")
  seed <- .Random.seed
  rm(".Random.seed", envir = globalenv())
  expect_null(ns$fast_purity_drop(1:p, X, NULL, FALSE, 3, 0.5, NULL))
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  assign(".Random.seed", seed, envir = globalenv())
  # null index
  res$null_index <- p
  expect_same_run(function() susie_get_cs(res, Xcorr = diag(p)))
})

test_that("a certified drop always has upstream min |corr| below the threshold", {
  ns <- asNamespace("susieR")
  old <- options(susieR.fast = "exact")
  on.exit(options(old))
  set.seed(9)
  for (trial in 1:60) {
    m <- sample(c(2, 3, 8, 40), 1)
    n <- sample(c(5, 30, 200), 1)
    X <- matrix(rnorm(n * m), n, m) %*% matrix(rnorm(m * m, sd = 0.3), m) +
         matrix(rnorm(n), n, m)
    if (trial %% 4 == 0) X <- X * 1e150
    thr <- runif(1, 0.05, 0.95)
    R <- cor(X)
    XtX <- crossprod(X) * runif(1, 1e-5, 1e5)
    for (args in list(list(X = X, Xcorr = NULL),
                      list(X = NULL, Xcorr = R),
                      list(X = NULL, Xcorr = structure(list(XtX = XtX),
                                                       class = "scaled_XtX")))) {
      v <- ns$fast_purity_drop(seq_len(m), args$X, args$Xcorr, FALSE, "auto", thr, NULL)
      up <- tryCatch(ns$get_purity(seq_len(m), args$X, args$Xcorr, FALSE, "auto"),
                     error = function(e) NULL)
      if (!is.null(v)) {
        expect_false(is.null(up))
        expect_lt(up[1], thr)
      }
    }
  }
})

test_that("z-scores are computed once per fit with identical results", {
  d <- post_test_data(b = c(0.6, -0.5, 0.4))
  ns <- asNamespace("susieR")
  cnt <- new.env()
  trace("calc_z", bquote(assign("n", get("n", envir = .(cnt)) + 1L, envir = .(cnt))),
        where = ns, print = FALSE)
  on.exit(untrace("calc_z", where = ns))
  ncalls <- function(fun, mode) {
    cnt$n <- 0L
    out <- post_run(fun, mode)
    list(out = out, n = cnt$n)
  }
  cases <- list(
    function() susie(d$X, d$y, refine = TRUE, compute_univariate_zscore = TRUE),
    function() susie(d$X, d$y, refine = TRUE, null_weight = 0.1,
                     compute_univariate_zscore = TRUE),
    function() susie(d$X, d$y, L = 6, L_greedy = 2, compute_univariate_zscore = TRUE),
    function() susie(d$X, d$y, standardize = FALSE, intercept = FALSE, refine = TRUE,
                     compute_univariate_zscore = TRUE),
    function() susie(Matrix::Matrix(d$X, sparse = TRUE), d$y, refine = TRUE,
                     compute_univariate_zscore = TRUE))
  n <- NULL
  for (fun in cases) {
    ref <- ncalls(fun, "off")
    out <- ncalls(fun, "exact")
    expect_identical(out$out, ref$out)
    expect_lte(out$n, ref$n)
    n <- rbind(n, c(ref$n, out$n))
  }
  # refine: once instead of once per fit
  expect_gt(n[1, 1], 1L)
  expect_equal(n[1, 2], 1L)
  # greedy rounds share the z-scores too
  expect_equal(n[3, 2], 1L)
})

test_that("the z-score memo checks X and y on every lookup", {
  ns <- asNamespace("susieR")
  old <- options(susieR.fast = "exact")
  on.exit(options(old))
  set.seed(5)
  X <- matrix(rnorm(200), 40, 5)
  data <- structure(list(X = X, y = rnorm(40)), class = "individual")
  params <- list(compute_univariate_zscore = TRUE, intercept = TRUE, standardize = TRUE)
  data <- ns$fast_zscore_memo_attach(data, params)
  expect_true(is.environment(data$fast_zscore_memo))
  model <- list()
  z1 <- ns$get_zscore.individual(data, params, model)
  expect_identical(ns$get_zscore.individual(data, params, model), z1)
  expect_identical(data$fast_zscore_memo$hits, 1L)
  d2 <- data
  d2$y <- rnorm(40)
  expect_identical(ns$get_zscore.individual(d2, params, model), calc_z(X, d2$y, center = TRUE, scale = TRUE))
  d3 <- data
  d3$X[1, 1] <- 0
  expect_identical(ns$get_zscore.individual(d3, params, model), calc_z(d3$X, d3$y, center = TRUE, scale = TRUE))
  p2 <- params; p2$standardize <- FALSE
  expect_identical(ns$get_zscore.individual(data, p2, model), calc_z(X, data$y, center = TRUE, scale = FALSE))
  options(susieR.fast = "off")
  fresh <- structure(list(X = X, y = data$y), class = "individual")
  expect_null(ns$fast_zscore_memo_attach(fresh, params)$fast_zscore_memo)
  options(susieR.fast = "exact")
  sub <- structure(list(X = X, y = data$y), class = c("mysub", "individual"))
  expect_null(ns$fast_zscore_memo_attach(sub, params)$fast_zscore_memo)
})

test_that("parallel refine is opt-in and gives the serial fits", {
  d <- post_test_data(b = c(0.6, -0.5, 0.4))
  ns <- asNamespace("susieR")
  cases <- list(
    function() susie_rss(d$z, d$R, n = d$n, refine = TRUE),
    function() susie(d$X, d$y, refine = TRUE),
    function() susie(d$X, d$y, refine = TRUE, compute_univariate_zscore = TRUE),
    function() susie(d$X, d$y, refine = TRUE, n_purity = 5),
    function() susie(d$X, d$y, refine = TRUE, max_iter = 3),
    function() susie(d$X, d$y, refine = TRUE, null_weight = 0.1))
  old <- options(susieR.refine_cores = 2)
  on.exit(options(old))
  if (.Platform$OS.type != "unix") {
    expect_null(ns$fast_refine_parallel(list(sets = list(cs = list(1, 2))),
                                        structure(list(), class = "ss"), list(), c(0.5, 0.5)))
  }
  for (fun in cases) {
    ref <- post_run(fun, "off")
    for (m in c("exact", "full")) {
      out <- post_run(fun, m)
      if (m == "full") { out$value$elbo <- NULL; ref$value$elbo <- NULL }
      expect_identical(out, ref)
    }
  }
  # serial fast run equals the parallel fast run
  for (fun in cases) {
    options(susieR.refine_cores = 1)
    a <- post_run(fun, "exact")
    options(susieR.refine_cores = 2)
    b <- post_run(fun, "exact")
    expect_identical(a, b)
  }
  # on Linux/macOS the candidates really ran in forked processes
  if (.Platform$OS.type == "unix")
    expect_gt(count_hits("fast_refine_parallel", post_run(cases[[1]], "exact")), 0)
})
