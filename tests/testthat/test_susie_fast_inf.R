context("susie-fast: SuSiE-inf code paths")

# Each fit is run with the upstream code path (options(susieR.fast = "off"))
# and in the two fast modes; fits (or errors), warnings and messages must be
# identical.
inf_fit_modes <- function(fun) {
  old <- options(susieR.fast = "off")
  on.exit(options(old))
  run <- function(mode) {
    options(susieR.fast = mode)
    conds <- character(0)
    set.seed(1)
    fit <- tryCatch(withCallingHandlers(fun(),
      warning = function(w) {
        conds <<- c(conds, paste("W:", conditionMessage(w)))
        invokeRestart("muffleWarning")
      },
      message = function(m) {
        conds <<- c(conds, paste("M:", conditionMessage(m)))
        invokeRestart("muffleMessage")
      }),
      error = function(e) paste("E:", conditionMessage(e)))
    list(fit = fit, conds = conds)
  }
  list(ref = run("off"), exact = run("exact"), full = run("full"))
}

inf_test_data <- function() {
  set.seed(42)
  n <- 400; p <- 120
  X <- matrix(as.double(rbinom(n * p, 2, 0.3)), n, p)
  for (j in 2:p) X[, j] <- ifelse(runif(n) < 0.7, X[, j - 1], X[, j])
  X <- X[, apply(X, 2, sd) > 0]
  b <- rep(0, ncol(X)); b[c(10, 60)] <- c(0.5, -0.4)
  y <- as.vector(X %*% b + 0.02 * rowSums(X) + rnorm(n))
  R <- cor(X)
  r <- as.vector(cor(X, y))
  z <- r * sqrt((n - 2) / (1 - r^2))
  Xs <- scale(X)
  yc <- y - mean(y)
  list(X = X, y = y, R = R, z = z, n = n,
       XtX = crossprod(Xs), Xty = as.vector(crossprod(Xs, yc)),
       yty = sum(yc^2))
}

expect_modes_identical <- function(fun) {
  m <- inf_fit_modes(fun)
  expect_true(identical(m$exact$fit, m$ref$fit, num.eq = FALSE))
  expect_true(identical(m$full$fit, m$ref$fit, num.eq = FALSE))
  expect_identical(m$exact$conds, m$ref$conds)
  expect_identical(m$full$conds, m$ref$conds)
}

test_that("SuSiE-inf fits on ss data are bit-identical in every mode", {
  d <- inf_test_data()
  cases <- list(
    function() susie_rss(d$z, d$R, n = d$n, unmappable_effects = "inf"),
    function() susie_rss(d$z, d$R, n = d$n, L = 1, unmappable_effects = "inf"),
    function() susie_rss(d$z, d$R, n = d$n, max_iter = 2, unmappable_effects = "inf"),
    function() susie_rss(d$z, d$R, n = d$n, track_fit = TRUE, verbose = TRUE,
                         unmappable_effects = "inf"),
    function() susie_rss(d$z, d$R, n = d$n, L = 9, L_greedy = 3,
                         unmappable_effects = "inf"),
    function() susie_rss(d$z, d$R, n = d$n, estimate_prior_variance = FALSE,
                         unmappable_effects = "inf"),
    function() susie_rss(d$z, d$R, n = d$n, R_finite = 300,
                         unmappable_effects = "inf"),
    function() susie_rss(d$z, d$R, n = d$n, R_finite = 300, R_mismatch = "eb",
                         unmappable_effects = "inf"),
    function() susie_rss(d$z, d$R, n = d$n, slot_prior = slot_prior_betabinom(),
                         unmappable_effects = "inf"),
    function() susie_ss(d$XtX, d$Xty, d$yty, d$n, unmappable_effects = "inf"),
    # NIG with inf on ss data stops inside the SER (upstream behaviour).
    function() susie_ss(d$XtX, d$Xty, d$yty, d$n,
                        estimate_residual_method = "NIG",
                        unmappable_effects = "inf"))
  for (fun in cases) expect_modes_identical(fun)
})

test_that("SuSiE-inf fits on individual data are bit-identical in every mode", {
  d <- inf_test_data()
  cases <- list(
    function() susie(d$X, d$y, unmappable_effects = "inf"),
    function() susie(d$X, d$y, standardize = FALSE, intercept = FALSE,
                     unmappable_effects = "inf"),
    function() susie(d$X, d$y, track_fit = TRUE, verbose = TRUE,
                     unmappable_effects = "inf"),
    function() susie(d$X, d$y, L = 8, L_greedy = 4, unmappable_effects = "inf"),
    function() susie(d$X, d$y, slot_prior = slot_prior_poisson(C = 3),
                     unmappable_effects = "inf"),
    function() susie(d$X, d$y, max_iter = 2, unmappable_effects = "inf"))
  for (fun in cases) expect_modes_identical(fun)
})

test_that("the transpose-free diagonal term equals rowSums(sweep(t(V^2)))", {
  set.seed(5)
  for (p in c(1, 2, 17, 300)) {
    r <- max(1, p - 3)
    V2 <- matrix(rnorm(p * r)^2 * 10^runif(p * r, -300, 300), p, r)
    tmpD <- rnorm(p) * 10^runif(p, -20, 20)
    if (p > 2) {
      tmpD[1] <- -0; V2[2, 1] <- NaN; V2[3, r] <- Inf; tmpD[p] <- NA
    }
    colnames(V2) <- paste0("k", seq_len(r))
    ref <- rowSums(sweep(t(V2), 2, tmpD, `*`))
    expect_true(identical(susieR:::fast_inf_diag_term(V2, tmpD), ref,
                          num.eq = FALSE))
  }
  expect_null(susieR:::fast_inf_diag_term(matrix(1, 3, 2), c(1, 2)))
})

test_that("a null effect's eigenspace product is skipped only when V is finite", {
  cache <- new.env(parent = emptyenv())
  V <- matrix(c(1, 2, 3, 4), 2, 2)
  expect_true(susieR:::fast_inf_null(cache, V, c(0, -0)))
  expect_false(susieR:::fast_inf_null(cache, V, c(0, 1e-300)))
  expect_false(susieR:::fast_inf_null(cache, V, c(0, NaN)))
  expect_true(all(t(V) %*% c(0, -0) == 0))
  cache <- new.env(parent = emptyenv())
  expect_false(susieR:::fast_inf_null(cache, matrix(c(1, NaN, 3, 4), 2, 2),
                                      c(0, 0)))
})

test_that("SuSiE-inf helpers fall back to upstream code outside the workhorse", {
  data  <- structure(list(), class = c("ss", "list"))
  model <- list(alpha = matrix(1 / 3, 3, 3))
  expect_null(susieR:::fast_inf_cache(model))
  expect_false(susieR:::fast_inf_skip_fitted(data, list(), model, 1L))
  expect_null(susieR:::fast_residuals_inf(data, model, rep(0, 3)))
  # A runtime cache is ignored in "off" mode.
  model$runtime <- list(fast_cache = new.env())
  old <- options(susieR.fast = "off")
  on.exit(options(old))
  expect_null(susieR:::fast_inf_cache(model))
})
