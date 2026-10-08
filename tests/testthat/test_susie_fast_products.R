context("susie-fast matrix-vector products")

# Fits with the products code compared with the upstream code path
# (options(susieR.fast = "off")) on the same input.
products_fit_modes <- function(fun, threads = 1) {
  old <- options(susieR.fast = "off", susieR.threads = threads)
  on.exit(options(old))
  set.seed(1); ref <- suppressWarnings(suppressMessages(fun()))
  options(susieR.fast = "exact")
  set.seed(1); exact <- suppressWarnings(suppressMessages(fun()))
  options(susieR.fast = "full")
  set.seed(1); full <- suppressWarnings(suppressMessages(fun()))
  list(ref = ref, exact = exact, full = full)
}

products_test_data <- function() {
  set.seed(43)
  n <- 300; p <- 90
  X <- matrix(as.double(rbinom(n * p, 2, 0.3)), n, p)
  for (j in 2:p) X[, j] <- ifelse(runif(n) < 0.7, X[, j - 1], X[, j])
  X <- X[, apply(X, 2, sd) > 0]
  b <- rep(0, ncol(X)); b[c(10, 50)] <- c(0.5, -0.4)
  y <- as.vector(X %*% b + rnorm(n))
  R <- cor(X)
  r <- as.vector(cor(X, y))
  z <- r * sqrt((n - 2) / (1 - r^2))
  list(X = X, y = y, R = R, z = z, n = n)
}

test_that("the products self-tests pass", {
  expect_true(susieR:::fast_products_self_test(FALSE, 1L))
  expect_true(susieR:::fast_products_ok())
  blas <- extSoftVersion()[["BLAS"]]
  if (blas == "" || grepl("^(lib)?Rblas", basename(blas))) {
    expect_true(susieR:::fast_products_self_test(TRUE, 1L))
    expect_true(susieR:::fast_products_self_test(TRUE, 3L))
  }
})

test_that("the NaN/Inf scan reproduces R's pairwise test", {
  scan <- susieR:::fast_mayhave_cpp
  expect_false(scan(c(1, 2, 3)))
  expect_true(scan(c(Inf, 1, 2)))         # odd length: first element alone
  expect_true(scan(c(1, NA)))
  expect_true(scan(c(Inf, -Inf)))
  expect_true(scan(c(1e308, 1e308)))      # finite values whose sum overflows
  expect_false(scan(c(1e308, -1e308)))
  expect_false(scan(c(1, 1e308, 1e308 / 2)))  # pairs are (1e308, 5e307): finite
  expect_false(scan(c(-0, 0)))
  x <- rnorm(20001)
  expect_false(scan(x))
  expect_true(scan(replace(x, 20001, NaN)))
  expect_true(scan(replace(x, 12345, -Inf)))
})

test_that("direct products equal R's products and refuse what R would not send to BLAS", {
  gemv <- susieR:::fast_gemv_cpp
  set.seed(3)
  for (ref in c(FALSE, TRUE)) {
    if (ref && !susieR:::fast_refblas_ok(1L)) next
    for (nt in if (ref) c(1L, 2L, 3L) else 1L) {
      for (s in list(c(2, 2), c(5, 3), c(3, 9), c(130, 77), c(1100, 13))) {
        X <- matrix(rnorm(prod(s)), s[1], s[2])
        v <- rnorm(s[2]); v[1] <- -0
        y <- rnorm(s[1])
        expect_identical(gemv(X, v, FALSE, ref, nt), as.vector(X %*% v))
        expect_identical(gemv(X, y, TRUE, ref, nt), as.vector(crossprod(y, X)))
      }
    }
  }
  X <- matrix(rnorm(12), 4, 3)
  expect_null(gemv(X, c(1, NaN, 2), FALSE, FALSE, 1L))
  expect_null(gemv(X, c(2, 1e308, 1e308), FALSE, FALSE, 1L))  # pair sum overflows
  expect_null(gemv(X, c(1, 2), FALSE, FALSE, 1L))
  expect_null(gemv(X, matrix(1, 3, 1), FALSE, FALSE, 1L))
  expect_null(gemv(X, c(1L, 2L, 3L), FALSE, FALSE, 1L))
  expect_identical(susieR:::fast_symmetric_cpp(X[1:3, ] + t(X[1:3, ])), TRUE)
  expect_identical(susieR:::fast_symmetric_cpp(X[1:3, ]), FALSE)
})

test_that("a product that disagrees with R switches the fit back to R", {
  st <- new.env()
  st$ok <- TRUE; st$checked <- character(0)
  out <- susieR:::fast_checked(st, "Rv", 1, 2, function() 3)
  expect_identical(out, 3)
  expect_false(st$ok)
  st$ok <- TRUE
  expect_identical(susieR:::fast_checked(st, "Rv", 0, 3, function() 3), 3)
  expect_identical(st$checked, character(0))      # a zero input does not count
  expect_identical(susieR:::fast_checked(st, "Rv", 1, 3, function() 3), 3)
  expect_identical(st$checked, "Rv")
  expect_identical(susieR:::fast_checked(st, "Rv", 1, 4, function() stop("not called")), 4)
})

test_that("the product state applies only to plain finite matrices", {
  d <- products_test_data()
  old <- options(susieR.fast = "exact")
  on.exit(options(old))
  ss <- list(XtX = d$R, p = ncol(d$R))
  class(ss) <- "ss"
  expect_true(is.environment(susieR:::fast_products_new(ss)))
  ss$XtX[2, 2] <- Inf
  expect_null(susieR:::fast_products_new(ss))
  ss$XtX <- Matrix::Matrix(d$R, sparse = TRUE)
  expect_null(susieR:::fast_products_new(ss))
  ss$XtX <- d$R
  class(ss) <- c("ss_other", "ss")
  expect_null(susieR:::fast_products_new(ss))
  class(ss) <- "ss"
  options(susieR.fast = "off")
  expect_null(susieR:::fast_products_new(ss))
})

test_that("summary-statistics fits are bit-identical to the upstream code path", {
  d <- products_test_data()
  Xu <- d$X[, 1:60]
  XtX <- crossprod(Xu); Xty <- as.vector(crossprod(Xu, d$y))
  init <- susie_rss(d$z, d$R, n = d$n, L = 3)
  cases <- list(
    function() susie_rss(d$z, d$R, n = d$n),
    function() susie_rss(d$z, d$R, n = d$n, estimate_residual_variance = TRUE),
    function() susie_rss(d$z, d$R, n = d$n, model_init = init),
    function() susie_rss(d$z, d$R, n = d$n, null_weight = 0.1),
    function() susie_rss(d$z, X = d$X[1:80, ], n = d$n),
    function() susie_ss(XtX, Xty, sum(d$y^2), d$n),
    function() susie_rss(d$z, d$R, n = d$n, unmappable_effects = "inf"),
    function() { old <- options(matprod = "internal"); on.exit(options(old))
                 susie_rss(d$z, d$R, n = d$n) })
  for (fun in cases) for (nt in c(1, 2)) {
    m <- products_fit_modes(fun, threads = nt)
    expect_identical(m$exact, m$ref)
    m$full$elbo <- m$ref$elbo <- NULL
    expect_identical(m$full, m$ref)
  }
})

test_that("individual-data fits are bit-identical to the upstream code path", {
  d <- products_test_data()
  init <- susie(d$X, d$y, L = 3)
  Xn <- d$X; colnames(Xn) <- paste0("v", seq_len(ncol(Xn)))
  cases <- list(
    function() susie(d$X, d$y),
    function() susie(Xn, d$y, standardize = FALSE),
    function() susie(d$X, d$y, model_init = init),
    function() susie(d$X, d$y, null_weight = 0.2, refine = TRUE),
    function() susie(d$X[, 1:2], d$y, L = 1),
    function() susie(d$X, d$y, unmappable_effects = "inf"),
    function() susie(Matrix::Matrix(d$X, sparse = TRUE), d$y))
  for (fun in cases) for (nt in c(1, 3)) {
    m <- products_fit_modes(fun, threads = nt)
    expect_identical(m$exact, m$ref)
    m$full$elbo <- m$ref$elbo <- NULL
    expect_identical(m$full, m$ref)
  }
})
