# Usage: Rscript run_suite.R <lib> <fast-mode|-> <out.rds> [filter-regex]
# Fits the configurations in fast/cases/*.R on data/ made by make_data.R
# (seeds 1 and 2) and saves the fits; compare two runs with compare.R.
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

# Case files: every fast/cases/*.R appends named functions to `cases`.
# Set SUSIE_SUITE_CASES to a comma-separated list of file stems to run a
# subset (e.g. "core,ser"); default: all files.
cases <- list()
case_dir <- file.path(dirname(normalizePath(sub("--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "cases")
case_files <- sort(list.files(case_dir, pattern = "[.]R$", full.names = TRUE))
sel <- Sys.getenv("SUSIE_SUITE_CASES", "")
if (nzchar(sel)) case_files <- case_files[sub("[.]R$", "", basename(case_files)) %in% strsplit(sel, ",")[[1]]]
for (cf in case_files) source(cf, local = TRUE)
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
