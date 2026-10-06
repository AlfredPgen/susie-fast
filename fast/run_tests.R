# Usage: Rscript run_tests.R <lib> [tests/testthat dir]
args <- commandArgs(TRUE); lib <- normalizePath(args[1])
.libPaths(c(lib, .libPaths()))
library(testthat)
res <- testthat::test_dir(if (length(args) > 1) args[2] else "tests/testthat", package = "susieR", load_package = "installed",
                          reporter = ProgressReporter$new(show_praise = FALSE), stop_on_failure = FALSE)
df <- as.data.frame(res)
cat("\nSUMMARY files", length(unique(df$file)), "tests", nrow(df), "failed", sum(df$failed), "errors", sum(df$error), "skipped", sum(df$skipped), "warnings", sum(df$warning), "\n")
bad <- df[df$failed > 0 | df$error, c("file", "test")]
if (nrow(bad)) print(bad)
