# Usage: Rscript bench_one.R <lib> <case> [fast-mode]
# Times one fit of <case> (e.g. rss_p3000, ind_p1000) on data/ made by make_data.R.
args <- commandArgs(TRUE); lib <- args[1]; case <- args[2]
if (length(args) > 2) options(susieR.fast = args[3])
suppressPackageStartupMessages(library(susieR, lib.loc = lib))
f <- switch(sub("_.*", "", sub("^[^_]+_", "", case)),
  p1000 = "data/region_n3000_p1000_s1.rds", p3000 = "data/region_n3000_p3000_s2.rds",
  p6000 = "data/region_n5000_p6000_s3.rds")
d <- readRDS(f)
run <- switch(sub("_p[0-9]+", "", case),
  rss        = function() susie_rss(d$z, d$R, n = d$n),
  rsserv     = function() susie_rss(d$z, d$R, n = d$n, estimate_residual_variance = TRUE),
  rssL20     = function() susie_rss(d$z, d$R, n = d$n, L = 20),
  rssrefine  = function() susie_rss(d$z, d$R, n = d$n, refine = TRUE),
  ind        = function() susie(d$X, d$y),
  indrefine  = function() susie(d$X, d$y, refine = TRUE))
gc()
t <- system.time(fit <- suppressWarnings(suppressMessages(run())))[["elapsed"]]
cat(sprintf("RESULT %s %s %.3f %d\n", case, lib, t, fit$niter))
