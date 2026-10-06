b <- read.csv("out/bench.csv")
w <- reshape(b[, c("round", "case", "build", "seconds")], idvar = c("round", "case"), timevar = "build", direction = "wide")
w$ratio <- w$seconds.base / w$seconds.fast
s <- do.call(rbind, lapply(split(w, w$case), function(x) data.frame(case = x$case[1],
  base_median = median(x$seconds.base), fast_median = median(x$seconds.fast),
  speedup_median = round(median(x$ratio), 2), speedup_min = round(min(x$ratio), 2), speedup_max = round(max(x$ratio), 2))))
s <- s[order(match(s$case, unique(b$case))), ]
niter <- tapply(b$niter, b$case, function(x) paste(unique(x), collapse = "/"))
s$niter <- niter[s$case]
print(s, row.names = FALSE)
write.csv(s, "out/bench_summary.csv", row.names = FALSE)
