# Usage: Rscript compare.R base.rds other.rds
args <- commandArgs(TRUE)
A <- readRDS(args[1]); B <- readRDS(args[2])
maxdiff <- function(a, b) {
  if (is.list(a) && !is.data.frame(a)) {
    if (!identical(names(a), names(b))) return(Inf)
    if (length(a) == 0) return(0)
    return(max(vapply(seq_along(a), function(i) maxdiff(a[[i]], b[[i]]), 0)))
  }
  if (is.numeric(a) && is.numeric(b) && length(a) == length(b)) {
    if (length(a) == 0) return(0)
    na <- is.na(a) | is.na(b)
    if (any(is.na(a) != is.na(b))) return(Inf)
    a <- a[!na]; b <- b[!na]; if (!length(a)) return(0)
    fin <- is.finite(a) & is.finite(b)
    if (any(a[!fin] != b[!fin])) return(Inf)
    return(max(c(0, abs(a[fin] - b[fin]) / pmax(1, abs(a[fin])))))
  }
  if (is.environment(a) && is.environment(b)) return(maxdiff(as.list(a, sorted = TRUE), as.list(b, sorted = TRUE)))
  if (identical(a, b)) 0 else Inf
}
cat(sprintf("%-16s %-9s %-10s %-10s %-10s %s\n", "case", "identical", "max_rel", "pip", "elbo", "sets/niter"))
for (nm in names(A)) {
  a <- A[[nm]]$fit; b <- B[[nm]]$fit
  if (is.null(b)) { cat(nm, "missing\n"); next }
  if (inherits(a, "err") || inherits(b, "err")) {
    cat(sprintf("%-16s %s\n", nm, if (identical(a, b)) "both error (same)" else "ERROR MISMATCH")); next }
  strip_env <- function(x) { x[vapply(x, is.environment, TRUE)] <- NULL; x }
  id <- identical(strip_env(unclass(a)), strip_env(unclass(b)))
  md <- maxdiff(unclass(a), unclass(b))
  md_pip <- maxdiff(a$pip, b$pip); md_elbo <- maxdiff(a$elbo, b$elbo)
  same_sets <- identical(a$sets$cs, b$sets$cs) && identical(a$niter, b$niter)
  cat(sprintf("%-16s %-9s %-10.2e %-10.2e %-10.2e %s\n", nm, id, md, md_pip, md_elbo, same_sets))
}
