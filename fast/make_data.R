# Simulated fine-mapping regions with realistic LD (AR(1) latent haplotypes -> 0/1/2 genotypes)
# Usage: Rscript make_data.R <n> <p> <seed>  -> data/region_n<n>_p<p>_s<seed>.rds
# The suite and benchmarks use: 3000 1000 1, 3000 3000 2, 5000 6000 3.
args <- commandArgs(TRUE)
n <- as.integer(args[1]); p <- as.integer(args[2]); seed <- as.integer(args[3])
set.seed(seed)
hap <- function(n, p, rho) {
  H <- matrix(0, n, p); H[, 1] <- rnorm(n)
  e <- sqrt(1 - rho^2)
  for (j in 2:p) H[, j] <- rho * H[, j - 1] + e * rnorm(n)
  thr <- qnorm(1 - runif(p, 0.05, 0.5))   # MAF 5-50%
  sweep(H, 2, thr, ">") * 1
}
rho <- 0.97
X <- hap(n, p, rho) + hap(n, p, rho)
storage.mode(X) <- "double"
X <- X[, apply(X, 2, sd) > 0]
p <- ncol(X)
b <- rep(0, p); causal <- sort(sample(p, 3)); b[causal] <- c(0.30, -0.22, 0.16)
y <- as.vector(X %*% b + rnorm(n))
Xs <- scale(X)
ys <- y - mean(y)
bhat <- as.vector(crossprod(Xs, ys)) / (n - 1)
# z from simple regression
r <- as.vector(cor(X, y))
z <- r * sqrt((n - 2) / (1 - r^2))
R <- cor(X)
saveRDS(list(X = X, y = y, z = z, R = R, n = n, causal = causal),
        sprintf("data/region_n%d_p%d_s%d.rds", n, p, seed), compress = FALSE)
cat("saved p =", p, " max|z| =", max(abs(z)), " causal:", causal, "\n")
