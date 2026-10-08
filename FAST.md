# susie-fast

**This repository is a fork of [stephenslab/susieR](https://github.com/stephenslab/susieR)**,
based on upstream commit
[`8e56a8e`](https://github.com/stephenslab/susieR/commit/8e56a8e038e989856d106d9ca5175cc664fea9d2)
(susieR 0.16.6); fork version 0.16.6.1. The SuSiE methods, the code and the credit belong to the
susieR authors (Wang, Zou, McCreight, Zhang, Denault, Carbonetto, Stephens).
Please cite their papers, not this fork.

The fork changes **how fast** susieR computes its fits, not **what** it
computes. It is a drop-in replacement: the package is still called `susieR`,
the functions and arguments are the same, and the fits are bit-identical to
upstream (see [Verification](#verification)).

## Install

```r
remotes::install_github("AlfredPgen/susie-fast")
```

This replaces an installed upstream `susieR`. To go back:
`install.packages("susieR")` or `remotes::install_github("stephenslab/susieR")`.

## What changed

The time in a SuSiE fit goes to matrix-vector products with the p x p
correlation matrix R (or the n x p genotype matrix X) and to the
one-dimensional optimisation of each effect's prior variance. Upstream repeats
some of that work. The fork computes each quantity once.

| | upstream, per effect update | fork |
|---|---|---|
| `susie_rss` / `susie_ss` | 2 products with R | 1 product, or none for a null effect |
| `susie` (individual data) | 3 products with X, 1 with X' | 1 with X, 1 with X'; none with X for a null effect |
| ELBO, per iteration | L + 1 products with R (or L with X) | none when the residual variance is fixed (default for `susie_rss`) |
| prior-variance optimiser | 20-40 R-level evaluations of the SER Bayes factor | same evaluations in one C++ pass each; repeats looked up |

In detail:

1. **Per-effect product cache.** To update effect l, IBSS removes the effect's
   old contribution R b_l from the fit and adds the new one. Upstream
   recomputes the old product, although it was computed when effect l was last
   updated. The fork keeps it and reuses it only when b_l is bit-for-bit the
   vector it was computed from. For individual data the KL term needs X b_l
   for the new b_l, which the fitted-value update also needs, so that product
   is now computed once as well.
2. **Null effects.** An effect with zero estimated prior variance has b_l = 0,
   so R b_l = 0. Upstream still multiplies (BLAS does not skip zeros). With the
   default L = 10 and a few real signals, most effects are null.
3. **Single-effect Bayes factor kernel.** Brent's method evaluates the SER log
   Bayes factor about 20-40 times per effect update. The terms that do not
   depend on the prior variance V are computed once per update and the rest in
   one C++ pass, with the same operations in the same order as the R code. The
   kernel is checked against the R code at first use in each session and is
   not used if a platform ever disagrees.
4. **Expected residual sum of squares.** The ELBO and the residual-variance
   update evaluate it on the same state in each iteration; it is now computed
   once. When it feeds only the ELBO, it is assembled from the cached
   per-effect products instead of L + 1 more products with R.
5. **Constructor.** XtX is standardised in one pass instead of three full-size
   copies (`t((1 / csd) * XtX) / csd`).

The new code is in [`R/susie_fast.R`](R/susie_fast.R) and
[`src/susie_fast.cpp`](src/susie_fast.cpp). The upstream files have one-line
changes at the call sites.

## Modes

```r
options(susieR.fast = "full")   # default
options(susieR.fast = "exact")
options(susieR.fast = "off")    # upstream code paths
```

- **full** (default): everything above. Fits are bit-identical to upstream.
  When the residual variance is fixed, the ELBO trace (`fit$elbo`) is summed
  in a different order and may differ in its last one or two digits
  (relative difference below 2e-15 in the tests).
- **exact**: everything except the reordered ELBO sum. Bit-identical to
  upstream, ELBO included.
- **off**: the upstream code.

The fast paths apply to the standard `ss` (`susie_rss`, `susie_ss`) and
`individual` (`susie`) data classes. Other classes (`susie_rss_lambda`,
multi-panel, and downstream packages' classes) use the upstream code.

### Post-processing ([`R/susie_fast_post.R`](R/susie_fast_post.R))

- **Purity early exit.** A diffuse effect (V > 0, flat alpha) gives a credible
  set of most variables, which the purity filter then drops. Its full purity
  costs O(m^2) (O(n m^2) from X). A few rows of the correlation matrix now
  show a pair below `min_abs_corr`, and the set is dropped without the rest.
  From R or XtX the probed values are the exact upstream elements; from X a
  rounding-error bound valid for any summation order certifies the decision.
  Applies when `median_abs_corr` is unset and `squared = FALSE`.
- **z-score memo.** With `compute_univariate_zscore = TRUE`, the univariate
  z-scores are computed once per top-level fit instead of once per refine
  candidate.
- **Parallel refinement (opt-in).** `options(susieR.refine_cores = k)` with
  k >= 2 fits the refine candidates of each step in k forked processes
  (Linux/macOS; the serial loop elsewhere, inside a forked worker, or with
  fewer than two credible sets). The fits, messages, warnings and RNG state
  are those of the serial loop: conditions are replayed in serial order,
  and if any candidate used the RNG (the X-path purity subsample) the step
  is redone serially. Use it with a single-threaded BLAS or OpenBLAS
  (pthreads): a BLAS or other code that started an OpenMP thread pool in
  the parent can hang the forked children. Apple's Accelerate, MKL and
  OpenMP builds are recognised by the BLAS library name and run serially;
  a BLAS whose name does not show its threading is not detected. MKL
  without CNR is not run-to-run reproducible even serially.

## Verification

Every claim was checked against upstream susieR built from the same commit
with the same toolchain, on two platforms: Windows 11 (R 4.3.2, Rtools43,
R's reference BLAS) and Linux (Ubuntu 26.04 under WSL2, conda-forge R 4.4.3,
GCC, OpenBLAS 0.3.34, one thread).

- **35 configurations**: `susie_rss` (defaults, residual variance estimated,
  L = 20 with prior weights, null weight, EM and simple prior methods, fixed
  prior, refine, PIP convergence, no n, bhat/shat input, greedy L, null
  threshold, track_fit, SuSiE-inf, SuSiE-ash, R mismatch, slot prior),
  `susie_ss` (MoM, MLE) and `susie` (defaults, unstandardised without
  intercept, EM, refine, fixed residual variance, L = 1, greedy L, NIG, inf,
  ash, slot prior, sparse X, the package's N3finemapping data). Results, on
  both platforms: `exact` bit-identical in all 35 (`identical()` on the whole
  fit; `off` too, checked on Windows); `full` bit-identical in all 35 except
  `fit$elbo` in 5 (Windows) or 7 (Linux) cases, largest relative difference
  1.8e-15. Tables in [`fast/results`](fast/results).
- **The package's own test suite** (1,240 tests): same result as upstream
  (8 failed expectations and 1 error in both builds on this machine, all
  pre-existing: 6 plotting tests that need ImageMagick and 2 R_mismatch tests
  that disagree with current upstream code).
- **New tests** in [`tests/testthat/test_susie_fast.R`](tests/testthat/test_susie_fast.R)
  compare the three modes on every run (pass on both platforms).

## Benchmarks

Simulated regions with realistic LD (n = 3,000 to 5,000 samples, p = 1,000 to
6,000 variants, 3 causal variants), L = 10 unless stated. Upstream and the fork
were run alternately, 5 rounds (3 for L = 20); speed-up is the median of the
paired per-round ratios. The machine had other load, so absolute times are
noisy; the pairing absorbs most of it.

Windows 11, R 4.3.2, reference BLAS, one thread (5 rounds):

| fit | upstream | fork | speed-up |
|---|---|---|---|
| `susie_rss`, p = 1,000 | 0.75 s | 0.40 s | 1.9x |
| `susie_rss`, p = 3,000 | 2.63 s | 1.02 s | 2.8x |
| `susie_rss`, p = 6,000 | 8.94 s | 2.98 s | 3.2x |
| `susie_rss`, p = 3,000, L = 20 | 8.31 s | 2.69 s | 2.8x |
| `susie_rss`, p = 3,000, residual variance estimated | 3.52 s | 1.50 s | 2.2x |
| `susie_rss`, p = 1,000, refine | 3.73 s | 1.69 s | 2.1x |
| `susie`, n = 3,000, p = 1,000 | 3.91 s | 1.96 s | 2.2x |
| `susie`, n = 3,000, p = 3,000 | 6.98 s | 3.92 s | 1.9x |
| `susie`, n = 3,000, p = 1,000, refine | 17.22 s | 8.22 s | 2.3x |

Linux (WSL2), R 4.4.3, OpenBLAS 0.3.34, one thread (3 rounds):

| fit | upstream | fork | speed-up |
|---|---|---|---|
| `susie_rss`, p = 1,000 | 0.73 s | 0.32 s | 2.2x |
| `susie_rss`, p = 3,000 | 4.15 s | 1.29 s | 2.7x |
| `susie_rss`, p = 6,000 | 9.23 s | 3.27 s | 2.4x |
| `susie_rss`, p = 3,000, L = 20 | 4.39 s | 1.65 s | 2.9x |
| `susie_rss`, p = 3,000, residual variance estimated | 2.71 s | 1.71 s | 1.6x |
| `susie`, n = 3,000, p = 1,000 | 2.22 s | 1.18 s | 1.9x |
| `susie`, n = 3,000, p = 3,000 | 5.29 s | 2.85 s | 1.9x |

The gain grows with p (the products dominate) and with the share of null
effects (L larger than the number of signals).

## Reproducing

[`fast/`](fast) holds the tools: `make_data.R` (simulated regions),
`run_suite.R` + `compare.R` (the 35-configuration comparison),
`bench_one.R` + `bench.sh` + `summarize_bench.R` (alternating benchmark) and
`run_tests.R` (the package test suite against an installed library). Install
upstream into one R library and this fork into another, then for example:

```sh
Rscript fast/make_data.R 3000 1000 1          # also 3000 3000 2, 5000 6000 3
Rscript fast/run_suite.R lib_base - out/base.rds
Rscript fast/run_suite.R lib_fast exact out/exact.rds
Rscript fast/compare.R out/base.rds out/exact.rds
LIB_BASE=lib_base LIB_FAST=lib_fast fast/bench.sh "rss_p1000 ind_p1000" 5
```

## Licence

Same as upstream: BSD 3-clause (see [LICENSE](LICENSE)).
