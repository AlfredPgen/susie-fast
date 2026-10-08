# susie-fast

**This repository is a fork of [stephenslab/susieR](https://github.com/stephenslab/susieR)**,
based on upstream commit
[`8e56a8e`](https://github.com/stephenslab/susieR/commit/8e56a8e038e989856d106d9ca5175cc664fea9d2)
(susieR 0.16.6); fork version 0.16.6.2. The SuSiE methods, the code and the credit belong to the
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

Round 2 (0.16.6.2) added, after an audit in which every proposed change was
checked by two independent reviewers (one for exactness, one for whether the
gain is real):

6. **Repeated single-effect regressions.** Consecutive null effects see the
   same residuals, so their SER (and, for individual data, the X'r product)
   is the same. A whole SER is reused when every input it reads is
   bit-for-bit the same as the previous one's.
7. **Fewer logarithms in the SER kernel.** In the standardised `susie_rss`
   path the prior-variance terms take very few distinct values across
   variants; `log(1 + V/s)` is computed once per distinct value. The final
   SER step for summary data (Bayes factors, posterior moments, KL) runs in
   C++ with the same operation order.
8. **Products without R's NaN scan.** Before each `%*%`, R scans the whole
   matrix for NaN/Inf, which costs about as much as the product. The matrix
   is now scanned once per fit and the same BLAS routine is called directly
   with the same arguments. With R's reference BLAS, C++ kernels that do the
   reference routine's operations in the same order for every output element
   are used instead (self-tested per session, and the first product of each
   fit is checked against R's own result).
9. **Constructors.** One pass checks XtX for symmetry and non-finite values;
   a correlation matrix with unit diagonal is not rescaled; the original-scale
   (`bhat`, `shat`, `var_y`) path builds its XtX in one pass; constant-column
   screening for individual data reads each column once.
10. **Post-processing.** See below.
11. **Other paths.** `unmappable_effects = "ash"` keeps its correlation
    matrices for the whole fit and uses a faster `mr.ash.rss` kernel;
    `unmappable_effects = "inf"` skips dead and repeated eigenspace products;
    `susie_rss_lambda` uses the round-1 techniques and reuses its p^3
    product across refine fits.

The new code is in `R/susie_fast*.R` and `src/susie_fast.cpp`,
`src/fast_*.cpp`. The upstream files have one- or two-line changes at the
call sites.

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

The fast paths apply to the standard `ss` (`susie_rss`, `susie_ss`),
`individual` (`susie`) and `rss_lambda` (`susie_rss_lambda`) data classes.
Other classes (multi-panel R lists and downstream packages' classes) use the
upstream code.

Two opt-in options, both with bit-identical results:

- `options(susieR.threads = n)`: threads for the matrix-vector kernels used
  with R's reference BLAS (default 1). Other BLAS libraries (OpenBLAS, MKL)
  use their own threading. Keep 1 when you run regions in parallel.
- `options(susieR.refine_cores = k)`: fit the candidates of each refine step
  in k forked processes (Linux/macOS; default 1). See below.

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

- **Round 2: 249 configurations** (the 35 below plus configurations for
  every changed path and its edge cases: NA/Inf inputs, null weight, prior
  weights, names, sparse matrices, model_init, strong signals, weak signals
  with dropped credible sets, ash, inf, rss_lambda, constructors,
  threads, parallel refine). On both platforms: `exact` bit-identical in all
  249 (`off` too, on Windows); `full` identical except `fit$elbo` in at most
  15 cases (largest relative difference 1.8e-15). Five deliberate error-path
  configurations stop with the same error in both builds.
- **Round 1: 35 configurations**: `susie_rss` (defaults, residual variance estimated,
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
- **The package's own test suite** (1,289 tests in round 2, including the
  new fast-path tests): same result as upstream
  (8 failed expectations and 1 error in both builds on this machine, all
  pre-existing: 6 plotting tests that need ImageMagick and 2 R_mismatch tests
  that disagree with current upstream code).
- **New tests** in `tests/testthat/test_susie_fast*.R` compare the fast
  paths with the upstream code paths on every run.

## Benchmarks

Simulated regions with realistic LD (n = 3,000 samples, p = 1,000 or 3,000
variants, 3 causal variants), L = 10 unless stated. Upstream and the fork
were run alternately, 3 rounds; speed-up is the median of the paired
per-round ratios, one thread. The machine had other load, so absolute times
are noisy; the pairing absorbs most of it. Round-1 figures (0.16.6.1) are in
[`fast/results`](fast/results).

| fit | Linux, OpenBLAS: upstream | fork | speed-up | Windows, reference BLAS: upstream | fork | speed-up |
|---|---|---|---|---|---|---|
| `susie_rss`, p = 1,000 | 0.54 s | 0.20 s | 2.6x | 0.98 s | 0.37 s | 2.6x |
| `susie_rss`, p = 3,000 | 2.06 s | 0.57 s | 3.6x | 3.26 s | 0.67 s | 4.6x |
| `susie_rss`, p = 3,000, L = 20 | 3.37 s | 0.77 s | 4.4x | | | |
| `susie_rss`, p = 3,000, residual variance estimated | 2.03 s | 0.64 s | 3.2x | | | |
| `susie_rss`, p = 1,000, refine | 2.42 s | 0.46 s | 5.2x | 4.30 s | 0.83 s | 5.0x |
| `susie`, n = 3,000, p = 1,000 | 1.62 s | 0.42 s | 3.6x | 3.86 s | 0.97 s | 4.0x |
| `susie`, n = 3,000, p = 3,000 | 3.03 s | 0.77 s | 3.9x | 6.99 s | 1.28 s | 5.5x |
| `susie`, n = 3,000, p = 1,000, refine | 7.03 s | 1.29 s | 5.4x | | | |

Round 1 alone gave 1.6-3.2x on the same kinds of fits. The work packages
for the other paths measured, against upstream: SuSiE-ash 1.5-2.1x,
SuSiE-inf 1.1-1.4x (the eigendecomposition is unchanged), `susie_rss_lambda`
1.1-2.3x (more with refine), and opt-in parallel refine a further 1.1-1.5x
with 3 cores.

## Reproducing

[`fast/`](fast) holds the tools: `make_data.R` (simulated regions),
`run_suite.R` + `compare.R` (the configuration comparison; the configurations
are in `fast/cases/*.R`),
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
