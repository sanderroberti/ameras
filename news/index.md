# Changelog

## ameras (development version)

### Bug fixes

- Fixed a bug causing the `inpar` check in
  [`ameras()`](https://ameras.sanderroberti.com/reference/ameras.md) to
  generate an error when it should not.
- Fixed an issue where FMA generated an error instead of returning
  `NULL` for generated samples when all individual fits were excluded.
- Fixed an issue where setting `keep.data=FALSE` and passing data to
  [`confint()`](https://rdrr.io/r/stats/confint.html) failed an internal
  check when effect modifiers were present.
- Fixed an issue where profile likelihood confidence intervals for the
  ERC method of the proportional hazards family were computed using the
  non-ERC likelihood, silently ignoring the measurement error
  correction.
- Removed the `isSymmetric` check for FMA variance matrices, which
  caused platform-dependent differences in included realizations due to
  numerical differences between Cholesky and `solve`-based computations.
  The Cholesky-based variance computation now used guarantees exact
  symmetry without an explicit check.

### Improvements

- Reduced memory usage and improved computation speed for large
  datasets: removed the use of N x N matrices for ERC for the
  proportional hazards family, and precomputed centered dose matrices
  for ERC are now stored and reused across likelihood evaluations rather
  than recomputed at each call.
- [`confint()`](https://rdrr.io/r/stats/confint.html) no longer
  recomputes confidence intervals by default if they have already been
  computed. Use `force=TRUE` to recompute with different settings.
  [`confint()`](https://rdrr.io/r/stats/confint.html) now also prints
  the computed confidence intervals.
- FMA and BMA output now includes a variance-covariance matrix `vcov`,
  computed from the model-averaged posterior samples.
- It is now possible to use a `dosevars` variable defined locally (e.g.,
  within a simulation script) through `all_of(dosevars)` in the formula
  passed to
  [`ameras()`](https://ameras.sanderroberti.com/reference/ameras.md).

### New methods and accessors

- [`residuals()`](https://rdrr.io/r/stats/residuals.html): computes
  Pearson, deviance, and response residuals for all supported families,
  and Schoenfeld residuals for `family="prophaz"`, supporting both raw
  and scaled versions following Grambsch and Therneau (1994).
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html): diagnostic
  plots including residuals versus fitted values and normal Q-Q plots.
  For `family="prophaz"`, Schoenfeld residual plots are produced to
  assess the proportional hazards assumption.
- [`vcov()`](https://rdrr.io/r/stats/vcov.html): extracts the
  variance-covariance matrix of parameter estimates for one or more
  estimation methods.
- [`included_realizations()`](https://ameras.sanderroberti.com/reference/included_realizations.md):
  returns the indices of realizations included in FMA and BMA model
  averaging.
- [`Rhat()`](https://ameras.sanderroberti.com/reference/Rhat.md):
  returns the Gelman-Rubin convergence diagnostics and effective sample
  sizes for BMA results.
- [`summary_table()`](https://ameras.sanderroberti.com/reference/summary.md):
  extracts the summary table from a `summary.amerasfit` object as a data
  frame, for programmatic access to parameter estimates, standard
  errors, and confidence intervals.

## ameras 0.3.0

CRAN release: 2026-05-07

### Breaking changes

- Replaced all occurrences of `replicate` with `realization`, including
  in names of arguments and output. As a result, the argument
  `included.replicates.BMA` is now deprecated and will be removed in
  version 1.0.0. Use `included.realizations.BMA` instead.

### New features

- [`ecdfplot()`](https://ameras.sanderroberti.com/reference/ecdfplot.md)
  now has an argument `show.mean` (default `TRUE`) to toggle whether to
  plot curves for the distribution of the mean doses across realizations
  and individuals.

### Minor improvements and fixes

- Shortened column names `CI.lowerbound` and `CI.upperbound` in the
  [`summary()`](https://rdrr.io/r/base/summary.html) output to
  `CI.lower` and `CI.upper`, respectively.

- Substantially increased ERC computation speed for the `clogit` and
  `prophaz` families.

## ameras 0.2.0

CRAN release: 2026-04-26

### Breaking changes

- Confidence intervals are no longer computed inside
  [`ameras()`](https://ameras.sanderroberti.com/reference/ameras.md).
  The arguments `CI`, `params.profCI`, `maxit.profCI`, and `tol.profCI`
  are deprecated and will be removed in version 1.0.0. Use the new
  [`confint()`](https://rdrr.io/r/stats/confint.html) method instead.
  See
  [`?confint.amerasfit`](https://ameras.sanderroberti.com/reference/confint.md)
  for details.

- The direct argument interface to
  [`ameras()`](https://ameras.sanderroberti.com/reference/ameras.md) is
  deprecated and will be removed in version 1.0.0. The arguments `Y`,
  `dosevars`, `doseRRmod`, `deg`, `M`, `X`, `offset`, `entry`, `exit`,
  and `setnr` are deprecated. Please use the new formula interface
  instead. See
  [`?ameras`](https://ameras.sanderroberti.com/reference/ameras.md) for
  details.

### New features

- Implemented a formula interface for
  [`ameras()`](https://ameras.sanderroberti.com/reference/ameras.md).
  The dose variable is specified using the special `dose()` term, which
  supports tidyselect syntax for selecting dose columns and allows
  specifying the dose-response model and effect modifiers directly in
  the formula. See
  [`?ameras`](https://ameras.sanderroberti.com/reference/ameras.md) for
  details and examples.

- Added
  [`confint.amerasfit()`](https://ameras.sanderroberti.com/reference/confint.md)
  for computing confidence intervals separately from model fitting. See
  [`?confint.amerasfit`](https://ameras.sanderroberti.com/reference/confint.md)
  for details.

- Added function
  [`ecdfplot()`](https://ameras.sanderroberti.com/reference/ecdfplot.md)
  for exploratory visualization of the dose realizations before model
  fitting.

### General improvements

- Reduced memory usage for large datasets.

  - Removed the use of an N x N matrix for ERC for the Poisson family,
    improving both memory and computation speed.
  - Removed internal duplication of data for RC and ERC for all
    families.

- [`summary.amerasfit()`](https://ameras.sanderroberti.com/reference/summary.md)
  now only includes confidence interval columns after they have been
  computed via [`confint()`](https://rdrr.io/r/stats/confint.html).
  Before calling [`confint()`](https://rdrr.io/r/stats/confint.html), a
  note is printed directing the user to compute confidence intervals.

- Profile likelihood confidence interval bounds now include p-values in
  the summary table, making it easier to assess the accuracy of the
  root-finding algorithm.

- Removed memoization of the profile likelihood function as it was
  causing issues and likely not providing much benefit.

### New arguments

- `keep.data` added to
  [`ameras()`](https://ameras.sanderroberti.com/reference/ameras.md)
  (default `TRUE`). When `TRUE`, the data are stored on the returned
  `amerasfit` object, which is required for profile likelihood
  confidence interval computation via
  [`confint()`](https://rdrr.io/r/stats/confint.html) without
  re-supplying the data. Set to `FALSE` to reduce memory usage for large
  datasets, in which case the data must be supplied to
  [`confint()`](https://rdrr.io/r/stats/confint.html) explicitly. See
  [`?ameras`](https://ameras.sanderroberti.com/reference/ameras.md) and
  [`?confint.amerasfit`](https://ameras.sanderroberti.com/reference/confint.md)
  for details.

## ameras 0.1.1

CRAN release: 2026-03-29

- Initial CRAN submission.
