# ameras (development version)

## Bug fixes

* Fixed a bug causing the `inpar` check in `ameras()` to generate an
  error when it should not.
* Fixed an issue where FMA generated an error instead of returning
  `NULL` for generated samples when all individual fits were excluded.
* Fixed an issue where setting `keep.data=FALSE` and passing data to
  `confint()` failed an internal check when effect modifiers were
  present.
* Fixed an issue where profile likelihood confidence intervals for the
  ERC method of the proportional hazards family were computed using
  the non-ERC likelihood, silently ignoring the measurement error
  correction.
* Removed the `isSymmetric` check for FMA variance matrices, which
  caused platform-dependent differences in included realizations due
  to numerical differences between Cholesky and `solve`-based
  computations. The Cholesky-based variance computation now used
  guarantees exact symmetry without an explicit check.

## Improvements

* Reduced memory usage and improved computation speed for large
  datasets: removed the use of N x N matrices for ERC for the
  proportional hazards family, and precomputed centered dose matrices
  for ERC are now stored and reused across likelihood evaluations
  rather than recomputed at each call.
* `confint()` no longer recomputes confidence intervals by default if
  they have already been computed. Use `force=TRUE` to recompute with
  different settings. `confint()` now also prints the computed
  confidence intervals.
* FMA and BMA output now includes a variance-covariance matrix
  `vcov`, computed from the model-averaged posterior samples.
* It is now possible to use a `dosevars` variable defined locally (e.g.,
  within a simulation script) through `all_of(dosevars)` in the 
  formula passed to `ameras()`.

## New methods and accessors

* `residuals()`: computes Pearson, deviance, and response residuals
  for all supported families, and Schoenfeld residuals for
  `family="prophaz"`, supporting both raw and scaled versions
  following Grambsch and Therneau (1994).
* `plot()`: diagnostic plots including residuals versus fitted values
  and normal Q-Q plots. For `family="prophaz"`, Schoenfeld residual
  plots are produced to assess the proportional hazards assumption.
* `vcov()`: extracts the variance-covariance matrix of parameter
  estimates for one or more estimation methods.
* `included_realizations()`: returns the indices of realizations
  included in FMA and BMA model averaging.
* `Rhat()`: returns the Gelman-Rubin convergence diagnostics and
  effective sample sizes for BMA results.
* `summary_table()`: extracts the summary table from a
  `summary.amerasfit` object as a data frame, for programmatic
  access to parameter estimates, standard errors, and confidence
  intervals.


# ameras 0.3.0

## Breaking changes

* Replaced all occurrences of `replicate` with `realization`, including in names of arguments and output. As a result, the argument `included.replicates.BMA` is now deprecated and will be removed in version 1.0.0. Use `included.realizations.BMA` instead.

## New features

* `ecdfplot()` now has an argument `show.mean` (default `TRUE`) to toggle whether to plot curves for the distribution of the mean doses across realizations and individuals.

## Minor improvements and fixes

* Shortened column names `CI.lowerbound` and `CI.upperbound` in the `summary()` output to `CI.lower` and `CI.upper`, respectively.

* Substantially increased ERC computation speed for the `clogit` and `prophaz` families.




# ameras 0.2.0

## Breaking changes

* Confidence intervals are no longer computed inside `ameras()`. The
  arguments `CI`, `params.profCI`, `maxit.profCI`, and `tol.profCI` are
  deprecated and will be removed in version 1.0.0. Use the new
  `confint()` method instead. See `?confint.amerasfit` for details.

* The direct argument interface to `ameras()` is deprecated and will be
  removed in version 1.0.0. The arguments `Y`, `dosevars`, `doseRRmod`,
  `deg`, `M`, `X`, `offset`, `entry`, `exit`, and `setnr` are
  deprecated. Please use the new formula interface instead.
  See `?ameras` for details.


## New features

* Implemented a formula interface for `ameras()`. The dose variable is
  specified using the special `dose()` term, which supports tidyselect
  syntax for selecting dose columns and allows specifying the
  dose-response model and effect modifiers directly in the formula.
  See `?ameras` for details and examples.
  
* Added `confint.amerasfit()` for computing confidence intervals
  separately from model fitting. See `?confint.amerasfit` for details.
  
* Added function `ecdfplot()` for exploratory visualization of the dose 
  realizations before model fitting.


## General improvements

* Reduced memory usage for large datasets.
    - Removed the use of an N x N matrix for ERC for the Poisson family, improving both memory and computation speed.
    - Removed internal duplication of data for RC and ERC for all families.

* `summary.amerasfit()` now only includes confidence interval columns
  after they have been computed via `confint()`. Before calling
  `confint()`, a note is printed directing the user to compute
  confidence intervals.

* Profile likelihood confidence interval bounds now include p-values
  in the summary table, making it easier to assess the accuracy of
  the root-finding algorithm.

* Removed memoization of the profile likelihood function as it was causing issues and likely not providing much benefit.
  
## New arguments

* `keep.data` added to `ameras()` (default `TRUE`). When `TRUE`, the
  data are stored on the returned `amerasfit` object, which is required
  for profile likelihood confidence interval computation via `confint()`
  without re-supplying the data. Set to `FALSE` to reduce memory usage
  for large datasets, in which case the data must be supplied to
  `confint()` explicitly. See `?ameras` and `?confint.amerasfit` for details.
  

# ameras 0.1.1

* Initial CRAN submission.
