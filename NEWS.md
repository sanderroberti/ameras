# ameras (development version)

## Bug fixes

* Fixed a BMA dimensionality error for models with a single formula-based effect modifier.
* Fixed a numerical issue for subgroup-coded dose modifiers where likelihood evaluations could create invalid internal contrast values during optimization, even though subgroup-specific dose effects were properly transformed.
* Fixed an issue where `optim.method` was not forwarded to FMA fits.

## Improvements

* Removed the redundant `runtime` field from newly created method results. Structured fitting, confidence interval, and total computation times remain available through `timing`; older runtime-only fitted objects remain supported.

# ameras 0.5.1

## Bug fixes

* Fixed an issue where large differences in AIC values between dose realizations could cause FMA weights to become `NaN`, resulting in an error.
* Fixed an issue where a subsequent call to `confint()` would print intervals for dose-related parameters only. `confint()` now prints intervals for all parameters.
* Fixed validation for proportional hazards models specified with `Surv(entry, exit, status)`. `ameras()` now checks the observed entry and exit time values and errors when any subject has `entry > exit`.
* Fixed conditional logistic regression fits with `strata()` terms that use a matched set column name other than `setnr`.
* Fixed BMA handling for models with a single sampled model parameter, including during sample-based confidence interval construction.
* Fixed an issue where an FMA realization with a finite fitted Hessian but non-finite sampling covariance could trigger a low-level `rmvnorm()` error. Such realizations are now excluded with a warning before FMA weights and samples are computed.
* BMA now errors clearly if `included.realizations.BMA` leaves fewer than two dose realizations, since one-realization analyses should use RC.
* Conditional logistic regression now explicitly excludes matched sets of size 1 and matched sets with no cases, and errors for matched sets with more than one case.
* Proportional hazards regression now explicitly excludes rows with `entry == exit`, which have zero follow-up time under the `entry < t <= exit` risk-set convention.
* Profile likelihood confidence intervals now handle unusable Hessian matrices and difficult bounds more robustly. The search uses explicit adaptive bracketing and returns `NA` with diagnostic warnings when a bound cannot be reliably solved.

## New features

* Added formula-based effect modifier syntax. Use `modifier = ~ M` for reference-plus-contrast coding and `modifier = ~ 0 + M` or `modifier = ~ M - 1` for subgroup-specific dose-effect parameters. Formula modifiers support binary numeric/logical variables and factors with treatment coding based on the factor level order. The older `modifier = M1 + M2` syntax is deprecated.
* BMA now supports subgroup-specific formula modifiers. Priors and MCMC sampling remain on the internal reference-plus-contrast scale, while reported summaries, stored samples, diagnostics, and sample-based intervals use the subgroup-specific scale.
* `ameras()` now supports a `na.action` argument for model-input missing values after formula terms have been expanded. By default it follows `getOption("na.action")` (typically `na.omit`), stores the omitted-row action on the fitted object, and reapplies the same policy when data are supplied later for objects fitted with `keep.data=FALSE`. `na.fail` and `na.pass` are also supported. `na.exclude` is accepted for fitting, but residuals and diagnostic plots are returned for the fitted rows rather than padded back to the originally supplied row count.
* FMA realization-specific fits can now use the `future` framework via `future.apply` when available. Users can enable parallel execution by setting a `future` plan before calling `ameras()`, and `future.chunk.size.FMA` controls the chunk size passed to `future.apply::future_lapply()`.
* Right-hand-side covariate terms now support namespace-qualified model-matrix basis functions such as `splines::ns()` and `splines::bs()`. For `keep.data=FALSE` fits, later calls that use supplied data apply the same expanded covariate columns as the original fit.
* Added `convergence()` for `amerasfit` objects to extract or recompute optimizer gradient diagnostics for RC, ERC, and MCML fits.
* Added `dose_lrt()` for global and individual likelihood-ratio tests of dose-related parameters in RC, ERC, and MCML fits.
* Added vignettes for standard analyses with one dose realization, effect modification, manual FMA from realization-specific RC fits, and parallel FMA with the `future` framework.

## Improvements

* FMA now always reports the number of excluded realizations, rather than warning only when more than 20% of realizations are excluded.
* FMA no longer assigns a minimum sample size of 1 to every realization. Realizations with model-averaging weights that round to 0 samples are excluded, and the number excluded for this reason is reported with a message.
* FMA now warns when results are based on only 1 or 2-5 realizations after exclusions.
* Added numerical diagnostics that warn when right-hand-side covariates appear poorly scaled or ill-conditioned, and when `optim()` reports convergence but optimizer diagnostics suggest the solution may not be fully stationary. When the Hessian is usable, the optimizer warning uses the approximate remaining objective improvement on both absolute and relative scales.
* `summary()` now reports row counts in the order applied during fitting: supplied rows, nonzero omissions by `na.action`, nonzero family-specific exclusions, and rows used for fitting when this differs from the supplied row count. For conditional logistic regression, it reports rows excluded because they belong to uninformative matched sets. For proportional hazards regression with entry times, it reports rows excluded because they have zero follow-up time.
* Streamlined information printed by `summary()` and `confint()`: columns `pval.lower` and `pval.upper` for profile likelihood intervals are no longer printed. They are still accessible within the fit object, and warnings are still printed in case an inaccurate profile likelihood bound is suspected.
* Added structured timing information to fitted method results, separating fitting time, confidence interval computation time, and total time. Printed runtime summaries now use CPU time, so time spent while the computer is asleep is not counted. The existing `runtime` field is retained as a compatibility summary and is updated when `confint()` adds confidence interval timing.
* Added validation to reject input data containing the reserved ameras column name `rcdose_ameras`.

# ameras 0.4.0

## Bug fixes

* Fixed a bug causing the `inpar` check in `ameras()` to generate an
  error when it should not.
* Fixed an issue where FMA generated an error instead of returning
  `NULL` for generated samples when all individual fits were excluded.
* Fixed an issue where setting `keep.data=FALSE` and passing data to
  `confint()` failed a validation check when effect modifiers were
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
  confidence intervals when called with `print=TRUE` (default).
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
    - Removed duplicate data storage for RC and ERC for all families.

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
