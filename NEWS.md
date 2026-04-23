# ameras (development version)

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
