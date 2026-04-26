# Changelog

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
