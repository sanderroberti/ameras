## R CMD check results

0 errors | 0 warnings | 0 notes


## revdepcheck results

We checked 0 reverse dependencies, comparing R CMD check results across CRAN
and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages


## Changes in this version

This is a minor release relative to version 0.4.0. It adds formula-based
effect modifier syntax, subgroup-specific dose-effect parameterization,
subgroup-coded BMA output, likelihood-ratio tests for dose-related parameters,
optimizer convergence diagnostics, optional `future`-based FMA execution,
support for expanded model-matrix terms such as `splines::ns()`, and new
vignettes.

It also includes several robustness fixes for FMA, BMA, conditional logistic
regression, proportional hazards regression, profile likelihood confidence
intervals, missing-value handling, and printed summaries. Full details are
provided in `NEWS.md`.


## Deprecations and compatibility notes

The old effect-modifier syntax, for example `modifier = M1 + M2` inside
`dose()`, is deprecated in favor of formula syntax such as
`modifier = ~ M1 + M2`. The old syntax remains supported in this release and
emits a lifecycle deprecation warning.

No previously deprecated arguments were removed in this release.
