## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Sander Roberti <sander.roberti@nih.gov>'

Possibly misspelled words in DESCRIPTION:
  Kwon (29:91, 30:35)
  Stayner (29:17)
  al (26:192, 27:68, 28:39, 29:28, 29:99, 30:43)
  dosimetry (26:165)
  et (26:189, 27:65, 28:36, 29:25, 29:96, 30:40)

  These are not true misspellings.

## Breaking changes and deprecations

This release introduces deprecations that will become breaking changes
in version 1.0.0:

* The direct argument interface to `ameras()` has been deprecated in
  favor of a new formula interface. The deprecated arguments (`Y`,
  `dosevars`, `doseRRmod`, `deg`, `M`, `X`, `offset`, `entry`, `exit`, 
  `setnr`) still work in this release but produce deprecation
  warnings via the `lifecycle` package directing users to the new
  interface.

* The arguments `CI`, `params.profCI`, `maxit.profCI`, and `tol.profCI`
  have been deprecated in `ameras()` in favor of the new
  `confint.amerasfit()` method. These arguments still work in this
  release but produce deprecation warnings.

All deprecated functionality remains fully operational in this release.
No existing user code will break. The breaking changes will be
introduced in version 1.0.0 with sufficient advance notice.

## Updated dependencies

* `lifecycle` has been added to Imports for managing deprecation
  warnings.
* `tidyselect` has been added to Imports for the formula interface.
* `memoize` has been removed from Imports as it is no longer being used.

## New suggested packages

* `tidyr`, `dplyr`, `patchwork`, and `scales` have been added to Suggests for
  `ecdfplot()`. They are loaded conditionally via `requireNamespace()`
  and are not required for standard use.

## Comments

* Building the vignettes modelfitting and confidenceintervals takes a very long time.
  Code chunks of these vignettes are not evaluated on CRAN, and should therefore 
  not be re-built.

* Examples for `ameras()` `traceplot()` are wrapped in \donttest as they take longer than 5 seconds


