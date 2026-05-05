## R CMD check results

0 errors | 0 warnings | 0 notes


## revdepcheck results

We checked 0 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages

## Breaking changes and deprecations

This release introduces deprecations that will become breaking changes
in version 1.0.0:

* The argument `included.replicates.BMA` has been deprecated in `ameras()` 
  as the terminology of 'replicate' has been replaced with 'realization'.
  The arguments still work in this release but produces deprecation warnings.

All deprecated functionality remains fully operational in this release.
No existing user code will break. The breaking changes will be
introduced in version 1.0.0 with sufficient advance notice.



## Comments

* Building the vignettes modelfitting and confidenceintervals takes a very long time.
  Code chunks of these vignettes are not evaluated on CRAN, and should therefore 
  not be re-built.

* Examples for `ameras()` `traceplot()` are wrapped in \donttest as they take longer than 5 seconds


