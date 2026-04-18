# ameras (development version)

### New features

* Added function ecdfplot for descriptive data visualization.

### General improvements

* Removed the use of an N x N matrix for ERC for the Poisson family to better handle large data and speed up likelihood computation.
* Removed duplication of data for RC and ERC for all families to better handle large data.

### User interface changes

* Implemented a formula interface for ameras.
* Added a print method for the amerasfit class.
* Moved the computation of confidence intervals to a new method confint. The ameras function now returns an amerasfit object without confidence intervals. The method confint takes this object and returns the same object with confidence intervals added. Importantly, this means ameras no longer uses arguments related to confidence intervals.
* Added additional components to the output of ameras for use by the new confint method.
* Added the argument keep.data to ameras (default TRUE). If TRUE, the data object is included in the amerasfit object for future profile likelihood confidence interval computations.
* The summary method for amerasfit now only prints confidence intervals after they have been computed using confint.

# ameras 0.1.1

* Initial CRAN submission.
