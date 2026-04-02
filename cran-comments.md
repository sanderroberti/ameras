## Resubmission
This is a resubmission. In this version I have:

* Wrapped the example for the ameras function in \donttest as it was taking slightly too long

## Comments

* This is a new release.
* Building the vignettes modelfitting and confidenceintervals takes a very long time. Code chunks of these vignettes are not evaluated on CRAN, and should therefore not be re-built.
* The examples for ameras and traceplot take longer than 5 seconds, and are therefore wrapped in \donttest

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Sander Roberti <sander.roberti@nih.gov>'

New submission

Possibly misspelled words in DESCRIPTION:
  Kwon (29:91, 30:35)
  Stayner (29:17)
  al (26:192, 27:68, 28:39, 29:28, 29:99, 30:43)
  dosimetry (26:165)
  et (26:189, 27:65, 28:36, 29:25, 29:96, 30:40)

  These are not true misspellings.