# Extract Variance-Covariance Matrices from an amerasfit Object

Extracts the variance-covariance matrices of the parameter estimates
from a fitted `amerasfit` object.

## Usage

``` r
# S3 method for class 'amerasfit'
vcov(object, methods = c("RC", "ERC", "MCML", "FMA", "BMA"), ...)
```

## Arguments

- object:

  A fitted model object of class `amerasfit`, as returned by
  [`ameras`](https://ameras.sanderroberti.com/reference/ameras.md).

- methods:

  Character vector specifying which estimation methods to extract
  variance-covariance matrices for. One or more of `"RC"`, `"ERC"`,
  `"MCML"`, `"FMA"`, and `"BMA"`. Defaults to all. Methods that were not
  run are skipped; an error is raised if none of the requested methods
  are present.

- ...:

  Additional arguments, currently unused.

## Value

If a single method is requested or only one method was run, a named
numeric matrix with rows and columns corresponding to model parameters.
If multiple methods are requested and multiple methods were run, a named
list of such matrices, one per method.

Note that for FMA and BMA the variance-covariance matrix is based on the
generated samples rather than a Hessian-based approximation, and may
differ in interpretation from the RC, ERC, and MCML variance-covariance
matrices. For BMA, the reliability of the variance-covariance matrix
depends on MCMC convergence as indicated by `Rhat`.

## See also

[`ameras`](https://ameras.sanderroberti.com/reference/ameras.md) for
model fitting,
[`coef`](https://ameras.sanderroberti.com/reference/coef.md) for
extracting coefficients,
[`summary.amerasfit`](https://ameras.sanderroberti.com/reference/summary.md)
for a summary including standard errors,
[`vcov`](https://rdrr.io/r/stats/vcov.html) for the generic function.

## Examples

``` r
# \donttest{
data("data", package="ameras")
dosevars <- paste0("V", 1:10)

fit <- ameras(Y.binomial ~ dose(all_of(dosevars), model="ERR"),
              data=data, family="binomial", methods=c("RC", "ERC"))
#> Fitting RC
#> Fitting ERC

## Extract vcov for all available methods
vcov(fit)
#> $RC
#>              (Intercept)         dose
#> (Intercept)  0.005443953 -0.008716237
#> dose        -0.008716237  0.018908197
#> 
#> $ERC
#>              (Intercept)         dose
#> (Intercept)  0.005589966 -0.009297482
#> dose        -0.009297482  0.020459538
#> 

## Extract vcov for a single method, returns a matrix
vcov(fit, methods="RC")
#>              (Intercept)         dose
#> (Intercept)  0.005443953 -0.008716237
#> dose        -0.008716237  0.018908197

## Extract vcov for multiple methods, returns a named list
vcov(fit, methods=c("RC", "ERC"))
#> $RC
#>              (Intercept)         dose
#> (Intercept)  0.005443953 -0.008716237
#> dose        -0.008716237  0.018908197
#> 
#> $ERC
#>              (Intercept)         dose
#> (Intercept)  0.005589966 -0.009297482
#> dose        -0.009297482  0.020459538
#> 
# }
```
