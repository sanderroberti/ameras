# Simple summary for an amerasfit object

Prints a simple summary of a fitted `amerasfit` object.

## Usage

``` r
# S3 method for class 'amerasfit'
print(x, digits = max(3, getOption("digits") - 3), ...)
```

## Arguments

- x:

  A fitted model object of class `amerasfit`, as returned by
  [`ameras`](https://ameras.sanderroberti.com/reference/ameras.md).

- digits:

  Number of significant digits to be printed. Default is \`max(3,
  getOption("digits") - 3)\`

- ...:

  Additional arguments, currently unused.

## Value

Prints the \`ameras\` call, number of rows and dose replicates in the
data, runtime, and model coefficients.

## See also

[`ameras`](https://ameras.sanderroberti.com/reference/ameras.md) for
model fitting,
[`summary`](https://ameras.sanderroberti.com/reference/summary.md) for a
more detailed summary of the fitted models including confidence
intervals if computed.

## Examples

``` r
# \donttest{
data("data", package = "ameras")

## Fit the model
fit <- ameras(Y.binomial~dose(V1:V10, model="ERR"), data = data, family = "binomial", 
              methods = c("RC","ERC"))
#> Fitting RC
#> Fitting ERC

## Default print
fit
#> Call:
#> ameras(formula = Y.binomial ~ dose(V1:V10, model = "ERR"), data = data, 
#>     family = "binomial", methods = c("RC", "ERC"))
#> 
#> Number of rows: 3000
#> Number of dose replicates: 10
#> 
#> Total run time: 14 seconds
#> 
#> Runtime in seconds by method:
#> 
#>  Method Runtime
#>      RC     0.0
#>     ERC    14.0
#> 
#> Estimated model parameters:
#> 
#>                  RC     ERC
#> (Intercept) -0.8847 -0.8849
#> dose         0.8020  0.8214
print(fit)
#> Call:
#> ameras(formula = Y.binomial ~ dose(V1:V10, model = "ERR"), data = data, 
#>     family = "binomial", methods = c("RC", "ERC"))
#> 
#> Number of rows: 3000
#> Number of dose replicates: 10
#> 
#> Total run time: 14 seconds
#> 
#> Runtime in seconds by method:
#> 
#>  Method Runtime
#>      RC     0.0
#>     ERC    14.0
#> 
#> Estimated model parameters:
#> 
#>                  RC     ERC
#> (Intercept) -0.8847 -0.8849
#> dose         0.8020  0.8214

## More digits
print(fit, digits=5)
#> Call:
#> ameras(formula = Y.binomial ~ dose(V1:V10, model = "ERR"), data = data, 
#>     family = "binomial", methods = c("RC", "ERC"))
#> 
#> Number of rows: 3000
#> Number of dose replicates: 10
#> 
#> Total run time: 14 seconds
#> 
#> Runtime in seconds by method:
#> 
#>  Method Runtime
#>      RC     0.0
#>     ERC    14.0
#> 
#> Estimated model parameters:
#> 
#>                   RC      ERC
#> (Intercept) -0.88472 -0.88491
#> dose         0.80195  0.82142
# }
```
