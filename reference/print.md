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
dosevars <- paste0("V", 1:10)

## Fit the model
fit <- ameras(data = data, family = "binomial", Y = "Y.binomial",
              dosevars = dosevars, methods = c("RC","ERC"), doseRRmod = "ERR")
#> Fitting RC
#> Fitting ERC

## Default print
fit
#> Call:
#> ameras(data = data, family = "binomial", Y = "Y.binomial", dosevars = dosevars, 
#>     methods = c("RC", "ERC"), doseRRmod = "ERR")
#> 
#> Number of rows: 3000
#> Number of dose replicates: 10
#> 
#> Total run time: 17.5 seconds
#> 
#> Runtime in seconds by method:
#> 
#>  Method Runtime
#>      RC     0.0
#>     ERC    17.5
#> 
#> Table of estimated model parameters:
#> 
#>                  RC     ERC
#> (Intercept) -0.8847 -0.8849
#> dose         0.8020  0.8214
print(fit)
#> Call:
#> ameras(data = data, family = "binomial", Y = "Y.binomial", dosevars = dosevars, 
#>     methods = c("RC", "ERC"), doseRRmod = "ERR")
#> 
#> Number of rows: 3000
#> Number of dose replicates: 10
#> 
#> Total run time: 17.5 seconds
#> 
#> Runtime in seconds by method:
#> 
#>  Method Runtime
#>      RC     0.0
#>     ERC    17.5
#> 
#> Table of estimated model parameters:
#> 
#>                  RC     ERC
#> (Intercept) -0.8847 -0.8849
#> dose         0.8020  0.8214

## More digits
print(fit, digits=5)
#> Call:
#> ameras(data = data, family = "binomial", Y = "Y.binomial", dosevars = dosevars, 
#>     methods = c("RC", "ERC"), doseRRmod = "ERR")
#> 
#> Number of rows: 3000
#> Number of dose replicates: 10
#> 
#> Total run time: 17.5 seconds
#> 
#> Runtime in seconds by method:
#> 
#>  Method Runtime
#>      RC     0.0
#>     ERC    17.5
#> 
#> Table of estimated model parameters:
#> 
#>                   RC      ERC
#> (Intercept) -0.88472 -0.88491
#> dose         0.80195  0.82142
# }
```
