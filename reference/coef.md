# Estimated Coefficients for an amerasfit Object

Returns the estimated parameters of a fitted `amerasfit` object. For an
object containing one fitted method, the result is a named numeric
vector. For an object containing multiple methods, the result is a data
frame with one column per method and rows corresponding to parameters.

## Usage

``` r
# S3 method for class 'amerasfit'
coef(object, ...)
```

## Arguments

- object:

  A fitted model object of class `amerasfit`, as returned by
  [`ameras`](https://ameras.sanderroberti.com/reference/ameras.md).

- ...:

  Additional arguments, currently unused.

## Value

A named numeric vector when one method was fitted. When multiple methods
were fitted, a data frame whose column names correspond to the methods
and whose row names correspond to parameter names.

## See also

[`ameras`](https://ameras.sanderroberti.com/reference/ameras.md) for
model fitting,
[`summary`](https://ameras.sanderroberti.com/reference/summary.md) for a
summary of the fitted model including confidence intervals if computed.

## Examples

``` r

dosevars <- paste0("V", 1:10)

## Fit the model
fit <- ameras(Y.binomial~dose(V1:V10, model="ERR"), data = ameras_simulated, family = "binomial",
methods = c("RC","ERC"))
#> Fitting RC
#> Fitting ERC
## Data frame for a fit containing multiple methods
coef(fit)
#>                     RC        ERC
#> (Intercept) -0.8847176 -0.8849093
#> dose         0.8019528  0.8214156

## Vector with RC parameters
coef(fit)$RC
#> [1] -0.8847176  0.8019528
```
