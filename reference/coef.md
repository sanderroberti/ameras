# Estimated coefficients for an amerasfit object

Returns a data frame with all the parameters of a fitted `amerasfit`
object. The resulting object has a column for every method supplied to
\`methods\` when calling \`ameras\`, with rows corresponding to
parameters.

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

Data frame with estimated model parameters. Column names correspond to
the \`methods\` used in the \`ameras\` call, and row names correspond to
parameter names.

## See also

[`ameras`](https://ameras.sanderroberti.com/reference/ameras.md) for
model fitting,
[`summary`](https://ameras.sanderroberti.com/reference/summary.md) for a
summary of the fitted model including confidence intervals if computed.

## Examples

``` r
# \donttest{
data("data", package = "ameras")
dosevars <- paste0("V", 1:10)

## Fit the model
fit <- ameras(Y.binomial~dose(V1:V10, model="ERR"), data = data, family = "binomial", 
methods = c("RC","ERC"))
#> Fitting RC
#> Fitting ERC
## Full matrix
coef(fit)
#>                     RC        ERC
#> (Intercept) -0.8847176 -0.8849093
#> dose         0.8019528  0.8214156

## Vector with RC parameters
coef(fit)$RC
#> [1] -0.8847176  0.8019528
# }
```
