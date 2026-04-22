# Summarize an amerasfit object

Produces a summary of a fitted `amerasfit` object, including parameter
estimates, standard errors, and confidence intervals if computed via
[`confint`](https://ameras.sanderroberti.com/reference/confint.md).

## Usage

``` r
# S3 method for class 'amerasfit'
summary(object, ...)

# S3 method for class 'summary.amerasfit'
print(x, digits = max(3, getOption("digits") - 3), ...)
```

## Arguments

- object:

  A fitted model object of class `amerasfit`, as returned by
  [`ameras`](https://ameras.sanderroberti.com/reference/ameras.md).

- x:

  An object of class `summary.amerasfit`, as returned by
  `summary.amerasfit`.

- digits:

  The number of significant digits to use. Defaults to
  `max(3, getOption("digits") - 3)`.

- ...:

  Additional arguments, currently unused.

## Value

`summary.amerasfit` returns an object of class `summary.amerasfit`,
which is a list containing the following elements:

- `call`:

  The matched call from the original
  [`ameras`](https://ameras.sanderroberti.com/reference/ameras.md)
  invocation.

- `summary_table`:

  A data frame with one row per parameter per method, containing
  columns:

  `Method`

  :   The estimation method (RC, ERC, MCML, FMA, or BMA).

  `Term`

  :   The parameter name.

  `Estimate`

  :   The parameter estimate.

  `SE`

  :   The standard error.

  `CI.lower`

  :   The lower confidence bound, if confidence intervals have been
      computed via
      [`confint`](https://ameras.sanderroberti.com/reference/confint.md).
      `NA` otherwise.

  `CI.upper`

  :   The upper confidence bound, if confidence intervals have been
      computed via
      [`confint`](https://ameras.sanderroberti.com/reference/confint.md).
      `NA` otherwise.

  `Rhat`

  :   The Gelman-Rubin convergence diagnostic, included only when BMA
      results are present. Values substantially above 1.05 indicate
      potential convergence problems.

  `n.eff`

  :   The effective sample size, included only when BMA results are
      present.

- `runtime_table`:

  A data frame with columns `Method` and `Runtime`, reporting the
  computation time in seconds for each method.

- `total_runtime_seconds`:

  The total computation time in seconds across all methods.

- `has_CI`:

  Logical. `TRUE` if confidence intervals have been computed via
  [`confint`](https://ameras.sanderroberti.com/reference/confint.md),
  `FALSE` otherwise.

## Details

`summary.amerasfit` collects results from all estimation methods present
in the fitted object into a single summary table. Columns for confidence
intervals are only printed if they have been computed by
[`confint`](https://ameras.sanderroberti.com/reference/confint.md). When
BMA results are present in the fitted object, the summary table includes
columns `Rhat` and `n.eff`, with `NA` values for all other methods.

## See also

[`ameras`](https://ameras.sanderroberti.com/reference/ameras.md) for
model fitting,
[`confint`](https://ameras.sanderroberti.com/reference/confint.md) for
computing confidence intervals,
[`print`](https://ameras.sanderroberti.com/reference/print.md) for a
shorter printed summary,
[`coef`](https://ameras.sanderroberti.com/reference/coef.md) for
extracting coefficients.

## Examples

``` r
data("data", package = "ameras")

## Fit the model
fit <- ameras(Y.binomial~dose(V1:V10, model="ERR"), data = data, family = "binomial",
              methods = "RC")
#> Fitting RC

## Summary without confidence intervals
summary(fit)
#> Call:
#> ameras(formula = Y.binomial ~ dose(V1:V10, model = "ERR"), data = data, 
#>     family = "binomial", methods = "RC")
#> 
#> Total run time: 0 seconds
#> 
#> Runtime in seconds by method:
#> 
#>  Method Runtime
#>      RC     0.0
#> 
#> Summary of coefficients by method:
#> 
#>  Method        Term Estimate      SE
#>      RC (Intercept)  -0.8847 0.07378
#>      RC        dose   0.8020 0.13751
#> 
#> Note: confidence intervals not yet computed. Use confint() to add them.

## Summary with confidence intervals
fit <- confint(fit, method = "wald.orig")
#> Obtaining profile likelihood CI for dose
summary(fit)
#> Call:
#> ameras(formula = Y.binomial ~ dose(V1:V10, model = "ERR"), data = data, 
#>     family = "binomial", methods = "RC")
#> 
#> Total run time: 0 seconds
#> 
#> Runtime in seconds by method:
#> 
#>  Method Runtime
#>      RC     0.0
#> 
#> Summary of coefficients by method:
#> 
#>  Method        Term Estimate      SE CI.lowerbound CI.upperbound pval.lower
#>      RC (Intercept)  -0.8847 0.07378            NA            NA         NA
#>      RC        dose   0.8020 0.13751        0.5648         1.112    0.05024
#>  pval.upper
#>          NA
#>     0.04959

## Access the summary table directly
s <- summary(fit)
s$summary_table
#>   Method        Term   Estimate         SE CI.lowerbound CI.upperbound
#> 1     RC (Intercept) -0.8847176 0.07378315            NA            NA
#> 2     RC        dose  0.8019528 0.13750708     0.5647607      1.111609
#>   pval.lower pval.upper
#> 1         NA         NA
#> 2 0.05024383 0.04958787

## Multiple methods
if (FALSE) { # \dontrun{
fit2 <- ameras(Y.binomial~dose(V1:V10, model="ERR"), data = data, family = "binomial",
              methods = c("RC", "ERC", "MCML"))
fit2 <- confint(fit2, method = "wald.orig")
summary(fit2)
} # }
```
