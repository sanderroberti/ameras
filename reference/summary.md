# Summarize an amerasfit Object

Produces a summary of a fitted `amerasfit` object, including parameter
estimates, standard errors, and confidence intervals if computed via
[`confint`](https://ameras.sanderroberti.com/reference/confint.md).

## Usage

``` r
# S3 method for class 'amerasfit'
summary(object, ...)

# S3 method for class 'amerasfit'
summary_table(object, ...)

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

  The matched call from the original call to
  [`ameras`](https://ameras.sanderroberti.com/reference/ameras.md).

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

  `CI.upper`

  :   The upper confidence bound, if confidence intervals have been
      computed via
      [`confint`](https://ameras.sanderroberti.com/reference/confint.md).

  `pval.lower`

  :   p-value associated with the lower bound of the profile likelihood
      CI, the assess validity of the obtained bound. Only included if
      profile likelihood confidence intervals were computed via
      [`confint`](https://ameras.sanderroberti.com/reference/confint.md).

  `pval.upper`

  :   p-value associated with the upper bound of the profile likelihood
      CI, the assess validity of the obtained bound. Only included if
      profile likelihood confidence intervals were computed via
      [`confint`](https://ameras.sanderroberti.com/reference/confint.md).

  `Rhat`

  :   The Gelman-Rubin convergence diagnostic, included only when BMA
      results are present. Values above 1.05 indicate potential
      convergence problems.

  `n.eff`

  :   The effective sample size, included only when BMA results are
      present.

- `runtime_table`:

  A data frame with columns `Method` and `Runtime`, reporting the
  computation time in seconds for each method.

- `total_runtime_seconds`:

  The total computation time in seconds across all methods.

- `CI.computed`:

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

Printing the summary prints the original call to `ameras`, the runtime
(total and by method), and the table described above. This table can
also be accessed directly (i.e., to retrieve confidence intervals) using
`summary_table`.

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
#> 

## Summary with confidence intervals
fit <- confint(fit, method = "wald.orig")
#> Obtaining profile likelihood CI for dose
#> RC confidence intervals:
#> 
#>       lower upper pval.lower pval.upper iter.lower iter.upper
#> dose 0.5648 1.112    0.05024    0.04959          5          4
#> 
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
#>  Method        Term Estimate      SE CI.lower CI.upper pval.lower pval.upper
#>      RC (Intercept)  -0.8847 0.07378       NA       NA         NA         NA
#>      RC        dose   0.8020 0.13751   0.5648    1.112    0.05024    0.04959
#> 

## Access the summary table directly
s <- summary_table(fit)

## Multiple methods
# \donttest{
fit2 <- ameras(Y.binomial~dose(V1:V10, model="ERR"), data = data, family = "binomial",
              methods = c("RC", "ERC", "MCML"))
#> Fitting RC
#> Fitting ERC
#> Fitting MCML
fit2 <- confint(fit2, method = "wald.orig")
#> Obtaining profile likelihood CI for dose
#> Obtaining profile likelihood CI for dose
#> Obtaining profile likelihood CI for dose
#> RC confidence intervals:
#> 
#>       lower upper pval.lower pval.upper iter.lower iter.upper
#> dose 0.5648 1.112    0.05024    0.04959          5          4
#> 
#> ERC confidence intervals:
#> 
#>       lower upper pval.lower pval.upper iter.lower iter.upper
#> dose 0.5728 1.144    0.04953    0.04819          6          4
#> 
#> MCML confidence intervals:
#> 
#>       lower upper pval.lower pval.upper iter.lower iter.upper
#> dose 0.5554 1.098    0.05008    0.04957          5          4
#> 
summary(fit2)
#> Call:
#> ameras(formula = Y.binomial ~ dose(V1:V10, model = "ERR"), data = data, 
#>     family = "binomial", methods = c("RC", "ERC", "MCML"))
#> 
#> Total run time: 13.7 seconds
#> 
#> Runtime in seconds by method:
#> 
#>  Method Runtime
#>      RC     0.0
#>     ERC    13.5
#>    MCML     0.2
#> 
#> Summary of coefficients by method:
#> 
#>  Method        Term Estimate      SE CI.lower CI.upper pval.lower pval.upper
#>      RC (Intercept)  -0.8847 0.07378       NA       NA         NA         NA
#>      RC        dose   0.8020 0.13751   0.5648    1.112    0.05024    0.04959
#>     ERC (Intercept)  -0.8849 0.07477       NA       NA         NA         NA
#>     ERC        dose   0.8214 0.14304   0.5728    1.144    0.04953    0.04819
#>    MCML (Intercept)  -0.8758 0.07323       NA       NA         NA         NA
#>    MCML        dose   0.7910 0.13644   0.5554    1.098    0.05008    0.04957
#> 
# }
```
