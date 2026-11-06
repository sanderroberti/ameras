# Summarize an amerasfit Object

Produces a summary of a fitted `amerasfit` object, including parameter
estimates, standard errors, and confidence intervals if attached via
[`add_confint`](https://ameras.sanderroberti.com/reference/add_confint.md).

## Usage

``` r
# S3 method for class 'amerasfit'
summary(object, ...)

# S3 method for class 'summary.amerasfit'
coef(object, ...)

# S3 method for class 'summary.amerasfit'
print(x, digits = max(3, getOption("digits") - 3), ...)
```

## Arguments

- object:

  For [`summary()`](https://rdrr.io/r/base/summary.html), a fitted model
  object of class `amerasfit`, as returned by
  [`ameras`](https://ameras.sanderroberti.com/reference/ameras.md). For
  [`coef()`](https://rdrr.io/r/stats/coef.html), an object of class
  `summary.amerasfit`.

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

- `row_info`:

  A list describing row handling during model fitting, including the
  number of rows supplied, omitted by `na.action`, and used for fitting.
  For conditional logistic regression, this also includes the number of
  rows excluded because they belong to uninformative matched sets. For
  proportional hazards regression with entry times, this also includes
  the number of rows excluded because they have zero follow-up time
  (`entry == exit`).

- `coefficients`:

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
      attached via
      [`add_confint`](https://ameras.sanderroberti.com/reference/add_confint.md).

  `CI.upper`

  :   The upper confidence bound, if confidence intervals have been
      attached via
      [`add_confint`](https://ameras.sanderroberti.com/reference/add_confint.md).

  `Rhat`

  :   The Gelman-Rubin convergence diagnostic, included only when BMA
      results are present. Values above 1.05 indicate potential
      convergence problems.

  `n.eff`

  :   The effective sample size, included only when BMA results are
      present.

- `runtime_table`:

  A data frame reporting computation time in seconds for each method.
  When structured timing information is available this has columns
  `Method`, `Fit`, `CI`, and `Total`, with CPU time reported for each
  timing component. Older fitted objects that only contain the former
  `runtime` field are still supported and produce columns `Method` and
  `Runtime`.

- `total_runtime_seconds`:

  The total CPU computation time in seconds across all methods.

- `CI.computed`:

  Logical. `TRUE` if confidence intervals have been attached via
  [`add_confint`](https://ameras.sanderroberti.com/reference/add_confint.md),
  `FALSE` otherwise.

`coef.summary.amerasfit` returns the `coefficients` data frame from the
summary object.

## Details

`summary.amerasfit` collects results from all estimation methods present
in the fitted object into a single summary table. Columns for confidence
intervals are only printed if they have been attached by
[`add_confint`](https://ameras.sanderroberti.com/reference/add_confint.md).
When BMA results are present in the fitted object, the summary table
includes columns `Rhat` and `n.eff`, with `NA` values for all other
methods. Profile likelihood diagnostic columns such as `pval.lower` and
`pval.upper` remain available in the stored `CI` component but are not
printed in the summary table.

Printing the summary prints the original call to `ameras`, row counts,
the CPU runtime (total and by method), and the table described above.
Row counts reflect omissions due to missing values and model-specific
exclusions, when applicable. For conditional logistic regression,
model-specific exclusions are rows in uninformative matched sets, i.e.,
matched sets of size 1 or with no cases. For proportional hazards
regression with entry times, model-specific exclusions are rows with
zero follow-up time (`entry == exit`). When confidence intervals have
been computed, CI runtime is included in the total runtime. The
coefficient table can be accessed directly using `coef(summary(fit))` or
`summary(fit)$coefficients`.

## See also

[`ameras`](https://ameras.sanderroberti.com/reference/ameras.md) for
model fitting,
[`add_confint`](https://ameras.sanderroberti.com/reference/add_confint.md)
for attaching confidence intervals,
[`confint`](https://ameras.sanderroberti.com/reference/confint.md) for
returning interval bounds,
[`print`](https://ameras.sanderroberti.com/reference/print.md) for a
shorter printed summary,
[`coef`](https://ameras.sanderroberti.com/reference/coef.md) for
extracting coefficients.

## Examples

``` r
## Fit the model
fit <- ameras(Y.binomial~dose(V1:V10, model="ERR"), data = ameras_simulated, family = "binomial",
              methods = "RC")
#> Fitting RC

## Summary without confidence intervals
summary(fit)
#> Call:
#> ameras(formula = Y.binomial ~ dose(V1:V10, model = "ERR"), data = ameras_simulated, 
#>     family = "binomial", methods = "RC")
#> 
#> Rows: 3000
#> 
#> Total CPU runtime: 0.1 seconds
#> 
#> CPU runtime in seconds by method:
#> 
#>  Method   Fit  CI Total
#>      RC 0.058 0.0 0.058
#> 
#> Summary of coefficients by method:
#> 
#>  Method        Term Estimate      SE
#>      RC (Intercept)  -0.8847 0.07378
#>      RC        dose   0.8020 0.13751
#> 
#> Note: confidence intervals not yet computed. Use add_confint() to add them.
#> 

## Summary with confidence intervals
fit <- add_confint(fit, type = "wald.orig")
#> RC confidence intervals:
#> 
#>               lower   upper
#> (Intercept) -1.0293 -0.7401
#> dose         0.5324  1.0715
#> 
summary(fit)
#> Call:
#> ameras(formula = Y.binomial ~ dose(V1:V10, model = "ERR"), data = ameras_simulated, 
#>     family = "binomial", methods = "RC")
#> 
#> Rows: 3000
#> 
#> Total CPU runtime: 0.1 seconds
#> 
#> CPU runtime in seconds by method:
#> 
#>  Method   Fit  CI Total
#>      RC 0.058 0.0 0.058
#> 
#> Summary of coefficients by method:
#> 
#>  Method        Term Estimate      SE CI.lower CI.upper
#>      RC (Intercept)  -0.8847 0.07378  -1.0293  -0.7401
#>      RC        dose   0.8020 0.13751   0.5324   1.0715
#> 

## Access the coefficient summary directly
s <- summary(fit)
coefficient_table <- coef(s)

## Multiple methods

fit2 <- ameras(Y.binomial~dose(V1:V10, model="ERR"), data = ameras_simulated, family = "binomial",
              methods = c("RC", "ERC", "MCML"))
#> Fitting RC
#> Fitting ERC
#> Fitting MCML
fit2 <- add_confint(fit2, type = "wald.orig")
#> RC confidence intervals:
#> 
#>               lower   upper
#> (Intercept) -1.0293 -0.7401
#> dose         0.5324  1.0715
#> 
#> ERC confidence intervals:
#> 
#>               lower   upper
#> (Intercept) -1.0314 -0.7384
#> dose         0.5411  1.1018
#> 
#> MCML confidence intervals:
#> 
#>               lower   upper
#> (Intercept) -1.0193 -0.7323
#> dose         0.5236  1.0584
#> 
summary(fit2)
#> Call:
#> ameras(formula = Y.binomial ~ dose(V1:V10, model = "ERR"), data = ameras_simulated, 
#>     family = "binomial", methods = c("RC", "ERC", "MCML"))
#> 
#> Rows: 3000
#> 
#> Total CPU runtime: 94.7 seconds
#> 
#> CPU runtime in seconds by method:
#> 
#>  Method    Fit    CI  Total
#>      RC  0.059 0.000  0.059
#>     ERC 94.347 0.001 94.348
#>    MCML  0.274 0.000  0.274
#> 
#> Summary of coefficients by method:
#> 
#>  Method        Term Estimate      SE CI.lower CI.upper
#>      RC (Intercept)  -0.8847 0.07378  -1.0293  -0.7401
#>      RC        dose   0.8020 0.13751   0.5324   1.0715
#>     ERC (Intercept)  -0.8849 0.07477  -1.0314  -0.7384
#>     ERC        dose   0.8214 0.14304   0.5411   1.1018
#>    MCML (Intercept)  -0.8758 0.07323  -1.0193  -0.7323
#>    MCML        dose   0.7910 0.13644   0.5236   1.0584
#> 
```
