# Residuals for an amerasfit Object

Computes residuals for a fitted `amerasfit` object.

## Usage

``` r
# S3 method for class 'amerasfit'
residuals(object,
          method = "RC", type = c("pearson", "deviance", "response", "schoenfeld"),
          data = NULL, dose_col = NULL, scaled_schoenfeld = TRUE, ...)
```

## Arguments

- object:

  A fitted model object of class `amerasfit`, as returned by
  [`ameras`](https://ameras.sanderroberti.com/reference/ameras.md).

- method:

  Character string specifying which estimation method to compute
  residuals for. One of `"RC"`, `"ERC"`, `"MCML"`, `"FMA"`, or `"BMA"`.
  Defaults to `"RC"`.

- type:

  The type of residuals to compute. Let \\\mu_i\\ denote the fitted
  value for individual \\i\\, defined as the predicted mean response
  given the covariates and specified `dose_col`. For conditional
  logistic regression, these are the predicted conditional probabilities
  within matched sets. One of:

  `"pearson"`

  :   Pearson residuals. For `"gaussian"`: \\(Y_i - \mu_i) / \sigma\\
      where \\\sigma\\ is the estimated residual standard deviation
      (note that this differs from
      [`glm`](https://rdrr.io/r/stats/glm.html) which returns the raw
      (response) residuals). For `"binomial"` or `"clogit"`: \\(Y_i -
      \mu_i) / \sqrt{\mu_i(1 - \mu_i)}\\. For `"poisson"`: \\(Y_i -
      \mu_i) / \sqrt{\mu_i}\\. For `"multinomial"`: per-category Pearson
      residuals \\(Y\_{iz} - \mu\_{iz}) / \sqrt{\mu\_{iz}(1 -
      \mu\_{iz})}\\ returned as an \\N \times (Z-1)\\ matrix, where
      \\\mu\_{iz}\\ is the fitted probability of category \\z\\ for
      individual \\i\\.

  `"deviance"`

  :   Deviance residuals, defined as the signed square root of the
      individual contribution to the deviance. For `"gaussian"`: the raw
      residual \\Y_i - \mu_i\\. For `"binomial"` or `"clogit"`:
      \\\pm\sqrt{-2\log(\hat{p}\_i)}\\ where \\\hat{p}\_i = \mu_i\\ if
      \\Y_i = 1\\ and \\\hat{p}\_i = 1 - \mu_i\\ if \\Y_i = 0\\, with
      sign equal to \\\text{sign}(Y_i - \mu_i)\\. For `"poisson"`:
      \\\pm\sqrt{2(Y_i \log(Y_i/\mu_i) - (Y_i - \mu_i))}\\ for \\Y_i \>
      0\\ and \\\pm\sqrt{2\mu_i}\\ for \\Y_i = 0\\, with sign equal to
      \\\text{sign}(Y_i - \mu_i)\\. For `"multinomial"`: per-category
      deviance residuals returned as an \\N \times (Z-1)\\ matrix.

  `"response"`

  :   Raw residuals \\Y_i - \mu_i\\. For `"multinomial"` the matrix
      \\\bm{Y} - \hat{\bm{P}}\\ where \\\bm{Y}\\ is the \\N \times Z\\
      indicator matrix of observed categories and \\\hat{\bm{P}}\\ is
      the \\N \times Z\\ matrix of fitted probabilities.

  `"schoenfeld"`

  :   Schoenfeld residuals for proportional hazards-type fits. This is a
      matrix with one column for each covariate in the model, and a row
      for each event time. For each event time, the vector of raw
      Schoenfeld residuals is the observed covariate vector for the
      individual experiencing the event minus the risk-set weighted mean
      of the covariates at that event time. If
      `scaled_schoenfeld = TRUE`, the residuals are scaled using the
      full estimated coefficient covariance matrix, otherwise the raw
      Schoenfeld residuals are returned. This residual type is only used
      for `family = "prophaz"`.

  Defaults to `"pearson"` for all families other than `"prophaz"`, for
  which the default is `"schoenfeld"`.

- data:

  The original data frame used for fitting. Only required when the model
  was fitted with `keep.data=FALSE`.

- dose_col:

  Character string specifying the dose column to use when computing
  fitted values. If `NULL` (the default), the dose column is selected
  automatically based on the estimation method: the mean dose across
  realizations for RC and ERC, the realization with the highest
  likelihood for MCML and FMA, and the most frequently selected
  realization for BMA. Can be set to any dose column present in the data
  to override the automatic selection.

- scaled_schoenfeld:

  Logical; whether to compute scaled (`scaled_schoenfeld=TRUE`) or raw
  (`scaled_schoenfeld=FALSE`) Schoenfeld residuals. Not used for other
  residual types.

- ...:

  Additional arguments, currently unused.

## Value

For `family` in `"gaussian"`, `"binomial"`, and `"poisson"`, a numeric
vector containing the residuals. For `family="multinomial"`, a numeric
matrix of dimension \\N \times (Z-1)\\ where \\Z\\ is the number of
outcome categories and the reference category (the last factor level) is
excluded. Column names correspond to the non-reference factor levels.

## Details

Fitted values are computed using the mean dose across realizations
(`rcdose_ameras`) for RC and ERC, or the single "best fitting"
realization for MCML, FMA, and BMA. For `family="multinomial"`,
residuals are computed per outcome category.

## See also

[`plot.amerasfit`](https://ameras.sanderroberti.com/reference/plot.md)
for diagnostic plots using these residuals,
[`ameras`](https://ameras.sanderroberti.com/reference/ameras.md) for
model fitting, [`residuals`](https://rdrr.io/r/stats/residuals.html) for
the generic function.

## Examples

``` r
data("data", package = "ameras")
dosevars <- paste0("V", 1:10)

fit <- ameras(Y.binomial ~ dose(all_of(dosevars), model = "ERR"),
              data = data, family = "binomial", methods = "RC")
#> Error in resolve_dose_selection(sel_args, data): ℹ In argument: `all_of(dosevars)`.
#> Caused by error:
#> ! object 'dosevars' not found

## Pearson residuals (default)
res <- residuals(fit)
#> Error: object 'fit' not found
summary(res)
#> Error: object 'res' not found

## Deviance residuals
res_dev <- residuals(fit, type = "deviance")
#> Error: object 'fit' not found
summary(res_dev)
#> Error: object 'res_dev' not found

## Response residuals
res_raw <- residuals(fit, type = "response")
#> Error: object 'fit' not found
summary(res_raw)
#> Error: object 'res_raw' not found

## Specific method
# \donttest{
fit2 <- ameras(Y.binomial ~ dose(all_of(dosevars), model = "ERR"),
               data = data, family = "binomial", methods = c("RC", "ERC"))
#> Error in resolve_dose_selection(sel_args, data): ℹ In argument: `all_of(dosevars)`.
#> Caused by error:
#> ! object 'dosevars' not found
res_erc <- residuals(fit2, method = "ERC")
#> Error: object 'fit2' not found
# }

## Specific dose column
res_v1 <- residuals(fit, dose_col = "V1")
#> Error: object 'fit' not found

## With keep.data = FALSE
# \donttest{
fit3 <- ameras(Y.binomial ~ dose(all_of(dosevars), model = "ERR"),
               data = data, family = "binomial", methods = "RC",
               keep.data = FALSE)
#> Error in resolve_dose_selection(sel_args, data): ℹ In argument: `all_of(dosevars)`.
#> Caused by error:
#> ! object 'dosevars' not found
res <- residuals(fit3, data = data)
#> Error: object 'fit3' not found
# }
## Schoenfeld residuals for proportional hazards fits
# \donttest{
fit4 <- ameras(Surv(time, status) ~ dose(all_of(dosevars), model = "ERR"),
              data = data, family = "prophaz", methods = "RC")
#> Error in resolve_dose_selection(sel_args, data): ℹ In argument: `all_of(dosevars)`.
#> Caused by error:
#> ! object 'dosevars' not found
res_sch <- residuals(fit4, type = "schoenfeld", scaled_schoenfeld = TRUE)
#> Error: object 'fit4' not found
# }
```
