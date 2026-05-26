# Residuals for an amerasfit Object

Computes residuals for a fitted `amerasfit` object.

## Usage

``` r
# S3 method for class 'amerasfit'
residuals(object,
          method = "RC",
          type = NULL,
          data = NULL,
          dose.col = NULL,
          scaled.schoenfeld = TRUE,
          ...)
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

  The type of residuals to compute. Defaults to `"pearson"` for all
  families except `"prophaz"`, for which the default is `"schoenfeld"`.
  For types other than `"schoenfeld"`, let \\Y_i\\ and \\\mu_i\\ denote
  the observed response and fitted mean response for individual \\i\\,
  respectively. Then the residuals are defined as follows:

  `"pearson"`

  :   Pearson residuals. For `"gaussian"`: \\(Y_i - \mu_i) / \sigma\\
      where \\\sigma\\ is the estimated residual standard deviation
      (note: [`glm`](https://rdrr.io/r/stats/glm.html) returns raw
      residuals for Gaussian). For `"binomial"` and `"clogit"`: \\(Y_i -
      \mu_i) / \sqrt{\mu_i(1 - \mu_i)}\\. For `"poisson"`: \\(Y_i -
      \mu_i) / \sqrt{\mu_i}\\. For `"multinomial"`: per-category Pearson
      residuals \\(Y\_{iz} - \mu\_{iz}) / \sqrt{\mu\_{iz}(1 -
      \mu\_{iz})}\\ where \\\mu\_{iz}\\ is the fitted probability of
      category \\z\\ for individual \\i\\.

  `"deviance"`

  :   Deviance residuals, defined as the signed square root of the
      individual contribution to the deviance. For `"gaussian"`: the raw
      residual \\Y_i - \mu_i\\. For `"binomial"` and `"clogit"`:
      \\\pm\sqrt{-2\log(\hat{p}\_i)}\\ where \\\hat{p}\_i = \mu_i\\ if
      \\Y_i = 1\\ and \\\hat{p}\_i = 1 - \mu_i\\ if \\Y_i = 0\\, with
      sign equal to \\\text{sign}(Y_i - \mu_i)\\. For `"poisson"`:
      \\\pm\sqrt{2(Y_i \log(Y_i/\mu_i) - (Y_i - \mu_i))}\\ for \\Y_i \>
      0\\ and \\\pm\sqrt{2\mu_i}\\ for \\Y_i = 0\\, with sign equal to
      \\\text{sign}(Y_i - \mu_i)\\. For `"multinomial"`: per-category
      deviance residuals as for `"binomial"`.

  `"response"`

  :   Raw residuals \\Y_i - \mu_i\\. For `"multinomial"`: the matrix
      \\\bm{Y} - \hat{\bm{P}}\\ where \\\bm{Y}\\ is the \\N \times Z\\
      indicator matrix of observed categories and \\\hat{\bm{P}}\\ is
      the \\N \times Z\\ matrix of fitted probabilities.

  `"schoenfeld"`

  :   Schoenfeld residuals for `family="prophaz"` only. For each event
      time, the unscaled residual for the individual experiencing the
      event is the observed covariate vector minus the risk-set weighted
      mean covariate vector at that event time. See Details.

- data:

  The original data frame used for fitting. Only required when the model
  was fitted with `keep.data=FALSE`.

- dose.col:

  Character string specifying the dose column to use when computing
  fitted values. If `NULL` (the default), the dose column is selected
  automatically: the mean dose across realizations for RC and ERC, the
  realization with the highest likelihood for MCML and FMA, and the most
  frequently selected realization for BMA. Can be set to any dose column
  present in the data to override the automatic selection.

- scaled.schoenfeld:

  Logical. If `TRUE` (the default), scaled Schoenfeld residuals are
  returned, obtained by multiplying the raw Schoenfeld residuals by the
  estimated coefficient variance-covariance matrix following Grambsch
  and Therneau (1994). If `FALSE`, raw Schoenfeld residuals are
  returned. Only used when `type="schoenfeld"`.

- ...:

  Additional arguments, currently unused.

## Value

For families `"gaussian"`, `"binomial"`, `"poisson"`, and `"clogit"`, a
numeric vector of length \\N\\ containing the residuals.

For `family="multinomial"`, a numeric matrix of dimension \\N \times Z\\
where \\Z\\ is the number of outcome categories. Column names correspond
to the factor levels. Note that when plotting via
[`plot`](https://ameras.sanderroberti.com/reference/plot.md), the
reference category is excluded since its residuals are a linear
combination of the other categories.

For `family="prophaz"` with `type="schoenfeld"`, a data frame with
columns `id` (row index of the event), `time` (event time), and one
additional column per model parameter containing the Schoenfeld
residuals. Only rows corresponding to events (`status=1`) are included.

## Details

Fitted values are computed using the dose column specified by
`dose.col`.

For `family="clogit"`, fitted values \\\mu_i\\ are the conditional
probabilities of being a case within each matched set, computed as \\R_i
/ \sum\_{j \in \text{set}(i)} R_j\\ where \\R_i\\ is the relative risk
for individual \\i\\ and the sum is over all individuals in the same
matched set.

For `family="prophaz"` and `type="schoenfeld"`, the Schoenfeld residual
for the individual failing at time \\t\\ is: \$\$\bm{r} = \bm{x} -
\bar{\bm{x}}(t)\$\$ where \\\bm{x}\\ is the observed covariate vector
for the failing individual and \$\$\bar{\bm{x}}(t) = \frac{\sum\_{j \in
\mathcal{R}(t)} \bm{x}\_j R_j} {\sum\_{j \in \mathcal{R}(t)} R_j}\$\$ is
the weighted mean covariate vector in the risk set \\\mathcal{R}(t)\\ at
time \\t\\, with weights equal to the relative risks \\R_j\\. Ties in
event times are handled using the Breslow approximation.

## References

Schoenfeld, D. (1982). Partial residuals for the proportional hazards
regression model. *Biometrika*, **69**(1), 239–241.
[doi:10.1093/biomet/69.1.239](https://doi.org/10.1093/biomet/69.1.239)

Grambsch, P. M. and Therneau, T. M. (1994). Proportional hazards tests
and diagnostics based on weighted residuals. *Biometrika*, **81**(3),
515–526.
[doi:10.1093/biomet/81.3.515](https://doi.org/10.1093/biomet/81.3.515)

## See also

[`plot`](https://ameras.sanderroberti.com/reference/plot.md) for
diagnostic plots using these residuals,
[`ameras`](https://ameras.sanderroberti.com/reference/ameras.md) for
model fitting, [`residuals`](https://rdrr.io/r/stats/residuals.html) for
the generic function.

## Examples

``` r
data("data", package="ameras")
dosevars <- paste0("V", 1:10)

fit <- ameras(Y.binomial ~ dose(all_of(dosevars), model="ERR"),
              data=data, family="binomial", methods="RC")
#> Fitting RC

## Pearson residuals (default)
res <- residuals(fit)
summary(res)
#>     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#> -1.45550 -0.80586 -0.67823  0.00103  1.10709  1.55575 

## Deviance residuals
res_dev <- residuals(fit, type="deviance")

## Response residuals
res_raw <- residuals(fit, type="response")

## Specific dose column
res_v1 <- residuals(fit, dose.col="V1")

## Multiple methods
# \donttest{
fit2 <- ameras(Y.binomial ~ dose(all_of(dosevars), model="ERR"),
               data=data, family="binomial", methods=c("RC", "ERC"))
#> Fitting RC
#> Fitting ERC
res_erc <- residuals(fit2, method="ERC")
# }

## With keep.data=FALSE
# \donttest{
fit3 <- ameras(Y.binomial ~ dose(all_of(dosevars), model="ERR"),
               data=data, family="binomial", methods="RC",
               keep.data=FALSE)
#> Fitting RC
res <- residuals(fit3, data=data)
# }

## Schoenfeld residuals for proportional hazards model
# \donttest{
fit4 <- ameras(Surv(time, status) ~ dose(all_of(dosevars), model="ERR"),
               data=data, family="prophaz", methods="RC")
#> Fitting RC
res_sch <- residuals(fit4)
res_raw_sch <- residuals(fit4, scaled.schoenfeld=FALSE)
# }
```
