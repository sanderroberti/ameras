# Confidence intervals for an amerasfit object

Computes confidence intervals for the parameters of a fitted `amerasfit`
object. This is a separate step from model fitting, i.e.,
[`ameras`](https://ameras.sanderroberti.com/reference/ameras.md) fits
the model and `confint` computes intervals and attaches them to the
fitted object.

## Usage

``` r
# S3 method for class 'amerasfit'
confint(object, parm="dose", level=0.95, 
        type=c("proflik","percentile"), maxit.profCI=20, 
        tol.profCI=1e-2, data=NULL, ...)
```

## Arguments

- object:

  A fitted model object of class `amerasfit`, as returned by
  [`ameras`](https://ameras.sanderroberti.com/reference/ameras.md).

- parm:

  Either `"dose"` to compute intervals for dose-related parameters only,
  `"all"` for all parameters, or a character vector of specific
  parameter names. Only used when `type = "proflik"` since Wald
  intervals are cheap to compute for all parameters simultaneously.
  Defaults to `"dose"` since profile likelihood computation can be
  extensive.

- level:

  The confidence level (default `0.95`).

- type:

  The type(s) of confidence intervals to determine. For RC, ERC, and
  MCML, this can be one of:

  `"wald.orig"`

  :   Wald intervals on the original parameter scale using the delta
      method variance-covariance matrix.

  `"wald.transformed"`

  :   Wald intervals computed on the transformed (reparametrized) scale
      and then back-transformed. Only available when a transformation
      was used during fitting.

  `"proflik"`

  :   Profile likelihood intervals based on the chi-squared
      approximation. More accurate than Wald intervals but
      computationally intensive, especially for large datasets or
      complex models.

  For FMA and BMA, confidence intervals are based on the generated
  samples and possible confidence interval types are:

  `"percentile"`

  :   Equal-tailed percentile intervals.

  `"hpd"`

  :   Highest posterior density intervals via
      [`HPDinterval`](https://rdrr.io/pkg/coda/man/HPDinterval.html)
      from the coda package.

  If `object` contains results for at least one of RC, ERC, and MCML and
  at least one of FMA and BMA, `type` must be length 2 and specify one
  method for RC, ERC, and MCML, and one for FMA and BMA.

- data:

  The original data frame used for fitting. Only required when
  `type = "proflik"` and the model was fitted with `keep.data = FALSE`

- maxit.profCI:

  Maximum number of iterations for the root-finding algorithm used to
  locate profile likelihood interval bounds. Only used when
  `type = "proflik"`. Defaults to `20`.

- tol.profCI:

  Tolerance for the root-finding algorithm. Only used when
  `type = "proflik"`. Defaults to `1e-2`. Reduce for more precise bounds
  at the cost of additional computation.

- ...:

  Additional arguments, currently unused.

## Value

The original `amerasfit` object with a `CI` element added to each fitted
method result. For RC, ERC, and MCML the `CI` element is a data frame
with columns:

- `lower`:

  Lower confidence bound.

- `upper`:

  Upper confidence bound.

When `type = "proflik"`, four additional columns are included:

- `pval.lower`:

  P-value at the lower bound, should be close to \\1 -\\ `level`.

- `pval.upper`:

  P-value at the upper bound, should be close to \\1 -\\ `level`.

- `iter.lower`:

  Number of iterations used by the root-finding algorithm for the lower
  bound.

- `iter.upper`:

  Number of iterations used by the root-finding algorithm for the upper
  bound.

For FMA and BMA the `CI` element is a data frame with columns `lower`
and `upper`.

## Details

For (extended) regression calibration and Monte Carlo maximum
likelihood, Wald and profile likelihood intervals can be obtained. When
a parameter transformation \\\bm\theta = h(\bm\eta)\\ is used,
`type="wald.transformed"` yields the CI at significance level \\\alpha\\
of \\h(\bm\eta \pm z\_{1-\alpha/2} \bm V)\\ where \\z\_{1-\alpha/2}\\ is
the \\1-\alpha/2\\-quantile of the standard normal distribution and
\\\bm V\\ is the vector of standard deviations estimated using the
inverse Hessian matrix, and `type="wald.orig"` uses the delta method to
obtain the CI \\h(\bm\eta)\pm z\_{1-\alpha/2} \bm V\_\*\\ where \\\bm
V\_\*\\ is the vector of standard deviations estimated using \\J H^{-1}
J^T\\ with \\J\\ the Jacobian of the transformation and \\H\\ is the
Hessian. When no transformation is used, `type="wald.orig"` should be
used. The third option is `proflik`, which uses the profile likelihood
to compute confidence bounds.

For FMA and BMA, the options for confidence/credible intervals are
`type="percentile"` which uses percentiles, and `type="hpd"` which
computes highest posterior density intervals using `HPDinterval` from
the `coda` package, both using the FMA samples or Bayesian posterior
samples.

Profile likelihood intervals (`type="proflik"`) require re-evaluating
the likelihood repeatedly and can be time-consuming. The `parm` argument
can be used to restrict computation to dose parameters only (the
default) when intervals for the other parameters are not of interest.

When the model was fitted with `keep.data=FALSE` and `type="proflik"` is
used for `confint`, the original data must be supplied via the `data`
argument. Wald intervals do not require the data and can always be
computed from the stored Hessian and parameter estimates alone.

## See also

[`ameras`](https://ameras.sanderroberti.com/reference/ameras.md) for
model fitting,
[`summary`](https://ameras.sanderroberti.com/reference/summary.md) for a
summary of the fitted model including confidence intervals if computed,
[`confint`](https://rdrr.io/r/stats/confint.html) for the generic
function.

## Examples

``` r
data("data", package = "ameras")

## Fit the model
fit <- ameras(Y.binomial~dose(V1:V10, model="ERR"), data = data, family = "binomial",
               methods = "RC")
#> Fitting RC

## Wald intervals (fast)
fit <- confint(fit, type = "wald.orig")
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
#>  Method        Term Estimate      SE CI.lowerbound CI.upperbound
#>      RC (Intercept)  -0.8847 0.07378       -1.0293       -0.7401
#>      RC        dose   0.8020 0.13751        0.5324        1.0715

## Profile likelihood intervals for dose parameters only (slower)
# \donttest{
fit <- confint(fit, type = "proflik", parm = "dose")
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
# }

## With keep.data = FALSE, supply data explicitly for proflik
# \donttest{
fit2 <- ameras(Y.binomial~dose(V1:V10, model="ERR"), data = data, family = "binomial",
               methods = "RC", keep.data = FALSE)
#> Fitting RC
fit2 <- confint(fit2, type = "proflik", data = data)
#> Obtaining profile likelihood CI for dose
# }

## FMA and BMA with percentile intervals
# \donttest{
fit3 <- ameras(Y.binomial~dose(V1:V10, model="ERR"), data = data, family = "binomial", 
               methods = c("FMA", "BMA"))
#> Note: BMA may require extensive computation time
#> Fitting FMA
#> Fitting BMA
#> Defining model
#> Building model
#> Setting data and initial values
#> Running calculate on model
#>   [Note] Any error reports that follow may simply reflect missing values in model variables.
#> Checking model sizes and dimensions
#>   [Note] This model is not fully initialized. This is not an error.
#>          To see which variables are not initialized, use model$initializeInfo().
#>          For more information on model initialization, see help(modelInitialization).
#> Compiling
#>   [Note] This may take a minute.
#>   [Note] Use 'showCompilerOutput = TRUE' to see C++ compilation details.
#> Compiling
#>   [Note] This may take a minute.
#>   [Note] Use 'showCompilerOutput = TRUE' to see C++ compilation details.
#> running chain 1...
#> |-------------|-------------|-------------|-------------|
#> |-------------------------------------------------------|
#> running chain 2...
#> |-------------|-------------|-------------|-------------|
#> |-------------------------------------------------------|
fit3 <- confint(fit3, type = "percentile")
summary(fit3)
#> Call:
#> ameras(formula = Y.binomial ~ dose(V1:V10, model = "ERR"), data = data, 
#>     family = "binomial", methods = c("FMA", "BMA"))
#> 
#> Total run time: 103.4 seconds
#> 
#> Runtime in seconds by method:
#> 
#>  Method Runtime
#>     FMA     0.4
#>     BMA   103.0
#> 
#> Summary of coefficients by method:
#> 
#>  Method        Term Estimate      SE CI.lowerbound CI.upperbound Rhat   n.eff
#>     FMA (Intercept)  -0.8762 0.07336       -1.0194       -0.7315   NA      NA
#>     FMA        dose   0.7916 0.13661        0.5241        1.0595   NA      NA
#>     BMA (Intercept)  -0.8733 0.07569       -1.0266       -0.7309 1.00 1024.00
#>     BMA        dose   0.7923 0.14167        0.5401        1.0977 1.00 1001.00
# }
```
