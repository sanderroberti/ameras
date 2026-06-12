# Confidence Intervals for an amerasfit Object

Computes and/or prints confidence intervals for the parameters of a
fitted `amerasfit` object. This is a separate step from model fitting,
i.e., [`ameras`](https://ameras.sanderroberti.com/reference/ameras.md)
fits the model and `confint` computes intervals, attaches them to the
fitted object, and prints them. If `confint` is called a second time on
the same object, intervals are only printed and not recomputed, unless
`force=TRUE`.

## Usage

``` r
# S3 method for class 'amerasfit'
confint(object, parm="dose", level=0.95, 
        type=c("proflik","percentile"), maxit.profCI=20, 
        tol.profCI=1e-2, data=NULL, force=FALSE, print=TRUE,
        digits = max(3, getOption("digits") - 3), ...)
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

- maxit.profCI:

  Maximum number of iterations for the root-finding algorithm used to
  locate profile likelihood interval bounds. Only used when
  `type = "proflik"`. Defaults to `20`.

- tol.profCI:

  Tolerance for the root-finding algorithm. Only used when
  `type = "proflik"`. Defaults to `1e-2`. Reduce for more precise bounds
  at the cost of additional computation.

- data:

  The original data frame used for fitting. Only required when
  `type = "proflik"` and the model was fitted with `keep.data = FALSE`

- force:

  Logical. If `TRUE`, confidence intervals are recomputed even if they
  have already been computed for this object. Defaults to `FALSE`, in
  which case a warning is issued and the object is returned unchanged if
  confidence intervals are already present.

- print:

  Logical. If `TRUE` (default), confidence intervals are printed to the
  console.

- digits:

  Number of significant digits to be printed. Default is \`max(3,
  getOption("digits") - 3)\`

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
#> RC confidence intervals:
#> 
#>               lower   upper
#> (Intercept) -1.0293 -0.7401
#> dose         0.5324  1.0715
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
#>  Method        Term Estimate      SE CI.lower CI.upper
#>      RC (Intercept)  -0.8847 0.07378  -1.0293  -0.7401
#>      RC        dose   0.8020 0.13751   0.5324   1.0715
#> 

## Profile likelihood intervals for dose parameters only (slower)
# \donttest{
fit <- confint(fit, type = "proflik", parm = "dose")
#> Confidence intervals have already been computed for this object. Returning the object unchanged. Use force=TRUE to recompute.
#> RC confidence intervals:
#> 
#>               lower   upper
#> (Intercept) -1.0293 -0.7401
#> dose         0.5324  1.0715
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
#>  Method        Term Estimate      SE CI.lower CI.upper
#>      RC (Intercept)  -0.8847 0.07378  -1.0293  -0.7401
#>      RC        dose   0.8020 0.13751   0.5324   1.0715
#> 
# }

## With keep.data = FALSE, supply data explicitly for proflik
# \donttest{
fit2 <- ameras(Y.binomial~dose(V1:V10, model="ERR"), data = data, family = "binomial",
               methods = "RC", keep.data = FALSE)
#> Fitting RC
fit2 <- confint(fit2, type = "proflik", data = data)
#> Obtaining profile likelihood CI for dose
#> RC confidence intervals:
#> 
#>       lower upper
#> dose 0.5648 1.112
#> 
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
#> FMA confidence intervals:
#> 
#>               lower   upper
#> (Intercept) -1.0203 -0.7321
#> dose         0.5243  1.0583
#> 
#> BMA confidence intervals:
#> 
#>               lower   upper
#> (Intercept) -1.0083 -0.7337
#> dose         0.5554  1.0806
#> 
summary(fit3)
#> Call:
#> ameras(formula = Y.binomial ~ dose(V1:V10, model = "ERR"), data = data, 
#>     family = "binomial", methods = c("FMA", "BMA"))
#> 
#> Total run time: 95.2 seconds
#> 
#> Runtime in seconds by method:
#> 
#>  Method Runtime
#>     FMA     0.5
#>     BMA    94.7
#> 
#> Summary of coefficients by method:
#> 
#>  Method        Term Estimate      SE CI.lower CI.upper Rhat   n.eff
#>     FMA (Intercept)  -0.8759 0.07349  -1.0203  -0.7321   NA      NA
#>     FMA        dose   0.7913 0.13664   0.5243   1.0583   NA      NA
#>     BMA (Intercept)  -0.8780 0.07121  -1.0083  -0.7337 1.00 1073.00
#>     BMA        dose   0.8003 0.13344   0.5554   1.0806 1.00 1068.00
#> 
# }
```
