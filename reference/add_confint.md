# Attach Confidence Intervals to an amerasfit Object

Computes confidence or credible intervals, attaches them and any
available diagnostics to a fitted `amerasfit` object, and optionally
prints the bounds. Use
[`confint`](https://ameras.sanderroberti.com/reference/confint.md) when
only the interval bounds are required.

## Usage

``` r
add_confint(object, ...)

# S3 method for class 'amerasfit'
add_confint(object, parm="dose", level=0.95,
        type=c("proflik","percentile"), maxit.profCI=20,
        tol.profCI=1e-2, data=NULL, print=TRUE,
        digits=max(3, getOption("digits") - 3), ...)
```

## Arguments

- object:

  A fitted model object of class `amerasfit`, as returned by
  [`ameras`](https://ameras.sanderroberti.com/reference/ameras.md).

- parm:

  Either `"dose"` to compute profile likelihood intervals for
  dose-related parameters only, `"all"` for all parameters, or a
  character vector of parameter names. Only used when
  `type = "proflik"`.

- level:

  Confidence level, as a single number between zero and one. Defaults to
  `0.95`.

- type:

  Interval type. For RC, ERC, and MCML, use `"wald.orig"`,
  `"wald.transformed"`, or `"proflik"`. For FMA and BMA, use
  `"percentile"` or `"hpd"`. When both groups of methods are present,
  supply one type for each group.

- maxit.profCI:

  Maximum number of root-finding iterations for profile likelihood
  intervals. Defaults to `20`.

- tol.profCI:

  Root-finding tolerance for profile likelihood intervals. Defaults to
  `1e-2`.

- data:

  Original fitting data. Required for profile likelihood intervals when
  the model was fitted with `keep.data = FALSE`.

- print:

  Whether to print the attached interval bounds. Defaults to `TRUE`.

- digits:

  Number of significant digits used when printing.

- ...:

  Additional arguments, currently unused.

## Value

A copy of the `amerasfit` object, invisibly, with a `CI` element added
to each fitted method and `CI.computed` set to `TRUE`. Method-specific
timing is updated to include interval computation. The input object is
not modified.

Each `CI` element contains `lower` and `upper`. Profile likelihood
results additionally contain boundary p-values and iteration counts in
`pval.lower`, `pval.upper`, `iter.lower`, and `iter.upper`.

## Details

The requested intervals are always computed. Any intervals already
attached to `object` are replaced in the returned object; the input
object is left unchanged. This function is useful when the fitted object
will be kept for reporting or further inspection:
[`summary()`](https://rdrr.io/r/base/summary.html) and
`coef(summary(object))` include the attached bounds, profile-likelihood
diagnostics are stored when available, and confidence-interval
computation time is recorded in the method-specific timing information.
Use [`confint`](https://ameras.sanderroberti.com/reference/confint.md)
when only standalone interval bounds are needed. Attaching intervals is
particularly useful for computationally expensive profile likelihood
intervals because subsequent summaries use the attached bounds and
diagnostics without repeating the profile calculation. Computation and
failure handling are otherwise identical.

For RC, ERC, and MCML, Wald intervals on the original parameter scale
use the delta-method variance-covariance matrix. Wald intervals on the
transformed scale are computed in the optimization parametrization and
then back-transformed. More specifically, when \\\boldsymbol\theta =
h(\boldsymbol\eta)\\, transformed-scale intervals apply \\h\\ to the
endpoints on the \\\boldsymbol\eta\\ scale, whereas original-scale
intervals use the covariance matrix obtained from \\J H^{-1} J^T\\,
where \\J\\ is the transformation Jacobian and \\H\\ is the optimization
Hessian. Transformed-scale intervals can therefore be asymmetric on the
reported parameter scale.

For FMA and BMA, percentile intervals use equal-tailed sample quantiles
and HPD intervals use
[`HPDinterval`](https://rdrr.io/pkg/coda/man/HPDinterval.html).

## See also

[`confint`](https://ameras.sanderroberti.com/reference/confint.md) for
returning interval bounds,
[`ameras`](https://ameras.sanderroberti.com/reference/ameras.md) for
model fitting, and
[`summary`](https://ameras.sanderroberti.com/reference/summary.md) for
summaries containing attached intervals.

## Examples

``` r
fit <- ameras(
  Y.binomial ~ dose(V1:V3, model = "ERR"),
  data = ameras_simulated,
  family = "binomial",
  methods = "RC"
)
#> Fitting RC

fit_with_ci <- add_confint(fit, type = "wald.orig", print = FALSE)
summary(fit_with_ci)
#> Call:
#> ameras(formula = Y.binomial ~ dose(V1:V3, model = "ERR"), data = ameras_simulated, 
#>     family = "binomial", methods = "RC")
#> 
#> Rows: 3000
#> 
#> Total CPU runtime: 0.1 seconds
#> 
#> CPU runtime in seconds by method:
#> 
#>  Method   Fit    CI Total
#>      RC 0.051 0.001 0.052
#> 
#> Summary of coefficients by method:
#> 
#>  Method        Term Estimate     SE CI.lower CI.upper
#>      RC (Intercept)  -0.8851 0.0738  -1.0297  -0.7404
#>      RC        dose   0.8044 0.1380   0.5339   1.0748
#> 

## Profile likelihood intervals are slower

fit_with_profile_ci <- add_confint(
  fit,
  type = "proflik",
  parm = "dose",
  print = FALSE
)
#> Obtaining profile likelihood CI for dose
summary(fit_with_profile_ci)
#> Call:
#> ameras(formula = Y.binomial ~ dose(V1:V3, model = "ERR"), data = ameras_simulated, 
#>     family = "binomial", methods = "RC")
#> 
#> Rows: 3000
#> 
#> Total CPU runtime: 0.1 seconds
#> 
#> CPU runtime in seconds by method:
#> 
#>  Method   Fit    CI Total
#>      RC 0.051 0.065 0.116
#> 
#> Summary of coefficients by method:
#> 
#>  Method        Term Estimate     SE CI.lower CI.upper
#>      RC (Intercept)  -0.8851 0.0738       NA       NA
#>      RC        dose   0.8044 0.1380   0.5664    1.115
#> 
```
