# Compute or Attach Confidence Intervals for an amerasfit Object

Compute confidence or credible intervals for the parameters of a fitted
`amerasfit` object. [`confint()`](https://rdrr.io/r/stats/confint.html)
returns standalone interval bounds. `add_confint()` attaches bounds and
any available diagnostics to a copy of the fitted object for inclusion
in [`summary()`](https://rdrr.io/r/base/summary.html) output, and
records confidence-interval timing.

## Usage

``` r
# S3 method for class 'amerasfit'
confint(object, parm="dose", level=0.95,
        type=c("proflik","percentile"), maxit.profCI=100,
        tol.profCI=1e-6, data=NULL, ...)

add_confint(object, ...)

# S3 method for class 'amerasfit'
add_confint(object, parm="dose", level=0.95,
        type=c("proflik","percentile"), maxit.profCI=100,
        tol.profCI=1e-6, data=NULL, ...)
```

## Arguments

- object:

  A fitted model object of class `amerasfit`, as returned by
  [`ameras`](https://ameras.sanderroberti.com/reference/ameras.md).

- parm:

  Either `"dose"` to compute profile likelihood intervals for
  dose-related parameters only, `"all"` for all parameters, or a
  character vector of parameter names. Only used when
  `type = "proflik"`. Defaults to `"dose"` because profile likelihood
  intervals can be computationally intensive.

- level:

  Confidence level, as a single number between zero and one. Defaults to
  `0.95`.

- type:

  Interval type. For RC, ERC, and MCML, use `"wald.orig"`,
  `"wald.transformed"`, or `"proflik"`. For FMA and BMA, use
  `"percentile"` or `"hpd"`. When both groups of methods are present,
  supply one type for each group.

- maxit.profCI:

  Maximum number of root-finding iterations per profile likelihood
  bound. Defaults to `100`. This does not control the separate
  bracketing search or the nuisance-parameter optimizer.

- tol.profCI:

  Root-finding tolerance on the internal optimizer parameter scale.
  Defaults to `1e-6`. Larger values can reduce computation time,
  particularly for ERC and MCML, but may produce inaccurate bounds. A
  bound is returned as `NA` with a warning if root finding warns or its
  p-value differs from `1 - level` by more than `0.005`, even when a
  larger tolerance is explicitly requested.

- data:

  Original fitting data. Required for profile likelihood intervals when
  the model was fitted with `keep.data = FALSE`.

- ...:

  Additional arguments, currently unused.

## Value

- [`confint()`](https://rdrr.io/r/stats/confint.html):

  For an object containing one fitted method, a numeric matrix with
  parameters in rows and the lower and upper confidence limits in two
  percentage-labelled columns. For an object containing multiple fitted
  methods, a named list of such matrices, one per method. The result is
  returned visibly and contains only interval bounds.

- `add_confint()`:

  A copy of the `amerasfit` object, returned invisibly, with a `CI`
  element added to each fitted method. Method-specific timing is updated
  to include interval computation. Each `CI` element contains `lower`
  and `upper`. Profile likelihood results additionally contain boundary
  p-values and iteration counts in `pval.lower`, `pval.upper`,
  `iter.lower`, and `iter.upper`.

## Details

Both functions use the same interval calculations and always compute the
requested intervals, even if `object` already contains attached
intervals. The input object is left unchanged.

Use [`confint()`](https://rdrr.io/r/stats/confint.html) when only
standalone bounds are needed, for example for comparison or export. Use
`add_confint()` when the fit will be retained for reporting or further
inspection: [`summary()`](https://rdrr.io/r/base/summary.html) and
`coef(summary(object))` include the attached bounds, profile-likelihood
diagnostics are stored when available, and computation time is recorded.
Assign the invisibly returned object to retain this information. Any
previously attached intervals are replaced in that returned copy. Use
[`summary()`](https://rdrr.io/r/base/summary.html) to display the stored
intervals and `print(summary(fit), digits = 3)`, for example, to control
the printed precision. Attaching intervals is particularly useful for
expensive profile likelihood calculations because subsequent summaries
use the stored results without repeating the calculation.

## Interval types

For RC, ERC, and MCML, Wald intervals on the original parameter scale
use the delta-method variance-covariance matrix. Wald intervals on the
transformed scale are computed in the optimization parametrization and
then back-transformed. More specifically, when \\\boldsymbol\theta =
h(\boldsymbol\eta)\\, transformed-scale intervals apply \\h\\ to the
endpoints on the \\\boldsymbol\eta\\ scale, whereas original-scale
intervals use the covariance matrix obtained from \\J H^{-1} J^T\\,
where \\J\\ is the transformation Jacobian and \\H\\ is the optimization
Hessian. Transformed-scale intervals can therefore be asymmetric on the
reported parameter scale. Profile likelihood intervals repeatedly
re-evaluate the likelihood and can therefore be substantially more
expensive.

Transformations must be componentwise and strictly monotone, as required
by [`ameras`](https://ameras.sanderroberti.com/reference/ameras.md).
Decreasing transformations reverse endpoint order; bounds and profile
diagnostics are returned in reported lower/upper order. Transformed Wald
and profile likelihood intervals check the transformation at the fitted
parameters and determine its directions locally. Original-scale Wald
intervals and sample-based intervals do not require this check.

For FMA and BMA, percentile intervals use equal-tailed sample quantiles
and HPD intervals use
[`HPDinterval`](https://rdrr.io/pkg/coda/man/HPDinterval.html).

## Profile likelihood bracketing

Before determining bounds, a separate adaptive search first finds a
bracket for each lower and upper bound. A bracket is an interval within
which the likelihood-ratio p-value crosses the target `1 - level`. The
search proceeds outward from the fitted parameter value on the internal
optimizer scale. Initial search distances are four approximate standard
errors when the stored Hessian is usable; otherwise, fallback starting
points are used. If necessary, the search distance is doubled up to 12
times. Each profile evaluation may require re-optimizing the nuisance
parameters, so bracketing can itself be computationally expensive,
particularly for ERC and MCML.

Once a bracket is found,
[`uniroot`](https://rdrr.io/r/stats/uniroot.html) locates the bound
within it using `tol.profCI` and `maxit.profCI`. These arguments control
root finding only; they do not affect the bracketing search. If a bound
cannot be bracketed or reliably solved, it is returned as `NA` with a
warning reporting the last attempted parameter value and its p-value. A
bracketing failure can reflect a limited search range, parameter
scaling, numerical problems in the profile evaluations, or the absence
of a finite bound.

Nuisance optimizations use the effective fitting controls. Missing or
`NULL` controls on older objects default to `list(reltol = 1e-10)`. When
only one nuisance parameter remains,
[`optimize()`](https://rdrr.io/r/stats/optimize.html) uses its own
tolerance; `reltol` does not control that optimization.

When a model was fitted with `keep.data = FALSE`, the original data must
be supplied through `data` for profile likelihood intervals. The
missing-value policy used during fitting is reapplied before likelihood
evaluation. Wald and sample-based intervals do not require the original
data. For models with formula effect modifiers, subgroup-specific
profile intervals can be obtained by fitting subgroup-specific dose
coefficients with `modifier = ~ 0 + M` or `modifier = ~ M - 1`.

## See also

[`ameras`](https://ameras.sanderroberti.com/reference/ameras.md) for
model fitting,
[`summary`](https://ameras.sanderroberti.com/reference/summary.md) for
fit summaries, and [`confint`](https://rdrr.io/r/stats/confint.html) for
the generic function. See also
[`vignette("confidenceintervals", package = "ameras")`](https://ameras.sanderroberti.com/articles/confidenceintervals.md).

## Examples

``` r
fit <- ameras(
  Y.binomial ~ dose(V1:V3, model = "ERR"),
  data = ameras_simulated,
  family = "binomial",
  methods = "RC"
)
#> Fitting RC

## Compute standalone Wald bounds
confint(fit, type = "wald.orig")
#>                  2.5 %     97.5 %
#> (Intercept) -1.0297063 -0.7404073
#> dose         0.5338981  1.0748159

## Attach Wald intervals for use in summaries
fit_with_ci <- add_confint(fit, type = "wald.orig")
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
#>  Method   Fit  CI Total
#>      RC 0.078 0.0 0.078
#> 
#> Summary of coefficients by method:
#> 
#>  Method        Term Estimate     SE CI.lower CI.upper
#>      RC (Intercept)  -0.8851 0.0738  -1.0297  -0.7404
#>      RC        dose   0.8044 0.1380   0.5339   1.0748
#> 

## Attach expensive profile intervals for subsequent reporting

fit_with_profile_ci <- add_confint(
  fit,
  type = "proflik",
  parm = "dose"
)
#> Obtaining profile likelihood CI for dose
summary(fit_with_profile_ci)
#> Call:
#> ameras(formula = Y.binomial ~ dose(V1:V3, model = "ERR"), data = ameras_simulated, 
#>     family = "binomial", methods = "RC")
#> 
#> Rows: 3000
#> 
#> Total CPU runtime: 0.2 seconds
#> 
#> CPU runtime in seconds by method:
#> 
#>  Method   Fit    CI Total
#>      RC 0.078 0.076 0.154
#> 
#> Summary of coefficients by method:
#> 
#>  Method        Term Estimate     SE CI.lower CI.upper
#>      RC (Intercept)  -0.8851 0.0738       NA       NA
#>      RC        dose   0.8044 0.1380   0.5661    1.114
#> 
```
