# Confidence Intervals for an amerasfit Object

Computes confidence or credible intervals for the parameters of a fitted
`amerasfit` object and returns the interval bounds. The fitted object is
not amended. Use
[`add_confint`](https://ameras.sanderroberti.com/reference/add_confint.md)
to attach intervals and any available profile diagnostics to an
`amerasfit` object for inclusion in
[`summary()`](https://rdrr.io/r/base/summary.html) output and recording
of confidence-interval timing.

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

  Maximum number of root-finding iterations for profile likelihood
  intervals. Defaults to `20`.

- tol.profCI:

  Root-finding tolerance for profile likelihood intervals. Defaults to
  `1e-2`.

- data:

  Original fitting data. Required for profile likelihood intervals when
  the model was fitted with `keep.data = FALSE`.

- ...:

  Additional arguments, currently unused. To attach intervals to a
  fitted object or control printed output, use
  [`add_confint`](https://ameras.sanderroberti.com/reference/add_confint.md).

## Value

For an object containing one fitted method, a numeric matrix with
parameters in rows and the lower and upper confidence limits in two
percentage-labelled columns. For an object containing multiple fitted
methods, a named list of such matrices, one per method.

## Details

For RC, ERC, and MCML, Wald intervals on the original parameter scale
use the delta-method variance-covariance matrix. Wald intervals on the
transformed scale are computed in the optimization parametrization and
then back-transformed. More specifically, when \\\boldsymbol\theta =
h(\boldsymbol\eta)\\, transformed-scale intervals apply \\h\\ to the
endpoints on the \\\boldsymbol\eta\\ scale, whereas original-scale
intervals use the covariance matrix obtained from \\J H^{-1} J^T\\,
where \\J\\ is the transformation Jacobian and \\H\\ is the optimization
Hessian. Profile likelihood intervals repeatedly re-evaluate the
likelihood and can therefore be substantially more expensive.

For FMA and BMA, percentile intervals use equal-tailed sample quantiles
and HPD intervals use
[`HPDinterval`](https://rdrr.io/pkg/coda/man/HPDinterval.html).

Profile likelihood diagnostics such as boundary p-values and iteration
counts are not included in the returned matrices. Use
[`add_confint`](https://ameras.sanderroberti.com/reference/add_confint.md)
to store these diagnostics in the method-specific `CI` components of a
fitted object.

The requested intervals are always computed, even if `object` already
contains attached intervals. The input object is left unchanged. If a
profile likelihood bound cannot be reliably bracketed or solved, that
bound is returned as `NA` with a diagnostic warning.

When a model was fitted with `keep.data = FALSE`, the original data must
be supplied through `data` for profile likelihood intervals. The
missing-value policy used during fitting is reapplied before likelihood
evaluation. Wald and sample-based intervals do not require the original
data. For models with formula effect modifiers, subgroup-specific
profile intervals can be obtained by fitting subgroup-specific dose
coefficients with `modifier = ~ 0 + M` or `modifier = ~ M - 1`.

## See also

[`add_confint`](https://ameras.sanderroberti.com/reference/add_confint.md)
for attaching intervals and diagnostics to a fit,
[`ameras`](https://ameras.sanderroberti.com/reference/ameras.md) for
model fitting,
[`summary`](https://ameras.sanderroberti.com/reference/summary.md) for
fit summaries, and [`confint`](https://rdrr.io/r/stats/confint.html) for
the generic function.

## Examples

``` r
fit <- ameras(
  Y.binomial ~ dose(V1:V3, model = "ERR"),
  data = ameras_simulated,
  family = "binomial",
  methods = "RC"
)
#> Fitting RC

## A two-column matrix of Wald interval bounds
confint(fit, type = "wald.orig")
#>                  2.5 %     97.5 %
#> (Intercept) -1.0297063 -0.7404073
#> dose         0.5338981  1.0748159

## Profile likelihood intervals are slower

confint(fit, type = "proflik", parm = "dose")
#> Obtaining profile likelihood CI for dose
#>          2.5 %   97.5 %
#> dose 0.5663534 1.115075
```
