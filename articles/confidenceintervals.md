# Confidence intervals

``` r

library(ameras)
#> Loading required package: nimble
#> nimble version 1.4.3 is loaded.
#> For more information on NIMBLE and a User Manual,
#> please visit https://R-nimble.org.
#> 
#> Attaching package: 'nimble'
#> The following object is masked from 'package:stats':
#> 
#>     simulate
#> The following object is masked from 'package:base':
#> 
#>     declare
```

## Introduction

To compute confidence intervals, first fit the model using `ameras`,
then use [`confint()`](https://rdrr.io/r/stats/confint.html) when the
interval bounds are needed as a standalone result, for example for
immediate comparison or export. Use
[`add_confint()`](https://ameras.sanderroberti.com/reference/confint.md)
when the fit will be retained for subsequent reporting: it stores the
bounds with each fitted method so that
[`summary()`](https://rdrr.io/r/base/summary.html) and
`coef(summary(fit))` include them, stores profile-likelihood diagnostics
when available, and records confidence-interval computation time in the
fit. Because
[`add_confint()`](https://ameras.sanderroberti.com/reference/confint.md)
returns an amended copy invisibly, assign its result to retain this
information. This is particularly useful for computationally expensive
profile likelihood intervals: subsequent summaries use the attached
bounds and diagnostics without rerunning the profile calculation. The
interval table is displayed by
[`summary()`](https://rdrr.io/r/base/summary.html), rather than
automatically by
[`add_confint()`](https://ameras.sanderroberti.com/reference/confint.md);
use `print(summary(fit), digits = 3)`, for example, to control the
printed precision. Several interval types are supported through the
`type` argument. When the fitted object contains at least one of `RC`,
`ERC`, and `MCML` and at least one of `FMA` and `BMA`, `type` should
contain one interval type for each group.

## Regression calibration, extended regression calibration, and Monte Carlo maximum likelihood

For (extended) regression calibration and Monte Carlo maximum
likelihood, Wald and profile likelihood intervals can be obtained. When
a parameter transformation \mathbf\theta = h(\mathbf\eta) is used,
`type="wald.transformed"` yields the CI at significance level \alpha of
h(\mathbf\eta \pm z\_{1-\alpha/2} \mathbf V) where z\_{1-\alpha/2} is
the 1-\alpha/2-quantile of the standard normal distribution and \mathbf
V is the vector of standard deviations estimated using the inverse
Hessian matrix (returned by `optim`), and `type="wald.orig"` uses the
delta method to obtain the CI h(\mathbf\eta)\pm z\_{1-\alpha/2} \mathbf
V\_\* where \mathbf V\_\* is the vector of standard deviations estimated
using J H^{-1} J^T with J the Jacobian of the transformation (obtained
with `transform.jacobian`) and H is the Hessian. When no transformation
is used, `type="wald.orig"` should be used. The third option is
`type="proflik"`, which uses the profile likelihood to compute
confidence bounds. For this, the profile log (partial) likelihood for
parameter \theta_p is defined as PL_p (\theta_p^\*) =
\max\_{\mathbf\theta: \theta_p = \theta_p^\*} \ell (\mathbf \theta),
where \ell is the log (partial) likelihood. Next, profile confidence
intervals (\theta_p^l, \theta_p^h) are obtained for parameter \theta_p
at significance level \alpha by solving -2 \\PL_p(\theta_p^\*) -
\ell(\hat{\mathbf{\theta}})\\=\chi^2\_{1,1-\alpha} using
[`uniroot()`](https://rdrr.io/r/stats/uniroot.html) after finding a
bracket, with \hat{\mathbf{\theta}} the maximum likelihood estimate.
Note that profile likelihoods are more computationally intensive to
obtain. For this reason, `confint` offers the option to only determine
them for the exposure-related parameters, which is the default setting.
To obtain profile likelihood intervals for all parameters, use
`parm = "all"`.

### Bracketing and numerical controls

Before solving for the bounds, both
[`confint()`](https://rdrr.io/r/stats/confint.html) and
[`add_confint()`](https://ameras.sanderroberti.com/reference/confint.md)
first search for a bracket separately for the lower and upper bounds. A
bracket is an interval within which the likelihood-ratio p-value crosses
the target `1 - level` (for example, `0.05` for a 95% interval).

The search starts from the fitted parameter value on the internal
optimizer scale. When the stored Hessian is usable, the initial
candidate is four approximate standard errors away in the direction of
the bound; otherwise, fallback starting points are used. If the target
is not reached, the search distance is doubled, with up to 12 doublings.
Each profile evaluation may require re-optimizing the nuisance
parameters, so bracketing itself can be costly, especially for ERC and
MCML.

Once a bracket has been found,
[`uniroot()`](https://rdrr.io/r/stats/uniroot.html) locates the bound
within it. Root finding defaults to `tol.profCI = 1e-6` and
`maxit.profCI = 100`. The tolerance applies to the internal optimizer
parameter scale, rather than the p-value scale. A larger tolerance can
reduce computation time when less precision is acceptable, but a bound
is still returned as `NA` with a warning if root finding warns or its
p-value differs from the target `1 - level` by more than `0.005`.
Warnings distinguish inaccurate roots from iteration limits and
bracketing failures; changing `tol.profCI` or `maxit.profCI` does not
affect the separate bracketing search or the nuisance-parameter
optimizer.

If no suitable bracket is found, the corresponding bound is returned as
`NA` with a warning showing the last attempted value and its p-value.
This can reflect a limited search range, parameter scaling, numerical
instability in the constrained optimizations, or the absence of a finite
bound. Check the original fit’s convergence diagnostics and parameter
scaling when investigating such a failure.

### Example

To illustrate, we determine the three types of confidence intervals for
a regression calibration analysis using the example data, using a linear
excess relative risk model with the default exponential transformation
(see [Parameter
transformations](https://ameras.sanderroberti.com/articles/transformations.md)).
For subgroup-specific intervals in models with effect modifiers, fit
subgroup-specific dose coefficients with `modifier = ~ 0 + M` or
`modifier = ~ M - 1`; see [Effect
modification](https://ameras.sanderroberti.com/articles/effectmodification.md).

``` r

dosevars <- paste0("V", 1:10)
fit.ameras <- ameras(Y.binomial~dose(V1:V10, model="ERR")+X1+X2, data = ameras_simulated,
                            family="binomial", methods=c("RC"))
#> Fitting RC

ci.ameras.waldorig <- confint(fit.ameras, type="wald.orig")
ci.ameras.waldtransformed <- confint(fit.ameras, type="wald.transformed")
ci.ameras.proflik <- confint(fit.ameras, type="proflik", parm="all")
#> Obtaining profile likelihood CI for (Intercept)
#> Obtaining profile likelihood CI for X1
#> Obtaining profile likelihood CI for X2
#> Obtaining profile likelihood CI for dose

ci.ameras.waldorig
#>                  2.5 %     97.5 %
#> (Intercept) -1.2363205 -0.8918389
#> X1           0.2913841  0.5904082
#> X2          -0.5230346 -0.1488986
#> dose         0.5662850  1.1353411
ci.ameras.waldtransformed
#>                  2.5 %     97.5 %
#> (Intercept) -1.2363205 -0.8918389
#> X1           0.2913841  0.5904082
#> X2          -0.5230346 -0.1488986
#> dose         0.6050480  1.1827031
ci.ameras.proflik
#>                  2.5 %     97.5 %
#> (Intercept) -1.2410617 -0.8960193
#> X1           0.2916105  0.5906751
#> X2          -0.5241009 -0.1498246
#> dose         0.6005363  1.1776549

# Attach intervals when they should be included in the fit summary
fit.ameras.waldorig <- add_confint(
  fit.ameras,
  type="wald.orig"
)
summary(fit.ameras.waldorig)
#> Call:
#> ameras(formula = Y.binomial ~ dose(V1:V10, model = "ERR") + X1 + 
#>     X2, data = ameras_simulated, family = "binomial", methods = c("RC"))
#> 
#> Rows: 3000
#> 
#> Total CPU runtime: 0.3 seconds
#> 
#> CPU runtime in seconds by method:
#> 
#>  Method   Fit  CI Total
#>      RC 0.325 0.0 0.325
#> 
#> Summary of coefficients by method:
#> 
#>  Method        Term Estimate      SE CI.lower CI.upper
#>      RC (Intercept)  -1.0641 0.08788  -1.2363  -0.8918
#>      RC          X1   0.4409 0.07628   0.2914   0.5904
#>      RC          X2  -0.3360 0.09544  -0.5230  -0.1489
#>      RC        dose   0.8508 0.14517   0.5663   1.1353
```

## Frequentist and Bayesian model averaging

For frequentist and Bayesian model averaging methods, the options are
`percentile` which uses equal-tailed quantiles, and `hpd` which computes
highest posterior density intervals using `HPDinterval` from the `coda`
package, using either the FMA samples or Bayesian posterior samples.

Again, we use the example data to illustrate.

``` r

set.seed(123)
fit.ameras2 <- ameras(Y.binomial~dose(V1:V10, model="ERR")+X1+X2, data = ameras_simulated,
                            family="binomial", methods=c("FMA"))
#> Fitting FMA

ci.ameras.hpd <- confint(fit.ameras2, type="hpd")
ci.ameras.percentile <- confint(fit.ameras2, type="percentile")

ci.ameras.hpd
#>                  2.5 %     97.5 %
#> (Intercept) -1.2283922 -0.8856075
#> X1           0.2939251  0.5936532
#> X2          -0.5260442 -0.1530053
#> dose         0.5639982  1.1317257
ci.ameras.percentile
#>                  2.5 %     97.5 %
#> (Intercept) -1.2289745 -0.8860001
#> X1           0.2937597  0.5934911
#> X2          -0.5250004 -0.1517645
#> dose         0.5615876  1.1297517
```
