# Confidence intervals

``` r

library(ameras)
#> Loading required package: nimble
#> nimble version 1.4.2 is loaded.
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
then use the `confint` method to attach confidence intervals. Several
types of confidence intervals are supported, which should be supplied to
the `type` argument of `confint`, see below. When `confint` is called
with `methods` containing at least one of `RC`, `ERC`, and `MCML` and at
least one of `FMA` and `BMA`, `type` should be a vector of length 2 with
one method for `RC`, `ERC` and `MCML` and one for `FMA` and `BMA`.

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
\ell(\hat{\mathbf{\theta}})\\=\chi^2\_{1,1-\alpha} using the bisection
method, with \hat{\mathbf{\theta}} the maximum likelihood estimate. Note
that profile likelihoods are more computationally intensive to obtain.
For this reason, `confint` offers the option to only determine them for
the exposure-related parameters, which is the default setting. To obtain
profile likelihood intervals for all parameters, use `parm = "all"`.

To illustrate, we determine the three types of confidence intervals for
a regression calibration analysis using the example data, using a linear
excess relative risk model with the default exponential transformation
(see [Parameter
transformations](https://ameras.sanderroberti.com/articles/transformations.md)).

``` r

data(data, package="ameras")
dosevars <- paste0("V", 1:10)
fit.ameras <- ameras(Y.binomial~dose(V1:V10, model="ERR")+X1+X2, data=data, 
                            family="binomial", methods=c("RC"))
#> Fitting RC

fit.ameras.waldorig <- confint(fit.ameras, type="wald.orig")
#> RC confidence intervals:
#> 
#>               lower   upper
#> (Intercept) -1.2363 -0.8918
#> X1           0.2914  0.5904
#> X2          -0.5230 -0.1489
#> dose         0.5663  1.1353
fit.ameras.waldtransformed <- confint(fit.ameras, type="wald.transformed")
#> RC confidence intervals:
#> 
#>               lower   upper
#> (Intercept) -1.2363 -0.8918
#> X1           0.2914  0.5904
#> X2          -0.5230 -0.1489
#> dose         0.6050  1.1827
fit.ameras.proflik <- confint(fit.ameras, type="proflik", parm="all")
#> Obtaining profile likelihood CI for (Intercept)
#> Obtaining profile likelihood CI for X1
#> Obtaining profile likelihood CI for X2
#> Obtaining profile likelihood CI for dose
#> RC confidence intervals:
#> 
#>               lower   upper
#> (Intercept) -1.2391 -0.8971
#> X1           0.2908  0.5917
#> X2          -0.5256 -0.1490
#> dose         0.6009  1.1784

summary(fit.ameras.waldorig)
#> Call:
#> ameras(formula = Y.binomial ~ dose(V1:V10, model = "ERR") + X1 + 
#>     X2, data = data, family = "binomial", methods = c("RC"))
#> 
#> Total run time: 0.5 seconds
#> 
#> Runtime in seconds by method:
#> 
#>  Method Runtime
#>      RC     0.5
#> 
#> Summary of coefficients by method:
#> 
#>  Method        Term Estimate      SE CI.lower CI.upper
#>      RC (Intercept)  -1.0641 0.08788  -1.2363  -0.8918
#>      RC          X1   0.4409 0.07628   0.2914   0.5904
#>      RC          X2  -0.3360 0.09544  -0.5230  -0.1489
#>      RC        dose   0.8508 0.14517   0.5663   1.1353
summary(fit.ameras.waldtransformed)
#> Call:
#> ameras(formula = Y.binomial ~ dose(V1:V10, model = "ERR") + X1 + 
#>     X2, data = data, family = "binomial", methods = c("RC"))
#> 
#> Total run time: 0.5 seconds
#> 
#> Runtime in seconds by method:
#> 
#>  Method Runtime
#>      RC     0.5
#> 
#> Summary of coefficients by method:
#> 
#>  Method        Term Estimate      SE CI.lower CI.upper
#>      RC (Intercept)  -1.0641 0.08788  -1.2363  -0.8918
#>      RC          X1   0.4409 0.07628   0.2914   0.5904
#>      RC          X2  -0.3360 0.09544  -0.5230  -0.1489
#>      RC        dose   0.8508 0.14517   0.6050   1.1827
summary(fit.ameras.proflik)
#> Call:
#> ameras(formula = Y.binomial ~ dose(V1:V10, model = "ERR") + X1 + 
#>     X2, data = data, family = "binomial", methods = c("RC"))
#> 
#> Total run time: 0.5 seconds
#> 
#> Runtime in seconds by method:
#> 
#>  Method Runtime
#>      RC     0.5
#> 
#> Summary of coefficients by method:
#> 
#>  Method        Term Estimate      SE CI.lower CI.upper
#>      RC (Intercept)  -1.0641 0.08788  -1.2391  -0.8971
#>      RC          X1   0.4409 0.07628   0.2908   0.5917
#>      RC          X2  -0.3360 0.09544  -0.5256  -0.1490
#>      RC        dose   0.8508 0.14517   0.6009   1.1784
```

## Frequentist and Bayesian model averaging

For frequentist and Bayesian model averaging methods, the options are
`percentile` which uses equal-tailed quantiles, and `hpd` which computes
highest posterior density intervals using `HPDinterval` from the `coda`
package, using either the FMA samples or Bayesian posterior samples.

Again, we use the example data to illustrate.

``` r

set.seed(123)
fit.ameras2 <- ameras(Y.binomial~dose(V1:V10, model="ERR")+X1+X2, data=data, 
                            family="binomial", methods=c("FMA"))
#> Fitting FMA

fit.ameras.hpd <- confint(fit.ameras2, type="hpd")
#> FMA confidence intervals:
#> 
#>               lower   upper
#> (Intercept) -1.2284 -0.8856
#> X1           0.2939  0.5937
#> X2          -0.5260 -0.1530
#> dose         0.5640  1.1317
fit.ameras.percentile <- confint(fit.ameras2, type="percentile")
#> FMA confidence intervals:
#> 
#>               lower   upper
#> (Intercept) -1.2290 -0.8860
#> X1           0.2938  0.5935
#> X2          -0.5250 -0.1518
#> dose         0.5616  1.1298

summary(fit.ameras.hpd)
#> Call:
#> ameras(formula = Y.binomial ~ dose(V1:V10, model = "ERR") + X1 + 
#>     X2, data = data, family = "binomial", methods = c("FMA"))
#> 
#> Total run time: 1.6 seconds
#> 
#> Runtime in seconds by method:
#> 
#>  Method Runtime
#>     FMA     1.6
#> 
#> Summary of coefficients by method:
#> 
#>  Method        Term Estimate      SE CI.lower CI.upper
#>     FMA (Intercept)  -1.0576 0.08770  -1.2284  -0.8856
#>     FMA          X1   0.4429 0.07652   0.2939   0.5937
#>     FMA          X2  -0.3378 0.09572  -0.5260  -0.1530
#>     FMA        dose   0.8447 0.14494   0.5640   1.1317
summary(fit.ameras.percentile)
#> Call:
#> ameras(formula = Y.binomial ~ dose(V1:V10, model = "ERR") + X1 + 
#>     X2, data = data, family = "binomial", methods = c("FMA"))
#> 
#> Total run time: 1.6 seconds
#> 
#> Runtime in seconds by method:
#> 
#>  Method Runtime
#>     FMA     1.6
#> 
#> Summary of coefficients by method:
#> 
#>  Method        Term Estimate      SE CI.lower CI.upper
#>     FMA (Intercept)  -1.0576 0.08770  -1.2290  -0.8860
#>     FMA          X1   0.4429 0.07652   0.2938   0.5935
#>     FMA          X2  -0.3378 0.09572  -0.5250  -0.1518
#>     FMA        dose   0.8447 0.14494   0.5616   1.1298
```
