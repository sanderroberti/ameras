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
likelihood, there are two types of Wald intervals, obtained either
before or after transformation. If no transformation is specified,
`wald.orig` should be used to obtain the standard Wald intervals. When a
transformation is used, `wald.transformed` is determined before
transforming, and `wald.orig` is obtained after transforming using the
delta method (using `transform.jacobian`). The third option is
`proflik`, which uses the profile likelihood to compute confidence
bounds. For this, the profile log (partial) likelihood for parameter
$\theta_{p}$ is defined as
$$PL_{p}\left( \theta_{p}^{*} \right) = \max\limits_{{\mathbf{θ}}:\theta_{p} = \theta_{p}^{*}}\ell({\mathbf{θ}}),$$
where $\ell$ is the log (partial) likelihood. Next, profile confidence
intervals $\left( \theta_{p}^{l},\theta_{p}^{h} \right)$ are obtained
for parameter $\theta_{p}$ at significance level $\alpha = 0.95$ by
solving
$- 2\{ PL_{p}\left( \theta_{p}^{*} \right) - \ell\left( \widehat{\mathbf{θ}} \right)\} = \chi_{1,1 - \alpha}^{2}$
using the bisection method, with $\widehat{\mathbf{θ}}$ the maximum
likelihood estimate. Note that profile likelihoods are more
computationally intensive to obtain. For this reason, `confint` offers
the option to only determine them for the exposure-related parameters,
which is the default setting. To obtain profile likelihood intervals for
all parameters, use `parm = "all"`.

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
fit.ameras.waldtransformed <- confint(fit.ameras, type="wald.transformed")
fit.ameras.proflik <- confint(fit.ameras, type="proflik", parm="all")
#> Obtaining profile likelihood CI for (Intercept)
#> Obtaining profile likelihood CI for X1
#> Obtaining profile likelihood CI for X2
#> Obtaining profile likelihood CI for dose

summary(fit.ameras.waldorig)
#> Call:
#> ameras(formula = Y.binomial ~ dose(V1:V10, model = "ERR") + X1 + 
#>     X2, data = data, family = "binomial", methods = c("RC"))
#> 
#> Total run time: 0.4 seconds
#> 
#> Runtime in seconds by method:
#> 
#>  Method Runtime
#>      RC     0.4
#> 
#> Summary of coefficients by method:
#> 
#>  Method        Term Estimate      SE CI.lowerbound CI.upperbound
#>      RC (Intercept)  -1.0641 0.08788       -1.2363       -0.8918
#>      RC          X1   0.4409 0.07628        0.2914        0.5904
#>      RC          X2  -0.3360 0.09544       -0.5230       -0.1489
#>      RC        dose   0.8508 0.14517        0.5663        1.1353
summary(fit.ameras.waldtransformed)
#> Call:
#> ameras(formula = Y.binomial ~ dose(V1:V10, model = "ERR") + X1 + 
#>     X2, data = data, family = "binomial", methods = c("RC"))
#> 
#> Total run time: 0.4 seconds
#> 
#> Runtime in seconds by method:
#> 
#>  Method Runtime
#>      RC     0.4
#> 
#> Summary of coefficients by method:
#> 
#>  Method        Term Estimate      SE CI.lowerbound CI.upperbound
#>      RC (Intercept)  -1.0641 0.08788       -1.2363       -0.8918
#>      RC          X1   0.4409 0.07628        0.2914        0.5904
#>      RC          X2  -0.3360 0.09544       -0.5230       -0.1489
#>      RC        dose   0.8508 0.14517        0.6050        1.1827
summary(fit.ameras.proflik)
#> Call:
#> ameras(formula = Y.binomial ~ dose(V1:V10, model = "ERR") + X1 + 
#>     X2, data = data, family = "binomial", methods = c("RC"))
#> 
#> Total run time: 0.4 seconds
#> 
#> Runtime in seconds by method:
#> 
#>  Method Runtime
#>      RC     0.4
#> 
#> Summary of coefficients by method:
#> 
#>  Method        Term Estimate      SE CI.lowerbound CI.upperbound
#>      RC (Intercept)  -1.0641 0.08788       -1.2391       -0.8971
#>      RC          X1   0.4409 0.07628        0.2908        0.5917
#>      RC          X2  -0.3360 0.09544       -0.5256       -0.1490
#>      RC        dose   0.8508 0.14517        0.6009        1.1784
```

## Frequentist and Bayesian model averaging

For frequentist and Bayesian model averaging methods, the options are
`percentile` which uses 2.5% and 97.5% percentiles, and `hpd` which
computes highest posterior density intervals using `HPDinterval` from
the `coda` package, using either the FMA samples or Bayesian posterior
samples.

Again, we use the example data to illustrate.

``` r
fit.ameras2 <- ameras(Y.binomial~dose(V1:V10, model="ERR")+X1+X2, data=data, 
                            family="binomial", methods=c("FMA"))
#> Fitting FMA

fit.ameras.hpd <- confint(fit.ameras2, type="hpd")
fit.ameras.percentile <- confint(fit.ameras2, type="percentile")

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
#>  Method        Term Estimate      SE CI.lowerbound CI.upperbound
#>     FMA (Intercept)  -1.0571 0.08777       -1.2290       -0.8842
#>     FMA          X1   0.4426 0.07638        0.2922        0.5916
#>     FMA          X2  -0.3380 0.09564       -0.5255       -0.1492
#>     FMA        dose   0.8442 0.14511        0.5604        1.1301
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
#>  Method        Term Estimate      SE CI.lowerbound CI.upperbound
#>     FMA (Intercept)  -1.0571 0.08777       -1.2290       -0.8843
#>     FMA          X1   0.4426 0.07638        0.2929        0.5925
#>     FMA          X2  -0.3380 0.09564       -0.5270       -0.1503
#>     FMA        dose   0.8442 0.14511        0.5600        1.1299
```
