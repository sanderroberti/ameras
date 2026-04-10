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

There are several options for confidence intervals that can be supplied
to the `CI` argument, see below. When `ameras` is called with `methods`
containing at least one of `RC`, `ERC`, and `MCML` and at least one of
`FMA` and `BMA`, `CI` should be a vector of length 2 with one method for
`RC`, `ERC` and `MCML` and one for `FMA` and `BMA`.

## Regression calibration, extended regression calibration, and Monte Carlo maximum likelihood

For (extended) regression calibration and Monte Carlo maximum
likelihood, there are two types of Wald intervals, obtained either
before or after transformation. If no transformation is specified,
`wald.orig` should be used to obtain the standard Wald intervals. When a
transformation is used, `wald.transformed` is determined before
transforming, and `wald.orig` is obtained after transforming using the
delta method (using `transform.jacobian` is required). The third option
is `proflik`, which uses the profile likelihood to compute confidence
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
computationally intensive to obtain. For this reason, `ameras` offers
the option to only determine them for the exposure-related parameters,
which is the default setting. To obtain profile likelihood intervals for
all parameters, use `params.profCI = "all"`.

To illustrate, we determine the three types of confidence intervals for
a regression calibration analysis using the example data, using a linear
excess relative risk model with the default exponential transformation
(see [Parameter
transformations](https://ameras.sanderroberti.com/articles/transformations.md)).

``` r
data(data, package="ameras")
dosevars <- paste0("V", 1:10)
fit.ameras.waldorig <- ameras(Y="Y.binomial", dosevars=dosevars, X=c("X1","X2"), data=data, 
                            family="binomial", methods=c("RC"), CI="wald.orig", doseRRmod="ERR")
#> Fitting RC
fit.ameras.waldtransformed <- ameras(Y="Y.binomial", dosevars=dosevars, X=c("X1","X2"), 
                                     data=data, family="binomial", methods=c("RC"), 
                                     CI="wald.transformed", doseRRmod="ERR")
#> Fitting RC
fit.ameras.proflik <- ameras(Y="Y.binomial", dosevars=dosevars, X=c("X1","X2"), data=data, 
                            family="binomial", methods=c("RC"), CI="proflik", doseRRmod="ERR", 
                            params.profCI = "all")
#> Fitting RC
#> Obtaining profile likelihood CI for (Intercept)
#> Obtaining profile likelihood CI for X1
#> Obtaining profile likelihood CI for X2
#> Obtaining profile likelihood CI for dose
summary(fit.ameras.waldorig)
#> Call:
#> ameras(data = data, family = "binomial", Y = "Y.binomial", dosevars = dosevars, 
#>     X = c("X1", "X2"), methods = c("RC"), doseRRmod = "ERR", 
#>     CI = "wald.orig")
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
#>  Method        Term Estimate      SE CI.lowerbound CI.upperbound
#>      RC (Intercept)  -1.0641 0.08788       -1.2363       -0.8918
#>      RC          X1   0.4409 0.07628        0.2914        0.5904
#>      RC          X2  -0.3360 0.09544       -0.5230       -0.1489
#>      RC        dose   0.8508 0.14517        0.5663        1.1353
summary(fit.ameras.waldtransformed)
#> Call:
#> ameras(data = data, family = "binomial", Y = "Y.binomial", dosevars = dosevars, 
#>     X = c("X1", "X2"), methods = c("RC"), doseRRmod = "ERR", 
#>     CI = "wald.transformed")
#> 
#> Total run time: 0.2 seconds
#> 
#> Runtime in seconds by method:
#> 
#>  Method Runtime
#>      RC     0.2
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
#> ameras(data = data, family = "binomial", Y = "Y.binomial", dosevars = dosevars, 
#>     X = c("X1", "X2"), methods = c("RC"), doseRRmod = "ERR", 
#>     CI = "proflik", params.profCI = "all")
#> 
#> Total run time: 3.5 seconds
#> 
#> Runtime in seconds by method:
#> 
#>  Method Runtime
#>      RC     3.5
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
fit.ameras.hpd <- ameras(Y="Y.binomial", dosevars=dosevars, X=c("X1","X2"), data=data, 
                            family="binomial", methods=c("FMA"), CI="hpd", doseRRmod="ERR")
#> Fitting FMA
fit.ameras.percentile <- ameras(Y="Y.binomial", dosevars=dosevars, X=c("X1","X2"), data=data, 
                            family="binomial", methods=c("FMA"), CI="percentile", 
                            doseRRmod="ERR")
#> Fitting FMA

summary(fit.ameras.hpd)
#> Call:
#> ameras(data = data, family = "binomial", Y = "Y.binomial", dosevars = dosevars, 
#>     X = c("X1", "X2"), methods = c("FMA"), doseRRmod = "ERR", 
#>     CI = "hpd")
#> 
#> Total run time: 1.7 seconds
#> 
#> Runtime in seconds by method:
#> 
#>  Method Runtime
#>     FMA     1.7
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
#> ameras(data = data, family = "binomial", Y = "Y.binomial", dosevars = dosevars, 
#>     X = c("X1", "X2"), methods = c("FMA"), doseRRmod = "ERR", 
#>     CI = "percentile")
#> 
#> Total run time: 1.8 seconds
#> 
#> Runtime in seconds by method:
#> 
#>  Method Runtime
#>     FMA     1.8
#> 
#> Summary of coefficients by method:
#> 
#>  Method        Term Estimate      SE CI.lowerbound CI.upperbound
#>     FMA (Intercept)  -1.0576 0.08738       -1.2293       -0.8862
#>     FMA          X1   0.4428 0.07607        0.2940        0.5925
#>     FMA          X2  -0.3377 0.09565       -0.5256       -0.1516
#>     FMA        dose   0.8444 0.14491        0.5611        1.1297
```
