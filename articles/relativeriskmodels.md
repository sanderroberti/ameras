# Relative risk models

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
library(ggplot2)
data(data, package="ameras")
dosevars <- paste0("V", 1:10)
```

## Introduction

For non-Gaussian families, three relative risk models for the main
exposure are supported, the usual exponential model
$$RR_{i} = \exp\left( \beta_{1}D_{i} + \beta_{2}D_{i}^{2} + \mathbf{M}_{i}^{T}{\mathbf{β}}_{m1}D_{i} + \mathbf{M}_{i}^{T}{\mathbf{β}}_{m2}D_{i}^{2} \right),$$
the linear(-quadratic) excess relative risk (ERR) model
$$RR_{i} = 1 + \beta_{1}D_{i} + \beta_{2}D_{i}^{2} + \mathbf{M}_{i}^{T}{\mathbf{β}}_{\mathbf{m}\mathbf{1}}D_{i} + \mathbf{M}_{i}^{T}{\mathbf{β}}_{m2}D_{i}^{2},$$
and the linear-exponential model
$$RR_{i} = 1 + \left( \beta_{1} + \mathbf{M}_{i}^{T}{\mathbf{β}}_{m1} \right)D_{i}\exp\{\left( \beta_{2} + \mathbf{M}_{i}^{T}{\mathbf{β}}_{m2} \right)D_{i}\}.$$
This vignette illustrates fitting the three models using regression
calibration for logistic regression, but the same syntax applies to all
other settings.

## Exponential relative risk

The usual exponential relative risk model is given by
$RR_{i} = \exp\left( \beta_{1}D_{i} + \beta_{2}D_{i}^{2} + \mathbf{M}_{i}^{T}{\mathbf{β}}_{m1}D_{i} + \mathbf{M}_{i}^{T}{\mathbf{β}}_{m2}D_{i}^{2} \right)$,
where the quadratic and effect modification terms are optional (not fit
by setting `deg=1` and not passing anything to `M`, respectively). This
model is fit by setting `doseRRmod="EXP"` as follows:

``` r
fit.ameras.exp <- ameras(Y="Y.binomial", dosevars=dosevars, X=c("X1","X2"), data=data, 
                            family="binomial", deg=2, doseRRmod = "EXP", methods="RC")
#> Fitting RC
#> Obtaining profile likelihood CI for dose
#> Obtaining profile likelihood CI for dose_squared
#> Warning in ameras.rc(family = family, dosevars = dosevars, data = data, :
#> P-value for dose_squared upper bound more than 0.005 away from 0.05, reducing
#> tol.profCI and/or increasing maxit.profCI is recommended
#> Warning in ameras.rc(family = family, dosevars = dosevars, data = data, :
#> P-value for dose_squared lower bound more than 0.005 away from 0.05, reducing
#> tol.profCI and/or increasing maxit.profCI is recommended
summary(fit.ameras.exp)
#> Call:
#> ameras(data = data, family = "binomial", Y = "Y.binomial", dosevars = dosevars, 
#>     X = c("X1", "X2"), methods = "RC", deg = 2, doseRRmod = "EXP")
#> 
#> Total run time: 3.7 seconds
#> 
#> Runtime in seconds by method:
#> 
#>  Method Runtime
#>      RC     3.7
#> 
#> Summary of coefficients by method:
#> 
#>  Method         Term Estimate      SE CI.lowerbound CI.upperbound
#>      RC  (Intercept) -0.94461 0.08409            NA            NA
#>      RC           X1  0.44552 0.07667            NA            NA
#>      RC           X2 -0.33376 0.09601            NA            NA
#>      RC         dose  0.37904 0.10388       0.17336       0.57642
#>      RC dose_squared  0.01943 0.02750      -0.03282       0.07785
```

## Linear excess relative risk

The linear excess relative risk model is given by
$RR_{i} = 1 + \beta_{1}D_{i} + \beta_{2}D_{i}^{2} + \mathbf{M}_{i}^{T}{\mathbf{β}}_{m1}D_{i} + \mathbf{M}_{i}^{T}{\mathbf{β}}_{m2}D_{i}^{2}$,
where again the quadratic and effect modification terms are optional.
This model is fit by setting `doseRRmod="ERR"` as follows:

``` r
fit.ameras.err <- ameras(Y="Y.binomial", dosevars=dosevars, X=c("X1","X2"), data=data, 
                            family="binomial", deg=2, doseRRmod = "ERR", methods="RC")
#> Fitting RC
#> Obtaining profile likelihood CI for dose
#> Warning in ameras.rc(family = family, dosevars = dosevars, data = data, :
#> WARNING: Lower bound for dose is < 0 and may not exist if rescaling the
#> variable does not help
#> Obtaining profile likelihood CI for dose_squared
summary(fit.ameras.err)
#> Call:
#> ameras(data = data, family = "binomial", Y = "Y.binomial", dosevars = dosevars, 
#>     X = c("X1", "X2"), methods = "RC", deg = 2, doseRRmod = "ERR")
#> 
#> Total run time: 5.6 seconds
#> 
#> Runtime in seconds by method:
#> 
#>  Method Runtime
#>      RC     5.6
#> 
#> Summary of coefficients by method:
#> 
#>  Method         Term Estimate      SE CI.lowerbound CI.upperbound
#>      RC  (Intercept) -0.87359 0.09759            NA            NA
#>      RC           X1  0.44587 0.07672            NA            NA
#>      RC           X2 -0.33552 0.09610            NA            NA
#>      RC         dose  0.04878 0.21283        0.0000        0.5115
#>      RC dose_squared  0.28763 0.08100        0.1325        0.4108
```

## Linear-exponential relative risk

The linear-exponential relative risk model is given by
$RR_{i} = 1 + \left( \beta_{1} + \mathbf{M}_{i}^{T}{\mathbf{β}}_{m1} \right)D_{i}\exp\{\left( \beta_{2} + \mathbf{M}_{i}^{T}{\mathbf{β}}_{m2} \right)D_{i}\}$,
where the effect modification terms are optional. This model is fit by
setting `doseRRmod="LINEXP"` as follows:

``` r
fit.ameras.linexp <- ameras(Y="Y.binomial", dosevars=dosevars, X=c("X1","X2"), data=data, 
                            family="binomial", doseRRmod = "LINEXP", methods="RC")
#> Fitting RC
#> Obtaining profile likelihood CI for dose_linear
#> Obtaining profile likelihood CI for dose_exponential
summary(fit.ameras.linexp)
#> Call:
#> ameras(data = data, family = "binomial", Y = "Y.binomial", dosevars = dosevars, 
#>     X = c("X1", "X2"), methods = "RC", doseRRmod = "LINEXP")
#> 
#> Total run time: 5.2 seconds
#> 
#> Runtime in seconds by method:
#> 
#>  Method Runtime
#>      RC     5.2
#> 
#> Summary of coefficients by method:
#> 
#>  Method             Term Estimate      SE CI.lowerbound CI.upperbound
#>      RC      (Intercept)  -0.9326 0.08592            NA            NA
#>      RC               X1   0.4456 0.07668            NA            NA
#>      RC               X2  -0.3343 0.09603            NA            NA
#>      RC      dose_linear   0.3255 0.11919        0.1473        0.6339
#>      RC dose_exponential   0.3455 0.10814        0.1452        0.5770
```

## Comparison between models

To compare between models, it is easiest to do so visually:

``` r
ggplot(data.frame(x=c(0, 5)), aes(x))+
  theme_light(base_size=15)+
  xlab("Exposure")+
  ylab("Relative risk")+
  labs(col="Model", lty="Model") +
  theme(legend.position = "inside", 
        legend.position.inside = c(.2,.85),
        legend.box.background = element_rect(color = "black", fill = "white", linewidth = 1))+
  stat_function(aes(col="Linear-quadratic ERR", lty="Linear-quadratic ERR" ),fun=function(x){
    1+fit.ameras.err$RC$coefficients["dose"]*x + fit.ameras.err$RC$coefficients["dose_squared"]*x^2
  }, linewidth=1.2) + 
  stat_function(aes(col="Exponential", lty="Exponential"),fun=function(x){
    exp(fit.ameras.exp$RC$coefficients["dose"]*x + fit.ameras.exp$RC$coefficients["dose_squared"]*x^2)
  }, linewidth=1.2) +
  stat_function(aes(col="Linear-exponential", lty="Linear-exponential"),fun=function(x){
    1+fit.ameras.linexp$RC$coefficients["dose_linear"]*x * exp(fit.ameras.linexp$RC$coefficients["dose_exponential"]*x)
  }, linewidth=1.2)
```

![](relativeriskmodels_files/figure-html/comparison-1.png)
