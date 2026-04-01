
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ameras

<!-- badges: start -->

<!-- badges: end -->

The goal of ameras is to provide a user-friendly interface to analyze
association studies with multiple replicates of a noisy exposure using a
variety of methods. ameras supports continuous, count, binary,
multinomial, and right-censored time-to-event outcomes. For binary
outcomes, the nested case-control design is also accommodated. Besides
the common exponential relative risk model $RR=\exp(\beta D)$ for the
exposure-outcome association with noisy exposure $D$, linear excess
relative risk $RR=1+\beta D$ and linear-exponential excess relative risk
models $RR=1+\beta_1 D \exp(\beta_2 D)$ can be used.

## Installation

To install from CRAN:

``` r
# install.packages("ameras")
```

## Example

This is a basic example which shows you how to fit a simple logistic
regression model:

``` r
library(ameras)
#> Loading required package: nimble
#> nimble version 1.4.1 is loaded.
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
data(data, package="ameras")
dosevars <- paste0("V", 1:10)
fit <- ameras(data, family="binomial", Y="Y.binomial", methods=c("RC","ERC","MCML", "FMA", "BMA"), 
              dosevars=dosevars)
#> Fitting RC
#> Fitting ERC
#> Fitting MCML
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
summary(fit)
#> Call:
#> ameras(data = data, family = "binomial", Y = "Y.binomial", dosevars = dosevars, 
#>     methods = c("RC", "ERC", "MCML", "FMA", "BMA"))
#> 
#> Total run time: 31 seconds
#> 
#> Runtime in seconds by method:
#> 
#>  Method Runtime
#>      RC     0.0
#>     ERC     8.5
#>    MCML     0.1
#>     FMA     0.2
#>     BMA    22.2
#> 
#> Summary of coefficients by method:
#> 
#>  Method        Term Estimate      SE CI.lowerbound CI.upperbound Rhat  n.eff
#>      RC (Intercept)  -0.8847 0.07378       -1.0293       -0.7401   NA     NA
#>      RC        dose   0.8020 0.13751        0.5324        1.0715   NA     NA
#>     ERC (Intercept)  -0.8849 0.07477       -1.0315       -0.7384   NA     NA
#>     ERC        dose   0.8214 0.14304        0.5411        1.1018   NA     NA
#>    MCML (Intercept)  -0.8758 0.07323       -1.0193       -0.7323   NA     NA
#>    MCML        dose   0.7910 0.13644        0.5236        1.0584   NA     NA
#>     FMA (Intercept)  -0.8758 0.07321       -1.0197       -0.7329   NA     NA
#>     FMA        dose   0.7913 0.13635        0.5245        1.0580   NA     NA
#>     BMA (Intercept)  -0.8718 0.07342       -1.0201       -0.7342 1.00 291.00
#>     BMA        dose   0.7920 0.14133        0.5546        1.0999 1.00 281.00
```
