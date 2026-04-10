# ameras

The goal of ameras is to provide a user-friendly interface to analyze
association studies with multiple replicates of a noisy exposure using a
variety of methods. ameras supports continuous, count, binary,
multinomial, and right-censored time-to-event outcomes. For binary
outcomes, the nested case-control design is also accommodated. Besides
the common exponential relative risk model $RR = \exp(\beta D)$ for the
exposure-outcome association with noisy exposure $D$, linear excess
relative risk $RR = 1 + \beta D$ and linear-exponential excess relative
risk models $RR = 1 + \beta_{1}D\exp\left( \beta_{2}D \right)$ can be
used.

## Installation

To install from CRAN:

``` r
install.packages("ameras")
```

To install the development version from GitHub:

``` r
pak::pak("sanderroberti/ameras")
```

## Example

This is a basic example which shows you how to fit a simple logistic
regression model:

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
data(data, package="ameras")
dosevars <- paste0("V", 1:10)
fit <- ameras(data, family="binomial", Y="Y.binomial", methods=c("RC","ERC","MCML", "FMA", "BMA"), 
              dosevars=dosevars)
#> Note: BMA may require extensive computation time in the order of multiple hours
#> Fitting RC
#> Obtaining profile likelihood CI for dose
#> Fitting ERC
#> Note: computation times for profile likelihood intervals for ERC may be extensive with large datasets or complex models
#> Obtaining profile likelihood CI for dose
#> Fitting MCML
#> Note: computation times for profile likelihood intervals for MCML may be extensive with large datasets or complex models
#> Obtaining profile likelihood CI for dose
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
#> Total run time: 67.4 seconds
#> 
#> Runtime in seconds by method:
#> 
#>  Method Runtime
#>      RC     0.3
#>     ERC    21.6
#>    MCML     0.5
#>     FMA     0.2
#>     BMA    44.8
#> 
#> Summary of coefficients by method:
#> 
#>  Method        Term Estimate      SE CI.lowerbound CI.upperbound Rhat  n.eff
#>      RC (Intercept)  -0.8847 0.07378            NA            NA   NA     NA
#>      RC        dose   0.8020 0.13751        0.5648        1.1116   NA     NA
#>     ERC (Intercept)  -0.8849 0.07477            NA            NA   NA     NA
#>     ERC        dose   0.8214 0.14304        0.5728        1.1439   NA     NA
#>    MCML (Intercept)  -0.8758 0.07323            NA            NA   NA     NA
#>    MCML        dose   0.7910 0.13644        0.5554        1.0981   NA     NA
#>     FMA (Intercept)  -0.8761 0.07339       -1.0210       -0.7322   NA     NA
#>     FMA        dose   0.7918 0.13673        0.5257        1.0610   NA     NA
#>     BMA (Intercept)  -0.8736 0.07505       -1.0239       -0.7289 1.00 894.00
#>     BMA        dose   0.7951 0.14051        0.5438        1.1021 1.00 990.00
```
