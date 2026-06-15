# Standard analyses with one dose realization

``` r

library(ameras)
```

## Introduction

The main purpose of `ameras` is to fit association models when multiple
realizations of an exposure are available. The same interface can also
be used for a standard analysis with a single observed dose variable.

By default,
[`ameras()`](https://ameras.sanderroberti.com/reference/ameras.md) fits
regression calibration (`RC`). With a single dose realization, the
regression-calibrated dose is just the supplied dose column, so the
examples below omit the `methods` argument and use the default `RC` fit.

For non-Gaussian families,
[`ameras()`](https://ameras.sanderroberti.com/reference/ameras.md)
defaults to a linear excess relative risk model (`model = "ERR"`). To
compare with ordinary [`glm()`](https://rdrr.io/r/stats/glm.html) fits,
use the exponential relative risk model (`model = "EXP"`), which puts
the dose coefficient on the usual logit or log-link scale.

``` r

data(data, package = "ameras")

## Use one exposure realization as the single observed dose.
data$D <- data$V1
head(data[c("D", "Y.gaussian", "Y.binomial", "Y.poisson", "X1", "X2")])
#>            D  Y.gaussian Y.binomial Y.poisson X1 X2
#> 1 0.42868043 -0.32647093          0         0  0  0
#> 2 0.73321154 -0.18734036          1         0  1  0
#> 3 0.70369712  0.08404044          0         2  0  0
#> 4 0.01845324  0.22432504          0         0  0  0
#> 5 0.39389441 -0.46317255          0         0  1  0
#> 6 0.01493158 -1.41036573          0         0  1  1
```

This small helper prints matching coefficients from an `ameras` fit and
a standard R model. The dose coefficient is named `dose` by `ameras` and
`D` by the standard model.

``` r

compare_coefficients <- function(ameras_fit, standard_fit) {
  ameras_terms <- c("(Intercept)", "X1", "X2", "dose")
  standard_terms <- c("(Intercept)", "X1", "X2", "D")

  out <- data.frame(
    term = ameras_terms,
    ameras = unname(ameras_fit$RC$coefficients[ameras_terms]),
    standard = unname(stats::coef(standard_fit)[standard_terms])
  )

  out$ameras <- round(out$ameras, 4)
  out$standard <- round(out$standard, 4)
  out
}
```

## Gaussian outcome

For a Gaussian outcome, a single-dose
[`ameras()`](https://ameras.sanderroberti.com/reference/ameras.md) fit
gives the same regression coefficients as
[`lm()`](https://rdrr.io/r/stats/lm.html), with only a small numerical
difference. The `ameras` fit also estimates `sigma`, so the comparison
below only shows the regression coefficients.

``` r

fit_ameras_gaussian <- ameras(Y.gaussian ~ dose(D) + X1 + X2,
                              data = data,
                              family = "gaussian")

fit_lm <- lm(Y.gaussian ~ D + X1 + X2, data = data)

compare_coefficients(fit_ameras_gaussian, fit_lm)
#>          term  ameras standard
#> 1 (Intercept) -1.2903  -1.2902
#> 2          X1  0.4938   0.4937
#> 3          X2 -0.5152  -0.5152
#> 4        dose  1.0920   1.0920
```

The usual `amerasfit` methods are available, including
[`summary()`](https://rdrr.io/r/base/summary.html),
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html), and
[`confint()`](https://rdrr.io/r/stats/confint.html).

``` r

fit_ameras_gaussian <- confint(fit_ameras_gaussian,
                               type = "wald.orig",
                               print = FALSE)
summary(fit_ameras_gaussian)
#> Call:
#> ameras(formula = Y.gaussian ~ dose(D) + X1 + X2, data = data, 
#>     family = "gaussian")
#> 
#> Total CPU runtime: 0.4 seconds
#> 
#> CPU runtime in seconds by method:
#> 
#>  Method   Fit  CI Total
#>      RC 0.364 0.0 0.364
#> 
#> Summary of coefficients by method:
#> 
#>  Method        Term Estimate      SE CI.lower CI.upper
#>      RC (Intercept)  -1.2903 0.03818  -1.3651  -1.2154
#>      RC          X1   0.4938 0.04236   0.4107   0.5768
#>      RC          X2  -0.5152 0.05206  -0.6172  -0.4131
#>      RC        dose   1.0920 0.02095   1.0509   1.1330
#>      RC       sigma   1.1597 0.01497   1.1304   1.1891
```

## Binary outcome

For a binary outcome, specify `model = "EXP"` to match the usual
logistic regression parameterization used by
`glm(..., family = binomial())`.

``` r

fit_ameras_binomial <- ameras(Y.binomial ~ dose(D, model = "EXP") + X1 + X2,
                              data = data,
                              family = "binomial")

fit_glm_binomial <- glm(Y.binomial ~ D + X1 + X2,
                        data = data,
                        family = binomial())

compare_coefficients(fit_ameras_binomial, fit_glm_binomial)
#>          term  ameras standard
#> 1 (Intercept) -0.9661  -0.9661
#> 2          X1  0.4513   0.4513
#> 3          X2 -0.3356  -0.3355
#> 4        dose  0.4374   0.4374
```

## Count outcome

The same idea applies to Poisson regression. Again, `model = "EXP"`
makes the `ameras` dose term directly comparable to the coefficient for
`D` in `glm(..., family = poisson())`.

``` r

fit_ameras_poisson <- ameras(Y.poisson ~ dose(D, model = "EXP") + X1 + X2,
                             data = data,
                             family = "poisson")

fit_glm_poisson <- glm(Y.poisson ~ D + X1 + X2,
                       data = data,
                       family = poisson())

compare_coefficients(fit_ameras_poisson, fit_glm_poisson)
#>          term  ameras standard
#> 1 (Intercept) -0.9178  -0.9178
#> 2          X1  0.5140   0.5140
#> 3          X2 -0.3569  -0.3570
#> 4        dose  0.3821   0.3821
```

## A note on comparison with other software

For a single dose realization, or for `RC` with multiple realizations
after replacing the dose by its realization-specific mean, point
estimates should generally agree with other software when the same
likelihood, dose-response model, covariates, and parameterization are
used.

Standard errors and Wald confidence intervals may nevertheless differ.
For `RC`, `ERC`, and `MCML`, `ameras` computes variance-covariance
matrices from the observed information, using a numerical Hessian of the
negative log-likelihood at the optimum. Other software may use expected
Fisher information, a scoring-based approximation, a different treatment
of parameter constraints, or a different internal parameter scale. These
choices can lead to different standard errors even when coefficient
estimates agree closely.

## Moving to multiple realizations

Once multiple exposure realizations are available, the formula changes
only in the `dose()` term. For example, replacing `D` by `V1:V10` fits
the default `RC` analysis using all 10 realizations:

``` r

fit_ameras_rc <- ameras(Y.binomial ~ dose(V1:V10, model = "EXP") + X1 + X2,
                        data = data,
                        family = "binomial")
summary(fit_ameras_rc)
#> Call:
#> ameras(formula = Y.binomial ~ dose(V1:V10, model = "EXP") + X1 + 
#>     X2, data = data, family = "binomial")
#> 
#> Total CPU runtime: 0.2 seconds
#> 
#> CPU runtime in seconds by method:
#> 
#>  Method  Fit  CI Total
#>      RC 0.18 0.0  0.18
#> 
#> Summary of coefficients by method:
#> 
#>  Method        Term Estimate      SE
#>      RC (Intercept)  -0.9760 0.07190
#>      RC          X1   0.4460 0.07667
#>      RC          X2  -0.3359 0.09596
#>      RC        dose   0.4471 0.04101
#> 
#> Note: confidence intervals not yet computed. Use confint() to add them.
```

Other methods can then be added through the `methods` argument:

``` r

fit_ameras_all <- ameras(Y.binomial ~ dose(V1:V10, model = "EXP") + X1 + X2,
                         data = data,
                         family = "binomial",
                         methods = c("RC", "ERC", "MCML", "FMA", "BMA"))
```
