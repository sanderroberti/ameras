# Likelihood-ratio tests for dose-related parameters

Compute global and individual likelihood-ratio tests for dose-related
parameters in supported `amerasfit` methods.

## Usage

``` r
dose_lrt(object, ...)

# S3 method for class 'amerasfit'
dose_lrt(
  object,
  methods = c("RC", "ERC", "MCML"),
  type = c("global", "individual"),
  data = NULL,
  ...
)
```

## Arguments

- object:

  An object of class `amerasfit`.

- methods:

  Character vector of fitted methods to test. Supported methods are
  `"RC"`, `"ERC"`, and `"MCML"`.

- type:

  Character vector selecting which tests to return. `"global"` tests all
  dose-related parameters jointly. `"individual"` tests each
  dose-related parameter separately. By default, both are returned.

- data:

  Optional data frame required when `object` was fitted with
  `keep.data=FALSE`.

- ...:

  Currently unused.

## Value

A data frame with one row per likelihood-ratio test and columns:

- `method`:

  The ameras estimation method.

- `type`:

  Either `"global"` or `"individual"`.

- `term`:

  The tested dose term. For global tests this is `"all dose terms"`.

- `df`:

  Degrees of freedom for the likelihood-ratio test.

- `logLik`:

  The fitted model log likelihood.

- `logLik.null`:

  The constrained null-model log likelihood.

- `statistic`:

  The likelihood-ratio statistic.

- `p.value`:

  The chi-square reference p-value.

- `null.optim.convergence`:

  Optimizer convergence code for the constrained null fit. A warning is
  issued if this code is nonzero.

## Details

`dose_lrt()` computes likelihood-ratio tests for dose-related parameters
in fitted RC, ERC, and MCML models. The null hypothesis sets the tested
dose-related parameter(s) equal to zero. FMA and BMA are not supported
by this function.

For transformed fits, nuisance parameters are evaluated through the
fitted transformation. The tested dose parameters are then fixed to zero
before the likelihood is evaluated. This is most straightforward for
component-wise transformations, including the default transformations
used by ameras. Results from coupled custom transformations should be
interpreted with caution. For example, if a custom transformation
defines one reported parameter as a function of multiple optimizer-scale
parameters, or makes a non-dose reported parameter depend directly on a
dose parameter, fixing only the reported dose parameter to zero may not
correspond to the intended constrained model.

For speed, constrained null fits report only the optimizer convergence
code. They do not run the additional numerical-gradient diagnostics
stored for the original RC, ERC, and MCML fits.

## Examples

``` r
set.seed(1)
dat <- data.frame(
  Y = 1 + 0.8 * seq(-1, 1, length.out = 30) + rnorm(30, sd = 0.2),
  D = seq(-1, 1, length.out = 30)
)

fit <- ameras(Y ~ dose(D), data = dat, family = "gaussian",
              methods = "RC", print = FALSE)
#> Fitting RC

dose_lrt(fit)
#>   method       type           term df   logLik logLik.null statistic
#> 1     RC     global all dose terms  1 8.776784   -21.27384  60.10124
#> 2     RC individual           dose  1 8.776784   -21.27384  60.10124
#>        p.value null.optim.convergence
#> 1 9.010161e-15                      0
#> 2 9.010161e-15                      0
dose_lrt(fit, type = "global")
#>   method   type           term df   logLik logLik.null statistic      p.value
#> 1     RC global all dose terms  1 8.776784   -21.27384  60.10124 9.010161e-15
#>   null.optim.convergence
#> 1                      0
```
