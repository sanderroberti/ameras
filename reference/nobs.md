# Number of Observations Used to Fit an ameras Model

Returns the number of observations used to fit an `amerasfit` object
after missing-value handling and any family-specific row exclusions.

## Usage

``` r
# S3 method for class 'amerasfit'
nobs(object, ...)
```

## Arguments

- object:

  A fitted model object of class `amerasfit`, as returned by
  [`ameras`](https://ameras.sanderroberti.com/reference/ameras.md).

- ...:

  Additional arguments, currently unused.

## Value

The number of observations used for fitting. This excludes rows omitted
by `na.action`, uninformative matched-set rows for conditional logistic
regression, and zero-follow-up rows for proportional hazards regression
with entry times.

## See also

[`ameras`](https://ameras.sanderroberti.com/reference/ameras.md) for
model fitting,
[`summary`](https://ameras.sanderroberti.com/reference/summary.md) for
detailed row-exclusion counts, and
[`nobs`](https://rdrr.io/r/stats/nobs.html) for the generic function.

## Examples

``` r
example_data <- ameras_simulated[1:100, ]
example_data$X1[1] <- NA

fit <- ameras(
  Y.gaussian ~ dose(V1:V10) + X1,
  data = example_data,
  family = "gaussian",
  methods = "RC",
  na.action = na.omit
)
#> Fitting RC

nobs(fit)
#> [1] 99
```
