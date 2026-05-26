# Extract MCMC Convergence Diagnostics from an amerasfit Object

Extracts the Gelman-Rubin convergence diagnostics from the BMA results
of a fitted `amerasfit` object.

## Usage

``` r
Rhat(x, ...)

# S3 method for class 'amerasfit'
Rhat(x, ...)
```

## Arguments

- x:

  A fitted model object of class `amerasfit`, as returned by
  [`ameras`](https://ameras.sanderroberti.com/reference/ameras.md), with
  BMA results present.

- ...:

  Additional arguments, currently unused.

## Value

A data frame with columns `Rhat` and `n.eff`, containing the
Gelman-Rubin statistic and effective sample size for each parameter.
Values of `Rhat` substantially above 1.05 indicate potential convergence
problems, in which case longer chains via `niter.BMA` and `nburnin.BMA`
are recommended.

## See also

[`ameras`](https://ameras.sanderroberti.com/reference/ameras.md),
[`summary.amerasfit`](https://ameras.sanderroberti.com/reference/summary.md),
[`included_realizations.amerasfit`](https://ameras.sanderroberti.com/reference/included_realizations.md)

## Examples

``` r
# \donttest{
data("data", package="ameras")
dosevars <- paste0("V", 1:10)
fit <- ameras(Y.binomial ~ dose(all_of(dosevars), model="ERR"),
              data=data, family="binomial", methods="BMA")
#> Note: BMA may require extensive computation time
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
Rhat(fit)
#>             Rhat n.eff
#> (Intercept) 1.03  1236
#> dose        1.02   969
# }
```
