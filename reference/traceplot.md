# Traceplots for MCMC samples

Produce MCMC traceplots for `amerasfit` objects.

## Usage

``` r
traceplot(object, ...)

# S3 method for class 'amerasfit'
traceplot(object, iter = 5000, Rhat = TRUE, n.eff = TRUE, pdf = FALSE, ...)
```

## Arguments

- object:

  a `amerasfit` object containing BMA output to be plotted

- iter:

  number of iterations to include in the traceplot (defaults to last
  5000)

- Rhat:

  logical; whether to include R-hat diagnostics in the plot (default
  TRUE)

- n.eff:

  logical; whether to include effective sample size in the plot (default
  TRUE)

- pdf:

  logical; whether to save the output as a PDF (default FALSE)

- ...:

  additional arguments passed to `MCMCtrace`

## Details

Wrapper for
[`MCMCvis::MCMCtrace`](https://rdrr.io/pkg/MCMCvis/man/MCMCtrace.html)
to produce MCMC diagnostic plots. See `?MCMCtrace` for more plotting
options that can be provided through `...`.

## Value

Traceplots and posterior density plots.

## See also

[`MCMCtrace`](https://rdrr.io/pkg/MCMCvis/man/MCMCtrace.html)

## Examples

``` r
# \donttest{
   data(data, package="ameras")
   fit <- ameras(data, methods="BMA", Y="Y.gaussian", dosevars=paste0("V", 1:10))
#> Note: BMA may require extensive computation time in the order of multiple hours
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
   traceplot(fit)


# }
```
