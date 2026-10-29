# Traceplots for MCMC Samples

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

   fit <- ameras(Y.gaussian~dose(V1:V10), data, methods="BMA")
#> Error in check_df(data): data must be a data frame
   traceplot(fit)
#> Error: object 'fit' not found
```
