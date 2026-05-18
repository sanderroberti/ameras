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
#> Error in resolve_dose_selection(sel_args, data): ℹ In argument: `all_of(dosevars)`.
#> Caused by error:
#> ! object 'dosevars' not found
Rhat(fit)
#> Error: object 'fit' not found
# }
```
