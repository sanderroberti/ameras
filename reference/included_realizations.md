# Extract Included Realizations from an amerasfit Object

Extracts the indices of dose realizations included in model averaging
for FMA and/or BMA results of a fitted `amerasfit` object. Realizations
are excluded from FMA if the optimization did not converge or the
Hessian was not positive definite.

## Usage

``` r
included_realizations(x, ...)

# S3 method for class 'amerasfit'
included_realizations(x,
                              method = c("FMA", "BMA"),
                              ...)
```

## Arguments

- x:

  A fitted model object of class `amerasfit`, as returned by
  [`ameras`](https://ameras.sanderroberti.com/reference/ameras.md), with
  FMA and/or BMA results present.

- method:

  Character vector specifying which method(s) to extract included
  realizations for. One or both of `"FMA"` and `"BMA"`. Defaults to
  both. Only methods that were run are returned.

- ...:

  Additional arguments, currently unused.

## Value

If a single method is requested or only one method was run, an integer
vector of indices of the included dose realizations. If both FMA and BMA
results are present and both are requested, a named list with elements
`FMA` and `BMA`, each containing an integer vector of included
realization indices.

## Details

For FMA, realizations are excluded if the optimization did not converge
or the Hessian was not invertible or not positive definite. If more than
20% of realizations are excluded, a warning is issued during fitting.

For BMA, all realizations are included by default. A subset of
realizations can be specified via the `included.replicates.BMA` argument
to [`ameras`](https://ameras.sanderroberti.com/reference/ameras.md), for
example by passing the indices returned by
`included_realizations(fit, method="FMA")` from a prior FMA fit.

## See also

[`ameras`](https://ameras.sanderroberti.com/reference/ameras.md),
[`Rhat.amerasfit`](https://ameras.sanderroberti.com/reference/Rhat.md)

## Examples

``` r
# \donttest{
data("data", package="ameras")
dosevars <- paste0("V", 1:10)

## FMA only
fit <- ameras(Y.binomial ~ dose(all_of(dosevars), model="ERR"),
              data=data, family="binomial", methods="FMA")
#> Error in resolve_dose_selection(sel_args, data): ℹ In argument: `all_of(dosevars)`.
#> Caused by error:
#> ! object 'dosevars' not found
included_realizations(fit)
#> Error: object 'fit' not found
length(included_realizations(fit))
#> Error: object 'fit' not found

## Both FMA and BMA
fit2 <- ameras(Y.binomial ~ dose(all_of(dosevars), model="ERR"),
               data=data, family="binomial", methods=c("FMA", "BMA"))
#> Note: BMA may require extensive computation time
#> Error in resolve_dose_selection(sel_args, data): ℹ In argument: `all_of(dosevars)`.
#> Caused by error:
#> ! object 'dosevars' not found
included_realizations(fit2)
#> Error: object 'fit2' not found
included_realizations(fit2, method="FMA")
#> Error: object 'fit2' not found
included_realizations(fit2, method="BMA")
#> Error: object 'fit2' not found
# }
```
