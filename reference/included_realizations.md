# Extract Included Realizations from an amerasfit Object

Extracts the indices of dose realizations included in model averaging
for FMA and/or BMA results of a fitted `amerasfit` object. Realizations
are excluded from FMA if the optimization did not converge or the
Hessian was not positive definite, or if the realization receives no
samples after integer sample allocation.

## Usage

``` r
included_realizations(x, ...)

# S3 method for class 'amerasfit'
included_realizations(x,
                              methods = c("FMA", "BMA"),
                              ...)
```

## Arguments

- x:

  A fitted model object of class `amerasfit`, as returned by
  [`ameras`](https://ameras.sanderroberti.com/reference/ameras.md), with
  FMA and/or BMA results present.

- methods:

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
or the Hessian was not invertible or not positive definite. Realizations
with model-averaging weights too small to receive any samples for the
requested `MFMA` are also excluded from the final FMA result. Warnings
and messages are issued during fitting when realizations are excluded or
when the final FMA result is based on very few realizations.

For BMA, all realizations are included by default. A subset of
realizations can be specified via the `included.realizations.BMA`
argument to
[`ameras`](https://ameras.sanderroberti.com/reference/ameras.md), for
example by passing the indices returned by
`included_realizations(fit, methods="FMA")` from a prior FMA fit.

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
#> Fitting FMA
included_realizations(fit)
#>  [1]  1  2  3  4  5  6  7  8  9 10
length(included_realizations(fit))
#> [1] 10

## Both FMA and BMA
fit2 <- ameras(Y.binomial ~ dose(all_of(dosevars), model="ERR"),
               data=data, family="binomial", methods=c("FMA", "BMA"))
#> Note: BMA may require extensive computation time
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
included_realizations(fit2)
#> $FMA
#>  [1]  1  2  3  4  5  6  7  8  9 10
#> 
#> $BMA
#>  [1]  1  2  3  4  5  6  7  8  9 10
#> 
included_realizations(fit2, methods="FMA")
#>  [1]  1  2  3  4  5  6  7  8  9 10
included_realizations(fit2, methods="BMA")
#>  [1]  1  2  3  4  5  6  7  8  9 10
# }
```
