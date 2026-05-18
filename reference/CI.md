# Extract Confidence Intervals from an amerasfit Object

Extracts confidence intervals from a fitted `amerasfit` object after
they have been computed via
[`confint`](https://ameras.sanderroberti.com/reference/confint.md). This
accessor is primarily intended for programmatic use, such as computing
coverage in simulation studies, where printing is not desired.

## Usage

``` r
CI(x, ...)

# S3 method for class 'amerasfit'
CI(x, method = c("RC", "ERC", "MCML", "FMA", "BMA"),
   parm = NULL, ...)
```

## Arguments

- x:

  A fitted model object of class `amerasfit`, as returned by
  [`ameras`](https://ameras.sanderroberti.com/reference/ameras.md), with
  confidence intervals computed via
  [`confint`](https://ameras.sanderroberti.com/reference/confint.md).

- method:

  Character vector specifying which method(s) to extract confidence
  intervals for. One or more of `"RC"`, `"ERC"`, `"MCML"`, `"FMA"`, and
  `"BMA"`. Defaults to all. Only methods that were run are returned.

- parm:

  Either `NULL` to return all parameters, `"dose"` to return
  dose-related parameters only, `"all"` to explicitly return all
  parameters, or a character vector of specific parameter names.
  Defaults to `NULL`.

- ...:

  Additional arguments, currently unused.

## Value

If a single method is requested or only one method was run, a data frame
with columns `lower` and `upper` and rows corresponding to parameters is
returned. If multiple methods are requested and multiple methods were
run, a named list of such data frames is returned, one per method.

## Details

Confidence intervals must be computed via
[`confint`](https://ameras.sanderroberti.com/reference/confint.md)
before calling `CI`. If confidence intervals have not yet been computed,
`NULL` is returned for all methods.

## See also

[`confint`](https://ameras.sanderroberti.com/reference/confint.md) for
computing confidence intervals,
[`summary.amerasfit`](https://ameras.sanderroberti.com/reference/summary.md)
for a summary of the fitted model including confidence intervals,
[`Rhat.amerasfit`](https://ameras.sanderroberti.com/reference/Rhat.md),
[`included_realizations.amerasfit`](https://ameras.sanderroberti.com/reference/included_realizations.md)

## Examples

``` r
data("data", package="ameras")
dosevars <- paste0("V", 1:10)

fit <- ameras(Y.binomial ~ dose(all_of(dosevars), model="ERR"),
              data=data, family="binomial", methods=c("RC", "FMA"))
#> Error in resolve_dose_selection(sel_args, data): ℹ In argument: `all_of(dosevars)`.
#> Caused by error:
#> ! object 'dosevars' not found
fit <- confint(fit, method=c("wald.orig","percentile"))
#> Error: object 'fit' not found

## Extract all CIs for RC
CI(fit, method="RC")
#> Error: object 'fit' not found

## Extract dose parameters only for all methods
CI(fit, parm="dose")
#> Error: object 'fit' not found

## Extract specific parameter
CI(fit, method="RC", parm="dose")
#> Error: object 'fit' not found
```
