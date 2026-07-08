# Optimizer Convergence Diagnostics

Compute or extract numerical gradient diagnostics for fitted `ameras`
models.

## Usage

``` r
convergence(object, ...)

# S3 method for class 'amerasfit'
convergence(
  object,
  methods = c("RC", "ERC", "MCML"),
  data = NULL,
  recompute = FALSE,
  warn = FALSE,
  ...
)
```

## Arguments

- object:

  an object of class `amerasfit`.

- methods:

  character vector of fitted methods for which diagnostics should be
  returned. Gradient diagnostics are available for `"RC"`, `"ERC"`, and
  `"MCML"`.

- data:

  optional data frame required when diagnostics must be recomputed for a
  fit created with `keep.data=FALSE`.

- recompute:

  logical, whether to recompute the numerical gradient even when stored
  diagnostics are available.

- warn:

  logical, whether to emit the same warning used during fitting when
  [`optim()`](https://rdrr.io/r/stats/optim.html) reported convergence
  but the optimizer diagnostics suggest that the solution may not be
  fully stationary.

- ...:

  currently unused.

## Value

A data frame with one row per requested fitted method and columns:

- method:

  method name.

- optim.convergence:

  optimizer convergence code, where 0 indicates convergence.

- gradient.rms:

  root mean square of the numerical gradient at the fitted parameter
  vector.

- gradient.rms.scaled:

  `gradient.rms` divided by `max(1, abs(objective value))`.

- newton.improvement:

  approximate objective improvement still available from one Newton
  step, computed as \\0.5 g^T H^{-1} g\\.

- newton.improvement.relative:

  `newton.improvement` divided by `max(1, abs(objective value))`. This
  makes the Newton-step diagnostic less sensitive to the scale of the
  likelihood.

- convergence.warning:

  logical indicator for whether the ameras diagnostic threshold would
  warn about potentially unreliable optimizer convergence.

## Details

By default, `convergence()` returns stored diagnostics when available
and recomputes them only when needed or requested. If diagnostics must
be recomputed for a fit created with `keep.data=FALSE`, the original
analysis data must be supplied through `data`.

These diagnostics do not change the fitted object. They are intended to
help identify fits where [`optim()`](https://rdrr.io/r/stats/optim.html)
returned convergence code 0 but the objective is not locally flat enough
to be fully reassuring. The `convergence.warning` column is based on
`optim.convergence` and, when available, `newton.improvement` and
`newton.improvement.relative`. Both Newton-improvement diagnostics must
exceed their warning thresholds: `newton.improvement > 1e-2` and
`newton.improvement.relative > 1e-6`. This is less sensitive to
likelihood scale than the raw gradient norm, because raw gradients can
be large for likelihoods with many observations even when the fitted
coefficients are effectively at the optimum. If the Hessian cannot be
used for this calculation, ameras falls back to `gradient.rms > 1e-3`
and `gradient.rms.scaled > 1e-4`. Poorly scaled continuous covariates,
such as uncentered calendar year variables, can contribute to this kind
of numerical behavior.

Gradient diagnostics are not returned for FMA or BMA. For BMA, use MCMC
diagnostics such as
[`Rhat`](https://ameras.sanderroberti.com/reference/Rhat.md).

## See also

[`ameras`](https://ameras.sanderroberti.com/reference/ameras.md),
[`Rhat`](https://ameras.sanderroberti.com/reference/Rhat.md)

## Examples

``` r
# \donttest{
  data(data, package="ameras")
  fit <- ameras(Y.gaussian ~ dose(V1:V10) + X1 + X2,
                data=data, family="gaussian", methods="RC")
#> Fitting RC
  convergence(fit)
#>    method optim.convergence gradient.rms gradient.rms.scaled newton.improvement
#> RC     RC                 0   0.01218529        2.670265e-06       4.442178e-07
#>    newton.improvement.relative convergence.warning
#> RC                 9.73452e-11               FALSE
# }
```
