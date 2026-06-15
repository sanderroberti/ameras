# Manual FMA from regression calibration fits

``` r

library(ameras)
```

## Introduction

Frequentist model averaging (`FMA`) in
[`ameras()`](https://ameras.sanderroberti.com/reference/ameras.md) fits
the same association model once for each dose realization, computes AIC
weights across those realization specific fits, and samples from the
fitted normal approximation for each realization in proportion to its
weight.

Usually, the simplest approach is to let
[`ameras()`](https://ameras.sanderroberti.com/reference/ameras.md) do
this directly with `methods = "FMA"`. For large datasets or many dose
realizations, however, it can be useful to split the work into smaller
jobs. This vignette shows the basic manual workflow: fit separate `RC`
models for individual dose realizations, save the fitted summaries if
desired, then compute the FMA weights and samples across all
realization-specific fits.

If the direct FMA fit is feasible in one R session and the main goal is
to use multiple local workers, see the [parallel FMA
vignette](https://ameras.sanderroberti.com/articles/parallelfma.md)
instead. The manual workflow shown here is mainly for analyses that need
to be split into separate jobs or resumed from saved
realization-specific summaries.

The example uses a Gaussian model to keep the code compact and avoid any
extra parameter transformations. The same idea applies to other
families, provided that every realization is fit with the same formula
structure, family, dose-response model, starting values, optimizer
settings, and transformation arguments.

## Fit realization-specific RC models

``` r

data(data, package = "ameras")
dosevars <- paste0("V", 1:10)
```

The first helper mirrors the convergence and Hessian screen used by
built-in FMA. The second helper fits one standard regression calibration
model with a single dose realization. The returned list keeps only the
pieces needed for model averaging: the coefficient estimates,
variance-covariance matrix, maximized log-likelihood, and AIC.

``` r

passes_fma_screen <- function(rc) {
  hessian <- rc$optim$hessian

  if (is.null(hessian) || rc$optim$convergence != 0) {
    return(FALSE)
  }

  isTRUE(det(hessian) != 0 &&
    rcond(hessian) > .Machine$double.eps &&
    all(eigen(hessian)$values > 0))
}

fit_rc_realization <- function(dosevar, data) {
  formula <- stats::as.formula(
    paste0("Y.gaussian ~ dose(", dosevar, ") + X1 + X2")
  )

  fit <- ameras(formula, data = data, family = "gaussian")
  rc <- fit$RC

  list(
    dosevar = dosevar,
    coefficients = rc$coefficients,
    vcov = rc$vcov,
    loglik = rc$loglik,
    AIC = -2 * rc$loglik + 2 * length(rc$coefficients),
    include = passes_fma_screen(rc)
  )
}
```

Fit the model separately for each dose realization.

``` r

rc_summaries <- lapply(dosevars, fit_rc_realization, data = data)
```

For large analyses, the same call can be run in separate R sessions or
submitted as separate jobs. Save the summaries from each job and read
them back before model averaging.

``` r

dir.create("fma-fits", showWarnings = FALSE)

for (dosevar in dosevars) {
  fit_summary <- fit_rc_realization(dosevar, data = data)
  saveRDS(fit_summary, file = file.path("fma-fits", paste0(dosevar, ".rds")))
}

fit_files <- file.path("fma-fits", paste0(dosevars, ".rds"))
rc_summaries <- lapply(fit_files, readRDS)
```

## Assemble the FMA result

The next helper performs the model averaging step. It first drops any
realization that did not pass the FMA screen. It then computes
stabilized AIC weights, allocates `MFMA` samples across realizations,
draws from each realization-specific normal approximation, and
summarizes the combined samples. The interval calculation shown here is
the equal-tailed percentile interval, matching
`confint(..., type = "percentile")` for FMA. `ameras` also supports
highest posterior density intervals for FMA with `type = "hpd"`, but
that is not implemented in this manual helper.

``` r

assemble_manual_fma <- function(rc_summaries, MFMA = 100000) {
  n_total <- length(rc_summaries)
  valid <- vapply(rc_summaries, function(x) isTRUE(x$include), logical(1))
  n_screen_excluded <- sum(!valid)

  if (!any(valid)) {
    stop("No realization-specific fits are available for FMA.")
  }

  original_indices <- which(valid)
  rc_summaries <- rc_summaries[valid]

  AIC <- vapply(rc_summaries, `[[`, numeric(1), "AIC")
  weights <- exp(-0.5 * (AIC - min(AIC)))
  weights <- weights / sum(weights)

  allocated_samples <- round(weights * MFMA)
  keep <- allocated_samples > 0
  n_weight_excluded <- sum(!keep)

  if (!any(keep)) {
    stop("No realizations received FMA samples. Increase MFMA.")
  }

  rc_summaries <- rc_summaries[keep]
  original_indices <- original_indices[keep]
  weights <- weights[keep]
  allocated_samples <- allocated_samples[keep]
  names(weights) <- vapply(rc_summaries, `[[`, character(1), "dosevar")
  diagnostics <- data.frame(
    stage = c(
      "total realization-specific fits",
      "excluded by convergence/Hessian screen",
      "excluded after zero-sample allocation",
      "included in manual FMA"
    ),
    realizations = c(
      n_total,
      n_screen_excluded,
      n_weight_excluded,
      length(original_indices)
    ),
    row.names = NULL
  )

  samples <- do.call(
    "rbind",
    Map(
      function(fit, n) {
        mvtnorm::rmvnorm(n = n, mean = fit$coefficients, sigma = fit$vcov)
      },
      rc_summaries,
      allocated_samples
    )
  )

  samples <- as.data.frame(samples)
  names(samples) <- names(rc_summaries[[1]]$coefficients)

  coefficients <- colMeans(samples)
  se <- apply(samples, 2, stats::sd)
  percentile_CI <- t(
    apply(samples, 2, stats::quantile, probs = c(0.025, 0.975))
  )
  colnames(percentile_CI) <- c("lower", "upper")

  list(
    coefficients = coefficients,
    SE = se,
    percentile_CI = percentile_CI,
    vcov = stats::var(samples),
    weights = weights,
    samples = samples,
    included.realizations = original_indices,
    included.samples = nrow(samples),
    diagnostics = diagnostics
  )
}
```

The diagnostics table reports how many realization-specific fits were
removed by the convergence/Hessian screen and how many passed that
screen but received zero FMA samples for the chosen `MFMA`.

``` r

set.seed(100)
manual_fma <- assemble_manual_fma(rc_summaries, MFMA = 100000)
manual_fma$diagnostics
#>                                    stage realizations
#> 1        total realization-specific fits           10
#> 2 excluded by convergence/Hessian screen            0
#> 3  excluded after zero-sample allocation            6
#> 4                 included in manual FMA            4

manual_summary <- data.frame(
  term = names(manual_fma$coefficients),
  estimate = unname(manual_fma$coefficients),
  SE = unname(manual_fma$SE),
  percentile_lower = manual_fma$percentile_CI[, "lower"],
  percentile_upper = manual_fma$percentile_CI[, "upper"],
  row.names = NULL
)

manual_summary[-1] <- round(manual_summary[-1], 4)
manual_summary
#>          term estimate     SE percentile_lower percentile_upper
#> 1 (Intercept)  -1.2805 0.0370          -1.3531          -1.2083
#> 2          X1   0.4837 0.0412           0.4034           0.5644
#> 3          X2  -0.5168 0.0511          -0.6180          -0.4171
#> 4        dose   1.0792 0.0200           1.0401           1.1185
#> 5       sigma   1.1378 0.0147           1.1089           1.1665
```

The weights show which realization-specific fits contributed samples
after integer allocation of `MFMA`.

``` r

round(manual_fma$weights, 4)
#>     V3     V4     V8     V9 
#> 0.0000 0.0006 0.9942 0.0052
manual_fma$included.realizations
#> [1] 3 4 8 9
manual_fma$included.samples
#> [1] 100000
```

## Compare with built-in FMA

For this small example, the direct
[`ameras()`](https://ameras.sanderroberti.com/reference/ameras.md) FMA
fit is still fast. The manual calculation and built-in calculation use
the same likelihoods, AIC weights, and normal sampling strategy.

``` r

set.seed(100)
fit_builtin_fma <- suppressWarnings(
  ameras(Y.gaussian ~ dose(V1:V10) + X1 + X2,
         data = data,
         family = "gaussian",
         methods = "FMA",
         MFMA = 100000)
)

comparison <- data.frame(
  term = names(manual_fma$coefficients),
  manual = unname(manual_fma$coefficients),
  builtin = unname(fit_builtin_fma$FMA$coefficients[names(manual_fma$coefficients)]),
  row.names = NULL
)

comparison[-1] <- round(comparison[-1], 4)
comparison
#>          term  manual builtin
#> 1 (Intercept) -1.2805 -1.2805
#> 2          X1  0.4837  0.4837
#> 3          X2 -0.5168 -0.5168
#> 4        dose  1.0792  1.0792
#> 5       sigma  1.1378  1.1378
```

The corresponding model-averaging weights are also the same.

``` r

weight_comparison <- merge(
  data.frame(dose = names(manual_fma$weights),
             manual = unname(manual_fma$weights)),
  data.frame(dose = names(fit_builtin_fma$FMA$weights),
             builtin = unname(fit_builtin_fma$FMA$weights)),
  by = "dose",
  all = TRUE
)

weight_comparison[-1] <- round(weight_comparison[-1], 4)
weight_comparison
#>   dose manual builtin
#> 1   V3 0.0000  0.0000
#> 2   V4 0.0006  0.0006
#> 3   V8 0.9942  0.9942
#> 4   V9 0.0052  0.0052
```

## Adapting the example to another analysis

To adapt this template, start from the ordinary
[`ameras()`](https://ameras.sanderroberti.com/reference/ameras.md) call
that you would have used for built-in FMA. For example, if the direct
call would use `methods = "FMA"` and a dose term such as
`dose(V1:V100, model = "EXP", deg = 2)`, the manual version keeps the
same outcome, covariates, family, dose-response model, transformations,
starting values, and optimizer settings. The only structural change is
that `fit_rc_realization()` substitutes one dose column at a time for
the full set of dose realizations.

In practice:

- set `dosevars` to the dose realization columns from the original call;
- update the formula inside `fit_rc_realization()` to use the same
  outcome, covariates, modifiers, and `dose()` options, replacing only
  the multi-column dose selection with `dosevar`;
- update the `family` argument, and pass through any arguments such as
  `transform`, `transform.jacobian`, `inpar`, `optim.method`, or
  `control` that the direct FMA call would have used;
- leave `methods` out of the realization-specific fits, because `RC` is
  the default and each fit should use a single dose realization;
- pass the desired number of FMA samples to `assemble_manual_fma()`
  through `MFMA`.

The manual combining step does not re-apply transformations. It uses the
`RC$coefficients` and `RC$vcov` returned by
[`ameras()`](https://ameras.sanderroberti.com/reference/ameras.md),
which are already on the original output scale. The AIC weights must be
computed globally across all realizations, not separately within each
job.

## Practical notes

This manual workflow is mainly a computational workaround. If the direct
`ameras(..., methods = "FMA")` fit is feasible, it is less error-prone
and returns a regular `amerasfit` object that works directly with
methods such as [`summary()`](https://rdrr.io/r/base/summary.html) and
[`confint()`](https://rdrr.io/r/stats/confint.html).

When running the realization-specific fits separately, keep the
following points fixed across every `RC` fit:

- the model formula, except for the one dose realization being fit;
- the `family`, dose-response model, degree, modifiers, covariates,
  offsets, and survival or risk-set inputs;
- optimizer settings, starting values, and any transformation arguments;
- the scale on which coefficients and variance-covariance matrices are
  stored.

The AIC weights must be computed globally across all realizations, not
separately within separate jobs. For this reason, each job must save the
maximized log-likelihood, coefficients, and variance-covariance matrix
for its realization-specific fit.
