# Internal Helper Flow

This note describes how the main fitting and confidence interval helpers fit
together. It is intended for package maintainers; the functions named here are
internal implementation details, not user-facing API.

## Overview

```mermaid
flowchart TD
  A["ameras()"] --> B["parse/check inputs<br/>build transforms<br/>add rcdose_ameras"]
  B --> C["ameras_main()"]

  C -->|MCML| D["ameras.mcml()"]
  C -->|RC| E["ameras.rc(ERC = FALSE)"]
  C -->|ERC| F["ameras.rc(ERC = TRUE)"]
  C -->|FMA| G["ameras.fma()"]
  C -->|BMA| H["ameras.bma()"]

  D --> I["make_mcml_loglik_fn()"]
  I --> J["make_single_realization_loglik()"]
  J --> JA["family-specific likelihood closure"]
  I --> JB["MCML objective closure"]
  JA -->|wrapped by MCML objective| JB
  JB --> AB["mcml_average_neg_loglik()<br/>scalar objective value"]
  JB --> N["fit_objective_with_hessian()"]

  E --> K["family-specific RC objective"]
  F --> L["family-specific ERC objective"]
  K --> N
  L --> N

  G --> M["fit_fma_realizations()"]
  M --> ML["fma_realization_lapply()<br/>lapply() / future_lapply()"]
  ML --> J
  JA -->|FMA realization objective| N

  N -->|MCML, RC, ERC| O["assemble_frequentist_fit_result()<br/>MCML, RC, ERC"]
  N -->|FMA| P["summarize_fma_realization_fit()<br/>fit_passes_hessian_check()<br/>prepare_fma_sampling_inputs()"]

  G --> Q["assemble_fma_result()<br/>FMA weights and samples"]
  P --> Q
  H --> R["NIMBLE model and MCMC samples<br/>uses included.realizations.BMA"]

  O --> S["new_amerasfit()<br/>amerasfit object"]
  Q --> S
  R --> S

  S --> T["confint.amerasfit()"]
  T -->|wald.orig / wald.transformed| U["compute_wald_CI()"]
  T -->|proflik| V["compute_proflik_CI()"]
  T -->|percentile / hpd| W["compute_sample_CI()"]

  V --> X["make_loglik_fn()"]
  X --> Y["make_base_loglik_fn()"]
  Y --> Z["make_single_realization_loglik()<br/>ordinary RC/ERC/MCML cases"]
  Y --> AA["loglik.poisson.erc() / loglik.prophaz.erc()<br/>special ERC cases"]
  Z --> AE["base likelihood closure"]
  AA --> AE
  X --> AF["method-specific profile objective closure"]
  AE --> AF
  AF -->|MCML profile objective| AB
  AF --> AC["compute_proflik_ci_one()"]
  V --> AC
  AC --> AD["proflik()<br/>profile likelihood value"]
  U --> AG["CI table"]
  V --> AG
  W --> AG
  AG --> AH["confint stores CI on method result<br/>set_ci_timing()<br/>mark CI.computed"]
  AH --> AI["updated amerasfit object"]
```

## Fitting Path

`ameras()` handles the public interface. It parses the formula or legacy
arguments, validates inputs, constructs default transformations where needed,
adds `rcdose_ameras`, and stores enough model metadata for later methods such as
`confint()`, `residuals()`, and `plot()`.

`ameras_main()` dispatches to the method-specific fitting functions:

- `ameras.mcml()` fits MCML.
- `ameras.rc(ERC = FALSE)` fits regression calibration.
- `ameras.rc(ERC = TRUE)` fits extended regression calibration.
- `ameras.fma()` fits each dose realization and performs frequentist model
  averaging.
- `ameras.bma()` builds and samples the Bayesian model averaging model.

When FMA and BMA are requested together, `ameras_main()` runs FMA first so the
methods appear in a predictable order. BMA does not currently use FMA screening
results by default; `included.realizations.BMA` controls the BMA realization
set and defaults to all dose realizations.

The main likelihood construction helper is `make_single_realization_loglik()`.
It builds a family-specific likelihood closure from raw fitting arguments. It is
used directly by FMA and indirectly by MCML and profile-likelihood
reconstruction.
In the flowchart, arrows out of `make_single_realization_loglik()` point to the
likelihood closures it returns; those closures are then wrapped or passed to the
optimizer/profile-likelihood helpers.
Leaf helper nodes such as `mcml_average_neg_loglik()` and `proflik()` represent
values computed inside their caller rather than separate downstream package
steps.

`make_mcml_loglik_fn()` wraps `make_single_realization_loglik()` for MCML. The
family likelihoods return one negative log likelihood per dose realization, and
`mcml_average_neg_loglik()` combines those values with the stabilized
log-mean-exp calculation.

`fit_objective_with_hessian()` is the shared optimizer wrapper. It normalizes
the scalar `optimize()` path and the general `optim()` path into one fit object
shape, then attaches a numeric Hessian. MCML, RC/ERC, and FMA all use this
optimizer helper.

`assemble_frequentist_fit_result()` packages MCML and RC/ERC results. It applies
any transform, propagates variance through `transform.jacobian`, names
coefficients and covariance matrices, stores optimizer details, and preserves
method-specific fields such as `ERC`.

FMA uses `fit_fma_realizations()` to build and fit one likelihood per dose
realization. The realization loop goes through `fma_realization_lapply()`,
which uses `future.apply::future_lapply()` when available and otherwise falls
back to `lapply()`. `ameras()` never sets a `future::plan()` internally; users
control whether FMA runs sequentially or in parallel by setting the plan before
calling `ameras()`. Each fitted realization is summarized by
`summarize_fma_realization_fit()`, which calls `fit_passes_hessian_check()` to
screen optimizer convergence and Hessian usability before model averaging.
For Hessian-eligible realizations, `prepare_fma_sampling_inputs()` then creates
and validates the transformed or untransformed mean/covariance pair used for
FMA sampling. The summary object stores only these validated sampling inputs,
not the raw optimizer coefficients and Hessian.

`assemble_fma_result()` takes those realization summaries, reports the
realization-level exclusions, applies negligible-weight exclusions, computes
AIC weights, allocates `MFMA` samples, draws from the validated coefficient
sampling inputs, and packages the final FMA result. Keeping this assembly logic
separate makes it reusable for future manual or chunked FMA workflows without
adding exported API.

Method results store structured timing in a `timing` list with separate
`fit`, `ci`, and `total` phases. Printed runtime summaries use CPU time so they
do not count time spent while the machine is asleep. The legacy `runtime`
string is still kept for compatibility and mirrors the current total CPU time.

## Profile Likelihood Path

`confint.amerasfit()` dispatches interval calculations by type. For
profile-likelihood intervals, `compute_proflik_CI()` reconstructs the objective
function from the stored `amerasfit` object and the resolved data.

That reconstruction goes through two adapter helpers:

- `make_loglik_fn()` adds method-specific behavior. MCML averages over all dose
  realizations; RC and ERC evaluate at `rcdose_ameras`.
- `make_base_loglik_fn()` reconstructs a likelihood closure from the stored
  model metadata. Ordinary cases delegate to `make_single_realization_loglik()`.
  Poisson ERC and proportional hazards ERC stay explicit because they use
  precomputed centered dose-realization residuals.

The reconstructed likelihood is then passed into `proflik()` through
`compute_proflik_ci_one()`. For MCML profile likelihoods,
`make_loglik_fn()` reuses `mcml_average_neg_loglik()` to average over the dose
realizations.

After each method's interval calculation, `confint.amerasfit()` records the CI
timing for that method and updates the method's total runtime. When intervals
are recomputed with `force = TRUE`, the previous CI timing is replaced rather
than accumulated. The updated `amerasfit` object is returned invisibly.

## Diagnostics Data Path

For fitted objects that do not store data (`keep.data = FALSE`), diagnostic
methods reconstruct the analysis data with `resolve_data()`. That helper
validates user-supplied data, rebuilds expanded `X` columns when needed, and
adds the internal mean-dose column `rcdose_ameras`.

`residuals.amerasfit()` resolves data and then delegates residual calculation
to `compute_residuals()`, which assumes the data have already been resolved.
`plot.amerasfit()` also resolves data once, then calls `compute_fitted()` and
`compute_residuals()` directly. This avoids re-validating an internally
augmented data frame as if it were fresh user input, while preserving the
reserved-column check for data supplied through public entry points.

## Special Cases

Poisson ERC and proportional hazards ERC are intentionally special in
`make_base_loglik_fn()`. They call `loglik.poisson.erc()` and
`loglik.prophaz.erc()` with precomputed `Xc` and `Kmat_diag` values. This avoids
recomputing realization-level residual variance terms during profile likelihood
evaluation.

Conditional logistic likelihoods need risk-set structures. Fitted-object
reconstruction builds those from the supplied `data`, not from
`object$model$data`, so profile likelihoods also work for objects fit with
`keep.data = FALSE` when the user supplies data to `confint()`.

Transform and likelihood arguments supplied through `...` are stored on the
fit object as `other.args`. Fitted-object likelihood reconstruction forwards
`other.args` so profile-likelihood calculations use the same transform behavior
as the original fit.

## Compiled Helpers

The risk-set and ERC correction helpers in `src/ameras.cpp` are exposed to R
through Rcpp-generated `.Call` wrappers in `R/RcppExports.R` and
`src/RcppExports.cpp`. The package registers those `.Call` entry points in
`src/ameras_init.c`.

Do not add a parallel `.C` wrapper unless there is a specific reason to keep a
second native interface. If a new `// [[Rcpp::export]]` helper is added, run
`Rcpp::compileAttributes()` and make sure the generated symbol is registered in
`src/ameras_init.c`.

## Extension Checklist

When adding or changing a family or method:

- Add or update the low-level family likelihood first.
- Teach `make_single_realization_loglik()` how to build the ordinary
  single-realization closure.
- If MCML should support the family, confirm that `make_mcml_loglik_fn()` can
  use the same closure over all dose realizations.
- If profile likelihood intervals should support the method, confirm
  `make_base_loglik_fn()` can reconstruct the objective from a fitted object.
- If the method uses `optim()` or `optimize()`, route it through
  `fit_objective_with_hessian()`.
- If the method returns MCML-like or RC-like output, package it through
  `assemble_frequentist_fit_result()`.
- If a workflow combines fitted FMA realization summaries, use
  `assemble_fma_result()` rather than duplicating the FMA weighting and sampling
  logic.
- Add focused tests for parameter names, optimizer output shape, transform
  argument forwarding, and profile-likelihood reconstruction.
