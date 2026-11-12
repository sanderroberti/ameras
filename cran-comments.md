## R CMD check results

Final checks of the submission tarball are pending.

<!-- Before submission, replace this paragraph with the final --as-cran
results, the tested R versions/platforms, and explanations of any remaining
NOTEs. Do not reuse results from an earlier package version. -->

## Reverse dependencies

The [CRAN package page](https://CRAN.R-project.org/package=ameras) lists no
reverse dependencies as of 2026-11-11.


## Changes in this version

This major release updates CRAN version 0.5.2. It retires deprecated
interfaces and makes deliberate API changes to align coefficient and
confidence-interval access with standard R conventions. Full details and
migration guidance are provided in `NEWS.md` and the updated documentation.

### Deprecated interfaces made defunct

* The direct argument interface and confidence-interval arguments in
  `ameras()`, deprecated since version 0.2.0, now produce errors with
  migration guidance. Models must use the formula interface; intervals can
  be computed after fitting with `confint()` or attached with `add_confint()`.
* The `included.replicates.BMA` alias, deprecated since version 0.3.0, now
  directs users to `included.realizations.BMA`.
* Non-formula effect-modifier syntax, deprecated since version 0.5.0, now
  directs users to `modifier = ~ M1 + M2` inside `dose()`.

Removal of the direct interface, confidence-interval arguments, and BMA alias
was previously announced for version 1.0.0. Migration errors are retained so
that old calls fail explicitly rather than being silently reinterpreted.

### Additional breaking changes

* `confint()` now returns interval bounds as a matrix for one fitted method
  or a named list of matrices for multiple methods. Code that previously used
  `fit <- confint(fit, ...)` should use `fit <- add_confint(fit, ...)` to
  attach intervals, diagnostics, and computation timing. Both functions
  compute the requested intervals on every call. The `force` argument has
  been removed; `print` and `digits` are supported by `add_confint()` only.
* `coef(fit)` now returns a named numeric vector for a single fitted method.
  The data-frame return for multiple methods is unchanged.
* The coefficient table in `summary(fit)` is now named `coefficients` rather
  than `summary_table`. The `summary_table()` accessor has been removed;
  use `coef(summary(fit))` instead. Printed summaries are unchanged.
* The bundled simulated dataset has been renamed from `data` to
  `ameras_simulated` to give it a descriptive, package-specific name.

These additional changes are introduced directly in this major release,
rather than being presented as the removal of previously deprecated features.
The interval-interface refactoring uses a shared computation path and does
not change the underlying interval calculations.

The release also adds `nobs()`, strengthens formula validation and optimizer
diagnostics, and expands tests comparing equivalent RC models with `lm()`,
`glm()`, `coxph()`, `clogit()`, and `multinom()`. The README, help examples,
and vignette sources have been updated for the revised interfaces.
