# Parallel FMA with future

``` r

library(ameras)
```

``` r

data(data, package = "ameras")
```

## Introduction

Frequentist model averaging (`FMA`) fits the specified model once for
each dose realization, computes AIC weights across those
realization-specific fits, and then samples from the fitted normal
approximations. The realization-specific fits are independent of one
another, so they are a natural place to use parallel computation.

[`ameras()`](https://ameras.sanderroberti.com/reference/ameras.md)
supports optional parallel FMA through the
[`future`](https://future.futureverse.org/reference/plan.html) framework
and
[`future.apply`](https://future.apply.futureverse.org/reference/future_lapply.html).
The package does not set a future plan internally. Instead, the user
decides whether FMA should run sequentially or in parallel by setting
[`future::plan()`](https://future.futureverse.org/reference/plan.html)
before calling
[`ameras()`](https://ameras.sanderroberti.com/reference/ameras.md).

If `future.apply` is not installed,
[`ameras()`](https://ameras.sanderroberti.com/reference/ameras.md) falls
back to the ordinary sequential
[`lapply()`](https://rdrr.io/r/base/lapply.html) path.

For more detail on future plans and future-backed apply loops, see
[`?future::plan`](https://future.futureverse.org/reference/plan.html)
and
[`?future.apply::future_lapply`](https://future.apply.futureverse.org/reference/future_lapply.html),
or the linked Futureverse documentation.

## Choosing a future plan

The default future plan is sequential, so this call uses the ordinary
single-process FMA path.

``` r

future::plan(future::sequential)

fit <- ameras(Y.gaussian ~ dose(V1:V10) + X1 + X2,
              data = data,
              family = "gaussian",
              methods = "FMA")
```

For parallel FMA on a local machine, `multisession` is usually the most
portable choice. It starts separate background R sessions and works
across the main platforms supported by R, including Windows, macOS, and
Linux.

``` r

future::plan(future::multisession, workers = 2)

fit <- ameras(Y.gaussian ~ dose(V1:V10) + X1 + X2,
              data = data,
              family = "gaussian",
              methods = "FMA")

future::plan(future::sequential)
```

The `multicore` backend can be efficient on systems that support
forking, but it is not available on Windows. For code intended to run on
multiple platforms, prefer `multisession`.

For more specialized settings, the future ecosystem also supports
cluster and HPC-oriented backends. See
[`?future::plan`](https://future.futureverse.org/reference/plan.html) or
the future package documentation for details.

## Controlling chunk size

By default, `future.apply` chooses how to split dose realizations into
futures. For FMA, this can be adjusted with `future.chunk.size.FMA`, the
average number of dose realizations handled by each future.

Let `K` be the number of dose realizations, `W` the number of workers,
and `C` the value of `future.chunk.size.FMA`. `K` and `C` determine the
approximate number of future tasks, `ceiling(K / C)`, and at most `W` of
those tasks can run at the same time.

For example, with `K = 100` dose realizations and `W = 4` workers:

| `future.chunk.size.FMA` | Approximate chunks | Interpretation |
|---:|---:|----|
| 1 | 100 | many small tasks; more load balancing, more overhead |
| 5 | 20 | moderate load balancing and overhead |
| 25 | 4 | about one task per worker |
| 100 | 1 | one task; only one worker is used |

``` r

future::plan(future::multisession, workers = 2)

fit <- ameras(Y.gaussian ~ dose(V1:V10) + X1 + X2,
              data = data,
              family = "gaussian",
              methods = "FMA",
              future.chunk.size.FMA = 2)

future::plan(future::sequential)
```

Smaller chunks can improve load balancing when some realizations take
longer to fit than others. Larger chunks can reduce overhead and may be
preferable when each realization-specific fit is fast. If the
realization-specific fits have similar runtimes, a reasonable starting
point is a chunk size near `ceiling(K / W)` or somewhat smaller. If
runtimes vary substantially across realizations, smaller chunks may keep
workers busier.

## Example: sequential and parallel FMA

The example below compares the same Gaussian FMA fit with a sequential
plan and a two-worker multisession plan. The seed is reset before each
call so that the final FMA sampling step is comparable.

``` r

old_plan <- future::plan()
on.exit(future::plan(old_plan), add = TRUE)

future::plan(future::sequential)
set.seed(2024)
sequential_elapsed <- system.time(
  sequential_fit <- suppressWarnings(
    ameras(Y.gaussian ~ dose(V1:V10) + X1 + X2,
           data = data,
           family = "gaussian",
           methods = "FMA",
           MFMA = 10000)
  )
)[["elapsed"]]

future::plan(future::multisession, workers = 2)
set.seed(2024)
parallel_elapsed <- system.time(
  parallel_fit <- suppressWarnings(
    ameras(Y.gaussian ~ dose(V1:V10) + X1 + X2,
           data = data,
           family = "gaussian",
           methods = "FMA",
           MFMA = 10000,
           future.chunk.size.FMA = 2)
  )
)[["elapsed"]]

future::plan(old_plan)
```

Compare the coefficient estimates and elapsed runtimes.

``` r

coefficient_comparison <- data.frame(
  term = names(sequential_fit$FMA$coefficients),
  sequential = unname(sequential_fit$FMA$coefficients),
  parallel = unname(parallel_fit$FMA$coefficients),
  row.names = NULL
)

coefficient_comparison[-1] <- round(coefficient_comparison[-1], 4)
coefficient_comparison
#>          term sequential parallel
#> 1 (Intercept)    -1.2802  -1.2802
#> 2          X1     0.4829   0.4829
#> 3          X2    -0.5163  -0.5163
#> 4        dose     1.0793   1.0793
#> 5       sigma     1.1376   1.1376

data.frame(
  plan = c("sequential", "multisession"),
  elapsed_seconds = round(c(sequential_elapsed, parallel_elapsed), 2)
)
#>           plan elapsed_seconds
#> 1   sequential            1.78
#> 2 multisession            2.25
```

The estimates should agree up to ordinary numerical precision. Parallel
execution is not guaranteed to be faster for small examples because
starting workers and transferring data have overhead. It is most useful
when the realization-specific fits are slow enough for that overhead to
matter less.

## Reproducibility

The parallelized part of FMA is the realization-specific optimization
step, which is deterministic for ordinary model fits. The random FMA
sampling step is performed afterward in the main R process. For
reproducible FMA samples, set the seed before calling
[`ameras()`](https://ameras.sanderroberti.com/reference/ameras.md), just
as in the sequential case.

``` r

future::plan(future::multisession, workers = 2)

set.seed(2024)
fit <- ameras(Y.gaussian ~ dose(V1:V10) + X1 + X2,
              data = data,
              family = "gaussian",
              methods = "FMA")

future::plan(future::sequential)
```

User-supplied functions passed through `transform`,
`transform.jacobian`, or other runtime arguments should also be
deterministic. Random number generation inside those functions is not
recommended.

## Practical recommendations

Use parallel FMA when a direct `ameras(..., methods = "FMA")` call is
feasible in one R session, but the realization-specific fits are slow
enough that using multiple workers is worthwhile.

Start conservatively:

- use `future::plan(future::multisession, workers = 2)` for a portable
  local parallel plan;
- increase `workers` only after checking memory use;
- try `future.chunk.size.FMA = 2` or larger when each
  realization-specific fit is fast;
- reset the plan to
  [`future::sequential`](https://future.futureverse.org/reference/sequential.html)
  after the fit if the rest of the session should not use parallel
  futures.

With `multisession`, each worker is a separate R process. Large datasets
may therefore be copied to multiple workers. If memory is tight, use
fewer workers, increase the chunk size, or consider the manual FMA
workflow.

## Parallel FMA or manual FMA?

The parallel option shown here is best when the built-in FMA workflow
fits in memory and you want a regular `amerasfit` object returned by
[`ameras()`](https://ameras.sanderroberti.com/reference/ameras.md).

The [manual FMA
vignette](https://ameras.sanderroberti.com/articles/manualfma.md) is
intended for larger or more operationally complicated analyses where it
is useful to split the realization-specific fits across separate jobs,
save intermediate fit summaries, or resume work after partial
completion. Manual FMA is more flexible, but it requires more
bookkeeping and does not return a regular `amerasfit` object.

In short:

- use built-in parallel FMA for ordinary analyses that benefit from
  local parallel workers;
- use manual FMA when the analysis must be broken into independent
  chunks or run outside a single R session.
