# These tests cover the future-aware FMA realization loop. Most checks avoid
# parallel workers so they stay portable; optional smoke tests exercise
# multisession workers when the future packages are installed off CRAN.

test_that("FMA realization mapper preserves lapply output shape", {
  out <- ameras:::fma_realization_lapply(
    1:4,
    function(i, add) i + add,
    add = 10,
    future.chunk.size.FMA = 2
  )

  expect_equal(out, lapply(1:4, function(i) i + 10))
})

test_that("FMA realization mapper can use a multisession future plan", {
  # This is a lightweight parallel smoke test. It is skipped on CRAN and when
  # the optional future packages are unavailable, but it lets CI exercise the
  # actual worker path when Suggests dependencies are installed.
  skip_on_cran()
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multisession, workers = 2)

  out <- ameras:::fma_realization_lapply(
    1:4,
    function(i, add) i + add,
    add = 10,
    future.chunk.size.FMA = 1
  )

  expect_equal(out, lapply(1:4, function(i) i + 10))
})

test_that("FMA fitting works with a multisession future plan", {
  # This catches the main parallelization risk: worker processes must be able
  # to see ameras' internal likelihood and optimizer helpers.
  skip_on_cran()
  skip_if_not_installed("future")
  skip_if_not_installed("future.apply")

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multisession, workers = 2)

  set.seed(121)
  n <- 50
  dose <- seq(-1, 1, length.out = n)
  dat <- data.frame(
    Y = 1 + dose + rnorm(n, sd = 0.2),
    V1 = dose + rnorm(n, sd = 0.02),
    V2 = dose + rnorm(n, sd = 0.02)
  )

  fit <- suppressWarnings(
    suppressMessages(
      ameras(
        Y ~ dose(V1:V2),
        data = dat,
        family = "gaussian",
        methods = "FMA",
        MFMA = 40,
        future.chunk.size.FMA = 1
      )
    )
  )

  expect_true(length(fit$FMA$included.realizations) > 0)
  expect_named(fit$FMA$coefficients, c("(Intercept)", "dose", "sigma"))
})

test_that("FMA chunk size does not change sequential fit results", {
  # Use several well-behaved Gaussian dose realizations so all realization fits
  # are valid. Resetting the seed before each call makes the FMA sampling step
  # comparable after the realization-level fits are assembled.
  set.seed(12)
  n <- 80
  dose <- seq(-1, 1, length.out = n)
  dat <- data.frame(
    Y = 1 + 0.5 * dose + rnorm(n, sd = 0.2),
    V1 = dose + rnorm(n, sd = 0.02),
    V2 = dose + rnorm(n, sd = 0.02),
    V3 = dose + rnorm(n, sd = 0.02)
  )

  fit_fma <- function(chunk_size) {
    set.seed(120)
    suppressWarnings(
      suppressMessages(
        ameras(
          Y ~ dose(V1:V3),
          data = dat,
          family = "gaussian",
          methods = "FMA",
          unweightedFMA = TRUE,
          MFMA = 120,
          future.chunk.size.FMA = chunk_size
        )
      )
    )
  }

  default_fit <- fit_fma(NULL)
  chunked_fit <- fit_fma(1)

  expect_identical(
    chunked_fit$FMA$included.realizations,
    default_fit$FMA$included.realizations
  )
  expect_equal(chunked_fit$FMA$weights, default_fit$FMA$weights)
  expect_equal(chunked_fit$FMA$coefficients, default_fit$FMA$coefficients)
  expect_equal(chunked_fit$FMA$sd, default_fit$FMA$sd)
  expect_equal(chunked_fit$FMA$vcov, default_fit$FMA$vcov)
  expect_equal(chunked_fit$FMA$samples, default_fit$FMA$samples)
})
