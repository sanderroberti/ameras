# These tests target the second half of FMA: turning already fitted
# realization-level summaries into weights, sampled coefficients, and the final
# amerasfit-compatible FMA result. Keeping these tests on fake fits makes the
# assembly logic fast to exercise and independent of the optimizer.

make_fma_summary <- function(
  coef,
  AIC,
  hess = diag(length(coef)),
  include = TRUE,
  transform = NULL,
  transform.jacobian = NULL,
  ...
) {
  # Build the same realization-level summary produced by fit_fma_realizations().
  # AIC is easier to reason about in the tests, so convert it back to the
  # minimized objective value expected by summarize_fma_realization_fit().
  npar <- length(coef)
  fit <- list(
    par = coef,
    hessian = hess,
    value = (AIC - 2 * npar) / 2,
    convergence = if (include) 0 else 1
  )

  ameras:::summarize_fma_realization_fit(
    fit,
    npar,
    transform = transform,
    transform.jacobian = transform.jacobian,
    ...
  )
}

capture_fma_messages <- function(expr) {
  messages <- character()
  value <- withCallingHandlers(
    expr,
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  list(value = value, messages = messages)
}

capture_fma_warnings <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = warnings)
}

test_that("assemble_fma_result computes weights and sample allocations", {
  # The third fit has a much worse AIC. It should receive zero simulated
  # samples for this MFMA value and be excluded after weight calculation.
  fits <- list(
    make_fma_summary(c(0, 1), AIC = 10),
    make_fma_summary(c(1, 2), AIC = 12),
    make_fma_summary(c(2, 3), AIC = 100)
  )
  parnames <- c("alpha", "dose")
  expected_weights <- exp(-0.5 * c(0, 2, 90))
  expected_weights <- expected_weights / sum(expected_weights)

  set.seed(100)
  out <- NULL
  captured <- capture_fma_messages(
    expect_warning(
      out <- ameras:::assemble_fma_result(
        FMAfits = fits,
        dosevars = c("V1", "V2", "V3"),
        parnames = parnames,
        inpar = c(0, 0),
        MFMA = 100,
        unweighted = FALSE,
        t0 = proc.time()
      ),
      "FMA is based on only 2 realizations",
      fixed = TRUE
    )
  )
  expect_true(any(grepl(
    "excluded due to negligible model averaging weight",
    captured$messages,
    fixed = TRUE
  )))
  expect_identical(out$included.realizations, c(1L, 2L))
  expect_equal(out$included.samples, sum(round(expected_weights[1:2] * 100)))
  expect_equal(unname(out$weights), expected_weights[1:2], tolerance = 1e-12)
  expect_named(out$weights, c("V1", "V2"))
  expect_equal(nrow(out$samples), out$included.samples)
  expect_identical(names(out$samples), parnames)
  expect_identical(names(out$coefficients), parnames)
  expect_identical(rownames(out$vcov), parnames)
  expect_identical(colnames(out$vcov), parnames)
  expect_match(out$runtime, "seconds$")
})

test_that("assemble_fma_result excludes non-finite sampling covariance", {
  # A finite Hessian can still produce a non-finite FMA sampling covariance
  # after applying a transformation Jacobian. FMA should drop that realization
  # before rmvnorm() rather than exposing a low-level eigen() error from mvtnorm.
  bad_jacobian <- function(params, ...) matrix(Inf, length(params), length(params))
  fits <- list(
    make_fma_summary(
      c(dose = 0),
      AIC = 10,
      hess = matrix(1, 1, 1),
      transform = function(params, ...) params,
      transform.jacobian = bad_jacobian
    ),
    make_fma_summary(c(dose = 1), AIC = 12, hess = matrix(1, 1, 1))
  )

  expect_identical(fits[[1]]$exclude.reason, "sampling")
  expect_identical(fits[[2]]$exclude.reason, NULL)

  set.seed(102)
  captured <- capture_fma_warnings(
    ameras:::assemble_fma_result(
      FMAfits = fits,
      dosevars = c("V1", "V2"),
      parnames = "dose",
      inpar = 0,
      MFMA = 20,
      unweighted = FALSE,
      t0 = proc.time()
    )
  )

  out <- captured$value
  expect_true(any(grepl(
    "excluded because FMA sampling means or covariance matrices",
    captured$warnings,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "FMA is based on a single realization",
    captured$warnings,
    fixed = TRUE
  )))
  expect_identical(out$included.realizations, 2L)
  expect_equal(unname(out$weights), 1)
  expect_equal(out$included.samples, 20)
  expect_named(out$coefficients, "dose")
  expect_true(all(is.finite(out$samples$dose)))
})

test_that("assemble_fma_result handles zero allocated samples", {
  # With unweighted FMA and MFMA = 1, each of three valid fits receives
  # round(1 / 3) = 0 samples. The output should mirror the degenerate FMA
  # result currently returned by ameras.fma().
  fits <- list(
    make_fma_summary(c(0, 1), AIC = 10),
    make_fma_summary(c(1, 2), AIC = 10),
    make_fma_summary(c(2, 3), AIC = 10)
  )
  parnames <- c("alpha", "dose")

  out <- NULL
  captured <- capture_fma_messages(
    expect_warning(
      out <- ameras:::assemble_fma_result(
        FMAfits = fits,
        dosevars = c("V1", "V2", "V3"),
        parnames = parnames,
        inpar = c(0, 0),
        MFMA = 1,
        unweighted = TRUE,
        t0 = proc.time()
      ),
      "No realizations contributed to FMA",
      fixed = TRUE
    )
  )
  expect_true(any(grepl(
    "3 realization(s) excluded due to negligible model averaging weight",
    captured$messages,
    fixed = TRUE
  )))
  expect_identical(out$included.realizations, integer(0))
  expect_equal(out$included.samples, 0)
  expect_null(out$weights)
  expect_null(out$samples)
  expect_named(out$coefficients, parnames)
  expect_true(all(is.na(out$coefficients)))
  expect_true(all(is.na(out$sd)))
  expect_true(all(is.na(out$vcov)))
  expect_identical(rownames(out$vcov), parnames)
  expect_identical(colnames(out$vcov), parnames)
})

test_that("FMA summaries forward transform arguments for sampling", {
  # The marker arguments make sure the transform and Jacobian used to prepare
  # FMA sampling inputs receive the same extra arguments supplied through
  # ameras().
  transform_shift <- function(params, marker, shift, ...) {
    if (!identical(marker, "fma-assembler-ok")) {
      stop("marker was not forwarded")
    }
    params + shift
  }
  jacobian_scale <- function(params, marker, scale, ...) {
    if (!identical(marker, "fma-assembler-ok")) {
      stop("marker was not forwarded")
    }
    diag(scale, length(params))
  }
  fits <- list(
    make_fma_summary(
      c(0, 1),
      AIC = 10,
      hess = diag(c(4, 9)),
      transform = transform_shift,
      transform.jacobian = jacobian_scale,
      marker = "fma-assembler-ok",
      shift = c(10, 20),
      scale = c(2, 3)
    ),
    make_fma_summary(
      c(1, 2),
      AIC = 10,
      hess = diag(c(4, 9)),
      transform = transform_shift,
      transform.jacobian = jacobian_scale,
      marker = "fma-assembler-ok",
      shift = c(10, 20),
      scale = c(2, 3)
    )
  )
  parnames <- c("alpha", "dose")

  set.seed(101)
  out <- suppressWarnings(
    ameras:::assemble_fma_result(
      FMAfits = fits,
      dosevars = c("V1", "V2"),
      parnames = parnames,
      inpar = c(0, 0),
      MFMA = 20,
      unweighted = TRUE,
      t0 = proc.time()
    )
  )

  expect_identical(out$included.realizations, c(1L, 2L))
  expect_equal(out$included.samples, 20)
  expect_identical(names(out$samples), parnames)
  expect_true(all(is.finite(as.matrix(out$samples))))
  expect_true(all(is.finite(out$coefficients)))
  expect_true(all(is.finite(out$vcov)))
})
