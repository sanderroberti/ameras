# Fast tests for small amerasfit accessors that do not need model fitting.

test_that("public methods warn about unused dots", {
  timing <- ameras:::new_method_timing(
    fit = ameras:::new_runtime_phase(cpu = 0)
  )
  ci <- data.frame(lower = 0.8, upper = 1.2, row.names = "dose")
  method_fit <- list(
    coefficients = c(dose = 1),
    sd = c(dose = 0.1),
    vcov = matrix(0.01, dimnames = list("dose", "dose")),
    CI = ci,
    timing = timing
  )
  fit <- ameras:::new_amerasfit(list(
    call = quote(ameras(Y ~ dose(D), data = dat)),
    num.rows = 1L,
    num.rows.original = 1L,
    num.realizations = 1L,
    CI.computed = TRUE,
    model = list(family = "gaussian"),
    RC = method_fit,
    FMA = c(method_fit, list(included.realizations = 1L)),
    BMA = c(
      method_fit,
      list(
        included.realizations = 1L,
        Rhat = data.frame(Rhat = 1, n.eff = 100, row.names = "dose")
      )
    )
  ))
  count_unused_dots_warnings <- function(expr) {
    n <- 0L
    withCallingHandlers(
      force(expr),
      warning = function(cnd) {
        if (grepl("extra argument", conditionMessage(cnd), fixed = TRUE)) {
          n <<- n + 1L
          invokeRestart("muffleWarning")
        }
      }
    )
    n
  }
  expect_unused_dots_warning <- function(expr) {
    expect_warning(expr, "extra argument")
  }

  expect_unused_dots_warning(vcov(fit, typo = TRUE))
  expect_unused_dots_warning(coef(fit, typo = TRUE))
  expect_unused_dots_warning(summary(fit, typo = TRUE))
  expect_identical(
    count_unused_dots_warnings(summary_table(fit, typo = TRUE)),
    1L
  )
  expect_unused_dots_warning(Rhat(fit, typo = TRUE))
  expect_unused_dots_warning(included_realizations(fit, typo = TRUE))
  expect_unused_dots_warning(
    capture.output(suppressMessages(confint(fit, print = FALSE, typo = TRUE)))
  )
  expect_identical(
    count_unused_dots_warnings(capture.output(print(fit, typo = TRUE))),
    1L
  )
  expect_identical(
    count_unused_dots_warnings(capture.output(print(summary(fit), typo = TRUE))),
    1L
  )

  expect_unused_dots_warning(
    expect_error(
      convergence(fit, methods = "FMA", typo = TRUE),
      "available for RC, ERC, and MCML"
    )
  )
  expect_unused_dots_warning(
    expect_error(
      residuals(fit, type = "schoenfeld", typo = TRUE),
      "schoenfeld residuals not supported"
    )
  )
  expect_unused_dots_warning(
    expect_error(
      dose_lrt(fit, methods = "FMA", typo = TRUE),
      "supports RC, ERC, and MCML"
    )
  )
})

test_that("Rhat returns BMA convergence diagnostics", {
  diagnostics <- data.frame(
    Rhat = c(1.01, 1.02),
    n.eff = c(100, 80),
    row.names = c("dose", "dose_squared")
  )
  fit <- ameras:::new_amerasfit(list(
    BMA = list(Rhat = diagnostics)
  ))

  expect_identical(Rhat(fit), diagnostics)
  expect_error(
    Rhat(ameras:::new_amerasfit(list(RC = list()))),
    "BMA not present"
  )
})

test_that("included_realizations returns available FMA and BMA indices", {
  fit <- ameras:::new_amerasfit(list(
    FMA = list(included.realizations = c(1L, 3L)),
    BMA = list(included.realizations = c(2L, 4L))
  ))

  both <- included_realizations(fit)
  expect_named(both, c("FMA", "BMA"))
  expect_identical(both$FMA, c(1L, 3L))
  expect_identical(both$BMA, c(2L, 4L))

  expect_identical(
    included_realizations(fit, methods = "FMA"),
    c(1L, 3L)
  )
  expect_identical(
    included_realizations(fit, methods = "BMA"),
    c(2L, 4L)
  )
})

test_that("included_realizations errors when requested methods are absent", {
  fit <- ameras:::new_amerasfit(list(RC = list()))

  expect_error(
    included_realizations(fit, methods = "FMA"),
    "None of the requested methods were run"
  )
})

test_that("vcov reports when requested methods are absent", {
  fit <- ameras:::new_amerasfit(list(
    RC = list(vcov = matrix(1, dimnames = list("dose", "dose")))
  ))

  expect_error(
    vcov(fit, methods = "FMA"),
    "None of the requested methods were run"
  )
})

test_that("row-info printing omits zero-count exclusions", {
  row_info <- list(
    supplied = 10L,
    omitted.na = 0L,
    clogit.uninformative.rows = NULL,
    prophaz.zero.followup.rows = NULL,
    used = 10L
  )

  out <- capture.output(ameras:::print_row_info(row_info, detailed = TRUE))

  expect_equal(out[nzchar(out)], "Rows: 10")
  expect_false(any(grepl("Supplied:", out, fixed = TRUE)))
  expect_false(any(grepl("Omitted by na.action", out, fixed = TRUE)))
  expect_false(any(grepl("Excluded as", out, fixed = TRUE)))
  expect_false(any(grepl("Used for fitting", out, fixed = TRUE)))
})

test_that("traceplot requires BMA output", {
  fit <- ameras:::new_amerasfit(list(RC = list()))

  expect_error(
    traceplot(fit),
    "traceplot\\(\\) requires BMA output"
  )
})

test_that("traceplot forwards BMA samples and plotting arguments", {
  samples <- list(
    chain1 = matrix(1:6, ncol = 2, dimnames = list(NULL, c("a", "b")))
  )
  fit <- ameras:::new_amerasfit(list(
    BMA = list(samples = samples)
  ))
  captured <- NULL

  testthat::local_mocked_bindings(
    MCMCtrace = function(object, ...) {
      captured <<- c(list(object = object), list(...))
      invisible("trace-ok")
    },
    .package = "ameras"
  )

  out <- traceplot(
    fit,
    iter = 10,
    Rhat = FALSE,
    n.eff = FALSE,
    pdf = TRUE,
    ind = TRUE
  )

  expect_identical(out, "trace-ok")
  expect_identical(captured$object, samples)
  expect_identical(captured$iter, 10)
  expect_false(captured$Rhat)
  expect_false(captured$n.eff)
  expect_true(captured$pdf)
  expect_true(captured$ind)
})

test_that("confint validates method-specific inputs and empty sample output", {
  empty_fit <- ameras:::new_amerasfit(list(CI.computed = FALSE))
  fma_fit <- ameras:::new_amerasfit(list(
    CI.computed = FALSE,
    FMA = list(
      coefficients = c(dose = 1),
      samples = matrix(
        seq(0.8, 1.2, length.out = 10),
        ncol = 1,
        dimnames = list(NULL, "dose")
      ),
      timing = ameras:::new_method_timing()
    )
  ))
  rc_fit <- ameras:::new_amerasfit(list(
    CI.computed = FALSE,
    transform = NULL,
    other.args = list(),
    model = list(family = "gaussian"),
    RC = list(
      coefficients = c(dose = 1),
      sd = c(dose = 0.1),
      vcov = matrix(0.01, dimnames = list("dose", "dose")),
      optim = list(par = c(dose = 1), hessian = matrix(100)),
      loglik = 0,
      timing = ameras:::new_method_timing()
    )
  ))

  # confint() has separate allowed interval types for sampling-based methods
  # (FMA/BMA) and likelihood/optimization-based methods (RC/ERC/MCML).
  expect_error(
    confint(empty_fit, level = 1, print = FALSE),
    "level must be a single numeric value between 0 and 1"
  )
  expect_error(
    confint(fma_fit, type = c("percentile", "hpd"), print = FALSE),
    "Exactly one CI type .* for FMA"
  )
  expect_error(
    confint(rc_fit, type = c("wald.orig", "proflik"), print = FALSE),
    "Exactly one CI type .* for RC"
  )

  # If all FMA realizations were excluded upstream, confint() should warn and
  # return NA interval bounds rather than failing on the missing sample matrix.
  no_samples_fit <- ameras:::new_amerasfit(list(
    CI.computed = FALSE,
    FMA = list(
      coefficients = c(dose = 1),
      samples = NULL,
      timing = ameras:::new_method_timing()
    )
  ))
  fit_ci <- NULL
  expect_warning(
    fit_ci <- confint(no_samples_fit, type = "percentile", print = FALSE),
    "No samples available"
  )

  expect_true(fit_ci$CI.computed)
  expect_identical(rownames(fit_ci$FMA$CI), "dose")
  expect_true(all(is.na(fit_ci$FMA$CI$lower)))
  expect_true(all(is.na(fit_ci$FMA$CI$upper)))
})

test_that("residuals report incompatible residual requests before fitting work", {
  gaussian_fit <- ameras:::new_amerasfit(list(
    model = list(family = "gaussian"),
    RC = list()
  ))
  prophaz_fit <- ameras:::new_amerasfit(list(
    model = list(family = "prophaz"),
    RC = list()
  ))

  # These checks happen before data reconstruction or fitted-value calculation,
  # so small synthetic objects are enough to lock down the user-facing errors.
  expect_error(
    residuals(gaussian_fit, type = "schoenfeld"),
    "schoenfeld residuals not supported for family='gaussian'"
  )
  expect_error(
    residuals(prophaz_fit, type = "pearson"),
    "Only schoenfeld residuals are supported for family 'prophaz'"
  )
  expect_error(
    residuals(gaussian_fit, method = "FMA"),
    "Method 'FMA' not present in fitted object"
  )

  gaussian_data_fit <- ameras:::new_amerasfit(list(
    num.rows = 3L,
    model = list(
      family = "gaussian",
      dosevars = c("D1", "D2"),
      X = NULL,
      M = NULL,
      M_names = NULL,
      Y = "Y",
      offset = NULL,
      doseRRmod = NULL,
      deg = 1
    ),
    RC = list(coefficients = c("(Intercept)" = 0, dose = 1, sigma = 1))
  ))
  dat <- data.frame(
    Y = c(1, 2, 3),
    D1 = c(0.1, 0.2, 0.3),
    D2 = c(0.2, 0.3, 0.4)
  )

  expect_error(
    residuals(gaussian_data_fit, data = dat, dose.col = "missing_dose"),
    "dose.col 'missing_dose' not found in data"
  )
})

test_that("plot reports when requested methods are absent", {
  fit <- ameras:::new_amerasfit(list(
    model = list(family = "gaussian"),
    RC = list(coefficients = c(dose = 1))
  ))

  expect_error(
    plot(fit, methods = "FMA"),
    "None of the requested methods were run"
  )
})
