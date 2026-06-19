test_that("BMA requires at least two included dose realizations", {
  data("data", package = "ameras")

  # The formula has multiple dose realizations, but this explicit BMA subset
  # would leave BMA with only one realization. That case should use RC instead.
  expect_error(
    suppressMessages(
      ameras(
        Y.gaussian ~ dose(V1:V2, deg = 1),
        data = data[seq_len(20), ],
        family = "gaussian",
        methods = "BMA",
        included.realizations.BMA = 1
      )
    ),
    "BMA requires at least two included exposure realizations",
    fixed = TRUE
  )
})

test_that("one-chain BMA fits support amerasfit methods", {
  skip_on_cran()
  set.seed(123)

  data("data", package = "ameras")
  small_data <- data[seq_len(30), ]

  # A one-chain BMA fit exercises a different output shape than the usual
  # multi-chain fit: samples are returned as one matrix instead of a list of
  # chain matrices, and convergence diagnostics cannot be estimated.
  fit <- suppressWarnings(
    suppressMessages(
      ameras(
        Y.gaussian ~ dose(V1:V2, deg = 1),
        data = small_data,
        family = "gaussian",
        methods = "BMA",
        nchains.BMA = 1,
        niter.BMA = 30,
        nburnin.BMA = 10,
        thin.BMA = 1,
        print = FALSE
      )
    )
  )

  expect_s3_class(fit, "amerasfit")
  expect_true(all(c(
    "call", "formula", "num.rows", "num.realizations", "transform",
    "transform.jacobian", "other.args", "model", "CI.computed", "BMA"
  ) %in% names(fit)))
  expect_named(fit$BMA$coefficients, c("(Intercept)", "dose", "sigma"))
  expect_true(is.matrix(fit$BMA$samples))
  expect_identical(
    colnames(fit$BMA$samples),
    c("(Intercept)", "dose", "sigma", "col.ind")
  )

  # With one MCMC chain, Rhat and effective sample size are intentionally NA.
  diagnostics <- Rhat(fit)
  expect_s3_class(diagnostics, "data.frame")
  expect_named(diagnostics, c("Rhat", "n.eff"))
  expect_true(all(is.na(diagnostics$Rhat)))
  expect_true(all(is.na(diagnostics$n.eff)))

  expect_identical(included_realizations(fit, methods = "BMA"), 1:2)
  expect_identical(rownames(coef(fit)), names(fit$BMA$coefficients))
  expect_identical(rownames(vcov(fit, methods = "BMA")), names(fit$BMA$coefficients))

  summary_table <- summary(fit)$summary_table
  expect_identical(summary_table$Method, rep("BMA", 3))
  expect_identical(summary_table$Term, names(fit$BMA$coefficients))
  expect_true(all(is.na(summary_table$Rhat)))
  expect_true(all(is.na(summary_table$n.eff)))

  # Percentile intervals should work directly from the one-chain sample matrix.
  fit_ci <- suppressMessages(confint(fit, type = "percentile", print = FALSE))
  expect_true(fit_ci$CI.computed)
  expect_named(fit_ci$BMA$CI, c("lower", "upper"))
  expect_identical(rownames(fit_ci$BMA$CI), names(fit$BMA$coefficients))
  expect_true(all(is.finite(as.matrix(fit_ci$BMA$CI))))

  captured <- NULL
  testthat::local_mocked_bindings(
    MCMCtrace = function(object, ...) {
      captured <<- c(list(object = object), list(...))
      invisible("trace-ok")
    },
    .package = "ameras"
  )

  trace_result <- traceplot(
    fit,
    iter = 5,
    Rhat = FALSE,
    n.eff = FALSE,
    pdf = FALSE
  )

  expect_identical(trace_result, "trace-ok")
  expect_identical(captured$object, fit$BMA$samples)
  expect_identical(captured$iter, 5)
  expect_false(captured$Rhat)
  expect_false(captured$n.eff)
  expect_false(captured$pdf)
})

test_that("one-parameter BMA fits preserve sample matrix dimensions", {
  skip_on_cran()
  set.seed(123)

  data("data", package = "ameras")
  small_data <- data[seq_len(20), ]

  # With clogit, linear ERR, and no X or M terms, the only sampled model
  # parameter summarized by ameras is dose. Keep this real fit to guard against
  # R dropping one-column matrices before colMeans(), sd(), and var() summaries.
  fit <- suppressWarnings(
    suppressMessages(
      ameras(
        Y.clogit ~ dose(V1:V2, model = ERR, deg = 1) + strata(setnr),
        data = small_data,
        family = "clogit",
        methods = "BMA",
        niter.BMA = 20,
        nburnin.BMA = 5,
        thin.BMA = 1,
        nchains.BMA = 2,
        print = FALSE
      )
    )
  )

  expect_named(fit$BMA$coefficients, "dose")
  expect_named(fit$BMA$sd, "dose")
  expect_identical(dim(fit$BMA$vcov), c(1L, 1L))
  expect_identical(rownames(fit$BMA$vcov), "dose")
  expect_identical(colnames(fit$BMA$vcov), "dose")
  expect_true(all(vapply(fit$BMA$samples, is.matrix, logical(1))))
  expect_true(all(vapply(fit$BMA$samples, ncol, integer(1)) == 2L))
  expect_identical(colnames(fit$BMA$samples[[1]]), c("dose", "col.ind"))

  # CI computation should preserve the single model parameter name after
  # stacking chains and dropping the sampled realization indicator.
  fit_percentile <- suppressMessages(
    confint(fit, type = "percentile", print = FALSE)
  )
  expect_identical(rownames(fit_percentile$BMA$CI), "dose")

  fit_hpd <- suppressMessages(confint(fit, type = "hpd", print = FALSE))
  expect_identical(rownames(fit_hpd$BMA$CI), "dose")
})
