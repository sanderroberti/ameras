test_that("prophaz BMA returns the documented baseline-hazard cutpoints", {
  skip_on_cran()
  set.seed(123)

  small_data <- ameras_simulated[seq_len(30), ]
  numints <- 3L
  fit <- suppressWarnings(
    suppressMessages(
      ameras(
        Surv(time, status) ~ dose(V1:V2, model = "EXP"),
        data = small_data,
        family = "prophaz",
        methods = "BMA",
        prophaz.numints.BMA = numints,
        nchains.BMA = 1,
        niter.BMA = 30,
        nburnin.BMA = 10,
        thin.BMA = 1
      )
    )
  )

  timepoints <- fit$BMA$prophaz.timepoints
  expect_type(timepoints, "double")
  expect_length(timepoints, numints + 1L)
  expect_true(all(is.finite(timepoints)))
  expect_false("prophaz_timepoints" %in% names(fit$BMA))
})
