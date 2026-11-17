test_that("HPD intervals honor their level for FMA and BMA sample shapes", {
  samples <- cbind(dose = qnorm(seq(0.001, 0.999, length.out = 1001)))
  for (method in c("FMA", "BMA")) {
    for (multiple_chains in c(FALSE, TRUE)) {
      method_fit <- list(coefficients = c(dose = 0), samples = samples,
                         timing = ameras:::new_method_timing())
      if (method == "BMA" && multiple_chains) {
        method_fit$samples <- list(samples[1:500, , drop = FALSE], samples[501:1001, , drop = FALSE])
      }
      fit <- ameras:::new_amerasfit(setNames(list(method_fit), method))
      for (level in c(0.8, 0.9, 0.95)) {
        expected <- coda::HPDinterval(coda::as.mcmc(samples), prob = level)
        bounds <- confint(fit, type = "hpd", level = level)
        attached <- add_confint(fit, type = "hpd", level = level)
        expect_equal(unname(bounds), unname(expected), ignore_attr = TRUE)
        expect_equal(as.numeric(attached[[method]]$CI[1, ]), as.numeric(expected[1, ]))
      }
    }
  }
})

test_that("profile calculations preserve explicit and effective default controls", {
  dat <- transform_contract_data()
  default <- suppressMessages(ameras(Y ~ dose(D1:D2), data = dat))
  explicit <- suppressMessages(ameras(Y ~ dose(D1:D2), data = dat, control = list(reltol = 1e-10)))
  empty <- suppressMessages(ameras(Y ~ dose(D1:D2), data = dat, control = list()))
  expect_identical(default$model$control, list(reltol = 1e-10))
  expect_identical(empty$model$control, list())
  expected <- suppressMessages(confint(explicit, type = "proflik"))
  expect_equal(suppressMessages(confint(default, type = "proflik")), expected)
  for (missing in c(FALSE, TRUE)) {
    old <- default
    old$model$control <- NULL
    if (!missing) old$model["control"] <- list(NULL)
    expect_equal(suppressMessages(confint(old, type = "proflik")), expected)
  }

  seen <- list()
  original <- ameras:::proflik
  testthat::local_mocked_bindings(
    proflik = function(..., control) {
      seen <<- c(seen, list(control))
      original(..., control = control)
    }, .package = "ameras"
  )
  suppressMessages(confint(empty, type = "proflik"))
  expect_gt(length(seen), 0)
  expect_true(all(vapply(seen, identical, logical(1), list())))
})
