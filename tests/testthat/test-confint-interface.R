make_confint_test_fit <- function(include_fma = FALSE) {
  coefficients <- c(a = 1, b = 2)
  method_fit <- list(
    coefficients = coefficients,
    sd = c(a = 0.1, b = 0.2),
    vcov = matrix(
      c(0.01, 0, 0, 0.04),
      nrow = 2,
      dimnames = list(names(coefficients), names(coefficients))
    ),
    timing = ameras:::new_method_timing()
  )
  fit <- list(
    transform = NULL,
    other.args = list(),
    RC = method_fit
  )

  if (include_fma) {
    fit$FMA <- c(
      method_fit,
      list(samples = rbind(c(a = 0.8, b = 1.6), c(a = 1.2, b = 2.4)))
    )
  }

  ameras:::new_amerasfit(fit)
}


test_that("confint returns bounds without amending a one-method fit", {
  fit <- make_confint_test_fit()
  visibility <- withVisible(confint(fit, type = "wald.orig"))
  amended <- add_confint(fit, type = "wald.orig")

  expect_true(visibility$visible)
  expect_true(is.matrix(visibility$value))
  expect_identical(dimnames(visibility$value), list(c("a", "b"), c("2.5 %", "97.5 %")))
  expect_equal(
    unname(visibility$value),
    unname(as.matrix(amended$RC$CI[, c("lower", "upper")]))
  )
  expect_null(fit$RC$CI)
})


test_that("confint returns a named list for multiple fitted methods", {
  fit <- make_confint_test_fit(include_fma = TRUE)
  result <- confint(fit, type = c("wald.orig", "percentile"))
  amended <- add_confint(
    fit,
    type = c("wald.orig", "percentile")
  )

  expect_named(result, c("RC", "FMA"))
  expect_true(all(vapply(result, is.matrix, logical(1))))
  for (method in names(result)) {
    expect_equal(
      unname(result[[method]]),
      unname(as.matrix(amended[[method]]$CI[, c("lower", "upper")]))
    )
  }
})


test_that("confint formatting excludes stored profile diagnostics", {
  fit <- make_confint_test_fit()
  fit$RC$CI <- data.frame(
    lower = c(NA_real_, 1.5),
    upper = c(1.2, NA_real_),
    pval.lower = c(NA_real_, 0.05),
    pval.upper = c(0.05, NA_real_),
    iter.lower = c(NA_integer_, 4L),
    iter.upper = c(3L, NA_integer_),
    row.names = c("a", "b")
  )

  result <- ameras:::format_confint_bounds(fit, level = 0.95)

  expect_identical(colnames(result), c("2.5 %", "97.5 %"))
  expect_identical(rownames(result), c("a", "b"))
  expect_identical(
    unname(result),
    unname(as.matrix(fit$RC$CI[, c("lower", "upper")]))
  )
  expect_true(is.na(result["a", "2.5 %"]))
  expect_true(is.na(result["b", "97.5 %"]))
  expect_identical(
    colnames(ameras:::format_confint_bounds(fit, level = 0.90)),
    c("5 %", "95 %")
  )
})


test_that("summaries detect attached intervals from method results", {
  fit <- make_confint_test_fit(include_fma = TRUE)
  fit$call <- quote(ameras(Y ~ dose(D), data = dat))
  fit$num.rows <- 2L
  without_ci <- summary(fit)

  expect_false("CI.computed" %in% names(without_ci))
  expect_named(coef(without_ci), c("Method", "Term", "Estimate", "SE"))
  expect_output(print(without_ci), "confidence intervals not yet computed")

  # A stored CI counts even when every bound is NA. A missing parameter or
  # method receives NA bounds so the combined table keeps the same columns.
  for (bounds in list(c(NA_real_, NA_real_), c(1.5, 2.5))) {
    fit$RC$CI <- data.frame(
      lower = bounds[1],
      upper = bounds[2],
      row.names = "b"
    )
    with_ci <- summary(fit)

    expect_equal(coef(with_ci)$CI.lower, c(NA, bounds[1], NA, NA))
    expect_equal(coef(with_ci)$CI.upper, c(NA, bounds[2], NA, NA))
    expect_false(any(grepl(
      "confidence intervals not yet computed",
      capture.output(print(with_ci)),
      fixed = TRUE
    )))
  }

  fit$RC$CI <- NULL
  expect_identical(summary(fit), without_ci)
})


test_that("CI methods leave interval-table printing to the caller", {
  fit <- make_confint_test_fit()

  for (argument in c("print", "digits")) {
    args <- c(
      list(object = fit, type = "wald.orig"),
      setNames(list(TRUE), argument)
    )
    expect_error(
      do.call(confint, args),
      "not supported by `confint\\(\\)`.*Assign its result.*print\\(confint",
      info = argument
    )
    expect_warning(
      do.call(add_confint, args),
      "extra argument",
      info = argument
    )
  }

  expect_error(
    confint(fit, type = "wald.orig", force = TRUE),
    "`force` is no longer supported"
  )
  expect_error(
    add_confint(fit, type = "wald.orig", force = TRUE),
    "`force` is no longer supported"
  )

  expect_warning(
    confint(fit, type = "wald.orig", typo = TRUE),
    "extra argument"
  )
})


test_that("add_confint invisibly amends a copy with the requested intervals", {
  fit <- make_confint_test_fit()
  output <- capture.output(visibility <- withVisible(add_confint(
    fit,
    type = "wald.orig"
  )))

  expect_identical(output, character())
  expect_false(visibility$visible)
  expect_s3_class(visibility$value, "amerasfit")
  expect_false("CI.computed" %in% names(visibility$value))
  expect_named(visibility$value$RC$CI, c("lower", "upper"))
  expect_null(fit$RC$CI)

  replaced <- add_confint(
    visibility$value,
    type = "wald.orig",
    level = 0.90
  )
  expected <- confint(fit, type = "wald.orig", level = 0.90)

  expect_equal(
    unname(as.matrix(replaced$RC$CI[, c("lower", "upper")])),
    unname(expected)
  )
  expect_false(isTRUE(all.equal(replaced$RC$CI, visibility$value$RC$CI)))
})
