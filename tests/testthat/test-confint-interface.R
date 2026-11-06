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
    CI.computed = FALSE,
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
  amended <- add_confint(fit, type = "wald.orig", print = FALSE)

  expect_true(visibility$visible)
  expect_true(is.matrix(visibility$value))
  expect_identical(dimnames(visibility$value), list(c("a", "b"), c("2.5 %", "97.5 %")))
  expect_equal(
    unname(visibility$value),
    unname(as.matrix(amended$RC$CI[, c("lower", "upper")]))
  )
  expect_false(fit$CI.computed)
  expect_null(fit$RC$CI)
})


test_that("confint returns a named list for multiple fitted methods", {
  fit <- make_confint_test_fit(include_fma = TRUE)
  result <- confint(fit, type = c("wald.orig", "percentile"))
  amended <- add_confint(
    fit,
    type = c("wald.orig", "percentile"),
    print = FALSE
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
  fit$CI.computed <- TRUE

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


test_that("printing arguments direct confint users to add_confint", {
  fit <- make_confint_test_fit()

  for (argument in c("print", "digits")) {
    args <- c(
      list(object = fit, type = "wald.orig"),
      setNames(list(TRUE), argument)
    )
    expect_error(
      do.call(confint, args),
      "only supported by `add_confint\\(\\)`",
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
  visibility <- withVisible(add_confint(
    fit,
    type = "wald.orig",
    print = FALSE
  ))

  expect_false(visibility$visible)
  expect_s3_class(visibility$value, "amerasfit")
  expect_true(visibility$value$CI.computed)
  expect_named(visibility$value$RC$CI, c("lower", "upper"))
  expect_false(fit$CI.computed)

  replaced <- add_confint(
    visibility$value,
    type = "wald.orig",
    level = 0.90,
    print = FALSE
  )
  expected <- confint(fit, type = "wald.orig", level = 0.90)

  expect_equal(
    unname(as.matrix(replaced$RC$CI[, c("lower", "upper")])),
    unname(expected)
  )
  expect_false(isTRUE(all.equal(replaced$RC$CI, visibility$value$RC$CI)))
})
