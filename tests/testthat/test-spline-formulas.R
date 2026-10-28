make_spline_formula_data <- function(n = 40) {
  set.seed(462)
  age <- seq(20, 80, length.out = n)
  dose <- seq(0.1, 1.4, length.out = n)
  data.frame(
    Y = 1 + 0.5 * dose + sin(age / 12) + rnorm(n, sd = 0.15),
    D = dose,
    age = age,
    group = factor(rep(c("a", "b"), length.out = n))
  )
}

test_that("namespace-qualified spline basis functions are accepted in X terms", {
  dat <- make_spline_formula_data()

  # Spline support should be general RHS model-matrix support, not a special
  # case for unqualified ns(). Namespace-qualified calls are the safest user
  # spelling because they do not depend on splines being attached.
  parsed_ns <- ameras:::parse_ameras_formula(
    Y ~ dose(D) + splines::ns(age, df = 3),
    data = dat,
    family = "gaussian"
  )
  parsed_bs <- ameras:::parse_ameras_formula(
    Y ~ dose(D) + splines::bs(age, df = 4),
    data = dat,
    family = "gaussian"
  )

  expect_identical(parsed_ns$Y, "Y")
  expect_identical(parsed_ns$dosevars, "D")
  expect_match(
    paste(deparse(parsed_ns$X_formula), collapse = ""),
    "splines::ns",
    fixed = TRUE
  )
  expect_match(
    paste(deparse(parsed_bs$X_formula), collapse = ""),
    "splines::bs",
    fixed = TRUE
  )
  expect_equal(
    ncol(ameras:::build_X_design(parsed_ns$X_formula, dat)$matrix),
    3L
  )
  expect_equal(
    ncol(ameras:::build_X_design(parsed_bs$X_formula, dat)$matrix),
    4L
  )
})

test_that("namespace-qualified spline basis functions can be fitted", {
  dat <- make_spline_formula_data()

  fit_ns <- suppressWarnings(suppressMessages(
    ameras(
      Y ~ dose(D) + splines::ns(age, df = 3),
      data = dat,
      family = "gaussian",
      methods = "RC"
    )
  ))
  fit_bs <- suppressWarnings(suppressMessages(
    ameras(
      Y ~ dose(D) + splines::bs(age, df = 4),
      data = dat,
      family = "gaussian",
      methods = "RC"
    )
  ))

  expect_s3_class(fit_ns, "amerasfit")
  expect_s3_class(fit_bs, "amerasfit")
  expect_length(fit_ns$model$X, 3L)
  expect_length(fit_bs$model$X, 4L)
  expect_true(all(is.finite(fit_ns$RC$coefficients)))
  expect_true(all(is.finite(fit_bs$RC$coefficients)))
})

test_that("keep.data FALSE can rebuild spline design columns from supplied data", {
  dat <- make_spline_formula_data()

  # This test uses an unqualified ns() with splines attached because users often
  # write formulas that way interactively. Reconstructing the design matrix
  # should not depend on the stored formula having access to the caller's search
  # path at a later point.
  library(splines)
  fit <- suppressWarnings(suppressMessages(
    ameras(
      Y ~ dose(D) + ns(age, df = 3),
      data = dat,
      family = "gaussian",
      methods = "RC",
      keep.data = FALSE
    )
  ))

  expect_null(fit$model$data)
  expect_true(is.list(fit$model$X_design_info))
  res <- residuals(fit, data = dat, type = "response")
  expect_length(res, nrow(dat))
  expect_true(all(is.finite(res)))
})

test_that("rebuilt spline design columns match the fitted design columns", {
  dat <- make_spline_formula_data()

  # When data have to be supplied after fitting, the spline basis should be
  # reconstructed with the same model-matrix design information used at fit
  # time. This is important for all splines-package bases, including ns() and
  # bs().
  library(splines)
  fit_with_data <- suppressWarnings(suppressMessages(
    ameras(
      Y ~ dose(D) + bs(age, df = 4) + group,
      data = dat,
      family = "gaussian",
      methods = "RC",
      keep.data = TRUE
    )
  ))
  fit_without_data <- suppressWarnings(suppressMessages(
    ameras(
      Y ~ dose(D) + bs(age, df = 4) + group,
      data = dat,
      family = "gaussian",
      methods = "RC",
      keep.data = FALSE
    )
  ))

  resolved <- ameras:::resolve_data(fit_without_data, data = dat)

  expect_identical(fit_without_data$model$X, fit_with_data$model$X)
  expect_equal(
    resolved[, fit_without_data$model$X, drop = FALSE],
    fit_with_data$model$data[, fit_with_data$model$X, drop = FALSE],
    tolerance = 1e-12
  )
})

test_that("stored spline design information does not retain caller workspace", {
  dat <- make_spline_formula_data()

  fit <- local({
    library(splines)
    large_object_that_should_not_be_stored <- rnorm(100000)
    suppressWarnings(suppressMessages(
      ameras(
        Y ~ dose(D) + ns(age, df = 3),
        data = dat,
        family = "gaussian",
        methods = "RC",
        keep.data = FALSE
      )
    ))
  })

  design_env <- environment(fit$model$X_design_info$terms)

  # The terms object needs ns() to rebuild the model matrix, but it should not
  # keep the whole caller environment alive.
  expect_true(exists("ns", envir = design_env, inherits = FALSE))
  expect_false(
    exists(
      "large_object_that_should_not_be_stored",
      envir = design_env,
      inherits = FALSE
    )
  )
  expect_identical(parent.env(design_env), baseenv())
  res <- residuals(fit, data = dat, type = "response")
  expect_length(res, nrow(dat))
  expect_true(all(is.finite(res)))
})
