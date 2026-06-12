test_that("fit_objective_with_hessian handles one-parameter objectives", {
  fit <- ameras:::fit_objective_with_hessian(
    start = 0,
    fn = function(x) (x - 2)^2
  )

  expect_equal(fit$par, 2, tolerance = 1e-5)
  expect_equal(fit$value, 0, tolerance = 1e-8)
  expect_equal(fit$convergence, 0)
  expect_true(ameras:::fit_passes_hessian_check(fit))
})

test_that("fit_objective_with_hessian handles multi-parameter objectives", {
  fit <- ameras:::fit_objective_with_hessian(
    start = c(0, 0),
    fn = function(x) sum((x - c(1, -2))^2)
  )

  expect_equal(fit$par, c(1, -2), tolerance = 1e-5)
  expect_equal(fit$value, 0, tolerance = 1e-8)
  expect_equal(fit$convergence, 0)
  expect_true(ameras:::fit_passes_hessian_check(fit))
})

test_that("fit_objective_with_hessian forwards objective arguments to hessian", {
  fit <- ameras:::fit_objective_with_hessian(
    start = c(0, 0),
    fn = function(x, center, scale_factor) scale_factor * sum((x - center)^2),
    center = c(2, -3),
    scale_factor = 4
  )

  expect_equal(fit$par, c(2, -3), tolerance = 1e-5)
  expect_equal(fit$value, 0, tolerance = 1e-8)
  expect_equal(fit$hessian, diag(8, 2), tolerance = 1e-5)
})

test_that("fit_passes_hessian_check rejects unusable fits", {
  expect_false(ameras:::fit_passes_hessian_check(list(
    convergence = 1,
    hessian = diag(1)
  )))
  expect_false(ameras:::fit_passes_hessian_check(list(
    convergence = 0,
    hessian = matrix(0, nrow = 1)
  )))
  expect_false(ameras:::fit_passes_hessian_check(list(
    convergence = 0,
    hessian = matrix(-1, nrow = 1)
  )))
})

test_that("assemble_frequentist_fit_result handles untransformed fits", {
  fit <- list(
    par = c(1, 2),
    value = 5,
    convergence = 0,
    counts = c(`function` = 3, gradient = 1),
    hessian = diag(c(4, 9))
  )

  out <- ameras:::assemble_frequentist_fit_result(
    fit = fit,
    parnames = c("alpha", "beta"),
    t0 = proc.time()
  )

  expect_named(out$coefficients, c("alpha", "beta"))
  expect_equal(out$coefficients, c(alpha = 1, beta = 2))
  expect_equal(unname(out$vcov), diag(c(1 / 4, 1 / 9)))
  expect_identical(rownames(out$vcov), c("alpha", "beta"))
  expect_identical(colnames(out$vcov), c("alpha", "beta"))
  expect_equal(out$sd, sqrt(c(alpha = 1 / 4, beta = 1 / 9)))
  expect_identical(out$optim$counts, fit$counts)
  expect_equal(out$loglik, -5)
  expect_match(out$runtime, "seconds$")
})

test_that("assemble_frequentist_fit_result handles transformed fits", {
  fit <- list(
    par = c(1, 2),
    value = 7,
    convergence = 0,
    counts = NULL,
    hessian = diag(c(4, 9))
  )
  transform_shift <- function(params, shift, ...) params + shift
  jacobian_scale <- function(params, jac_scale, ...) diag(jac_scale, 2)

  out <- ameras:::assemble_frequentist_fit_result(
    fit = fit,
    parnames = c("theta1", "theta2"),
    t0 = proc.time(),
    transform = transform_shift,
    transform.jacobian = jacobian_scale,
    shift = c(10, 20),
    jac_scale = c(2, 3)
  )

  # The variance is jacobian %*% inverse(hessian) %*% t(jacobian).
  expect_equal(out$coefficients, c(theta1 = 11, theta2 = 22))
  expect_equal(unname(out$vcov), diag(c(1, 1)))
  expect_null(out$optim$counts)
})

test_that("assemble_frequentist_fit_result handles boundcheck transforms", {
  fit <- list(
    par = c(1, 2),
    value = 1,
    convergence = 0,
    counts = NULL,
    hessian = diag(2)
  )
  transform_boundcheck <- function(params, boundcheck = FALSE, ...) {
    if (!isTRUE(boundcheck)) {
      stop("boundcheck was not forwarded")
    }
    params
  }

  out <- ameras:::assemble_frequentist_fit_result(
    fit = fit,
    parnames = c("a", "b"),
    t0 = proc.time(),
    transform = transform_boundcheck,
    transform.jacobian = function(params, ...) diag(2)
  )

  expect_equal(out$coefficients, c(a = 1, b = 2))
})

test_that("assemble_frequentist_fit_result warns for unusable Hessians", {
  fit <- list(
    par = c(1, 2),
    value = 3,
    convergence = 0,
    counts = NULL,
    hessian = matrix(0, nrow = 2, ncol = 2)
  )

  expect_warning(
    out <- ameras:::assemble_frequentist_fit_result(
      fit = fit,
      parnames = c("a", "b"),
      t0 = proc.time()
    ),
    "Hessian was not invertible",
    fixed = TRUE
  )

  expect_equal(out$coefficients, c(a = 1, b = 2))
  expect_true(all(is.na(out$vcov)))
  expect_identical(dim(out$vcov), c(2L, 2L))
  expect_identical(rownames(out$vcov), c("a", "b"))
  expect_identical(colnames(out$vcov), c("a", "b"))
})

test_that("assemble_frequentist_fit_result validates transform functions", {
  fit <- list(
    par = 1,
    value = 3,
    convergence = 0,
    counts = NULL,
    hessian = matrix(1, nrow = 1, ncol = 1)
  )

  expect_error(
    ameras:::assemble_frequentist_fit_result(
      fit = fit,
      parnames = "dose",
      t0 = proc.time(),
      transform = "not a function",
      transform.jacobian = function(params, ...) matrix(1, 1, 1)
    ),
    "transform and transform.jacobian should be functions",
    fixed = TRUE
  )
})

test_that("BMA automatic realization screening is disabled explicitly", {
  expect_error(
    ameras:::ameras.bma(
      family = "gaussian",
      dosevars = "V1",
      data = data.frame(Y = 1:3, V1 = 1:3),
      deg = 1,
      Y = "Y",
      included.realizations = NULL
    ),
    "Automatic BMA realization screening is not currently enabled",
    fixed = TRUE
  )
})

test_that("ameras passes default BMA realization indices explicitly", {
  data("data", package = "ameras")
  captured <- NULL

  testthat::with_mocked_bindings(
    ameras.bma = function(...) {
      captured <<- list(...)
      list(
        coefficients = numeric(),
        sd = numeric(),
        vcov = matrix(numeric()),
        Rhat = NULL,
        samples = NULL,
        included.realizations = captured$included.realizations,
        runtime = "0 seconds"
      )
    },
    {
      fit <- suppressMessages(
        ameras(
          Y.gaussian ~ dose(V1:V2),
          data = data[1:20, ],
          family = "gaussian",
          methods = "BMA"
        )
      )
    }
  )

  expect_identical(captured$included.realizations, 1:2)
  expect_identical(fit$BMA$included.realizations, 1:2)
})
