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
