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
  expect_lt(fit$gradient.rms, 1e-5)
  expect_lt(fit$gradient.rms.scaled, 1e-5)
  expect_true(ameras:::fit_passes_hessian_check(fit))
})


test_that("fit_objective_with_hessian can skip expensive Hessian work", {
  fit <- ameras:::fit_objective_with_hessian(
    start = c(x = 1, y = -1),
    fn = function(par) sum((par - c(0.5, 0.25))^2),
    compute.hessian = FALSE,
    gradient.check = FALSE
  )

  expect_equal(unname(fit$par), c(0.5, 0.25), tolerance = 1e-6)
  expect_equal(fit$value, 0, tolerance = 1e-8)
  expect_null(fit$hessian)
  expect_null(fit$gradient)
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
  expect_false(ameras:::fit_passes_hessian_check(list(
    convergence = 0,
    hessian = matrix(Inf, nrow = 1)
  )))
  expect_false(ameras:::fit_passes_hessian_check(list(
    convergence = 0,
    hessian = matrix(NA_real_, nrow = 1)
  )))
  expect_false(ameras:::hessian_supports_vcov(matrix(Inf, nrow = 1)))
  expect_false(ameras:::hessian_supports_vcov(matrix(NA_real_, nrow = 1)))
  expect_false(ameras:::hessian_supports_vcov(matrix(numeric(0), nrow = 0)))
  expect_false(ameras:::hessian_supports_vcov(matrix(1:2, nrow = 1)))
  expect_false(ameras:::hessian_supports_vcov(matrix(c(1, 2, 2, 1), nrow = 2)))
})

test_that("optimizer gradient diagnostics warn only for suspicious convergence", {
  suspicious_fit <- list(
    convergence = 0,
    gradient.rms = 10,
    gradient.rms.scaled = 1e-2,
    gradient.max.abs = 12,
    newton.decrement = 0.2,
    newton.improvement = 0.02,
    newton.improvement.relative = 2e-5
  )
  acceptable_fit <- list(
    convergence = 0,
    gradient.rms = 10,
    gradient.rms.scaled = 1e-2,
    gradient.max.abs = 12,
    newton.decrement = 0.05531,
    newton.improvement = 0.5 * 0.05531^2,
    newton.improvement.relative = 1e-6
  )
  scaled_fit <- list(
    convergence = 0,
    gradient.rms = 10,
    gradient.rms.scaled = 1e-2,
    gradient.max.abs = 12,
    newton.decrement = 0.2,
    newton.improvement = 0.02,
    newton.improvement.relative = 1e-8
  )
  fallback_fit <- list(
    convergence = 0,
    gradient.rms = 10,
    gradient.rms.scaled = 1e-2,
    gradient.max.abs = 12
  )
  not_converged_fit <- suspicious_fit
  not_converged_fit$convergence <- 1

  # The diagnostic is intentionally warning-only. When available, it uses a
  # curvature-scaled approximate objective improvement so large raw gradients do
  # not warn if the remaining improvement is negligible in absolute terms or
  # relative to the objective scale.
  expect_warning(
    ameras:::warn_if_large_optimizer_gradient(suspicious_fit),
    "not be fully stationary"
  )
  expect_silent(ameras:::warn_if_large_optimizer_gradient(acceptable_fit))
  expect_silent(ameras:::warn_if_large_optimizer_gradient(scaled_fit))
  expect_warning(
    ameras:::warn_if_large_optimizer_gradient(fallback_fit),
    "not be fully stationary"
  )
  expect_silent(ameras:::warn_if_large_optimizer_gradient(not_converged_fit))
})

test_that("convergence extracts and recomputes optimizer gradients", {
  set.seed(20260708)
  D <- runif(40, 0.1, 1.5)
  X <- seq(-1, 1, length.out = 40)
  dat <- data.frame(
    Y = 1 + 0.7 * D + 0.4 * X + rnorm(40, sd = 0.25),
    D = D,
    X = X
  )
  fit <- suppressWarnings(suppressMessages(
    ameras(
      Y ~ dose(D) + X,
      data = dat,
      family = "gaussian",
      methods = "RC"
    )
  ))

  stored <- convergence(fit)
  expect_s3_class(stored, "data.frame")
  expect_identical(rownames(stored), "RC")
  expect_named(
    stored,
    c(
      "method",
      "optim.convergence",
      "gradient.rms",
      "gradient.rms.scaled",
      "newton.improvement",
      "newton.improvement.relative",
      "convergence.warning"
    )
  )
  expect_equal(stored$optim.convergence, 0)
  expect_true(is.finite(stored$gradient.rms))
  expect_true(is.finite(stored$newton.improvement))

  # Simulate an object fitted before gradient diagnostics were stored. The
  # method should reconstruct the likelihood and compute the gradient on demand.
  old_fit <- fit
  old_fit$RC$optim$gradient <- NULL
  old_fit$RC$optim$gradient.max.abs <- NULL
  old_fit$RC$optim$gradient.rms <- NULL
  old_fit$RC$optim$gradient.rms.scaled <- NULL

  recomputed <- convergence(old_fit)
  expect_equal(recomputed$gradient.rms, stored$gradient.rms, tolerance = 1e-6)
  expect_equal(
    recomputed$gradient.rms.scaled,
    stored$gradient.rms.scaled,
    tolerance = 1e-6
  )
  expect_equal(
    recomputed$newton.improvement,
    stored$newton.improvement,
    tolerance = 1e-6
  )
  expect_equal(
    recomputed$newton.improvement.relative,
    stored$newton.improvement.relative,
    tolerance = 1e-6
  )
})

test_that("convergence reports unsupported model-averaging methods", {
  fit <- list(FMA = list(samples = data.frame(dose = 1:3)))
  class(fit) <- "amerasfit"

  expect_error(
    convergence(fit, methods = "FMA"),
    "available for RC, ERC, and MCML"
  )
})

test_that("convergence can recompute gradients when data were not stored", {
  set.seed(20260708)
  D <- runif(40, 0.1, 1.5)
  X <- seq(-1, 1, length.out = 40)
  dat <- data.frame(
    Y = 1 + 0.7 * D + 0.4 * X + rnorm(40, sd = 0.25),
    D = D,
    X = X
  )
  fit <- suppressWarnings(suppressMessages(
    ameras(
      Y ~ dose(D) + X,
      data = dat,
      family = "gaussian",
      methods = "RC",
      keep.data = FALSE
    )
  ))

  fit$RC$optim$gradient <- NULL
  fit$RC$optim$gradient.max.abs <- NULL
  fit$RC$optim$gradient.rms <- NULL
  fit$RC$optim$gradient.rms.scaled <- NULL

  expect_error(
    convergence(fit),
    "Data not stored"
  )
  expect_no_error(diag <- convergence(fit, data = dat))
  expect_true(is.finite(diag$gradient.rms))
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
  expect_false("runtime" %in% names(out))
  expect_true(all(c("fit", "ci", "total") %in% names(out$timing)))
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
        timing = ameras:::new_method_timing()
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
