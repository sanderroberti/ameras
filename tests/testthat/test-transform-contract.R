test_that("coupled transformations are rejected before fitting for every method", {
  coupled <- function(params, ...) {
    params[2] <- params[1] + params[2]
    params[3] <- exp(params[3])
    params
  }
  jacobian <- function(params, ...) {
    matrix(c(1, 1, 0, 0, 1, 0, 0, 0, exp(params[3])), nrow = 3)
  }
  for (method in c("RC", "ERC", "MCML", "FMA", "BMA")) {
    expect_error(
      suppressMessages(ameras(
        Y ~ dose(D1:D2), data = transform_contract_data(), methods = method,
        transform = coupled, transform.jacobian = jacobian
      )),
      "Coupled transformations"
    )
  }
  expect_error(
    suppressMessages(ameras(
      Y ~ dose(D1:D2), data = transform_contract_data(), methods = "BMA",
      transform = coupled
    )),
    "Coupled transformations"
  )
})

test_that("validation detects coupling away from a stationary starting point", {
  coupled <- function(params, ...) c(params[1], params[2] + params[1]^2)
  jacobian <- function(params, ...) matrix(c(1, 2 * params[1], 0, 1), 2)
  expect_equal(jacobian(c(0, 0)), diag(2))
  expect_error(
    ameras:::check_transform(coupled, jacobian, c(0, 0)),
    "Coupled transformations"
  )
  expect_error(
    ameras:::check_transform(function(params) c(params, 1), NULL, c(0, 0)),
    "unchanged length"
  )
  expect_error(
    ameras:::check_transform(function(params) params^2, NULL, c(0, 0)),
    "strictly monotone"
  )
  expect_error(
    ameras:::check_transform(
      function(params) params, function(params) diag(1), c(0, 0)
    ),
    "square matrix"
  )
  expect_error(
    ameras:::check_transform(
      function(params) params + 0.2 * sin(10 * params),
      function(params) diag(1 + 2 * cos(10 * params)), c(0, 0)
    ),
    "strictly monotone"
  )
})

test_that("identity, exponential, and sigmoid componentwise transforms are accepted", {
  transforms <- list(
    function(params, ...) params,
    function(params, ...) exp(params),
    function(params, ...) plogis(params)
  )
  jacobians <- list(
    function(params, ...) diag(length(params)),
    function(params, ...) diag(exp(params)),
    function(params, ...) diag(plogis(params) * (1 - plogis(params)))
  )
  for (i in seq_along(transforms)) {
    direction <- ameras:::check_transform(
      transforms[[i]], jacobians[[i]], c(0.2, -0.1)
    )
    expect_equal(direction, c(1, 1))
  }
})

test_that("transformation directions are local to transformed interval calculations", {
  calls <- 0L
  transform <- function(params, ...) c(params[1], -params[2], exp(params[3]))
  jacobian <- function(params, ...) {
    calls <<- calls + 1L
    diag(c(1, -1, exp(params[3])))
  }
  fit <- suppressMessages(ameras(
    Y ~ dose(D1:D2), data = transform_contract_data(),
    transform = transform, transform.jacobian = jacobian,
    inpar = c(1, -0.5, log(0.2))
  ))
  expect_identical(fit$transform, transform)
  expect_identical(fit$transform.jacobian, jacobian)
  expect_false("transform.direction" %in% names(fit))
  # One initial evaluation and four probes per parameter, plus one Jacobian
  # evaluation to propagate the fitted covariance matrix.
  expect_identical(calls, 14L)
  confint(fit, type = "wald.orig")
  expect_identical(calls, 14L)
  confint(fit, type = "wald.transformed")
  expect_identical(calls, 27L)
  suppressMessages(confint(fit, type = "proflik"))
  expect_identical(calls, 40L)

  amended <- add_confint(fit, type = "wald.transformed")
  expect_identical(calls, 53L)
  expect_false("transform.direction" %in% names(amended))
  expect_identical(amended$transform, transform)
  confint(amended, type = "wald.transformed")
  expect_identical(calls, 66L)
})

test_that("sample intervals do not evaluate parameter transformations", {
  for (method in c("FMA", "BMA")) {
    fit <- ameras:::new_amerasfit(c(
      list(
        transform = function(...) stop("transformation should not be evaluated"),
        transform.jacobian = function(...) stop("Jacobian should not be evaluated")
      ),
      setNames(list(list(
        coefficients = c(dose = 0), samples = cbind(dose = seq(-1, 1, length.out = 100)),
        timing = ameras:::new_method_timing()
      )), method)
    ))
    for (type in c("percentile", "hpd")) {
      expect_true(all(is.finite(confint(fit, type = type))))
      expect_false("transform.direction" %in% names(add_confint(fit, type = type)))
    }
  }
})

test_that("decreasing transformations preserve fitted results and ordered intervals", {
  dat <- transform_contract_data()
  decreasing <- function(params, ...) {
    params[2] <- -params[2]
    params[3] <- exp(params[3])
    params
  }
  jacobian <- function(params, ...) diag(c(1, -1, exp(params[3])))
  plain <- suppressMessages(ameras(Y ~ dose(D1:D2), data = dat))
  fit <- suppressMessages(ameras(
    Y ~ dose(D1:D2), data = dat, transform = decreasing,
    transform.jacobian = jacobian, inpar = c(1, -0.5, log(0.2))
  ))
  expect_equal(coef(fit), coef(plain), tolerance = 1e-4)
  for (type in c("wald.orig", "wald.transformed", "proflik")) {
    ci <- suppressMessages(confint(fit, type = type))
    expected <- suppressMessages(confint(plain, type = type))
    expect_true(all(is.finite(ci["dose", ])))
    expect_lt(ci["dose", 1], ci["dose", 2])
    expect_equal(ci["dose", ], expected["dose", ], tolerance = 1e-4)
  }
})

test_that("decreasing profile bounds keep one-sided failures and their diagnostics", {
  transform <- function(params, ...) c(-params[1], params[2])
  method_fit <- list(
    coefficients = c(dose = 0, other = 0),
    optim = list(par = c(0, 0), hessian = diag(2)), loglik = 0
  )
  object <- list(
    transform = transform, transform.jacobian = function(params, ...) diag(c(-1, 1)),
    other.args = list(), RC = method_fit
  )
  testthat::local_mocked_bindings(
    make_loglik_fn = function(...) function(params) sum(params^2),
    compute_proflik_ci_one = function(..., decreasing) {
      expect_true(decreasing)
      list(lower = NA_real_, upper = 1, pval.lower = 0.3, pval.upper = 0.05,
           iter.lower = 2L, iter.upper = 7L)
    },
    .package = "ameras"
  )
  ci <- suppressMessages(ameras:::compute_proflik_CI(method_fit, object, "RC", data.frame()))
  expect_equal(ci$lower, -1)
  expect_true(is.na(ci$upper))
  expect_equal(ci$pval.lower, 0.05)
  expect_equal(ci$pval.upper, 0.3)
  expect_equal(ci$iter.lower, 7L)
  expect_equal(ci$iter.upper, 2L)
})

test_that("decreasing transformations use reported directions in profile warnings", {
  expect_warning(
    ci <- ameras:::compute_proflik_ci_one(
      index = 1, inpar = 0, optval = 0,
      loglik_fn = function(params) if (params[1] < 0) 0 else params[1]^2,
      lowlim = -1, uplim = 1, parname = "dose",
      transform = function(params, ...) -params, decreasing = TRUE
    ),
    "upper bound for dose could not be bracketed"
  )
  expect_true(is.na(ci$lower))
  expect_true(is.finite(ci$upper))
})

test_that("older coupled fits are rejected before recomputed inference", {
  fit <- suppressMessages(ameras(Y ~ dose(D1:D2), data = transform_contract_data()))
  fit$transform <- function(params, ...) c(params[1], params[1] + params[2], exp(params[3]))
  fit$transform.jacobian <- function(params, ...) {
    matrix(c(1, 1, 0, 0, 1, 0, 0, 0, exp(params[3])), 3)
  }
  for (type in c("wald.transformed", "proflik")) {
    expect_error(confint(fit, type = type), "Coupled transformations")
  }
  expect_error(dose_lrt(fit), "Coupled transformations")
  expect_error(convergence(fit, recompute = TRUE), "Coupled transformations")
})

test_that("BMA-only transformations do not require a Jacobian", {
  testthat::local_mocked_bindings(
    ameras.bma = function(...) list(coefficients = c("(Intercept)" = 0, dose = 0.5, sigma = 1)),
    .package = "ameras"
  )
  expect_s3_class(
    suppressMessages(ameras(
      Y ~ dose(D1:D2), data = transform_contract_data(), methods = "BMA",
      transform = function(params, ...) params
    )),
    "amerasfit"
  )
  expect_error(
    suppressMessages(ameras(
      Y ~ dose(D1:D2), data = transform_contract_data(), methods = c("RC", "BMA"),
      transform = function(params, ...) params
    )),
    "transform.jacobian is required"
  )
})
