test_that("profile likelihood CIs use optimize for two-parameter fits", {
  alpha <- 0.05

  # With this objective, the unconstrained optimum is at c(0, 0).
  # When one parameter is fixed at x, profiling over the other parameter leaves
  # an objective difference of x^2. compute_proflik_ci_one() therefore solves
  # 2 * x^2 = qchisq(1 - alpha, 1), giving an analytic target interval.
  quadratic_loglik <- function(params) {
    sum(params^2)
  }
  expected_limit <- sqrt(qchisq(1 - alpha, df = 1) / 2)

  # Cover both ways proflik() rebuilds the full parameter vector in the
  # length(inpar) == 2 branch: fixed parameter first and fixed parameter second.
  for (index in 1:2) {
    ci <- ameras:::compute_proflik_ci_one(
      index = index,
      inpar = c(theta1 = 0, theta2 = 0),
      optval = 0,
      loglik_fn = quadratic_loglik,
      lowlim = -3,
      uplim = 3,
      alpha = alpha,
      maxit.profCI = 100,
      tol.profCI = 1e-8,
      parname = paste0("theta", index)
    )

    expect_equal(ci$lower, -expected_limit, tolerance = 1e-4)
    expect_equal(ci$upper, expected_limit, tolerance = 1e-4)
    expect_equal(ci$pval.lower, alpha, tolerance = 0.005)
    expect_equal(ci$pval.upper, alpha, tolerance = 0.005)
    expect_true(ci$iter.lower > 0)
    expect_true(ci$iter.upper > 0)
  }
})


test_that("profile likelihood bracketing expands beyond initial limits", {
  alpha <- 0.05
  quadratic_loglik <- function(params) {
    params[1]^2
  }
  expected_limit <- sqrt(qchisq(1 - alpha, df = 1) / 2)

  # The initial limits are intentionally much narrower than the true profile
  # interval. The explicit bracketing step should expand efficiently and still
  # find the analytic roots.
  expect_no_warning(
    ci <- ameras:::compute_proflik_ci_one(
      index = 1,
      inpar = c(theta = 0),
      optval = 0,
      loglik_fn = quadratic_loglik,
      lowlim = -0.1,
      uplim = 0.1,
      alpha = alpha,
      maxit.profCI = 100,
      tol.profCI = 1e-8,
      parname = "theta"
    )
  )

  expect_equal(ci$lower, -expected_limit, tolerance = 1e-4)
  expect_equal(ci$upper, expected_limit, tolerance = 1e-4)
  expect_equal(ci$pval.lower, alpha, tolerance = 0.005)
  expect_equal(ci$pval.upper, alpha, tolerance = 0.005)
})


test_that("profile likelihood unresolved bounds return NA with diagnostics", {
  # This transform maps very small optimizer-scale values to a small negative
  # reported-scale value. The warning should preserve that value instead of
  # rounding it to 0.
  transform_lower_limit <- function(params, ...) {
    params[1] <- exp(params[1]) - 0.045
    params
  }

  one_sided_loglik <- function(params) {
    if (params[1] < 0) {
      1e-8 * params[1]^2
    } else {
      params[1]^2
    }
  }

  expect_warning(
    ci <- ameras:::compute_proflik_ci_one(
      index = 1,
      inpar = 0,
      optval = 0,
      loglik_fn = one_sided_loglik,
      lowlim = -3,
      uplim = 3,
      alpha = 0.05,
      maxit.profCI = 100,
      tol.profCI = 1e-8,
      parname = "dose",
      transform = transform_lower_limit
    ),
    paste0(
      "lower bound for dose could not be bracketed; returning NA.*",
      "Last attempted value: -0.045 ",
      "\\(internal optimizer-scale value: .*\\).*",
      "p-value at last attempted value: .*target p-value: 0.05"
    )
  )

  expect_true(is.na(ci$lower))
  expect_true(is.finite(ci$pval.lower))
  expect_true(is.finite(ci$upper))
})


test_that("profile likelihood inaccurate roots return NA with diagnostics", {
  quadratic_loglik <- function(params) {
    params[1]^2
  }

  expect_warning(
    expect_warning(
      ci <- ameras:::compute_proflik_ci_one(
        index = 1,
        inpar = c(theta = 0),
        optval = 0,
        loglik_fn = quadratic_loglik,
        lowlim = -3,
        uplim = 3,
        alpha = 0.05,
        maxit.profCI = 1,
        tol.profCI = 1e-12,
        parname = "theta"
      ),
      paste0(
        "lower bound for theta was not solved accurately; returning NA.*",
        "p-value at last attempted value: .*target p-value: 0.05"
      )
    ),
    paste0(
      "upper bound for theta was not solved accurately; returning NA.*",
      "p-value at last attempted value: .*target p-value: 0.05"
    )
  )

  expect_true(is.na(ci$lower))
  expect_true(is.na(ci$upper))
  expect_true(is.finite(ci$pval.lower))
  expect_true(is.finite(ci$pval.upper))
  expect_equal(ci$iter.lower, 1)
  expect_equal(ci$iter.upper, 1)
})
