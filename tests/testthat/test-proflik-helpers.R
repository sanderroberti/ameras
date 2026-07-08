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


test_that("profile likelihood boundary warnings use informative precision", {
  # This transform maps the lower profile-search boundary to a small negative
  # value. Using significant digits should preserve that value instead of
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
    "lower profile search boundary \\(corresponding parameter value: -0.04495\\)"
  )

  expect_equal(ci$lower, -Inf)
  expect_true(is.finite(ci$upper))
})
