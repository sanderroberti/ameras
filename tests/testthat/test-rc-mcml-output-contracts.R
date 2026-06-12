# These tests lock down behavior that will be easy to disturb when RC, ERC,
# and MCML start sharing likelihood, optimizer, and result-assembly helpers.
# They intentionally check names and argument forwarding rather than numerical
# estimates, because the refactor should preserve the shape of results exactly.

make_refactor_guard_data <- function(n = 80) {
  # Use a small synthetic dataset so these guard tests exercise the real fitting
  # code without paying the cost of the larger package example data.
  set.seed(42)
  dose <- seq(0, 1, length.out = n)
  X1 <- rnorm(n)
  M1 <- rbinom(n, size = 1, prob = 0.5)

  data.frame(
    Y.gaussian = 1 + 2 * dose + 0.4 * X1 + rnorm(n, sd = 0.25),
    Y.binomial = rbinom(n, size = 1, prob = plogis(-0.2 + 0.7 * dose)),
    Y.poisson = pmax(0, round(exp(0.2 + 0.6 * dose + 0.2 * X1))),
    Y.multinomial = factor(
      sample(c("low", "mid", "high"), size = n, replace = TRUE),
      levels = c("low", "mid", "high")
    ),
    Y.clogit = rep(c(1, 0), length.out = n),
    setnr = rep(seq_len(n / 2), each = 2),
    X1 = X1,
    M1 = M1,
    Mzero = 0,
    V1 = dose + rnorm(n, sd = 0.03),
    V2 = dose + rnorm(n, sd = 0.03)
  )
}

fit_refactor_guard <- function(formula, data, family, methods, ...) {
  # Suppress routine fit messages/warnings: these tests are checking structural
  # result contracts that should not depend on optimizer chatter.
  suppressWarnings(
    suppressMessages(
      ameras(
        formula,
        data = data,
        family = family,
        methods = methods,
        optim.method = "BFGS",
        ...
      )
    )
  )
}

expect_method_parameter_names <- function(fit, method, expected) {
  # RC, ERC, and MCML should name coefficients, standard errors, and both vcov
  # dimensions consistently. A shared naming helper must keep all four aligned.
  expect_identical(names(fit[[method]]$coefficients), expected)
  expect_identical(names(fit[[method]]$sd), expected)
  expect_identical(rownames(fit[[method]]$vcov), expected)
  expect_identical(colnames(fit[[method]]$vcov), expected)
}

expect_methods_parameter_names <- function(fit, methods, expected) {
  for (method in methods) {
    expect_method_parameter_names(fit, method, expected)
  }
}

test_that("RC, ERC, and MCML Gaussian names include modifiers and sigma", {
  guard_data <- make_refactor_guard_data()
  methods <- c("RC", "ERC", "MCML")

  fit <- fit_refactor_guard(
    Y.gaussian ~ dose(V1:V2, deg = 2, modifier = M1),
    data = guard_data,
    family = "gaussian",
    methods = methods
  )

  # Gaussian adds sigma after the regression and dose-response parameters.
  expect_methods_parameter_names(
    fit,
    methods,
    c(
      "(Intercept)",
      "dose",
      "dose_squared",
      "dose:M1",
      "dose_squared:M1",
      "sigma"
    )
  )
})

test_that("RC, ERC, and MCML LINEXP names are set with and without modifiers", {
  guard_data <- make_refactor_guard_data()
  methods <- c("RC", "ERC", "MCML")

  fit_without_modifier <- fit_refactor_guard(
    Y.binomial ~ dose(V1:V2, model = LINEXP),
    data = guard_data,
    family = "binomial",
    methods = methods
  )
  expect_methods_parameter_names(
    fit_without_modifier,
    methods,
    c("(Intercept)", "dose_linear", "dose_exponential")
  )

  fit_with_modifier <- fit_refactor_guard(
    Y.binomial ~ dose(V1:V2, model = LINEXP, modifier = Mzero),
    data = guard_data,
    family = "binomial",
    methods = methods
  )

  # Mzero is deliberately inert numerically, but it makes M non-null so the
  # LINEXP modifier-name branch is exercised for all three methods.
  expect_methods_parameter_names(
    fit_with_modifier,
    methods,
    c(
      "(Intercept)",
      "dose_linear",
      "dose_exponential",
      "dose_linear:Mzero",
      "dose_exponential:Mzero"
    )
  )
})

test_that("RC and ERC Poisson fits keep the optimizer output contract", {
  guard_data <- make_refactor_guard_data()

  fit <- fit_refactor_guard(
    Y.poisson ~ dose(V1:V2, model = EXP, deg = 2),
    data = guard_data,
    family = "poisson",
    methods = c("RC", "ERC")
  )

  # ERC Poisson precomputes realization-specific residual variance terms
  # before optimization. This guards that special setup while the optimizer
  # itself is moved to the shared helper.
  expect_methods_parameter_names(
    fit,
    c("RC", "ERC"),
    c("(Intercept)", "dose", "dose_squared")
  )
  for (method in c("RC", "ERC")) {
    expect_equal(fit[[method]]$optim$convergence, 0)
    expect_named(fit[[method]]$optim$counts, c("function", "gradient"))
    expect_identical(dim(fit[[method]]$optim$hessian), c(3L, 3L))
  }
})

test_that("RC and MCML no-intercept family names start with dose terms", {
  guard_data <- make_refactor_guard_data()
  methods <- c("RC", "MCML")

  fit <- fit_refactor_guard(
    Y.clogit ~ dose(V1:V2, model = EXP, deg = 2, modifier = M1) +
      strata(setnr),
    data = guard_data,
    family = "clogit",
    methods = methods
  )

  # Conditional logistic models do not include a global intercept. This guards
  # the family-specific branch that should not prepend "(Intercept)".
  expect_methods_parameter_names(
    fit,
    methods,
    c("dose", "dose_squared", "dose:M1", "dose_squared:M1")
  )
})

test_that("MCML one-parameter fits keep the scalar optimizer contract", {
  guard_data <- make_refactor_guard_data()

  fit <- fit_refactor_guard(
    Y.clogit ~ dose(V1:V2, model = EXP) + strata(setnr),
    data = guard_data,
    family = "clogit",
    methods = "MCML"
  )

  # Clogit has no intercept, and with deg = 1, no X, and no M there is only
  # one fitted parameter. MCML should therefore keep using the scalar optimizer
  # path, whose result has no `optim()` iteration counts.
  expect_method_parameter_names(fit, "MCML", "dose")
  expect_equal(fit$MCML$optim$convergence, 0)
  expect_null(fit$MCML$optim$counts)
  expect_identical(dim(fit$MCML$optim$hessian), c(1L, 1L))
})

test_that("RC one-parameter fits keep the scalar optimizer contract", {
  guard_data <- make_refactor_guard_data()

  fit <- fit_refactor_guard(
    Y.clogit ~ dose(V1:V2, model = EXP) + strata(setnr),
    data = guard_data,
    family = "clogit",
    methods = "RC"
  )

  # This is the same one-parameter conditional logistic setup as the MCML
  # guard, but it protects RC's scalar `optimize()` path during refactoring.
  expect_method_parameter_names(fit, "RC", "dose")
  expect_equal(fit$RC$optim$convergence, 0)
  expect_null(fit$RC$optim$counts)
  expect_identical(dim(fit$RC$optim$hessian), c(1L, 1L))
})

test_that("RC and MCML multinomial names include response-level prefixes", {
  guard_data <- make_refactor_guard_data()
  methods <- c("RC", "MCML")

  fit <- fit_refactor_guard(
    Y.multinomial ~ dose(V1:V2, model = LINEXP, modifier = Mzero),
    data = guard_data,
    family = "multinomial",
    methods = methods
  )

  # Multinomial fits repeat the same base parameter names for each non-reference
  # response level. The final factor level is the reference and is omitted.
  base_names <- c(
    "(Intercept)",
    "dose_linear",
    "dose_exponential",
    "dose_linear:Mzero",
    "dose_exponential:Mzero"
  )
  non_reference_levels <- levels(guard_data$Y.multinomial)
  non_reference_levels <- non_reference_levels[-length(non_reference_levels)]
  expected <- unlist(
    lapply(non_reference_levels, function(level) {
      paste0("(", level, ")_", base_names)
    }),
    use.names = FALSE
  )

  expect_methods_parameter_names(fit, methods, expected)
})

test_that("RC, ERC, and MCML forward transform arguments through fitting", {
  # The custom transform and Jacobian require marker = "rc-mcml-transform-ok".
  # If `...` is dropped during optimization, Hessian calculation, final
  # transformation, or Jacobian variance propagation, the fit errors here.
  require_marker <- function(marker) {
    if (!identical(marker, "rc-mcml-transform-ok")) {
      stop("transform marker was not forwarded")
    }
  }
  transform_sigma <- function(params, marker, ...) {
    require_marker(marker)
    transform1(params, index.t = 3)
  }
  transform_sigma_jacobian <- function(params, marker, ...) {
    require_marker(marker)
    transform1.jacobian(params, index.t = 3)
  }

  guard_data <- make_refactor_guard_data()

  fit <- fit_refactor_guard(
    Y.gaussian ~ dose(V1:V2),
    data = guard_data,
    family = "gaussian",
    methods = c("RC", "ERC", "MCML"),
    transform = transform_sigma,
    transform.jacobian = transform_sigma_jacobian,
    marker = "rc-mcml-transform-ok"
  )

  for (method in c("RC", "ERC", "MCML")) {
    expect_named(fit[[method]]$coefficients, c("(Intercept)", "dose", "sigma"))
    expect_true(all(is.finite(fit[[method]]$coefficients)))
    expect_true(all(is.finite(fit[[method]]$vcov)))
  }
})
