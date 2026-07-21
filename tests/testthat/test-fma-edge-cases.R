# These tests exercise FMA edge cases that are easy to miss in broad
# snapshot tests: extreme AIC differences, realizations receiving zero
# simulated samples after weight rounding, and warning paths for very few
# retained realizations.

capture_messages <- function(expr) {
  messages <- character()
  value <- withCallingHandlers(
    expr,
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  list(value = value, messages = messages)
}

make_extreme_aic_data <- function(two_good = TRUE) {
  # Build a small Gaussian dataset where one dose realization is almost
  # perfectly predictive and the remaining realizations are mostly noise.
  #
  # When two_good = TRUE, V1 and V2 both track the true dose. FMA should retain
  # two realizations and down-weight the noisy realizations to zero samples.
  #
  # When two_good = FALSE, only V1 tracks the true dose. FMA should retain one
  # realization and trigger the single-realization warning.
  set.seed(1)
  n <- 200
  dose <- seq(-1, 1, length.out = n)
  data.frame(
    Y = 1 + 2 * dose + rnorm(n, sd = 1e-4),
    V1 = dose,
    V2 = if (two_good) dose + rnorm(n, sd = 1e-5) else rnorm(n),
    V3 = rnorm(n),
    V4 = rnorm(n),
    V5 = rnorm(n)
  )
}

test_that("ameras forwards optim.method to FMA", {
  data("data", package = "ameras")
  captured <- NULL

  testthat::with_mocked_bindings(
    ameras.fma = function(...) {
      captured <<- list(...)
      list(
        coefficients = c("(Intercept)" = 0, dose = 0, sigma = 1),
        sd = c("(Intercept)" = 0, dose = 0, sigma = 0),
        vcov = diag(3),
        samples = NULL,
        included.realizations = 1:2,
        weights = c(V1 = 0.5, V2 = 0.5),
        runtime = "0 seconds"
      )
    },
    {
      fit <- suppressMessages(
        ameras(
          Y.gaussian ~ dose(V1:V2),
          data = data[1:20, ],
          family = "gaussian",
          methods = "FMA",
          optim.method = "BFGS"
        )
      )
    }
  )

  expect_identical(captured$optim.method, "BFGS")
  expect_identical(fit$FMA$included.realizations, 1:2)
})

test_that("FMA weights are finite with extreme AIC differences", {
  skip_on_cran()

  # This reproduces the numerical shape behind the original NaN-weight bug:
  # a few realizations have much smaller AIC values than the rest.
  edge_data <- make_extreme_aic_data(two_good = TRUE)

  fit <- NULL
  captured <- capture_messages(
    expect_warning(
      fit <- ameras(
        Y ~ dose(V1:V5),
        data = edge_data,
        family = "gaussian",
        methods = "FMA",
        MFMA = 1000
      ),
      "FMA is based on only 2 realizations"
    )
  )

  expect_true(any(grepl(
    "excluded due to negligible model averaging weight",
    captured$messages,
    fixed = TRUE
  )))
  expect_identical(fit$FMA$included.realizations, c(1L, 2L))
  expect_equal(fit$FMA$included.samples, 1000)
  expect_equal(nrow(fit$FMA$samples), 1000)
  expect_true(all(is.finite(fit$FMA$weights)))
  expect_false(any(is.nan(fit$FMA$weights)))
  expect_equal(sum(fit$FMA$weights), 1, tolerance = 1e-12)
})

test_that("FMA warns when weight filtering leaves a single realization", {
  skip_on_cran()

  # With only one good realization, the noisy realizations are allocated zero
  # samples and FMA should warn that model averaging uncertainty is not present.
  edge_data <- make_extreme_aic_data(two_good = FALSE)

  fit <- NULL
  captured <- capture_messages(
    expect_warning(
      fit <- ameras(
        Y ~ dose(V1:V5),
        data = edge_data,
        family = "gaussian",
        methods = "FMA",
        MFMA = 1000
      ),
      "FMA is based on a single realization"
    )
  )

  expect_true(any(grepl(
    "excluded due to negligible model averaging weight",
    captured$messages,
    fixed = TRUE
  )))
  expect_identical(fit$FMA$included.realizations, 1L)
  expect_equal(fit$FMA$included.samples, 1000)
  expect_named(fit$FMA$weights, "V1")
  expect_equal(unname(fit$FMA$weights), 1)
})

test_that("FMA returns NA estimates when all weight allocations round to zero", {
  skip_on_cran()

  # Unweighted FMA with MFMA = 1 assigns round(1 / 3) = 0 samples to each
  # realization. The fit should return empty FMA samples and NA summaries.
  set.seed(3)
  n <- 80
  dose <- seq(-1, 1, length.out = n)
  edge_data <- data.frame(
    Y = 1 + 2 * dose + rnorm(n, sd = 0.2),
    V1 = dose + rnorm(n, sd = 0.1),
    V2 = dose + rnorm(n, sd = 0.1),
    V3 = dose + rnorm(n, sd = 0.1)
  )

  fit <- NULL
  captured <- capture_messages(
    expect_warning(
      fit <- ameras(
        Y ~ dose(V1:V3),
        data = edge_data,
        family = "gaussian",
        methods = "FMA",
        unweightedFMA = TRUE,
        MFMA = 1
      ),
      "No realizations contributed to FMA"
    )
  )

  expect_true(any(grepl(
    "3 realization(s) excluded due to negligible model averaging weight",
    captured$messages,
    fixed = TRUE
  )))
  expect_identical(fit$FMA$included.realizations, integer(0))
  expect_equal(fit$FMA$included.samples, 0)
  expect_null(fit$FMA$weights)
  expect_null(fit$FMA$samples)
  expect_true(all(is.na(fit$FMA$coefficients)))
  expect_true(all(is.na(fit$FMA$sd)))
  expect_true(all(is.na(fit$FMA$vcov)))
})

test_that("FMA forwards transform arguments through fitting and sampling", {
  skip_on_cran()

  # This catches regressions where the likelihood closure receives `...` during
  # optimization but the Hessian or FMA sampling step does not.
  require_marker <- function(marker) {
    if (!identical(marker, "fma-transform-ok")) {
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

  set.seed(4)
  n <- 80
  dose <- seq(-1, 1, length.out = n)
  transform_data <- data.frame(
    Y = 1 + 2 * dose + rnorm(n, sd = 0.2),
    V1 = dose + rnorm(n, sd = 0.02),
    V2 = dose + rnorm(n, sd = 0.02)
  )

  fit <- suppressMessages(
    suppressWarnings(
      ameras(
        Y ~ dose(V1:V2),
        data = transform_data,
        family = "gaussian",
        methods = "FMA",
        MFMA = 100,
        transform = transform_sigma,
        transform.jacobian = transform_sigma_jacobian,
        marker = "fma-transform-ok"
      )
    )
  )

  expect_named(fit$FMA$coefficients, c("(Intercept)", "dose", "sigma"))
  expect_true(all(is.finite(fit$FMA$coefficients)))
  expect_true(all(is.finite(fit$FMA$vcov)))
})

expect_fma_parameter_names <- function(fit, expected) {
  # FMA names are copied into several parallel result objects. Checking all of
  # them catches cases where only one summary table was named correctly.
  expect_identical(names(fit$FMA$coefficients), expected)
  expect_identical(names(fit$FMA$sd), expected)
  expect_identical(rownames(fit$FMA$vcov), expected)
  expect_identical(colnames(fit$FMA$vcov), expected)

  # Degenerate FMA fits can have no sampled coefficients, but ordinary fits
  # should carry the same parameter names on the simulated coefficient draws.
  if (!is.null(fit$FMA$samples)) {
    expect_identical(names(fit$FMA$samples), expected)
  }
}

fit_fma_quietly <- function(formula, data, family) {
  # These tests only care about names. Suppressing routine FMA warnings/messages
  # keeps the checks focused on parname assignment.
  suppressWarnings(
    suppressMessages(
      ameras(
        formula,
        data = data,
        family = family,
        methods = "FMA",
        MFMA = 100
      )
    )
  )
}

test_that("FMA LINEXP parameter names are set with and without modifiers", {
  skip_on_cran()

  data("data", package = "ameras")
  linexp_data <- data
  # A zero modifier is intentionally inert: it makes M non-null so the modifier
  # parname branch is exercised without changing the fitted LINEXP dose effect.
  linexp_data$Mzero <- 0

  # GLM-style families include an intercept. Conditional logistic and
  # proportional hazards FMA branches do not add one to the dose parnames.
  expected_with_intercept <- c(
    "(Intercept)",
    "dose_linear",
    "dose_exponential"
  )
  expected_with_intercept_m <- c(
    expected_with_intercept,
    "dose_linear:Mzero",
    "dose_exponential:Mzero"
  )
  expected_no_intercept <- c("dose_linear", "dose_exponential")
  expected_no_intercept_m <- c(
    expected_no_intercept,
    "dose_linear:Mzero",
    "dose_exponential:Mzero"
  )

  # Binomial and Poisson share the same LINEXP parname shape.
  expect_fma_parameter_names(
    fit_fma_quietly(
      Y.binomial ~ dose(V1:V2, model = LINEXP),
      data = linexp_data,
      family = "binomial"
    ),
    expected_with_intercept
  )
  expect_fma_parameter_names(
    fit_fma_quietly(
      Y.binomial ~ dose(V1:V2, model = LINEXP, modifier = ~ Mzero),
      data = linexp_data,
      family = "binomial"
    ),
    expected_with_intercept_m
  )

  expect_fma_parameter_names(
    fit_fma_quietly(
      Y.poisson ~ dose(V1:V2, model = LINEXP),
      data = linexp_data,
      family = "poisson"
    ),
    expected_with_intercept
  )
  expect_fma_parameter_names(
    fit_fma_quietly(
      Y.poisson ~ dose(V1:V2, model = LINEXP, modifier = ~ Mzero),
      data = linexp_data,
      family = "poisson"
    ),
    expected_with_intercept_m
  )

  # Survival-style families exercise separate no-intercept LINEXP branches.
  expect_fma_parameter_names(
    fit_fma_quietly(
      Y.clogit ~ dose(V1:V2, model = LINEXP) + strata(setnr),
      data = linexp_data,
      family = "clogit"
    ),
    expected_no_intercept
  )
  expect_fma_parameter_names(
    fit_fma_quietly(
      Y.clogit ~ dose(V1:V2, model = LINEXP, modifier = ~ Mzero) + strata(setnr),
      data = linexp_data,
      family = "clogit"
    ),
    expected_no_intercept_m
  )

  expect_fma_parameter_names(
    fit_fma_quietly(
      Surv(time, status) ~ dose(V1:V2, model = LINEXP),
      data = linexp_data,
      family = "prophaz"
    ),
    expected_no_intercept
  )
  expect_fma_parameter_names(
    fit_fma_quietly(
      Surv(time, status) ~ dose(V1:V2, model = LINEXP, modifier = ~ Mzero),
      data = linexp_data,
      family = "prophaz"
    ),
    expected_no_intercept_m
  )
})

test_that("FMA non-LINEXP parameter names include modifier terms", {
  skip_on_cran()

  data("data", package = "ameras")

  # EXP reaches the doseRRmod != "LINEXP" parname branch without the lower-bound
  # fragility of ERR models. deg = 2 also covers the dose_squared modifier names.
  expected_with_intercept_m <- c(
    "(Intercept)",
    "dose",
    "dose_squared",
    "dose:M1",
    "dose_squared:M1"
  )
  expected_gaussian_m <- c(expected_with_intercept_m, "sigma")
  expected_no_intercept_m <- c(
    "dose",
    "dose_squared",
    "dose:M1",
    "dose_squared:M1"
  )

  expect_fma_parameter_names(
    fit_fma_quietly(
      Y.gaussian ~ dose(V1:V2, deg = 2, modifier = ~ M1),
      data = data,
      family = "gaussian"
    ),
    expected_gaussian_m
  )

  expect_fma_parameter_names(
    fit_fma_quietly(
      Y.binomial ~ dose(V1:V2, model = EXP, deg = 2, modifier = ~ M1),
      data = data,
      family = "binomial"
    ),
    expected_with_intercept_m
  )

  expect_fma_parameter_names(
    fit_fma_quietly(
      Y.poisson ~ dose(V1:V2, model = EXP, deg = 2, modifier = ~ M1),
      data = data,
      family = "poisson"
    ),
    expected_with_intercept_m
  )

  expect_fma_parameter_names(
    fit_fma_quietly(
      Y.clogit ~ dose(V1:V2, model = EXP, deg = 2, modifier = ~ M1) +
        strata(setnr),
      data = data,
      family = "clogit"
    ),
    expected_no_intercept_m
  )

  expect_fma_parameter_names(
    fit_fma_quietly(
      Surv(time, status) ~ dose(V1:V2, model = EXP, deg = 2, modifier = ~ M1),
      data = data,
      family = "prophaz"
    ),
    expected_no_intercept_m
  )
})

test_that("FMA multinomial LINEXP parameter names include level prefixes", {
  skip_on_cran()

  data("data", package = "ameras")
  linexp_data <- data
  # As above, this safely reaches the non-null M branch for LINEXP parnames.
  linexp_data$Mzero <- 0

  expected_base <- c(
    "(Intercept)",
    "dose_linear",
    "dose_exponential"
  )
  expected_base_m <- c(
    expected_base,
    "dose_linear:Mzero",
    "dose_exponential:Mzero"
  )
  # Multinomial FMA repeats the base parameter names for every modeled response
  # level, prefixing each block with the level label. The final factor level is
  # the reference level and is therefore omitted.
  non_reference_levels <- levels(linexp_data$Y.multinomial)
  non_reference_levels <- non_reference_levels[-length(non_reference_levels)]

  expect_fma_parameter_names(
    fit_fma_quietly(
      Y.multinomial ~ dose(V1:V2, model = LINEXP),
      data = linexp_data,
      family = "multinomial"
    ),
    unlist(
      lapply(
        non_reference_levels,
        function(level) paste0("(", level, ")_", expected_base)
      ),
      use.names = FALSE
    )
  )
  expect_fma_parameter_names(
    fit_fma_quietly(
      Y.multinomial ~ dose(V1:V2, model = LINEXP, modifier = ~ Mzero),
      data = linexp_data,
      family = "multinomial"
    ),
    unlist(
      lapply(
        non_reference_levels,
        function(level) paste0("(", level, ")_", expected_base_m)
      ),
      use.names = FALSE
    )
  )
})

test_that("FMA multinomial non-LINEXP names include modifiers and prefixes", {
  skip_on_cran()

  data("data", package = "ameras")

  expected_base_m <- c(
    "(Intercept)",
    "dose",
    "dose_squared",
    "dose:M1",
    "dose_squared:M1"
  )
  non_reference_levels <- levels(data$Y.multinomial)
  non_reference_levels <- non_reference_levels[-length(non_reference_levels)]

  expect_fma_parameter_names(
    fit_fma_quietly(
      Y.multinomial ~ dose(V1:V2, model = EXP, deg = 2, modifier = ~ M1),
      data = data,
      family = "multinomial"
    ),
    unlist(
      lapply(
        non_reference_levels,
        function(level) paste0("(", level, ")_", expected_base_m)
      ),
      use.names = FALSE
    )
  )
})
