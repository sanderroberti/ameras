make_dose_lrt_gaussian_data <- function(n = 50) {
  D <- seq(-1.5, 1.5, length.out = n)
  X <- rep(c(-0.5, 0.5), length.out = n)
  data.frame(
    Y = 1 + 0.4 * X + 0.7 * D - 0.25 * D^2 + 0.1 * sin(seq_len(n)),
    D = D,
    X = X
  )
}


normal_ml_loglik <- function(formula, data) {
  X <- model.matrix(formula, data)
  y <- model.response(model.frame(formula, data))
  fit <- lm.fit(X, y)
  rss <- sum(fit$residuals^2)
  n <- length(y)
  sigma2 <- rss / n
  -0.5 * n * (log(2 * pi * sigma2) + 1)
}


make_dose_lrt_binomial_data <- function(n = 80) {
  set.seed(1001)
  D <- seq(0.05, 2.4, length.out = n)
  X <- rep(c(0, 1), length.out = n)
  eta <- -1 + 0.6 * X + 0.7 * D
  data.frame(
    Y = rbinom(n, size = 1, prob = plogis(eta)),
    D = D,
    X = X
  )
}


test_that("dose parameter selection covers supported parameter names", {
  expect_identical(
    ameras:::dose_parameter_indices(c(
      "(Intercept)", "X1", "dose", "dose_squared", "dose:M1",
      "dose_squared:M1", "sigma", "h0[1]", "dose_age"
    )),
    3:6
  )

  expect_identical(
    ameras:::dose_parameter_indices(c(
      "X1", "dose_linear", "dose_exponential", "dose_linear:M1",
      "dose_exponential:M1"
    )),
    2:5
  )

  expect_identical(
    ameras:::dose_parameter_indices(c(
      "(low)_(Intercept)", "(low)_X1", "(low)_dose",
      "(low)_dose_squared", "(mid)_dose:M1", "(mid)_sigma"
    )),
    3:5
  )
})


test_that("reported-scale likelihood helper preserves fitted likelihoods", {
  dat <- make_dose_lrt_binomial_data()
  fit_untransformed <- suppressMessages(
    ameras(
      Y ~ dose(D, model = EXP) + X,
      data = dat,
      family = "binomial",
      methods = "RC",
      print = FALSE
    )
  )

  original_fn <- ameras:::make_loglik_fn(
    fit_untransformed,
    "RC",
    fit_untransformed$RC,
    fit_untransformed$model$data
  )
  reported_fn <- ameras:::make_reported_scale_loglik_fn(
    fit_untransformed,
    "RC",
    fit_untransformed$RC,
    fit_untransformed$model$data
  )
  expect_equal(
    reported_fn(fit_untransformed$RC$optim$par),
    original_fn(fit_untransformed$RC$optim$par),
    tolerance = 1e-10
  )

  dat_gaussian <- make_dose_lrt_gaussian_data()
  fit_transformed <- suppressMessages(
    ameras(
      Y ~ dose(D) + X,
      data = dat_gaussian,
      family = "gaussian",
      methods = "RC",
      print = FALSE
    )
  )
  original_fn <- ameras:::make_loglik_fn(
    fit_transformed,
    "RC",
    fit_transformed$RC,
    fit_transformed$model$data
  )
  reported_fn <- ameras:::make_reported_scale_loglik_fn(
    fit_transformed,
    "RC",
    fit_transformed$RC,
    fit_transformed$model$data
  )
  reported_par <- ameras:::apply_fitted_transform(
    fit_transformed,
    fit_transformed$RC$optim$par
  )

  expect_equal(
    reported_fn(reported_par),
    original_fn(fit_transformed$RC$optim$par),
    tolerance = 1e-8
  )
})


test_that("constrained objective fixes reported dose parameters", {
  dat <- make_dose_lrt_binomial_data()
  fit_untransformed <- suppressMessages(
    ameras(
      Y ~ dose(D, model = EXP) + X,
      data = dat,
      family = "binomial",
      methods = "RC",
      print = FALSE
    )
  )
  dose_idx <- ameras:::dose_parameter_indices(names(fit_untransformed$RC$coefficients))
  objective <- ameras:::make_constrained_lrt_objective(
    fit_untransformed,
    "RC",
    fit_untransformed$RC,
    fit_untransformed$model$data,
    fixed.idx = dose_idx
  )
  shifted_free <- objective$start + 0.1
  constrained_par <- objective$parameter_fn(shifted_free)

  expect_equal(unname(constrained_par[dose_idx]), 0)
  expect_equal(
    unname(constrained_par[objective$free.idx]),
    unname(shifted_free),
    tolerance = 1e-12
  )

  dat_gaussian <- make_dose_lrt_gaussian_data()
  fit_transformed <- suppressMessages(
    ameras(
      Y ~ dose(D) + X,
      data = dat_gaussian,
      family = "gaussian",
      methods = "RC",
      print = FALSE
    )
  )
  dose_idx <- ameras:::dose_parameter_indices(names(fit_transformed$RC$coefficients))
  objective <- ameras:::make_constrained_lrt_objective(
    fit_transformed,
    "RC",
    fit_transformed$RC,
    fit_transformed$model$data,
    fixed.idx = dose_idx
  )
  constrained_par <- objective$parameter_fn(objective$start)

  expect_equal(unname(constrained_par[dose_idx]), 0)
  expect_true(constrained_par["sigma"] > 0)

  fit_linexp <- suppressMessages(
    ameras(
      Y ~ dose(D, model = LINEXP) + X,
      data = dat,
      family = "binomial",
      methods = "RC",
      print = FALSE
    )
  )
  dose_idx <- ameras:::dose_parameter_indices(names(fit_linexp$RC$coefficients))
  objective <- ameras:::make_constrained_lrt_objective(
    fit_linexp,
    "RC",
    fit_linexp$RC,
    fit_linexp$model$data,
    fixed.idx = dose_idx[1]
  )

  expect_equal(unname(objective$parameter_fn(objective$start)[dose_idx[1]]), 0)
  expect_true(is.finite(objective$fn(objective$start)))
})


test_that("constrained LRT optimizer handles nuisance and no-nuisance cases", {
  objective <- list(
    start = c(x = 0, y = 0),
    fn = function(par) (par[1] - 1)^2 + (par[2] + 2)^2
  )
  fit <- ameras:::fit_constrained_lrt_objective(objective)

  expect_equal(unname(fit$par), c(1, -2), tolerance = 1e-6)
  expect_equal(fit$value, 0, tolerance = 1e-8)
  expect_null(fit$hessian)

  direct <- ameras:::fit_constrained_lrt_objective(list(
    start = numeric(0),
    fn = function(par) 3.5
  ))
  expect_identical(direct$par, numeric(0))
  expect_equal(direct$value, 3.5)
  expect_equal(direct$convergence, 0)
})


test_that("dose LRT rows report constrained-null convergence", {
  out <- NULL
  expect_warning(
    out <- ameras:::dose_lrt_row(
      method_name = "RC",
      type = "global",
      term = "all dose terms",
      df = 1,
      full_loglik = 0,
      null_fit = list(value = 1, convergence = 1L)
    ),
    "Constrained null fit"
  )

  expect_named(out, c(
    "method", "type", "term", "df", "logLik", "logLik.null",
    "statistic", "p.value", "null.optim.convergence"
  ))
  expect_equal(out$null.optim.convergence, 1L)
})


test_that("dose_lrt matches manual Gaussian RC likelihood-ratio tests", {
  dat <- make_dose_lrt_gaussian_data()
  fit <- suppressMessages(
    ameras(
      Y ~ dose(D, deg = 2) + X,
      data = dat,
      family = "gaussian",
      methods = "RC",
      print = FALSE
    )
  )

  out <- dose_lrt(fit)
  full_ll <- normal_ml_loglik(Y ~ X + D + I(D^2), dat)
  global_ll <- normal_ml_loglik(Y ~ X, dat)
  dose_ll <- normal_ml_loglik(Y ~ X + I(D^2), dat)
  dose_sq_ll <- normal_ml_loglik(Y ~ X + D, dat)

  expect_equal(out$type, c("global", "individual", "individual"))
  expect_equal(out$term, c("all dose terms", "dose", "dose_squared"))
  expect_equal(out$df, c(2, 1, 1))
  expect_equal(out$statistic[1], 2 * (full_ll - global_ll), tolerance = 1e-4)
  expect_equal(out$statistic[2], 2 * (full_ll - dose_ll), tolerance = 1e-4)
  expect_equal(out$statistic[3], 2 * (full_ll - dose_sq_ll), tolerance = 1e-4)
})


test_that("dose_lrt matches nested glm tests for binomial EXP RC", {
  dat <- make_dose_lrt_binomial_data()
  fit <- suppressMessages(
    ameras(
      Y ~ dose(D, model = EXP) + X,
      data = dat,
      family = "binomial",
      methods = "RC",
      print = FALSE
    )
  )

  out <- dose_lrt(fit)
  glm_full <- glm(Y ~ X + D, data = dat, family = binomial())
  glm_null <- glm(Y ~ X, data = dat, family = binomial())
  expected <- 2 * (as.numeric(logLik(glm_full)) - as.numeric(logLik(glm_null)))

  expect_equal(out$type, c("global", "individual"))
  expect_equal(out$term, c("all dose terms", "dose"))
  expect_equal(out$statistic, c(expected, expected), tolerance = 1e-4)
  expect_equal(out$p.value, pchisq(out$statistic, df = 1, lower.tail = FALSE))
})


test_that("dose_lrt handles transformed ERR and boundary LINEXP fits", {
  dat <- make_dose_lrt_binomial_data()

  fit_err <- suppressMessages(
    ameras(
      Y ~ dose(D, model = ERR) + X,
      data = dat,
      family = "binomial",
      methods = "RC",
      print = FALSE
    )
  )
  err_lrt <- dose_lrt(fit_err)
  expect_true(all(is.finite(err_lrt$statistic)))
  expect_true(all(is.finite(err_lrt$p.value)))

  fit_linexp <- suppressMessages(
    ameras(
      Y ~ dose(D, model = LINEXP) + X,
      data = dat,
      family = "binomial",
      methods = "RC",
      print = FALSE
    )
  )
  linexp_lrt <- dose_lrt(fit_linexp)
  expect_equal(
    linexp_lrt$term,
    c("all dose terms", "dose_linear", "dose_exponential")
  )
  expect_true(all(is.finite(linexp_lrt$statistic)))
  expect_true(all(is.finite(linexp_lrt$p.value)))
})


test_that("dose_lrt smoke-tests ERC and MCML", {
  n <- 24
  D <- seq(-1, 1, length.out = n)
  dat <- data.frame(
    Y = 1 + 0.6 * D + 0.1 * sin(seq_len(n)),
    V1 = D - 0.05,
    V2 = D,
    V3 = D + 0.05
  )

  fit <- suppressWarnings(suppressMessages(
    ameras(
      Y ~ dose(V1:V3),
      data = dat,
      family = "gaussian",
      methods = c("ERC", "MCML"),
      print = FALSE
    )
  ))
  out <- dose_lrt(fit, methods = c("ERC", "MCML"), type = "global")

  expect_equal(out$method, c("ERC", "MCML"))
  expect_equal(out$type, c("global", "global"))
  expect_true(all(is.finite(out$statistic)))
  expect_true(all(is.finite(out$p.value)))
})


test_that("dose_lrt supports type/method filtering and keep.data FALSE", {
  dat <- make_dose_lrt_gaussian_data()
  fit <- suppressMessages(
    ameras(
      Y ~ dose(D, deg = 2) + X,
      data = dat,
      family = "gaussian",
      methods = "RC",
      keep.data = FALSE,
      print = FALSE
    )
  )

  expect_error(dose_lrt(fit), "Data not stored")

  global <- dose_lrt(fit, data = dat, type = "global")
  individual <- dose_lrt(fit, data = dat, type = "individual")

  expect_equal(global$type, "global")
  expect_equal(individual$type, c("individual", "individual"))
  expect_equal(unique(global$method), "RC")
  expect_equal(unique(individual$method), "RC")
})


test_that("dose_lrt errors clearly for unsupported methods", {
  fake <- structure(
    list(FMA = list(), model = list(data = data.frame())),
    class = "amerasfit"
  )

  expect_error(dose_lrt(fake, methods = "FMA"), "supports RC, ERC, and MCML")
  expect_error(
    dose_lrt(fake),
    "None of the requested supported methods are present"
  )
})
