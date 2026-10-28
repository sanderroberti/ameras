make_rc_doses <- function(dose) {
  perturbation <- 0.08 * sin(seq_along(dose) * 1.7)
  data.frame(
    D = dose,
    D1 = dose - perturbation,
    D2 = dose + perturbation
  )
}


expect_matching_loglik <- function(
  ameras_fit,
  standard_fit,
  tolerance = 1e-4
) {
  expect_equal(
    ameras_fit$RC$loglik,
    as.numeric(logLik(standard_fit)),
    tolerance = tolerance
  )
}


test_that("Gaussian RC agrees with lm on the mean dose", {
  n <- 80L
  i <- seq_len(n)
  dat <- cbind(
    make_rc_doses(seq(0.1, 1.8, length.out = n)),
    X = cos(i / 7)
  )
  dat$Y <- 1.5 + 0.8 * dat$D - 0.4 * dat$X + 0.15 * sin(i * 0.9)

  fit <- suppressMessages(ameras(
    Y ~ dose(D1:D2) + X,
    data = dat,
    family = "gaussian",
    methods = "RC"
  ))
  reference <- lm(Y ~ D + X, data = dat)

  expect_equal(
    unname(fit$RC$coefficients[c("(Intercept)", "dose", "X")]),
    unname(coef(reference)[c("(Intercept)", "D", "X")]),
    tolerance = 1e-4
  )
  expect_matching_loglik(fit, reference)
})


test_that("binomial RC with model EXP agrees with logistic regression", {
  set.seed(20261029)
  n <- 180L
  dose <- runif(n, 0.05, 1.8)
  X <- rnorm(n)
  dat <- cbind(make_rc_doses(dose), X = X)
  dat$Y <- rbinom(n, 1, plogis(-0.6 + 0.9 * dose - 0.5 * X))

  fit <- suppressMessages(ameras(
    Y ~ dose(D1:D2, model = "EXP") + X,
    data = dat,
    family = "binomial",
    methods = "RC"
  ))
  reference <- glm(Y ~ D + X, data = dat, family = binomial())

  expect_equal(
    unname(fit$RC$coefficients[c("(Intercept)", "dose", "X")]),
    unname(coef(reference)[c("(Intercept)", "D", "X")]),
    tolerance = 1e-4
  )
  expect_matching_loglik(fit, reference)
})


test_that("Poisson RC with model EXP agrees with Poisson regression", {
  set.seed(20261030)
  n <- 160L
  dose <- runif(n, 0.05, 1.8)
  X <- rnorm(n)
  dat <- cbind(make_rc_doses(dose), X = X)
  dat$Y <- rpois(n, exp(0.2 + 0.45 * dose - 0.3 * X))

  fit <- suppressMessages(ameras(
    Y ~ dose(D1:D2, model = "EXP") + X,
    data = dat,
    family = "poisson",
    methods = "RC"
  ))
  reference <- glm(Y ~ D + X, data = dat, family = poisson())

  expect_equal(
    unname(fit$RC$coefficients[c("(Intercept)", "dose", "X")]),
    unname(coef(reference)[c("(Intercept)", "D", "X")]),
    tolerance = 1e-4
  )
  # ameras normalizes the Poisson likelihood against the saturated model,
  # making its stored value equal to minus half the model deviance.
  expect_equal(
    fit$RC$loglik,
    -deviance(reference) / 2,
    tolerance = 1e-4
  )
})


test_that("proportional-hazards RC agrees with Breslow coxph under ties", {
  set.seed(20261031)
  n <- 180L
  dose <- runif(n, 0.05, 1.8)
  X <- rnorm(n)
  event_time <- rexp(n, rate = 0.12 * exp(0.55 * dose - 0.4 * X))
  censor_time <- rexp(n, rate = 0.07)
  status <- as.integer(event_time <= censor_time)
  time <- pmax(round(pmin(event_time, censor_time), 1), 0.1)
  dat <- cbind(
    make_rc_doses(dose),
    X = X,
    time = time,
    status = status
  )

  expect_true(any(duplicated(dat$time[dat$status == 1])))

  fit <- suppressMessages(ameras(
    Surv(time, status) ~ dose(D1:D2, model = "EXP") + X,
    data = dat,
    family = "prophaz",
    methods = "RC"
  ))
  reference <- survival::coxph(
    survival::Surv(time, status) ~ D + X,
    data = dat,
    ties = "breslow"
  )

  expect_equal(
    unname(fit$RC$coefficients[c("dose", "X")]),
    unname(coef(reference)[c("D", "X")]),
    tolerance = 1e-4
  )
  expect_matching_loglik(fit, reference)
})


test_that("conditional-logistic RC agrees with clogit", {
  set.seed(20261101)
  n_sets <- 60L
  set <- rep(seq_len(n_sets), each = 3L)
  dose <- runif(length(set), 0.05, 1.8)
  X <- rnorm(length(set))
  linear_predictor <- 0.7 * dose - 0.45 * X
  Y <- integer(length(set))
  for (s in seq_len(n_sets)) {
    idx <- which(set == s)
    probabilities <- exp(linear_predictor[idx])
    Y[sample(idx, 1L, prob = probabilities)] <- 1L
  }
  dat <- cbind(make_rc_doses(dose), X = X, Y = Y, set = set)

  fit <- suppressMessages(ameras(
    Y ~ dose(D1:D2, model = "EXP") + X + strata(set),
    data = dat,
    family = "clogit",
    methods = "RC"
  ))
  reference <- local({
    # clogit() rewrites its call to an unqualified coxph() call evaluated in
    # the caller, so expose the namespace function without attaching survival.
    coxph <- survival::coxph
    Surv <- survival::Surv
    survival::clogit(
      Y ~ D + X + survival::strata(set),
      data = dat,
      method = "exact"
    )
  })

  expect_equal(
    unname(fit$RC$coefficients[c("dose", "X")]),
    unname(coef(reference)[c("D", "X")]),
    tolerance = 1e-4
  )
  expect_matching_loglik(fit, reference)
})


test_that("multinomial RC with model EXP agrees with multinom", {
  set.seed(20261102)
  n <- 240L
  dose <- runif(n, 0.05, 1.8)
  X <- rnorm(n)
  eta_a <- 0.4 + 0.55 * dose - 0.3 * X
  eta_b <- -0.2 - 0.35 * dose + 0.45 * X
  denominator <- 1 + exp(eta_a) + exp(eta_b)
  probabilities <- cbind(
    a = exp(eta_a) / denominator,
    b = exp(eta_b) / denominator,
    ref = 1 / denominator
  )
  Y <- apply(probabilities, 1L, function(prob) {
    sample(colnames(probabilities), 1L, prob = prob)
  })
  dat <- cbind(make_rc_doses(dose), X = X)
  dat$Y <- factor(Y, levels = c("a", "b", "ref"))
  dat$Y_reference <- relevel(dat$Y, ref = "ref")

  fit <- suppressMessages(ameras(
    Y ~ dose(D1:D2, model = "EXP") + X,
    data = dat,
    family = "multinomial",
    methods = "RC"
  ))
  reference <- nnet::multinom(
    Y_reference ~ D + X,
    data = dat,
    trace = FALSE,
    maxit = 1000
  )
  reference_coef <- coef(reference)
  expected <- c(
    reference_coef["a", "(Intercept)"],
    reference_coef["a", "D"],
    reference_coef["a", "X"],
    reference_coef["b", "(Intercept)"],
    reference_coef["b", "D"],
    reference_coef["b", "X"]
  )

  expect_equal(
    unname(fit$RC$coefficients[c(
      "(a)_(Intercept)", "(a)_dose", "(a)_X",
      "(b)_(Intercept)", "(b)_dose", "(b)_X"
    )]),
    unname(expected),
    tolerance = 1e-4
  )
  expect_matching_loglik(fit, reference, tolerance = 1e-4)
})
