test_that("all FMA family dispatches forward loglim", {
  seen <- numeric()
  testthat::local_mocked_bindings(
    fit_fma_realizations = function(..., loglim, inpar, dosevars) {
      seen <<- c(seen, loglim)
      lapply(dosevars, function(dose) list(
        include = TRUE, AIC = 0,
        sampling = list(mean = rep(0, length(inpar)), var = diag(length(inpar)))
      ))
    }, .package = "ameras"
  )
  formulas <- list(
    gaussian = Y.gaussian ~ dose(V1:V3),
    binomial = Y.binomial ~ dose(V1:V3, model = "EXP"),
    poisson = Y.poisson ~ dose(V1:V3, model = "EXP"),
    multinomial = Y.multinomial ~ dose(V1:V3, model = "EXP"),
    prophaz = Surv(time, status) ~ dose(V1:V3, model = "EXP"),
    clogit = Y.clogit ~ dose(V1:V3, model = "EXP") + strata(setnr)
  )
  for (family in names(formulas)) {
    fit <- suppressWarnings(suppressMessages(ameras(
      formulas[[family]], data = ameras_simulated[1:100, ], family = family,
      methods = "FMA", loglim = 0.1, MFMA = 30
    )))
    expect_identical(fit$model$loglim, 0.1)
  }
  expect_equal(seen, rep(0.1, length(formulas)))
})

test_that("FMA realization likelihoods receive the requested floor", {
  seen <- numeric()
  original <- ameras:::make_single_realization_loglik
  testthat::local_mocked_bindings(
    make_single_realization_loglik = function(..., loglim = 1e-30) {
      seen <<- c(seen, loglim)
      original(..., loglim = loglim)
    }, .package = "ameras"
  )
  suppressWarnings(suppressMessages(ameras(
    Y.poisson ~ dose(V1:V3, model = "EXP"), data = ameras_simulated[1:100, ],
    family = "poisson", methods = "FMA", loglim = 0.1, MFMA = 30
  )))
  expect_equal(seen, rep(0.1, 3))
  likelihood <- original(
    family = "poisson", dose.col = "D", data = data.frame(Y = rep(0, 30), D = 1:30),
    Y = "Y", doseRRmod = "EXP", loglim = 0.1
  )
  expect_equal(likelihood(c(-10, 0)), 3)
})
