test_that("MCML aggregation helper preserves clipped log-mean-exp behavior", {
  # MCML averages likelihoods across dose realizations on a stabilized scale.
  # This locks down the existing clipping rule so future refactors do not
  # accidentally replace it with a subtly different log-sum-exp expression.
  logliks <- c(-100, -2, 10, 120)
  shifted <- pmax(pmin(-1 * logliks - max(-1 * logliks), 7e1), -7e1)
  expected <- -1 * log(mean(exp(shifted))) - max(-1 * logliks)

  expect_equal(ameras:::mcml_average_neg_loglik(logliks), expected)

  # For moderate values, the stabilized expression is equivalent to the
  # ordinary negative log of the mean likelihood.
  moderate <- c(3, 4, 5)
  expect_equal(
    ameras:::mcml_average_neg_loglik(moderate),
    -log(mean(exp(-moderate))),
    tolerance = 1e-12
  )
})

test_that("MCML loglik helper matches direct Gaussian likelihood aggregation", {
  data("data", package = "ameras")
  params <- c(0, 0.1, 1)
  dosevars <- c("V1", "V2")

  helper <- ameras:::make_mcml_loglik_fn(
    family = "gaussian",
    dosevars = dosevars,
    data = data,
    Y = "Y.gaussian"
  )

  direct_logliks <- ameras:::loglik.gaussian(
    params = params,
    D = dosevars,
    X = NULL,
    Y = "Y.gaussian",
    M = NULL,
    data = data,
    deg = 1,
    ERC = FALSE,
    Kmat = NULL,
    loglim = 1e-30,
    transform = NULL
  )

  expect_equal(
    helper(params),
    ameras:::mcml_average_neg_loglik(direct_logliks),
    tolerance = 1e-12
  )
})

test_that("MCML loglik helper forwards runtime transform arguments", {
  data("data", package = "ameras")
  params <- c(0, 0.1, 1)
  shift <- c(0.2, 0.3, 0.4)
  dosevars <- c("V1", "V2")

  # The transform requires a runtime-only argument. If the MCML helper drops
  # `...`, this test errors before the likelihood can be compared.
  transform_shift <- function(params, shift, ...) params + shift
  helper <- ameras:::make_mcml_loglik_fn(
    family = "gaussian",
    dosevars = dosevars,
    data = data,
    Y = "Y.gaussian",
    transform = transform_shift
  )

  direct_logliks <- ameras:::loglik.gaussian(
    params = params + shift,
    D = dosevars,
    X = NULL,
    Y = "Y.gaussian",
    M = NULL,
    data = data,
    deg = 1,
    ERC = FALSE,
    Kmat = NULL,
    loglim = 1e-30,
    transform = NULL
  )

  expect_equal(
    helper(params, shift = shift),
    ameras:::mcml_average_neg_loglik(direct_logliks),
    tolerance = 1e-12
  )
})
