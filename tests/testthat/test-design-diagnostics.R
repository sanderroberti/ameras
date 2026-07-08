test_that("X conditioning diagnostic warns for calendar-year-like covariates", {
  year <- 2000 + seq(0, 5, length.out = 30)

  # Calendar-year variables can be statistically valid but numerically awkward:
  # their mean can be huge relative to their spread, especially with an
  # intercept. The diagnostic should guide users toward centering/scaling rather
  # than changing the fitted data internally.
  expect_warning(
    ameras:::warn_if_poorly_conditioned_X(
      X_matrix = cbind(year = year),
      family = "gaussian"
    ),
    "poorly scaled"
  )
})

test_that("X conditioning diagnostic is quiet for ordinary scaled covariates", {
  X <- cbind(
    age = seq(-2, 2, length.out = 30),
    group = rep(c(0, 1), length.out = 30)
  )

  expect_silent(
    ameras:::warn_if_poorly_conditioned_X(
      X_matrix = X,
      family = "gaussian"
    )
  )
})
