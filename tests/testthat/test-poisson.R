set.seed(123)
data("data", package = "ameras")


# Poisson, basic model, all methods

for (method in all_methods) {
  test_that(paste("Poisson snapshot:", method), {
    if (method %in% c("ERC", "MCML", "BMA")) {
      skip_on_cran()
    }

    fit <- fit_combination(
      family = "poisson",
      Y = "Y.poisson",
      deg = 2,
      doseRRmod = "EXP",
      X = NULL,
      M = NULL,
      data = data,
      methods = method,
      niter.BMA = 1000,
      nburnin.BMA = 500
    )
    fit <- confint(fit, type = c("wald.orig", "percentile"))
    expect_snapshot_value(
      fit[[method]]$coefficients,
      tolerance = 1e-4,
      style = "deparse"
    )
    expect_snapshot_value(fit[[method]]$sd, tolerance = 1e-4, style = "deparse")
    expect_snapshot_value(fit[[method]]$CI, tolerance = 1e-4, style = "deparse")
  })
}
