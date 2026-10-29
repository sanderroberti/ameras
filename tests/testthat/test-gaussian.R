set.seed(123)
data <- ameras_simulated


# Gaussian, basic model, all methods

for (method in all_methods) {
  test_that(paste("gaussian snapshot:", method), {
    skip_on_cran()

    if (method == "BMA") {
      skip_on_cran()
    }

    fit <- fit_combination(
      family = "gaussian",
      Y = "Y.gaussian",
      deg = 2,
      X = NULL,
      M = NULL,
      data = data,
      methods = method,
      niter.BMA = 5000,
      nburnin.BMA = 1000
    )
    fit <- confint(fit, type = c("wald.orig", "percentile"))
    expect_snapshot_value(
      coef(fit),
      tolerance = 1e-4,
      style = "deparse"
    )
    expect_snapshot_value(fit[[method]]$sd, tolerance = 1e-4, style = "deparse")
    expect_snapshot_value(fit[[method]]$CI, tolerance = 1e-4, style = "deparse")
    if (method == "BMA") {
      expect_snapshot(print_summary_snapshot(
        fit, max_rhat = 1.2, min_n_eff = 300
      ))
    } else {
      expect_snapshot(print_summary_snapshot(fit))
    }
  })
}
