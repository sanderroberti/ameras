set.seed(123)
data("data", package = "ameras")


# Clogit, basic model, all methods

for (method in all_methods) {
  test_that(paste("clogit snapshot:", method), {
    if (method %in% c("ERC", "MCML", "BMA")) {
      skip_on_cran()
    }

    fit <- fit_combination(
      family = "clogit",
      Y = "Y.clogit",
      doseRRmod = "EXP",
      deg = 2,
      X = NULL,
      M = NULL,
      data = data,
      methods = method,
      niter.BMA = 1000,
      nburnin.BMA = 500
    )
    fit <- confint(fit, type = c("wald.orig", "percentile"))
    expect_snapshot_value(
      coef(fit),
      tolerance = 1e-4,
      style = "deparse"
    )
    expect_snapshot_value(fit[[method]]$sd, tolerance = 1e-4, style = "deparse")
    expect_snapshot_value(fit[[method]]$CI, tolerance = 1e-4, style = "deparse")
    expect_snapshot(print(summary(fit)$summary_table, digits = 2))
  })
}


test_that(paste("clogit snapshot: FMA 1-par"), {
  skip_on_cran()

  fit <- fit_combination(
    family = "clogit",
    Y = "Y.clogit",
    doseRRmod = "ERR",
    deg = 1,
    X = NULL,
    M = NULL,
    data = data,
    methods = "FMA"
  )
  fit <- confint(fit)
  expect_snapshot_value(
    coef(fit),
    tolerance = 1e-4,
    style = "deparse"
  )
  expect_snapshot_value(fit[["FMA"]]$sd, tolerance = 1e-4, style = "deparse")
  expect_snapshot_value(fit[["FMA"]]$CI, tolerance = 1e-4, style = "deparse")
  expect_snapshot(print(summary(fit)$summary_table, digits = 2))
})
