set.seed(123)
data("data", package = "ameras")


# Clogit, basic model, all methods

test_that("clogit supports non-default strata column names", {
  dat <- data
  dat$Set <- dat$setnr

  # The formula parser stores the variable named inside strata(). This guards
  # against later clogit setup accidentally assuming the column is named setnr.
  fit_setnr <- ameras(
    Y.clogit ~ dose(V1:V10, model = EXP, deg = 1) + strata(setnr),
    data = dat,
    family = "clogit"
  )
  fit_Set <- ameras(
    Y.clogit ~ dose(V1:V10, model = EXP, deg = 1) + strata(Set),
    data = dat,
    family = "clogit"
  )

  expect_equal(coef(fit_Set), coef(fit_setnr), tolerance = 1e-8)
  expect_equal(fit_Set$model$setnr, "Set")
})

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
