set.seed(123)
data("data", package = "ameras")


for (method in c("RC", "ERC", "MCML")) {
  test_that(paste("proflik/wald.transformed snapshot:", method), {
    if (method %in% c("ERC", "MCML")) {
      skip_on_cran()
    }

    fit <- fit_combination(
      family = "binomial",
      Y = "Y.binomial",
      doseRRmod = "ERR",
      deg = 2,
      X = NULL,
      M = NULL,
      data = data,
      methods = method
    )
    fit1 <- confint(fit, type = c("proflik"))
    expect_snapshot_value(
      fit1[[method]]$coefficients,
      tolerance = 1e-4,
      style = "deparse"
    )
    expect_snapshot_value(
      fit1[[method]]$sd,
      tolerance = 1e-4,
      style = "deparse"
    )
    expect_snapshot_value(
      fit1[[method]]$CI,
      tolerance = 1e-4,
      style = "deparse"
    )

    fit2 <- confint(fit, type = c("wald.transformed"))
    expect_snapshot_value(
      fit2[[method]]$coefficients,
      tolerance = 1e-4,
      style = "deparse"
    )
    expect_snapshot_value(
      fit2[[method]]$sd,
      tolerance = 1e-4,
      style = "deparse"
    )
    expect_snapshot_value(
      fit2[[method]]$CI,
      tolerance = 1e-4,
      style = "deparse"
    )
  })
}


for (method in c("FMA", "BMA")) {
  test_that(paste("percentile/hpd snapshot:", method), {
    skip_on_cran()

    fit <- fit_combination(
      family = "binomial",
      Y = "Y.binomial",
      doseRRmod = "ERR",
      deg = 2,
      X = NULL,
      M = NULL,
      data = data,
      methods = method
    )
    fit1 <- confint(fit, type = c("percentile"))
    expect_snapshot_value(
      fit1[[method]]$coefficients,
      tolerance = 1e-4,
      style = "deparse"
    )
    expect_snapshot_value(
      fit1[[method]]$sd,
      tolerance = 1e-4,
      style = "deparse"
    )
    expect_snapshot_value(
      fit1[[method]]$CI,
      tolerance = 1e-4,
      style = "deparse"
    )

    fit2 <- confint(fit, type = c("hpd"))
    expect_snapshot_value(
      fit2[[method]]$coefficients,
      tolerance = 1e-4,
      style = "deparse"
    )
    expect_snapshot_value(
      fit2[[method]]$sd,
      tolerance = 1e-4,
      style = "deparse"
    )
    expect_snapshot_value(
      fit2[[method]]$CI,
      tolerance = 1e-4,
      style = "deparse"
    )
  })
}
