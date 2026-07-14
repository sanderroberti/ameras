set.seed(123)
data("data", package = "ameras")


test_that("public profile likelihood CIs hit target p-values for stable RC fits", {
  fit <- fit_combination(
    family = "gaussian",
    Y = "Y.gaussian",
    doseRRmod = NULL,
    deg = 1,
    X = NULL,
    M = NULL,
    data = data,
    methods = "RC"
  )
  fit <- suppressMessages(confint(
    fit,
    type = "proflik",
    parm = "dose",
    maxit.profCI = 100,
    tol.profCI = 1e-8,
    print = FALSE
  ))

  ci <- fit$RC$CI["dose", , drop = FALSE]
  expect_true(all(is.finite(as.matrix(ci[, c("lower", "upper")]))))
  expect_equal(ci$pval.lower, 0.05, tolerance = 0.005)
  expect_equal(ci$pval.upper, 0.05, tolerance = 0.005)
})


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
    if (method == "RC") {
      expect_warning(
        fit1 <- confint(fit, type = c("proflik")),
        "Profile likelihood lower bound for dose could not be bracketed"
      )
    } else {
      fit1 <- confint(fit, type = c("proflik"))
    }
    expect_snapshot_value(
      coef(fit1),
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
    expect_snapshot(print(summary(fit1)$summary_table, digits = 2))

    fit2 <- confint(fit, type = c("wald.transformed"))
    expect_snapshot_value(
      coef(fit2),
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
    expect_snapshot(print(summary(fit2)$summary_table, digits = 2))
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
      coef(fit1),
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
    expect_snapshot(print(summary(fit1)$summary_table, digits = 2))
    fit2 <- confint(fit, type = c("hpd"))
    expect_snapshot_value(
      coef(fit2),
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
    expect_snapshot(print(summary(fit2)$summary_table, digits = 2))
  })
}
