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

test_that("clogit filters uninformative matched sets before fitting", {
  base <- do.call(rbind, lapply(seq_len(8), function(s) {
    data.frame(
      Y = c(1, 0, 0),
      setnr = s,
      D1 = (3 * s + 0:2) / 10,
      D2 = (3 * s + 1:3) / 10
    )
  }))
  dat <- rbind(
    base,
    data.frame(Y = c(0, 0), setnr = 9, D1 = c(3, 3.1), D2 = c(3.2, 3.3)),
    data.frame(Y = 1, setnr = 10, D1 = 3.4, D2 = 3.5)
  )

  # The last two sets are uninformative: one has no case and one has size 1.
  # They should be removed once, before any clogit method-specific fit runs.
  expect_warning(
    fit <- suppressMessages(ameras(
      Y ~ dose(D1:D2, model = EXP, deg = 1) + strata(setnr),
      data = dat,
      family = "clogit",
      methods = "RC",
      print = FALSE
    )),
    "Excluding 1 matched set"
  )

  expect_equal(fit$num.rows, nrow(base))
  expect_false(any(fit$model$data$setnr %in% c(9, 10)))
  expect_true(all(is.finite(fit$RC$coefficients)))
})

test_that("clogit rejects unsupported matched set case patterns", {
  multi_case <- data.frame(
    Y = c(1, 1, 0),
    setnr = c(1, 1, 1),
    D1 = c(0.1, 0.2, 0.3),
    D2 = c(0.2, 0.3, 0.4)
  )
  expect_error(
    ameras(
      Y ~ dose(D1:D2, model = EXP, deg = 1) + strata(setnr),
      data = multi_case,
      family = "clogit"
    ),
    "more than one case"
  )

  no_informative <- data.frame(
    Y = c(0, 0, 1),
    setnr = c(1, 1, 2),
    D1 = c(0.1, 0.2, 0.3),
    D2 = c(0.2, 0.3, 0.4)
  )
  expect_error(
    suppressWarnings(ameras(
      Y ~ dose(D1:D2, model = EXP, deg = 1) + strata(setnr),
      data = no_informative,
      family = "clogit"
    )),
    "No informative matched sets remain"
  )
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
