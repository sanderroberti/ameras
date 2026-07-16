set.seed(123)
data("data", package = "ameras")


# Prophaz, basic model, all methods

test_that("prophaz filters zero-follow-up rows before fitting", {
  base <- data.frame(
    entry = rep(0, 8),
    exit = seq_len(8),
    status = rep(c(1, 0), 4),
    D = seq(0.1, 0.8, length.out = 8)
  )
  dat <- rbind(
    base,
    data.frame(entry = c(2, 5), exit = c(2, 5), status = c(1, 0), D = c(0.9, 1))
  )

  expect_warning(
    fit <- suppressMessages(ameras(
      Surv(entry, exit, status) ~ dose(D, model = EXP),
      data = dat,
      family = "prophaz",
      methods = "RC",
      print = FALSE
    )),
    "entry == exit"
  )

  expect_equal(fit$num.rows, nrow(base))
  expect_false(any(fit$model$data$entry == fit$model$data$exit))
  expect_true(all(is.finite(fit$RC$coefficients)))

  summ <- summary(fit)
  expect_equal(summ$row_info$supplied, nrow(dat))
  expect_equal(summ$row_info$omitted.na, 0L)
  expect_equal(summ$row_info$prophaz.zero.followup.rows, nrow(dat) - nrow(base))
  expect_equal(summ$row_info$used, nrow(base))

  out <- capture.output(print(summ))
  expect_true(any(grepl(
    paste0(
      "Excluded as zero-follow-up proportional hazards rows: ",
      nrow(dat) - nrow(base)
    ),
    out,
    fixed = TRUE
  )))
  expect_true(any(grepl("entry == exit", out, fixed = TRUE)))
  expect_false(any(grepl("Omitted by na.action: 0", out, fixed = TRUE)))
  expect_true(any(grepl(
    paste0("Used for fitting: ", nrow(base)),
    out,
    fixed = TRUE
  )))

  fit_nodata <- NULL
  expect_warning(
    fit_nodata <- suppressMessages(ameras(
      Surv(entry, exit, status) ~ dose(D, model = EXP),
      data = dat,
      family = "prophaz",
      methods = "RC",
      keep.data = FALSE,
      print = FALSE
    )),
    "entry == exit"
  )
  resolved <- NULL
  expect_warning(
    resolved <- ameras:::resolve_data(fit_nodata, data = dat),
    "entry == exit"
  )
  expect_equal(nrow(resolved), nrow(base))
  expect_false(any(resolved$entry == resolved$exit))
})

for (method in all_methods) {
  test_that(paste("prophaz snapshot:", method), {
    skip_on_cran()

    fit <- fit_combination(
      family = "prophaz",
      Y = "status",
      doseRRmod = "EXP",
      deg = 2,
      X = NULL,
      M = NULL,
      data = data,
      methods = method,
      niter.BMA = 5000,
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


test_that(paste("prophaz snapshot: FMA 1-par"), {
  skip_on_cran()

  fit <- fit_combination(
    family = "prophaz",
    Y = "status",
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
