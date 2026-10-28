set.seed(123)
data("data", package = "ameras")

test_that("binomial accessors and diagnostic plots return structured results", {
  skip_on_cran()

  fit <- ameras(
    Y.binomial ~ dose(V1:V10) + X1 * X2,
    data = data,
    family = "binomial",
    methods = c("RC", "FMA")
  )

  vcovs <- vcov(fit)
  expect_named(vcovs, c("RC", "FMA"))
  expect_identical(vcovs$RC, fit$RC$vcov)
  expect_identical(vcovs$FMA, fit$FMA$vcov)
  expect_identical(vcov(fit, methods = "FMA"), fit$FMA$vcov)

  summ <- summary(fit)
  table <- summary_table(fit)
  expect_s3_class(summ, "summary.amerasfit")
  expect_identical(table, summ$summary_table)
  expect_named(table, c("Method", "Term", "Estimate", "SE"))
  expect_setequal(unique(table$Method), c("RC", "FMA"))

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_identical(plot(fit, ask = FALSE), fit)
  expect_identical(plot(fit, methods = "RC", ask = FALSE), fit)
  expect_identical(
    plot(fit, which = "qq", ask = FALSE),
    fit
  )
  expect_identical(
    plot(fit, type = "deviance", ask = FALSE),
    fit
  )
})

test_that("plot works with supplied data when keep.data is FALSE", {
  dat <- data
  fit <- ameras(
    Y.binomial ~ dose(V1:V10),
    data = dat,
    family = "binomial",
    methods = "RC",
    keep.data = FALSE
  )

  # plot() resolves the supplied data once and reuses that resolved version
  # while computing residuals; the caller's data frame should stay untouched.
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_identical(plot(fit, data = dat, ask = FALSE), fit)
  expect_false("rcdose_ameras" %in% names(dat))
})


test_that("multinomial accessors and plots preserve category structure", {
  fit <- ameras(
    Y.multinomial ~ dose(V1:V10) + X1 * X2,
    data = data,
    family = "multinomial",
    methods = "RC"
  )

  expect_identical(vcov(fit), fit$RC$vcov)
  table <- summary_table(fit)
  expect_identical(unique(table$Method), "RC")
  expect_identical(table$Term, names(fit$RC$coefficients))

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_identical(plot(fit, ask = FALSE), fit)
  expect_identical(
    plot(fit, which = "qq", ask = FALSE),
    fit
  )
  expect_identical(
    plot(fit, type = "deviance", ask = FALSE),
    fit
  )
})


test_that("proportional-hazards accessors and plots preserve coefficient terms", {
  fit <- ameras(
    Surv(time, status) ~ dose(V1:V10) + X1 * X2,
    data = data,
    family = "prophaz",
    methods = "RC"
  )

  expect_identical(vcov(fit), fit$RC$vcov)
  table <- summary_table(fit)
  expect_identical(unique(table$Method), "RC")
  expect_identical(table$Term, names(fit$RC$coefficients))

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_identical(plot(fit, ask = FALSE), fit)
})
