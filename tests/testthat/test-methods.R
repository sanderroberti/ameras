set.seed(123)
data("data", package = "ameras")

test_that("vcov & plot & summary_table binomial", {
  fit <- ameras(
    Y.binomial ~ dose(V1:V10) + X1 * X2,
    data = data,
    family = "binomial",
    methods = c("RC", "FMA")
  )

  expect_no_error({
    vcov(fit)
    vcov(fit, methods = "FMA")
    summary(fit)
    summary_table(fit)
    plot(fit, ask = FALSE)
    plot(fit, methods = "RC", ask = FALSE)
  })

  expect_no_error(
    plot(fit, which = "qq", ask = FALSE)
  )

  expect_no_error(
    plot(fit, type = "deviance", ask = FALSE)
  )
})


test_that("vcov & plot & summary_table multinomial", {
  fit <- ameras(
    Y.multinomial ~ dose(V1:V10) + X1 * X2,
    data = data,
    family = "multinomial",
    methods = "RC"
  )

  expect_no_error({
    vcov(fit)
    summary_table(fit)
    plot(fit, ask = FALSE)
  })

  expect_no_error(
    plot(fit, which = "qq", ask = FALSE)
  )

  expect_no_error(
    plot(fit, type = "deviance", ask = FALSE)
  )
})


test_that("vcov & plot & summary_table prophaz", {
  fit <- ameras(
    Surv(time, status) ~ dose(V1:V10) + X1 * X2,
    data = data,
    family = "prophaz",
    methods = "RC"
  )

  expect_no_error({
    vcov(fit)
    summary_table(fit)
    plot(fit, ask = FALSE)
  })
})
