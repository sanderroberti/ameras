make_missing_x_prophaz_data <- function() {
  data.frame(
    entry = c(0, 0, 1, 2),
    exit = c(1, 2, 3, 4),
    status = c(1, 0, 1, 0),
    D = c(0.1, 0.2, 0.3, 0.4),
    X = c(1, NA, 2, 3)
  )
}


test_that("missing ordinary X values are reported by X validation", {
  dat <- make_missing_x_prophaz_data()

  expect_error(
    ameras(
      Surv(entry, exit, status) ~ dose(D, model = "ERR") + X,
      data = dat,
      family = "prophaz",
      methods = "RC"
    ),
    "X:X must contain finite values"
  )
})


test_that("X terms that drop missing rows fail before data assignment", {
  dat <- make_missing_x_prophaz_data()

  expect_error(
    ameras(
      Surv(entry, exit, status) ~ dose(D, model = "ERR") +
        stats::na.omit(X),
      data = dat,
      family = "prophaz",
      methods = "RC"
    ),
    "X formula terms must preserve one value per row"
  )
})
