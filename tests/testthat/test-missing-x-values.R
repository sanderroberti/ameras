make_missing_x_prophaz_data <- function() {
  data.frame(
    entry = rep(0, 12),
    exit = seq_len(12),
    status = rep(c(1, 0), 6),
    D = seq(0.1, 1.2, length.out = 12),
    X = c(0, 0.1, NA, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1)
  )
}


make_missing_x_gaussian_data <- function() {
  data.frame(
    Y = c(1.2, 1.8, 2.4, 2.7, 3.1, 3.8, 4.0, 4.4, 5.1, 5.3, 5.8, 6.2),
    D = seq(0.1, 1.2, length.out = 12),
    X = c(0, 1, NA, 0, 1, 0, 1, 0, 1, 0, 1, 0)
  )
}


test_that("na.omit drops missing ordinary X values before fitting", {
  dat <- make_missing_x_prophaz_data()

  fit <- suppressWarnings(suppressMessages(
    ameras(
      Surv(entry, exit, status) ~ dose(D, model = "ERR") + X,
      data = dat,
      family = "prophaz",
      methods = "RC",
      na.action = na.omit
    )
  ))

  expect_equal(fit$num.rows.original, nrow(dat))
  expect_equal(fit$num.rows, nrow(dat) - 1)
  expect_s3_class(fit$na.action, "omit")
  expect_equal(as.integer(fit$na.action), 3L)
  expect_false(anyNA(fit$model$data$X))

  summ <- summary(fit)
  expect_equal(summ$row_info$supplied, nrow(dat))
  expect_equal(summ$row_info$omitted.na, 1L)
  expect_null(summ$row_info$clogit.uninformative.rows)
  expect_equal(summ$row_info$prophaz.zero.followup.rows, 0L)
  expect_equal(summ$row_info$used, nrow(dat) - 1)

  out <- capture.output(print(summ))
  expect_true(any(grepl("Rows:", out, fixed = TRUE)))
  expect_true(any(grepl(paste0("Supplied: ", nrow(dat)), out, fixed = TRUE)))
  expect_true(any(grepl("Omitted by na.action: 1", out, fixed = TRUE)))
  expect_false(any(grepl("uninformative matched-set", out, fixed = TRUE)))
  expect_false(any(grepl("zero-follow-up", out, fixed = TRUE)))
  expect_true(any(grepl(
    paste0("Used for fitting: ", nrow(dat) - 1),
    out,
    fixed = TRUE
  )))
})


test_that("resolve_data reapplies fitted na.action for keep.data FALSE", {
  dat <- make_missing_x_prophaz_data()

  fit <- suppressWarnings(suppressMessages(
    ameras(
      Surv(entry, exit, status) ~ dose(D, model = "ERR") + X,
      data = dat,
      family = "prophaz",
      methods = "RC",
      keep.data = FALSE,
      na.action = na.omit
    )
  ))

  resolved <- ameras:::resolve_data(fit, data = dat)

  expect_equal(nrow(resolved), nrow(dat) - 1)
  expect_false(anyNA(resolved$X))
})


test_that("na.exclude does not pad residuals to original row count", {
  dat <- make_missing_x_gaussian_data()

  fit <- suppressWarnings(suppressMessages(
    ameras(
      Y ~ dose(D) + X,
      data = dat,
      family = "gaussian",
      methods = "RC",
      na.action = na.exclude
    )
  ))

  res <- residuals(fit, type = "response")

  expect_s3_class(fit$na.action, "exclude")
  expect_equal(fit$num.rows.original, nrow(dat))
  expect_equal(fit$num.rows, nrow(dat) - 1)
  expect_length(res, fit$num.rows)
})


test_that("na.pass leaves missing X values to existing validation", {
  dat <- make_missing_x_prophaz_data()

  expect_error(
    ameras(
      Surv(entry, exit, status) ~ dose(D, model = "ERR") + X,
      data = dat,
      family = "prophaz",
      methods = "RC",
      na.action = na.pass
    ),
    "X:X must contain finite values"
  )
})


test_that("na.fail rejects missing model inputs before fitting", {
  dat <- make_missing_x_prophaz_data()

  expect_error(
    ameras(
      Surv(entry, exit, status) ~ dose(D, model = "ERR") + X,
      data = dat,
      family = "prophaz",
      methods = "RC",
      na.action = na.fail
    ),
    "missing values"
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
