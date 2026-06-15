# Fast tests for small amerasfit accessors that do not need model fitting.

test_that("Rhat returns BMA convergence diagnostics", {
  diagnostics <- data.frame(
    Rhat = c(1.01, 1.02),
    n.eff = c(100, 80),
    row.names = c("dose", "dose_squared")
  )
  fit <- ameras:::new_amerasfit(list(
    BMA = list(Rhat = diagnostics)
  ))

  expect_identical(Rhat(fit), diagnostics)
  expect_error(
    Rhat(ameras:::new_amerasfit(list(RC = list()))),
    "BMA not present"
  )
})

test_that("included_realizations returns available FMA and BMA indices", {
  fit <- ameras:::new_amerasfit(list(
    FMA = list(included.realizations = c(1L, 3L)),
    BMA = list(included.realizations = c(2L, 4L))
  ))

  both <- included_realizations(fit)
  expect_named(both, c("FMA", "BMA"))
  expect_identical(both$FMA, c(1L, 3L))
  expect_identical(both$BMA, c(2L, 4L))

  expect_identical(
    included_realizations(fit, methods = "FMA"),
    c(1L, 3L)
  )
  expect_identical(
    included_realizations(fit, methods = "BMA"),
    c(2L, 4L)
  )
})

test_that("included_realizations errors when requested methods are absent", {
  fit <- ameras:::new_amerasfit(list(RC = list()))

  expect_error(
    included_realizations(fit, methods = "FMA"),
    "None of the requested methods were run"
  )
})

test_that("traceplot requires BMA output", {
  fit <- ameras:::new_amerasfit(list(RC = list()))

  expect_error(
    traceplot(fit),
    "traceplot\\(\\) requires BMA output"
  )
})

test_that("traceplot forwards BMA samples and plotting arguments", {
  samples <- list(
    chain1 = matrix(1:6, ncol = 2, dimnames = list(NULL, c("a", "b")))
  )
  fit <- ameras:::new_amerasfit(list(
    BMA = list(samples = samples)
  ))
  captured <- NULL

  testthat::local_mocked_bindings(
    MCMCtrace = function(object, ...) {
      captured <<- c(list(object = object), list(...))
      invisible("trace-ok")
    },
    .package = "ameras"
  )

  out <- traceplot(
    fit,
    iter = 10,
    Rhat = FALSE,
    n.eff = FALSE,
    pdf = TRUE,
    ind = TRUE
  )

  expect_identical(out, "trace-ok")
  expect_identical(captured$object, samples)
  expect_identical(captured$iter, 10)
  expect_false(captured$Rhat)
  expect_false(captured$n.eff)
  expect_true(captured$pdf)
  expect_true(captured$ind)
})
