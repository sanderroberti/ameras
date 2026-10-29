test_that("simulated example data have the documented structure", {
  expect_s3_class(ameras_simulated, "data.frame")
  expect_equal(dim(ameras_simulated), c(3000L, 22L))
  expect_named(
    ameras_simulated,
    c(
      "Y.gaussian", "Y.binomial", "Y.poisson", "time", "status", "setnr",
      "Y.clogit", "Y.multinomial", "X1", "X2", "M1", "M2",
      paste0("V", seq_len(10))
    )
  )
})
