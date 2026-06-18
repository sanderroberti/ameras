test_that("ameras rejects input data using reserved internal column names", {
  dat <- data.frame(
    Y = c(1, 2, 3),
    D = c(1, 2, 3),
    rcdose_ameras = 0
  )

  # ameras() creates rcdose_ameras internally for RC/ERC-style mean dose use.
  # Rejecting the name in user data avoids silently overwriting a real column.
  expect_error(
    ameras(Y ~ dose(D), data = dat),
    "reserved ameras column name"
  )
})
