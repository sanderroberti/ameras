test_that("dldd_clogit matches direct risk-set calculations", {
  designmat <- matrix(c(1, 0, 1, 1), nrow = 2)
  RRs <- c(2, 3)

  expected <- c(
    1 / (RRs[1] + RRs[2]),
    1 / (RRs[1] + RRs[2]) + 1 / RRs[2]
  )

  expect_equal(ameras:::dldd_clogit(designmat, RRs), expected)
})

test_that("dldd_prophaz matches direct risk-set calculations", {
  entry <- c(0, 0, 1)
  exit <- c(1, 2, 3)
  status <- c(1, 1, 0)
  RRs <- c(2, 3, 5)

  expected <- c(
    1 / (RRs[1] + RRs[2]),
    1 / (RRs[1] + RRs[2]) + 1 / (RRs[2] + RRs[3]),
    1 / (RRs[2] + RRs[3])
  )

  expect_equal(ameras:::dldd_prophaz(entry, exit, status, RRs), expected)
})
