# These tests cover data reconstruction for amerasfit objects created with
# keep.data = FALSE. Profile likelihood CIs, residuals, and plots all rely on
# resolve_data() to validate user-supplied data and rebuild helper columns.

make_resolve_data_object <- function(stored_data = NULL, num.rows = 3) {
  ameras:::new_amerasfit(list(
    num.rows = num.rows,
    model = list(
      data = stored_data,
      family = "poisson",
      dosevars = c("D1", "D2"),
      Y = "Y",
      M_names = "M",
      X = c("groupb", "Z"),
      X_formula = ~ group + Z,
      offset = "offset",
      entry = NULL,
      exit = NULL,
      status = NULL,
      setnr = NULL
    )
  ))
}

make_resolve_data_frame <- function() {
  data.frame(
    Y = c(1, 2, 3),
    D1 = c(1, 2, 3),
    D2 = c(2, 4, 6),
    M = c(0, 1, 0),
    offset = c(1, 1, 1),
    group = factor(c("a", "b", "a")),
    Z = c(0.1, 0.2, 0.3)
  )
}

test_that("resolve_data returns stored data when available", {
  stored <- make_resolve_data_frame()
  supplied <- stored
  supplied$D1 <- supplied$D1 + 100
  object <- make_resolve_data_object(stored_data = stored)

  # Stored data is authoritative for keep.data = TRUE objects; supplied data is
  # ignored because no reconstruction is needed.
  expect_identical(ameras:::resolve_data(object, data = supplied), stored)
})

test_that("resolve_data requires data for keep.data FALSE objects", {
  object <- make_resolve_data_object(stored_data = NULL)

  expect_error(
    ameras:::resolve_data(object),
    "Data not stored on object"
  )
  expect_error(
    ameras:::resolve_data(object, data = list(D1 = 1)),
    "data must be a data frame"
  )
})

test_that("resolve_data reports missing supplied data columns by role", {
  object <- make_resolve_data_object(stored_data = NULL)
  dat <- make_resolve_data_frame()

  expect_error(
    ameras:::resolve_data(object, data = dat[setdiff(names(dat), "D2")]),
    "missing dose columns.*D2"
  )
  expect_error(
    ameras:::resolve_data(object, data = dat[setdiff(names(dat), "offset")]),
    "missing columns present during fitting.*offset"
  )
  expect_error(
    ameras:::resolve_data(object, data = dat[setdiff(names(dat), "group")]),
    "missing X variables present during fitting.*group"
  )
})

test_that("resolve_data validates supplied dose columns and row count", {
  object <- make_resolve_data_object(stored_data = NULL)
  dat <- make_resolve_data_frame()

  nonnumeric_dose <- dat
  nonnumeric_dose$D1 <- as.character(nonnumeric_dose$D1)
  expect_error(
    ameras:::resolve_data(object, data = nonnumeric_dose),
    "dosevars:D1 must be numeric"
  )

  nonfinite_dose <- dat
  nonfinite_dose$D2[1] <- Inf
  expect_error(
    ameras:::resolve_data(object, data = nonfinite_dose),
    "dosevars:D2 must contain finite values"
  )

  extra_row <- rbind(dat, dat[1, ])
  expect_warning(
    resolved <- ameras:::resolve_data(object, data = extra_row),
    "Supplied data has 4 rows"
  )
  expect_equal(nrow(resolved), 4)
})

test_that("resolve_data rebuilds model matrix columns and mean dose", {
  object <- make_resolve_data_object(stored_data = NULL)
  dat <- make_resolve_data_frame()

  resolved <- ameras:::resolve_data(object, data = dat)

  # X_formula is re-expanded on the supplied data. In this example, the raw
  # factor column creates the expanded groupb column used by the original fit.
  expect_true("groupb" %in% names(resolved))
  expect_equal(as.numeric(resolved$groupb), c(0, 1, 0))

  # rcdose_ameras is always rebuilt from the supplied dose realizations for RC,
  # ERC, residuals, plots, and profile-likelihood reconstruction.
  expect_equal(resolved$rcdose_ameras, rowMeans(dat[c("D1", "D2")]))
})
