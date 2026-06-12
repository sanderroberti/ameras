test_that("string helpers accept and format single strings", {
  expect_true(ameras:::isString("RC"))
  expect_false(ameras:::isString(c("RC", "FMA")))
  expect_false(ameras:::isString(1))

  expect_identical(ameras:::getCharVecStr(c("RC", "FMA")), "'RC','FMA'")
  expect_identical(
    ameras:::getCharVecStr(c("RC", "FMA"), sep = " or "),
    "'RC' or 'FMA'"
  )

  # check_string() trims user input before matching against the allowed values.
  expect_identical(
    ameras:::check_string(" RC ", valid = c("RC", "FMA"), parm = "methods"),
    "RC"
  )
  expect_error(
    ameras:::check_string(c("RC", "FMA"), valid = c("RC", "FMA"), parm = "methods"),
    "invalid values 'RC','FMA'"
  )
  expect_error(
    ameras:::check_string("MCML", valid = c("RC", "FMA"), parm = "methods"),
    "invalid values 'RC','FMA'"
  )
})

test_that("offset checks cover empty, valid, and invalid offsets", {
  dat <- data.frame(
    y = 1:3,
    offset = c(0, 1, 2),
    bad_offset = c(0, -1, 2)
  )

  # Empty offset specifications are allowed and return early.
  expect_null(ameras:::check_offset(NULL, dat))
  expect_null(ameras:::check_offset(character(), dat))

  # A supplied offset must name a numeric, finite, non-negative column.
  expect_null(ameras:::check_offset("offset", dat))
  expect_error(
    ameras:::check_offset("bad_offset", dat),
    "offset must contain non-negative values"
  )
  expect_error(
    ameras:::check_offset("missing_offset", dat),
    "offset contains invalid values"
  )
})

test_that("entry and exit checks compare the observed times", {
  dat <- data.frame(
    entry = c(0, 1, 2),
    exit = c(1, 2, 3),
    bad_entry = c(0, 3, 2)
  )

  # The no-entry branch only checks that exit exists and is numeric.
  expect_null(ameras:::check_entry_exit(NULL, "exit", dat))

  # With entry supplied, the helper should compare entry/exit values row-wise.
  expect_null(ameras:::check_entry_exit("entry", "exit", dat))
  expect_error(
    ameras:::check_entry_exit("bad_entry", "exit", dat),
    "entry > exit"
  )
})

test_that("initial parameter checks use family-specific parameter counts", {
  X <- c("X1", "X2")
  M <- "M1"
  deg <- 2

  # Empty initial values are allowed; ameras() then chooses defaults later.
  expect_null(
    ameras:::check_inpar(NULL, family = "gaussian", M = M, X = X, deg = deg)
  )

  expected_lengths <- c(
    gaussian = 2 + length(X) + length(M) * deg + deg,
    binomial = 1 + length(X) + length(M) * deg + deg,
    poisson = 1 + length(X) + length(M) * deg + deg,
    prophaz = length(X) + length(M) * deg + deg,
    clogit = length(X) + length(M) * deg + deg,
    multinomial = (4 - 1) * (1 + length(X) + length(M) * deg + deg)
  )

  for (family in names(expected_lengths)) {
    multinom_levels <- if (family == "multinomial") 4 else 0
    valid <- seq_len(expected_lengths[[family]])

    expect_identical(
      ameras:::check_inpar(
        valid,
        family = family,
        M = M,
        X = X,
        deg = deg,
        multinom_levels = multinom_levels
      ),
      valid
    )

    # One missing value exercises the length check for the same family branch.
    expect_error(
      ameras:::check_inpar(
        valid[-1],
        family = family,
        M = M,
        X = X,
        deg = deg,
        multinom_levels = multinom_levels
      ),
      "inpar must be a numeric vector of length"
    )
  }

  expect_error(
    ameras:::check_inpar(1, family = "unknown", M = M, X = X, deg = deg),
    "ERROR"
  )
})

test_that("general check helpers cover representative validation branches", {
  dat <- data.frame(
    x = c(1, 2, 3),
    y = c(0, 1, 0),
    z = factor(c("a", "b", "c")),
    z_unused = factor(c("a", "b", "b"), levels = c("a", "b", "c"))
  )

  expect_null(ameras:::check_df(dat))
  expect_error(ameras:::check_df(list(x = 1:3)), "must be a data frame")
  expect_error(ameras:::check_df(dat[0, ]), "has no rows")
  expect_error(ameras:::check_df(dat[, 0]), "has no columns")

  expect_identical(ameras:::check_methods(c("RC", "RC", "FMA")), c("RC", "FMA"))
  expect_error(ameras:::check_methods("BAD"), "methods contains invalid values")

  expect_identical(ameras:::check_deg(NULL), 2)
  expect_identical(ameras:::check_deg(1), 1)
  expect_error(ameras:::check_deg(3), "deg must be <= 2")

  expect_null(ameras:::check_factor_vec(dat$z, "Y"))
  expect_error(ameras:::check_factor_vec(dat$x, "Y"), "Y must be numeric")
  expect_error(
    ameras:::check_factor_vec(dat$z_unused, "Y"),
    "Y contains unused levels"
  )

  expect_identical(
    ameras:::check_vars(dat, c("x", "y"), "vars", minlen = 1),
    c("x", "y")
  )
  expect_identical(
    ameras:::check_vars(dat, 1:2, "vars", minlen = 1),
    1:2
  )
  expect_error(
    ameras:::check_vars(dat, c("x", "x"), "vars", minlen = 1),
    "vars contains duplicated values"
  )
  expect_error(
    ameras:::check_vars(dat, "missing", "vars", minlen = 1),
    "vars contains invalid values"
  )
})
