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

test_that("FMA future chunk size checks accept positive scalar values", {
  # This argument is passed through to future.apply when available. Validate it
  # before fitting so invalid values fail with an ameras-specific error.
  expect_null(ameras:::check_future_chunk_size_FMA(NULL))
  expect_identical(ameras:::check_future_chunk_size_FMA(1), 1)
  expect_identical(ameras:::check_future_chunk_size_FMA(Inf), Inf)

  expect_error(
    ameras:::check_future_chunk_size_FMA(0),
    "future.chunk.size.FMA must be a positive number"
  )
  expect_error(
    ameras:::check_future_chunk_size_FMA(NA_real_),
    "future.chunk.size.FMA must be a positive number"
  )
  expect_error(
    ameras:::check_future_chunk_size_FMA(c(1, 2)),
    "future.chunk.size.FMA must be a positive number"
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
  expect_null(ameras:::check_reserved_names(dat))
  expect_error(
    ameras:::check_reserved_names(cbind(dat, rcdose_ameras = 1:3)),
    "reserved ameras column name"
  )

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

test_that("family and dose-response checks validate scalar choices", {
  # These wrappers are thin, but they are the public-facing messages for bad
  # family/model choices parsed from ameras() calls.
  expect_identical(ameras:::check_family("gaussian"), "gaussian")
  expect_identical(ameras:::check_doseRRmod("ERR"), "ERR")

  expect_error(
    ameras:::check_family(c("gaussian", "poisson")),
    "family must have length 1"
  )
  expect_error(
    ameras:::check_family("cox"),
    "family contains invalid values"
  )
  expect_error(
    ameras:::check_doseRRmod("linear"),
    "doseRRmod contains invalid values"
  )
})

test_that("outcome checks enforce family-specific response requirements", {
  dat <- data.frame(
    y_gaussian = c(0.2, 1.3, 2.1),
    y_bad_numeric = c(0, Inf, 1),
    y_binary = c(0, 1, 0),
    y_bad_binary = c(0, 2, 1),
    y_count = c(0, 1, 3),
    y_bad_count = c(0, 1.5, 3),
    y_factor = factor(c("a", "b", "c")),
    y_bad_factor = factor(c("a", "b", "b"))
  )

  # The same Y argument is checked differently depending on family.
  expect_null(ameras:::check_Y("y_gaussian", dat, family = "gaussian"))
  expect_null(ameras:::check_Y("y_binary", dat, family = "binomial"))
  expect_null(ameras:::check_Y("y_count", dat, family = "poisson"))
  expect_null(ameras:::check_Y("y_binary", dat, family = "prophaz"))
  expect_null(ameras:::check_Y("y_binary", dat, family = "clogit"))
  expect_null(ameras:::check_Y("y_factor", dat, family = "multinomial"))

  expect_error(
    ameras:::check_Y("y_bad_numeric", dat, family = "gaussian"),
    "Y must contain finite values"
  )
  expect_error(
    ameras:::check_Y("y_bad_binary", dat, family = "binomial"),
    "Y must contain binary"
  )
  expect_error(
    ameras:::check_Y("y_bad_count", dat, family = "poisson"),
    "Y must be integer"
  )
  expect_error(
    ameras:::check_Y("y_binary", dat, family = "multinomial"),
    "Y must be numeric"
  )
  expect_error(
    ameras:::check_Y("missing_y", dat, family = "gaussian"),
    "Y contains invalid values"
  )
})

test_that("dose, modifier, covariate, and set checks validate data columns", {
  dat <- data.frame(
    D1 = c(0.1, 0.2, 0.3),
    D2 = c(0.2, 0.4, 0.6),
    D_bad = c(0.1, NA, 0.3),
    M = c(0, 1, 0),
    M_bad = c(0, 2, 1),
    X = c(1.1, 2.2, 3.3),
    X_bad = c(1, NaN, 3),
    setnr = c(1, 1, 2),
    bad_setnr = c(1, 1.5, 2)
  )

  # Non-RC methods require more than one dose realization.
  expect_null(ameras:::check_D(c("D1", "D2"), dat, methods = "FMA"))
  expect_null(ameras:::check_D("D1", dat, methods = "RC"))
  expect_error(
    ameras:::check_D("D1", dat, methods = "FMA"),
    "Multiple exposure realizations required"
  )
  expect_error(
    ameras:::check_D("D_bad", dat, methods = "RC"),
    "dosevars:D_bad must contain finite values"
  )

  # Modifiers must be binary, while ordinary X covariates only need to be
  # numeric and finite. NULL X is allowed because a model can have no covariates.
  expect_null(ameras:::check_M("M", dat))
  expect_error(ameras:::check_M("M_bad", dat), "M:M_bad must contain binary")
  expect_null(ameras:::check_X(NULL, dat))
  expect_null(ameras:::check_X("X", dat))
  expect_error(ameras:::check_X("X_bad", dat), "X:X_bad must contain finite")

  # Matched sets of size 1 are allowed but warned about because they do not
  # contribute to conditional likelihood estimation.
  expect_warning(
    ameras:::check_setnr("setnr", dat),
    "matched sets of size 1"
  )
  expect_error(
    ameras:::check_setnr("bad_setnr", dat),
    "setnr must be integer"
  )
})

test_that("primitive vector checks report length, type, and bounds errors", {
  expect_null(ameras:::check_num_vec(c(1, 2), "x", len = 2))
  expect_error(ameras:::check_num_vec("1", "x"), "x must be numeric")
  expect_error(
    ameras:::check_num_vec(c(1, 2, 3), "x", len = 2),
    "x must be a numeric vector of length 2"
  )
  expect_error(ameras:::check_num_vec(c(1, Inf), "x"), "finite values")
  expect_error(
    ameras:::check_num_vec(c(0, 2), "x", binary = 1),
    "binary"
  )
  expect_error(
    ameras:::check_num_vec(c(0, -1), "x", nonneg = 1),
    "non-negative"
  )
  expect_error(
    ameras:::check_num_vec(c(1, 1.2), "x", integer = 1),
    "x must be integer"
  )

  expect_null(ameras:::check_integer(1:2, "idx", min = 1, max = 2))
  expect_error(ameras:::check_integer(integer(), "idx"), "length >= 1")
  expect_error(
    ameras:::check_integer(c(1, 2), "idx", minlen = 1, maxlen = 1),
    "length 1"
  )
  expect_error(ameras:::check_integer("1", "idx"), "idx must be integer")
  expect_error(ameras:::check_integer(Inf, "idx"), "idx must be integer")
  expect_error(ameras:::check_integer(1.5, "idx"), "idx must be integer")
  expect_error(ameras:::check_integer(0, "idx", min = 1), "idx must be >= 1")
  expect_error(ameras:::check_integer(3, "idx", max = 2), "idx must be <= 2")

  expect_identical(
    ameras:::check_char_vec(character(), "choice", def = "RC"),
    "RC"
  )
  expect_identical(ameras:::check_char_vec("RC", "choice"), "RC")
  expect_error(
    ameras:::check_char_vec(c("RC", "FMA"), "choice", len = 1),
    "choice must have length 1"
  )
  expect_error(ameras:::check_char_vec(1, "choice"), "choice must be character")
  expect_error(
    ameras:::check_char_vec("bad", "choice", valid = c("RC", "FMA")),
    "choice contains invalid values"
  )
})

test_that("variable utilities handle empty, numeric, character, and bad inputs", {
  dat <- data.frame(x = 1:3, y = 4:6)

  expect_null(ameras:::check_vars(dat, NULL, "vars"))
  expect_error(ameras:::check_vars(dat, NULL, "vars", minlen = 1), "length >= 1")
  expect_error(
    ameras:::check_vars(dat, NULL, "vars", minlen = 1, maxlen = 1),
    "length 1"
  )
  expect_error(
    ameras:::check_vars(dat, matrix(1, nrow = 1), "vars", minlen = 1),
    "must be a vector"
  )
  expect_error(
    ameras:::check_vars(dat, 0, "vars", minlen = 1),
    "vars must be >= 1"
  )
  expect_error(
    ameras:::check_vars(dat, 3, "vars", minlen = 1),
    "vars must be <= 2"
  )

  expect_null(ameras:::getVarNumbers(NULL, dat))
  expect_identical(ameras:::getVarNumbers(2, dat), 2)
  expect_identical(ameras:::getVarNumbers("y", dat), 2L)
})

test_that("required variable helper lists family-specific source columns", {
  base_model <- list(
    dosevars = c("D1", "D2"),
    X = "X1",
    M_names = "M1",
    Y = "Y",
    offset = "offset",
    status = "status",
    exit = "exit",
    entry = "entry",
    setnr = "setnr"
  )

  # required_vars() is used when data must be supplied again after the fitted
  # object no longer stores it, so each family needs the right reconstruction
  # columns.
  expect_identical(
    ameras:::required_vars(c(base_model, family = "gaussian")),
    c("D1", "D2", "X1", "M1", "Y")
  )
  expect_identical(
    ameras:::required_vars(c(base_model, family = "poisson")),
    c("D1", "D2", "X1", "M1", "Y", "offset")
  )
  expect_identical(
    ameras:::required_vars(c(base_model, family = "multinomial")),
    c("D1", "D2", "X1", "M1", "Y")
  )
  expect_identical(
    ameras:::required_vars(c(base_model, family = "prophaz")),
    c("D1", "D2", "X1", "M1", "status", "exit", "entry")
  )
  expect_identical(
    ameras:::required_vars(c(base_model, family = "clogit")),
    c("D1", "D2", "X1", "M1", "status", "setnr")
  )

  no_optional <- base_model
  no_optional$X <- NULL
  no_optional$M_names <- NULL
  no_optional$offset <- NULL
  no_optional$entry <- NULL
  expect_identical(
    ameras:::required_vars(c(no_optional, family = "prophaz")),
    c("D1", "D2", "status", "exit")
  )
})
