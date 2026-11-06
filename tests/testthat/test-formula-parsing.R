formula_data <- data.frame(
  Y = c(0, 1, 0, 1),
  D1 = 1:4,
  D2 = 2:5,
  X = c(0.2, 0.4, 0.1, 0.3),
  off = rep(1, 4),
  set = c(1, 1, 2, 2),
  time = 1:4,
  status = c(1, 0, 1, 0)
)


parse_test_formula <- function(formula, family = "gaussian") {
  ameras:::parse_ameras_formula(
    formula,
    formula_data,
    family,
    env = parent.frame()
  )
}


test_that("valid dose selections and arguments remain supported", {
  dose_names <- c("D1", "D2")

  parsed <- parse_test_formula(
    Y ~ dose(tidyselect::all_of(dose_names), model = EXP, deg = 2) + X
  )
  expect_identical(parsed$dosevars, dose_names)
  expect_identical(parsed$doseRRmod, "EXP")
  expect_identical(parsed$deg, 2L)

  parsed_multiple <- parse_test_formula(Y ~ dose(D1, D2))
  expect_identical(parsed_multiple$dosevars, dose_names)
})


test_that("dose argument names and values are validated", {
  for (unknown in c("degg", "anything")) {
    formula <- as.formula(paste0("Y ~ dose(D1, ", unknown, " = 2)"))
    expect_error(
      parse_test_formula(formula),
      paste0("Unknown named argument\\(s\\) in dose\\(\\): ", unknown),
      info = unknown
    )
  }

  expect_error(
    parse_test_formula(Y ~ dose(D1, mod = "EXP")),
    "Unknown named argument\\(s\\) in dose\\(\\): mod"
  )
  expect_error(
    parse_test_formula(Y ~ dose(D1, deg = 1, deg = 2)),
    "specified more than once in dose\\(\\): deg"
  )
  expect_error(
    parse_test_formula(Y ~ dose(model = EXP)),
    "at least one unnamed dose-selection expression"
  )
  expect_error(
    parse_test_formula(Y ~ dose(tidyselect::starts_with("absent"))),
    "selection must select at least one data column"
  )
  expect_error(
    parse_test_formula(Y ~ dose(D1, model = "typo")),
    "model argument in dose\\(\\) must be one of"
  )
  expect_error(
    parse_test_formula(Y ~ dose(D1, deg = 3)),
    "deg argument in dose\\(\\) must be the number 1 or 2"
  )
})


test_that("dose must occur once as a standalone additive term", {
  expect_error(
    parse_test_formula(Y ~ X),
    "exactly one standalone dose\\(\\) term"
  )
  expect_error(
    parse_test_formula(Y ~ dose(D1) + dose(D2)),
    "exactly one standalone dose\\(\\) term"
  )
  expect_error(
    parse_test_formula(Y ~ I(dose(D1))),
    "standalone additive term"
  )
  expect_error(
    parse_test_formula(Y ~ dose(D1) * X),
    "cannot be used in an interaction or transformation"
  )
  expect_error(
    parse_test_formula(Y ~ dose(D1)^2),
    "cannot be used in an interaction or transformation"
  )
})


test_that("response and family-specific special terms are validated", {
  expect_error(
    parse_test_formula(~ dose(D1)),
    "Model formula must be two-sided"
  )
  expect_error(
    parse_test_formula(log(Y) ~ dose(D1)),
    "response must be a single, untransformed data-column name"
  )
  expect_error(
    parse_test_formula(Y ~ dose(D1) + strata(set)),
    "strata\\(\\) is only supported for family='clogit'"
  )
  expect_error(
    parse_test_formula(Y ~ dose(D1) + offset(off)),
    "offset\\(\\) is only supported for family='poisson'"
  )
  expect_error(
    parse_test_formula(Y ~ dose(D1), family = "clogit"),
    "must contain exactly one strata\\(\\) term"
  )
  expect_error(
    parse_test_formula(
      Y ~ dose(D1) + strata(set) + strata(X),
      family = "clogit"
    ),
    "must contain exactly one strata\\(\\) term"
  )
  expect_error(
    parse_test_formula(
      Y ~ dose(D1) + strata(interaction(set, X)),
      family = "clogit"
    ),
    "exactly one untransformed data-column name"
  )
  expect_error(
    parse_test_formula(
      Y ~ dose(D1) + offset(off) + offset(X),
      family = "poisson"
    ),
    "can contain at most one offset\\(\\) term"
  )
  expect_error(
    parse_test_formula(
      Y ~ dose(D1) + offset(log(off)),
      family = "poisson"
    ),
    "exactly one untransformed data-column name"
  )
  expect_error(
    parse_test_formula(
      Surv(log(time), status) ~ dose(D1),
      family = "prophaz"
    ),
    "Every Surv\\(\\) argument must be"
  )

  expect_identical(
    parse_test_formula(
      Y ~ dose(D1) + strata(set),
      family = "clogit"
    )$setnr,
    "set"
  )
  expect_identical(
    parse_test_formula(
      Y ~ dose(D1) + offset(off),
      family = "poisson"
    )$offset,
    "off"
  )
  survival <- parse_test_formula(
    Surv(time, status) ~ dose(D1),
    family = "prophaz"
  )
  expect_identical(survival$exit, "time")
  expect_identical(survival$status, "status")
})
