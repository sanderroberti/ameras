make_parname_test_data <- function() {
  # Keep the fixture tiny: the helper only reads column names and multinomial
  # factor levels, so no realistic outcome values are needed.
  data.frame(
    Y.multinomial = factor(
      c("low", "mid", "high"),
      levels = c("low", "mid", "high")
    ),
    X1 = 1:3,
    X2 = 4:6,
    M1 = c(0, 1, 0),
    M2 = c(1, 0, 1)
  )
}

test_that("method parnames preserve Gaussian sigma and modifier names", {
  parname_data <- make_parname_test_data()

  # Gaussian methods append sigma after the mean-model and dose-response terms.
  expect_identical(
    ameras:::make_method_parnames(
      family = "gaussian",
      data = parname_data,
      X = c(2, 3),
      M = c(4, 5),
      deg = 2
    ),
    c(
      "(Intercept)",
      "X1",
      "X2",
      "dose",
      "dose_squared",
      "dose:M1",
      "dose:M2",
      "dose_squared:M1",
      "dose_squared:M2",
      "sigma"
    )
  )
})

test_that("method parnames distinguish LINEXP and no-intercept families", {
  parname_data <- make_parname_test_data()

  # GLM-like LINEXP models retain the intercept and use the two LINEXP-specific
  # dose labels, including separate modifier names for each dose component.
  expect_identical(
    ameras:::make_method_parnames(
      family = "binomial",
      data = parname_data,
      M = 4,
      deg = 2,
      doseRRmod = "LINEXP"
    ),
    c(
      "(Intercept)",
      "dose_linear",
      "dose_exponential",
      "dose_linear:M1",
      "dose_exponential:M1"
    )
  )

  # Conditional logistic and proportional hazards methods intentionally omit
  # the intercept-like parameter from otherwise identical dose-response names.
  expect_identical(
    ameras:::make_method_parnames(
      family = "clogit",
      data = parname_data,
      M = 4,
      deg = 2,
      doseRRmod = "EXP"
    ),
    c("dose", "dose_squared", "dose:M1", "dose_squared:M1")
  )
})

test_that("method parnames include multinomial prefixes and BMA hazard names", {
  parname_data <- make_parname_test_data()

  # Multinomial methods prefix each base parameter with the modeled response
  # level. The last factor level is the reference and is intentionally absent.
  expect_identical(
    ameras:::make_method_parnames(
      family = "multinomial",
      data = parname_data,
      Y = 1,
      M = 4,
      deg = 2,
      doseRRmod = "LINEXP"
    ),
    c(
      "(low)_(Intercept)",
      "(low)_dose_linear",
      "(low)_dose_exponential",
      "(low)_dose_linear:M1",
      "(low)_dose_exponential:M1",
      "(mid)_(Intercept)",
      "(mid)_dose_linear",
      "(mid)_dose_exponential",
      "(mid)_dose_linear:M1",
      "(mid)_dose_exponential:M1"
    )
  )

  # BMA proportional hazards summaries append sampled baseline hazards. The
  # frequentist proportional hazards methods leave prophaz_numints as NULL.
  expect_identical(
    ameras:::make_method_parnames(
      family = "prophaz",
      data = parname_data,
      deg = 1,
      doseRRmod = "EXP",
      prophaz_numints = 3
    ),
    c("dose", "h0[1]", "h0[2]", "h0[3]")
  )
})
