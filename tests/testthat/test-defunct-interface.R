test_that("the direct argument interface is defunct", {
  dat <- ameras_simulated[seq_len(10), ]

  expect_error(
    ameras(
      data = dat,
      family = "binomial",
      Y = "Y.binomial",
      dosevars = "V1"
    ),
    "Use the formula interface",
    class = "lifecycle_error_deprecated"
  )

  # Detect old arguments before R's partial argument matching can interpret M
  # as MFMA.
  expect_error(
    ameras(
      Y.binomial ~ dose(V1),
      data = dat,
      family = "binomial",
      M = "M1"
    ),
    "Use the formula interface",
    class = "lifecycle_error_deprecated"
  )

  # Positional calls to the former data-first interface should receive the
  # same migration guidance rather than failing during data validation.
  expect_error(
    ameras(
      dat,
      family = "binomial",
      Y = "Y.binomial",
      dosevars = "V1"
    ),
    "Use the formula interface",
    class = "lifecycle_error_deprecated"
  )
})


test_that("a missing formula is not mistaken for the direct interface", {
  error <- tryCatch(
    ameras(data = ameras_simulated[seq_len(10), ]),
    error = identity
  )

  expect_s3_class(error, "error")
  expect_false(inherits(error, "lifecycle_error_deprecated"))
  expect_match(conditionMessage(error), "`formula` is required", fixed = TRUE)
})


test_that("defunct arguments are absent from the public interface", {
  defunct <- c(
    "Y", "dosevars", "doseRRmod", "deg", "M", "X", "offset", "entry",
    "exit", "setnr", "CI", "params.profCI", "maxit.profCI", "tol.profCI",
    "included.replicates.BMA"
  )

  expect_false(any(defunct %in% names(formals(ameras))))
})


test_that("confidence-interval arguments in ameras are defunct", {
  dat <- ameras_simulated[seq_len(10), ]
  base_args <- list(
    formula = Y.binomial ~ dose(V1),
    data = dat,
    family = "binomial"
  )
  defunct_args <- list(
    CI = "wald.orig",
    params.profCI = "dose",
    maxit.profCI = 20,
    tol.profCI = 1e-2
  )

  for (arg in names(defunct_args)) {
    call_args <- c(base_args, setNames(list(defunct_args[[arg]]), arg))
    expect_error(
      do.call(ameras, call_args),
      "call add_confint",
      class = "lifecycle_error_deprecated",
      info = arg
    )
  }
})


test_that("included.replicates.BMA is defunct", {
  expect_error(
    ameras(
      Y.binomial ~ dose(V1:V2),
      data = ameras_simulated[seq_len(10), ],
      family = "binomial",
      methods = "BMA",
      included.replicates.BMA = 1:2
    ),
    "included.realizations.BMA",
    class = "lifecycle_error_deprecated"
  )
})
