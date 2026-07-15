test_that("formula modifiers build binary contrast and subgroup designs", {
  dat <- data.frame(
    Y = rnorm(4),
    D = 1:4,
    M_num = c(0, 1, 0, 1),
    M_log = c(FALSE, TRUE, FALSE, TRUE),
    M_fac = factor(c("low", "high", "low", "high"), levels = c("low", "high"))
  )

  # Formula modifiers are converted to numeric design columns before the
  # likelihood sees them. This lets factors/logicals work without changing the
  # likelihood code, which already expects a numeric modifier matrix.
  parsed_fac <- ameras:::parse_ameras_formula(
    Y ~ dose(D, modifier = ~ M_fac),
    dat,
    family = "gaussian"
  )
  inputs_fac <- ameras:::prepare_modifier_inputs(dat, parsed_fac$modifier_info)
  expect_identical(parsed_fac$modifier_info$coding, "contrast")
  expect_equal(inputs_fac$data[[inputs_fac$design_vars]], c(0, 1, 0, 1))
  expect_identical(inputs_fac$modifier_info$parameter_names, "M_fac")

  parsed_group <- ameras:::parse_ameras_formula(
    Y ~ dose(D, modifier = ~ 0 + M_log),
    dat,
    family = "gaussian"
  )
  inputs_group <- ameras:::prepare_modifier_inputs(
    dat,
    parsed_group$modifier_info
  )
  expect_identical(parsed_group$modifier_info$coding, "group")
  expect_equal(inputs_group$data[[inputs_group$design_vars]], c(0, 1, 0, 1))
  expect_identical(
    inputs_group$modifier_info$group_labels,
    c("M_log=FALSE", "M_log=TRUE")
  )

  parsed_group_alt <- ameras:::parse_ameras_formula(
    Y ~ dose(D, modifier = ~ M_num - 1),
    dat,
    family = "gaussian"
  )
  inputs_group_alt <- ameras:::prepare_modifier_inputs(
    dat,
    parsed_group_alt$modifier_info
  )
  expect_identical(parsed_group_alt$modifier_info$coding, "group")
  expect_equal(
    inputs_group_alt$data[[inputs_group_alt$design_vars]],
    dat$M_num
  )
})


test_that("formula modifier validation keeps the first release binary-only", {
  dat <- data.frame(
    Y = rnorm(5),
    D = 1:5,
    M = c(0, 1, 0, 1, 0),
    M2 = c(1, 0, 1, 0, 1),
    continuous = seq(0, 1, length.out = 5),
    multi = factor(c("a", "b", "c", "a", "b"))
  )

  expect_error(
    ameras:::parse_ameras_formula(
      Y ~ dose(D, modifier = ~ M:M2),
      dat,
      family = "gaussian"
    ),
    "Interactions and transformed terms"
  )
  expect_error(
    ameras:::parse_ameras_formula(
      Y ~ dose(D, modifier = ~ I(M)),
      dat,
      family = "gaussian"
    ),
    "Interactions and transformed terms"
  )
  expect_error(
    ameras:::parse_ameras_formula(
      Y ~ dose(D, modifier = ~ 0 + M + M2),
      dat,
      family = "gaussian"
    ),
    "exactly one binary"
  )

  parsed_cont <- ameras:::parse_ameras_formula(
    Y ~ dose(D, modifier = ~ continuous),
    dat,
    family = "gaussian"
  )
  expect_error(
    ameras:::prepare_modifier_inputs(dat, parsed_cont$modifier_info),
    "must be binary"
  )

  parsed_multi <- ameras:::parse_ameras_formula(
    Y ~ dose(D, modifier = ~ multi),
    dat,
    family = "gaussian"
  )
  expect_error(
    ameras:::prepare_modifier_inputs(dat, parsed_multi$modifier_info),
    "exactly two levels"
  )
})


test_that("legacy modifier syntax is deprecated but preserves RC output", {
  data("data", package = "ameras")

  lifecycle::expect_deprecated(
    legacy <- suppressMessages(ameras(
      Y.binomial ~ dose(V1:V2, model = "EXP", modifier = M1),
      data = data,
      family = "binomial",
      methods = "RC"
    ))
  )
  formula <- suppressMessages(ameras(
    Y.binomial ~ dose(V1:V2, model = "EXP", modifier = ~ M1),
    data = data,
    family = "binomial",
    methods = "RC"
  ))

  expect_equal(coef(legacy), coef(formula), tolerance = 1e-7)
  expect_equal(vcov(legacy), vcov(formula), tolerance = 1e-6)
  expect_equal(legacy$RC$loglik, formula$RC$loglik, tolerance = 1e-7)
})


test_that("subgroup-coded modifiers report subgroup dose coefficients", {
  data("data", package = "ameras")

  contrast <- suppressMessages(ameras(
    Y.binomial ~ dose(V1:V2, model = "EXP", modifier = ~ M1),
    data = data,
    family = "binomial",
    methods = "RC"
  ))
  group <- suppressMessages(ameras(
    Y.binomial ~ dose(V1:V2, model = "EXP", modifier = ~ 0 + M1),
    data = data,
    family = "binomial",
    methods = "RC"
  ))
  group_alt <- suppressMessages(ameras(
    Y.binomial ~ dose(V1:V2, model = "EXP", modifier = ~ M1 - 1),
    data = data,
    family = "binomial",
    methods = "RC"
  ))

  expect_identical(
    rownames(coef(group)),
    c("(Intercept)", "dose[M1=0]", "dose[M1=1]")
  )
  expect_identical(rownames(coef(group)), rownames(coef(group_alt)))
  expect_equal(coef(group), coef(group_alt), tolerance = 1e-7)

  contrast_coef <- setNames(coef(contrast)$RC, rownames(coef(contrast)))
  expected_group <- c(
    "(Intercept)" = unname(contrast_coef["(Intercept)"]),
    "dose[M1=0]" = unname(contrast_coef["dose"]),
    "dose[M1=1]" = unname(contrast_coef["dose"] + contrast_coef["dose:M1"])
  )
  expect_equal(
    setNames(coef(group)$RC, rownames(coef(group))),
    expected_group,
    tolerance = 1e-5
  )
  expect_equal(group$RC$loglik, contrast$RC$loglik, tolerance = 1e-6)
})


test_that("subgroup-coded multinomial modifiers keep response-level prefixes", {
  data("data", package = "ameras")

  fit <- suppressMessages(ameras(
    Y.multinomial ~ dose(V1:V2, model = "EXP", modifier = ~ 0 + M1),
    data = data,
    family = "multinomial",
    methods = "RC"
  ))

  expect_identical(
    rownames(coef(fit)),
    c(
      "(1)_(Intercept)", "(1)_dose[M1=0]", "(1)_dose[M1=1]",
      "(2)_(Intercept)", "(2)_dose[M1=0]", "(2)_dose[M1=1]"
    )
  )
})


test_that("subgroup-coded modifiers work for profile CIs with keep.data FALSE", {
  data("data", package = "ameras")

  fit <- suppressMessages(ameras(
    Y.binomial ~ dose(V1:V2, model = "EXP", modifier = ~ 0 + M1),
    data = data,
    family = "binomial",
    methods = "RC",
    keep.data = FALSE
  ))

  fit <- suppressMessages(confint(
    fit,
    type = "proflik",
    parm = "dose",
    data = data,
    maxit.profCI = 10,
    print = FALSE
  ))

  expect_identical(
    rownames(fit$RC$CI),
    c("dose[M1=0]", "dose[M1=1]")
  )
})


test_that("subgroup-coded modifiers smoke-test ERC, MCML, FMA, and BMA guard", {
  data("data", package = "ameras")

  fit_freq <- suppressMessages(ameras(
    Y.gaussian ~ dose(V1:V2, modifier = ~ 0 + M1),
    data = data,
    family = "gaussian",
    methods = c("RC", "ERC", "MCML")
  ))
  expect_true(all(c("RC", "ERC", "MCML") %in% names(fit_freq)))
  expect_true(all(grepl("^dose\\[M1=", rownames(coef(fit_freq))[2:3])))

  fit_fma <- suppressWarnings(suppressMessages(ameras(
    Y.gaussian ~ dose(V1:V2, modifier = ~ 0 + M1),
    data = data,
    family = "gaussian",
    methods = "FMA",
    MFMA = 1000
  )))
  expect_true("FMA" %in% names(fit_fma))
  expect_true(all(grepl("^dose\\[M1=", names(fit_fma$FMA$coefficients)[2:3])))

  expect_error(
    suppressMessages(ameras(
      Y.gaussian ~ dose(V1:V2, modifier = ~ 0 + M1),
      data = data,
      family = "gaussian",
      methods = "BMA",
      niter.BMA = 20,
      nburnin.BMA = 5,
      nchains.BMA = 1
    )),
    "BMA does not yet support subgroup-coded modifiers"
  )
})


test_that("formula contrast modifiers remain supported for BMA", {
  skip_on_cran()
  data("data", package = "ameras")

  # BMA still uses the existing reference-plus-contrast parameterization. This
  # keeps the new formula syntax compatible with the current BMA model code.
  fit <- suppressWarnings(suppressMessages(ameras(
    Y.gaussian ~ dose(V1:V2, modifier = ~ M1),
    data = data[seq_len(20), ],
    family = "gaussian",
    methods = "BMA",
    niter.BMA = 20,
    nburnin.BMA = 5,
    nchains.BMA = 1,
    print = FALSE
  )))

  expect_named(
    fit$BMA$coefficients,
    c("(Intercept)", "dose", "dose:M1", "sigma")
  )
})


test_that("dose_lrt treats subgroup-coded dose coefficients as dose terms", {
  data("data", package = "ameras")

  fit <- suppressMessages(ameras(
    Y.gaussian ~ dose(V1:V2, modifier = ~ 0 + M1),
    data = data,
    family = "gaussian",
    methods = "RC"
  ))

  lrt <- dose_lrt(fit, methods = "RC", type = "individual")
  expect_identical(lrt$term, c("dose[M1=0]", "dose[M1=1]"))
})
