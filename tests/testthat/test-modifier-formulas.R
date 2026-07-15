three_level_modifier_data <- function() {
  data("data", package = "ameras")
  data$F3 <- factor(
    rep(c("low", "mid", "high"), length.out = nrow(data)),
    levels = c("low", "mid", "high")
  )
  data
}


test_that("formula modifiers build contrast and subgroup designs", {
  dat <- data.frame(
    Y = rnorm(4),
    D = 1:4,
    M_num = c(0, 1, 0, 1),
    M_log = c(FALSE, TRUE, FALSE, TRUE),
    M_fac = factor(c("low", "high", "low", "high"), levels = c("low", "high")),
    F3 = factor(c("a", "b", "c", "a"), levels = c("a", "b", "c"))
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
  expect_equal(
    unname(as.vector(inputs_fac$data[[inputs_fac$design_vars]])),
    c(0, 1, 0, 1)
  )
  expect_identical(inputs_fac$modifier_info$parameter_names, "M_fac=high")

  parsed_f3 <- ameras:::parse_ameras_formula(
    Y ~ dose(D, modifier = ~ F3),
    dat,
    family = "gaussian"
  )
  inputs_f3 <- ameras:::prepare_modifier_inputs(dat, parsed_f3$modifier_info)
  expect_identical(inputs_f3$modifier_info$parameter_names, c("F3=b", "F3=c"))
  expect_equal(
    unname(as.matrix(inputs_f3$data[, inputs_f3$design_vars])),
    cbind(c(0, 1, 0, 0), c(0, 0, 1, 0))
  )

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
  expect_equal(
    unname(as.vector(inputs_group$data[[inputs_group$design_vars]])),
    c(0, 1, 0, 1)
  )
  expect_identical(
    inputs_group$modifier_info$group_labels,
    c("M_log=FALSE", "M_log=TRUE")
  )

  parsed_group_alt <- ameras:::parse_ameras_formula(
    Y ~ dose(D, modifier = ~ F3 - 1),
    dat,
    family = "gaussian"
  )
  inputs_group_alt <- ameras:::prepare_modifier_inputs(
    dat,
    parsed_group_alt$modifier_info
  )
  expect_identical(parsed_group_alt$modifier_info$coding, "group")
  expect_identical(
    inputs_group_alt$modifier_info$group_labels,
    c("F3=a", "F3=b", "F3=c")
  )
  expect_identical(
    inputs_group_alt$modifier_info$parameter_names,
    c("F3=b", "F3=c")
  )
})


test_that("formula modifier validation rejects unsupported terms", {
  dat <- data.frame(
    Y = rnorm(5),
    D = 1:5,
    M = c(0, 1, 0, 1, 0),
    M2 = c(1, 0, 1, 0, 1),
    continuous = seq(0, 1, length.out = 5),
    character = c("a", "b", "c", "a", "b"),
    one_level = factor(rep("a", 5))
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
    "exactly one"
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

  parsed_char <- ameras:::parse_ameras_formula(
    Y ~ dose(D, modifier = ~ character),
    dat,
    family = "gaussian"
  )
  expect_error(
    ameras:::prepare_modifier_inputs(dat, parsed_char$modifier_info),
    "must be binary"
  )

  parsed_one_level <- ameras:::parse_ameras_formula(
    Y ~ dose(D, modifier = ~ one_level),
    dat,
    family = "gaussian"
  )
  expect_error(
    ameras:::prepare_modifier_inputs(dat, parsed_one_level$modifier_info),
    "at least two levels"
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


test_that("subgroup-coded multi-level factors match contrast-coded fits", {
  data <- three_level_modifier_data()

  contrast <- suppressMessages(ameras(
    Y.gaussian ~ dose(V1:V2, modifier = ~ F3),
    data = data,
    family = "gaussian",
    methods = "RC"
  ))
  group <- suppressMessages(ameras(
    Y.gaussian ~ dose(V1:V2, modifier = ~ 0 + F3),
    data = data,
    family = "gaussian",
    methods = "RC"
  ))
  group_alt <- suppressMessages(ameras(
    Y.gaussian ~ dose(V1:V2, modifier = ~ F3 - 1),
    data = data,
    family = "gaussian",
    methods = "RC"
  ))

  expect_identical(
    rownames(coef(group)),
    c("(Intercept)", "dose[F3=low]", "dose[F3=mid]", "dose[F3=high]", "sigma")
  )
  expect_identical(rownames(coef(group)), rownames(coef(group_alt)))
  expect_equal(coef(group), coef(group_alt), tolerance = 1e-7)

  contrast_coef <- setNames(coef(contrast)$RC, rownames(coef(contrast)))
  expected_group <- c(
    "(Intercept)" = unname(contrast_coef["(Intercept)"]),
    "dose[F3=low]" = unname(contrast_coef["dose"]),
    "dose[F3=mid]" = unname(contrast_coef["dose"] + contrast_coef["dose:F3=mid"]),
    "dose[F3=high]" = unname(contrast_coef["dose"] + contrast_coef["dose:F3=high"]),
    "sigma" = unname(contrast_coef["sigma"])
  )
  expect_equal(
    setNames(coef(group)$RC, rownames(coef(group))),
    expected_group,
    tolerance = 1e-5
  )
  expect_equal(group$RC$loglik, contrast$RC$loglik, tolerance = 1e-6)
})


test_that("BMA subgroup sample conversion maps contrasts to subgroup effects", {
  modifier_info <- ameras:::new_modifier_info(
    coding = "group",
    design_vars = c(".ameras_modifier_1", ".ameras_modifier_2"),
    parameter_names = c("F3=mid", "F3=high"),
    group_labels = c("F3=low", "F3=mid", "F3=high")
  )

  # Internal BMA samples are ordered as base dose components followed by
  # contrast components. Reported subgroup samples should be group-major:
  # reference subgroup first, then reference + each contrast.
  samples <- matrix(
    c(
      0, 1, 2, 0.5, 1.5, 0.25, 0.75, 3,
      10, 2, 4, -1, 3, -2, 4, 5
    ),
    nrow = 2,
    byrow = TRUE
  )

  out <- ameras:::modifier_internal_to_reported_sample_matrix(
    samples = samples,
    family = "gaussian",
    M = 1:2,
    deg = 2,
    modifier_info = modifier_info
  )

  expect_equal(
    out[, 2:7],
    matrix(
      c(
        1, 2, 1.5, 2.25, 2.5, 2.75,
        2, 4, 1, 2, 5, 8
      ),
      nrow = 2,
      byrow = TRUE
    )
  )
  expect_equal(out[, c(1, 8)], samples[, c(1, 8)])
})


test_that("BMA subgroup sample conversion handles multinomial blocks", {
  modifier_info <- ameras:::new_modifier_info(
    coding = "group",
    design_vars = c(".ameras_modifier_1", ".ameras_modifier_2"),
    parameter_names = c("F3=mid", "F3=high"),
    group_labels = c("F3=low", "F3=mid", "F3=high")
  )
  dat <- data.frame(Y = factor(c("a", "b", "c")))

  # Two modeled multinomial response levels repeat the same parameter block.
  samples <- matrix(
    c(
      0, 1, 0.5, 1.5,
      10, 2, -1, 3
    ),
    nrow = 1
  )

  out <- ameras:::modifier_internal_to_reported_sample_matrix(
    samples = samples,
    family = "multinomial",
    M = 1:2,
    deg = 1,
    modifier_info = modifier_info,
    Y = "Y",
    data = dat
  )

  expect_equal(out, matrix(c(0, 1, 1.5, 2.5, 10, 2, 1, 5), nrow = 1))
})


test_that("subgroup-coded multinomial modifiers keep response-level prefixes", {
  data <- three_level_modifier_data()

  fit <- suppressMessages(ameras(
    Y.multinomial ~ dose(V1:V2, model = "EXP", modifier = ~ 0 + F3),
    data = data,
    family = "multinomial",
    methods = "RC"
  ))

  expect_identical(
    rownames(coef(fit)),
    c(
      "(1)_(Intercept)",
      "(1)_dose[F3=low]", "(1)_dose[F3=mid]", "(1)_dose[F3=high]",
      "(2)_(Intercept)",
      "(2)_dose[F3=low]", "(2)_dose[F3=mid]", "(2)_dose[F3=high]"
    )
  )
})


test_that("subgroup-coded modifiers work for profile CIs with keep.data FALSE", {
  data <- three_level_modifier_data()

  fit <- suppressMessages(ameras(
    Y.binomial ~ dose(V1:V2, model = "EXP", modifier = ~ 0 + F3),
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
    c("dose[F3=low]", "dose[F3=mid]", "dose[F3=high]")
  )
})


test_that("subgroup-coded modifiers smoke-test ERC, MCML, and FMA", {
  data <- three_level_modifier_data()

  fit_freq <- suppressWarnings(suppressMessages(ameras(
    Y.gaussian ~ dose(V1:V2, modifier = ~ 0 + F3),
    data = data,
    family = "gaussian",
    methods = c("RC", "ERC", "MCML")
  )))
  expect_true(all(c("RC", "ERC", "MCML") %in% names(fit_freq)))
  expect_true(all(grepl("^dose\\[F3=", rownames(coef(fit_freq))[2:4])))

  fit_fma <- suppressWarnings(suppressMessages(ameras(
    Y.gaussian ~ dose(V1:V2, modifier = ~ 0 + F3),
    data = data,
    family = "gaussian",
    methods = "FMA",
    MFMA = 1000
  )))
  expect_true("FMA" %in% names(fit_fma))
  expect_true(all(grepl("^dose\\[F3=", names(fit_fma$FMA$coefficients)[2:4])))
})


test_that("subgroup-coded modifiers report subgroup BMA samples", {
  skip_on_cran()
  set.seed(123)
  data <- three_level_modifier_data()

  # BMA still samples the NIMBLE model on the internal contrast scale. The
  # stored samples and summaries should be post-processed back to the subgroup
  # scale used by coef(), vcov(), summary(), confint(), and traceplot().
  fit <- suppressWarnings(suppressMessages(ameras(
    Y.gaussian ~ dose(V1:V2, modifier = ~ 0 + F3),
    data = data[seq_len(24), ],
    family = "gaussian",
    methods = "BMA",
    niter.BMA = 20,
    nburnin.BMA = 5,
    nchains.BMA = 1,
    thin.BMA = 1,
    print = FALSE
  )))

  expected <- c(
    "(Intercept)",
    "dose[F3=low]", "dose[F3=mid]", "dose[F3=high]",
    "sigma"
  )
  expect_named(fit$BMA$coefficients, expected)
  expect_named(fit$BMA$sd, expected)
  expect_identical(rownames(fit$BMA$vcov), expected)
  expect_identical(colnames(fit$BMA$samples), c(expected, "col.ind"))

  fit_ci <- suppressMessages(confint(fit, type = "percentile", print = FALSE))
  expect_identical(rownames(fit_ci$BMA$CI), expected)
})


test_that("formula contrast modifiers remain supported for BMA", {
  skip_on_cran()
  data <- three_level_modifier_data()

  # BMA still uses the existing reference-plus-contrast parameterization. This
  # keeps the new formula syntax compatible with the current BMA model code.
  fit <- suppressWarnings(suppressMessages(ameras(
    Y.gaussian ~ dose(V1:V2, modifier = ~ F3),
    data = data[seq_len(24), ],
    family = "gaussian",
    methods = "BMA",
    niter.BMA = 20,
    nburnin.BMA = 5,
    nchains.BMA = 1,
    print = FALSE
  )))

  expect_named(
    fit$BMA$coefficients,
    c("(Intercept)", "dose", "dose:F3=mid", "dose:F3=high", "sigma")
  )
})


test_that("dose_lrt treats subgroup-coded dose coefficients as dose terms", {
  data <- three_level_modifier_data()

  fit <- suppressMessages(ameras(
    Y.gaussian ~ dose(V1:V2, modifier = ~ 0 + F3),
    data = data,
    family = "gaussian",
    methods = "RC"
  ))

  lrt <- dose_lrt(fit, methods = "RC", type = "individual")
  expect_identical(
    lrt$term,
    c("dose[F3=low]", "dose[F3=mid]", "dose[F3=high]")
  )
})
