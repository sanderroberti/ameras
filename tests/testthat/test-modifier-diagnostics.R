test_that("diagnostic dose designs follow fitted modifier names and ordering", {
  dat <- ameras_simulated[seq_len(120), ]
  dat$F3 <- factor(rep(c("a", "b", "c"), length.out = nrow(dat)))
  modifiers <- list(~ M1, ~ M1 + M2, ~ F3, ~ 0 + M1, ~ 0 + F3)
  for (modifier in modifiers) {
    for (model in c("EXP", "ERR", "LINEXP")) {
      for (degree in 1:2) {
        formula <- substitute(
          Surv(time, status) ~ dose(V1:V3, model = MODEL, deg = DEG, modifier = MOD) + X1,
          list(MODEL = model, DEG = degree, MOD = modifier)
        )
        parsed <- ameras:::parse_ameras_formula(as.formula(formula), dat, "prophaz")
        prepared <- ameras:::prepare_modifier_inputs(dat, parsed$modifier_info)
        specification <- list(
          M = prepared$design_vars, modifier_info = prepared$modifier_info,
          doseRRmod = model, deg = if (model == "LINEXP") 2 else degree
        )
        design <- ameras:::modifier_dose_design(prepared$data, "V1", specification)
        expected_names <- ameras:::make_method_parnames(
          "prophaz", prepared$data, X = "X1", M = specification$M,
          deg = specification$deg, doseRRmod = model,
          modifier_info = specification$modifier_info
        )[-1]
        expect_identical(colnames(design), expected_names)
        expect_equal(nrow(design), nrow(dat))
        if (identical(specification$modifier_info$coding, "group")) {
          expect_equal(rowSums(design[, seq(1, ncol(design), by = specification$deg), drop = FALSE]), dat$V1)
        } else {
          expect_equal(unname(design[, 1]), dat$V1)
        }
      }
    }
  }
})

test_that("Schoenfeld residuals and plots work for factor and subgroup modifiers", {
  dat <- ameras_simulated[seq_len(300), ]
  dat$F2 <- factor(dat$M1, levels = c(0, 1), labels = c("reference", "other"))
  dat$F3 <- factor(rep(c("a", "b", "c"), length.out = nrow(dat)))
  dat$X1[1] <- NA
  formulas <- list(
    Surv(time, status) ~ dose(V1:V3, model = "EXP", modifier = ~ F2) + X1,
    Surv(time, status) ~ dose(V1:V3, model = "EXP", deg = 2, modifier = ~ F3) + X1,
    Surv(time, status) ~ dose(V1:V3, model = "EXP", modifier = ~ 0 + M1) + X1,
    Surv(time, status) ~ dose(V1:V3, model = "EXP", deg = 2, modifier = ~ 0 + F3) + X1
  )
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (formula in formulas) {
    fit <- suppressWarnings(suppressMessages(ameras(
      formula, data = dat, family = "prophaz", methods = "RC",
      optim.method = "BFGS", keep.data = FALSE
    )))
    expect_identical(fit$model$M_names, all.vars(formula[[3]])[grepl("^(M|F)", all.vars(formula[[3]]))])
    for (scaled in c(FALSE, TRUE)) {
      result <- residuals(fit, data = dat, scaled.schoenfeld = scaled)
      expect_setequal(names(result), c("id", "time", names(coef(fit))))
      expect_true(all(is.finite(as.matrix(result))))
      expect_equal(nrow(result), sum(dat$status[-1] == 1))
    }
    expect_identical(plot(fit, data = dat, ask = FALSE), fit)
  }
})

test_that("subgroup raw Schoenfeld residuals agree with recoded contrast residuals", {
  dat <- ameras_simulated[seq_len(300), ]
  fit_contrast <- suppressMessages(ameras(
    Surv(time, status) ~ dose(V1:V3, model = "EXP", modifier = ~ M1),
    data = dat, family = "prophaz", optim.method = "BFGS"
  ))
  fit_group <- suppressMessages(ameras(
    Surv(time, status) ~ dose(V1:V3, model = "EXP", modifier = ~ 0 + M1),
    data = dat, family = "prophaz", optim.method = "BFGS"
  ))
  contrast <- residuals(fit_contrast, scaled.schoenfeld = FALSE)
  group <- residuals(fit_group, scaled.schoenfeld = FALSE)
  expect_equal(group[["dose[M1=0]"]], contrast$dose - contrast[["dose:M1"]], tolerance = 1e-5)
  expect_equal(group[["dose[M1=1]"]], contrast[["dose:M1"]], tolerance = 1e-5)
})
