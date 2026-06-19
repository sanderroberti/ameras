## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE
)

## ----setup--------------------------------------------------------------------
library(ameras)

## ----data---------------------------------------------------------------------
data(data, package = "ameras")
dosevars <- paste0("V", 1:10)

## ----fit-helper---------------------------------------------------------------
passes_fma_screen <- function(rc) {
  hessian <- rc$optim$hessian

  if (is.null(hessian) || rc$optim$convergence != 0) {
    return(FALSE)
  }

  isTRUE(det(hessian) != 0 &&
    rcond(hessian) > .Machine$double.eps &&
    all(eigen(hessian)$values > 0))
}

fit_rc_realization <- function(dosevar, data) {
  formula <- stats::as.formula(
    paste0("Y.gaussian ~ dose(", dosevar, ") + X1 + X2")
  )

  fit <- ameras(formula, data = data, family = "gaussian")
  rc <- fit$RC

  list(
    dosevar = dosevar,
    coefficients = rc$coefficients,
    vcov = rc$vcov,
    loglik = rc$loglik,
    AIC = -2 * rc$loglik + 2 * length(rc$coefficients),
    include = passes_fma_screen(rc)
  )
}

## ----fit-realizations---------------------------------------------------------
rc_summaries <- lapply(dosevars, fit_rc_realization, data = data)

## ----save-fits, eval = FALSE--------------------------------------------------
# dir.create("fma-fits", showWarnings = FALSE)
# 
# for (dosevar in dosevars) {
#   fit_summary <- fit_rc_realization(dosevar, data = data)
#   saveRDS(fit_summary, file = file.path("fma-fits", paste0(dosevar, ".rds")))
# }
# 
# fit_files <- file.path("fma-fits", paste0(dosevars, ".rds"))
# rc_summaries <- lapply(fit_files, readRDS)

## ----assemble-helper----------------------------------------------------------
assemble_manual_fma <- function(rc_summaries, MFMA = 100000) {
  n_total <- length(rc_summaries)
  valid <- vapply(rc_summaries, function(x) isTRUE(x$include), logical(1))
  n_screen_excluded <- sum(!valid)

  if (!any(valid)) {
    stop("No realization-specific fits are available for FMA.")
  }

  original_indices <- which(valid)
  rc_summaries <- rc_summaries[valid]

  AIC <- vapply(rc_summaries, `[[`, numeric(1), "AIC")
  weights <- exp(-0.5 * (AIC - min(AIC)))
  weights <- weights / sum(weights)

  allocated_samples <- round(weights * MFMA)
  keep <- allocated_samples > 0
  n_weight_excluded <- sum(!keep)

  if (!any(keep)) {
    stop("No realizations received FMA samples. Increase MFMA.")
  }

  rc_summaries <- rc_summaries[keep]
  original_indices <- original_indices[keep]
  weights <- weights[keep]
  allocated_samples <- allocated_samples[keep]
  names(weights) <- vapply(rc_summaries, `[[`, character(1), "dosevar")
  diagnostics <- data.frame(
    stage = c(
      "total realization-specific fits",
      "excluded by convergence/Hessian screen",
      "excluded after zero-sample allocation",
      "included in manual FMA"
    ),
    realizations = c(
      n_total,
      n_screen_excluded,
      n_weight_excluded,
      length(original_indices)
    ),
    row.names = NULL
  )

  samples <- do.call(
    "rbind",
    Map(
      function(fit, n) {
        mvtnorm::rmvnorm(n = n, mean = fit$coefficients, sigma = fit$vcov)
      },
      rc_summaries,
      allocated_samples
    )
  )

  samples <- as.data.frame(samples)
  names(samples) <- names(rc_summaries[[1]]$coefficients)

  coefficients <- colMeans(samples)
  se <- apply(samples, 2, stats::sd)
  percentile_CI <- t(
    apply(samples, 2, stats::quantile, probs = c(0.025, 0.975))
  )
  colnames(percentile_CI) <- c("lower", "upper")

  list(
    coefficients = coefficients,
    SE = se,
    percentile_CI = percentile_CI,
    vcov = stats::var(samples),
    weights = weights,
    samples = samples,
    included.realizations = original_indices,
    included.samples = nrow(samples),
    diagnostics = diagnostics
  )
}

## ----manual-fma---------------------------------------------------------------
set.seed(100)
manual_fma <- assemble_manual_fma(rc_summaries, MFMA = 100000)
manual_fma$diagnostics

manual_summary <- data.frame(
  term = names(manual_fma$coefficients),
  estimate = unname(manual_fma$coefficients),
  SE = unname(manual_fma$SE),
  percentile_lower = manual_fma$percentile_CI[, "lower"],
  percentile_upper = manual_fma$percentile_CI[, "upper"],
  row.names = NULL
)

manual_summary[-1] <- round(manual_summary[-1], 4)
manual_summary

## ----manual-weights-----------------------------------------------------------
round(manual_fma$weights, 4)
manual_fma$included.realizations
manual_fma$included.samples

## ----builtin-fma--------------------------------------------------------------
set.seed(100)
fit_builtin_fma <- suppressWarnings(
  ameras(Y.gaussian ~ dose(V1:V10) + X1 + X2,
         data = data,
         family = "gaussian",
         methods = "FMA",
         MFMA = 100000)
)

comparison <- data.frame(
  term = names(manual_fma$coefficients),
  manual = unname(manual_fma$coefficients),
  builtin = unname(fit_builtin_fma$FMA$coefficients[names(manual_fma$coefficients)]),
  row.names = NULL
)

comparison[-1] <- round(comparison[-1], 4)
comparison

## ----compare-weights----------------------------------------------------------
weight_comparison <- merge(
  data.frame(dose = names(manual_fma$weights),
             manual = unname(manual_fma$weights)),
  data.frame(dose = names(fit_builtin_fma$FMA$weights),
             builtin = unname(fit_builtin_fma$FMA$weights)),
  by = "dose",
  all = TRUE
)

weight_comparison[-1] <- round(weight_comparison[-1], 4)
weight_comparison

