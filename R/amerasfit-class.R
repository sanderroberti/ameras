new_amerasfit <- function(x = list()) {
  stopifnot(is.list(x))
  structure(
    x,
    class = "amerasfit"
  )
}


print.amerasfit <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  object0 <- x[intersect(names(x), c("RC", "ERC", "MCML", "FMA", "BMA"))]

  coefs <- coef.amerasfit(x, ...)

  runtime_table <- do.call(
    "rbind",
    lapply(1:length(object0), function(i) {
      y <- object0[[i]]
      method <- names(object0)[i]

      runtime <- as.numeric(strsplit(y$runtime, " seconds")[[1]])

      res <- data.frame(Method = method, Runtime = runtime)
      rownames(res) <- NULL
      res
    })
  )

  total_runtime_seconds <- sum(sapply(object0, function(x) {
    as.numeric(strsplit(x$runtime, " seconds")[[1]])
  }))

  cat("Call:\n")
  print(x$call)

  cat(paste0("\nNumber of rows: ", x$num.rows, "\n"))
  cat(paste0("Number of dose realizations: ", x$num.realizations, "\n"))

  cat(paste0("\nTotal run time: ", total_runtime_seconds, " seconds\n\n"))

  cat("Runtime in seconds by method:\n\n")
  print(format(runtime_table, digits = digits, nsmall = 1), row.names = FALSE)

  cat("\nEstimated model parameters:\n\n")
  print(format(coefs, digits = digits, nsmall = 1), row.names = TRUE)

  cat("\n")

  invisible(x)
}


coef.amerasfit <- function(object, ...) {
  object <- object[intersect(
    names(object),
    c("RC", "ERC", "MCML", "FMA", "BMA")
  )]

  bma <- ("BMA" %in% names(object))

  res <- as.data.frame(do.call(
    "cbind",
    lapply(1:length(object), function(i) {
      y <- object[[i]]

      if (bma) {
        tmp <- c(
          y$coefficients,
          object$BMA$coefficients[-(1:length(y$coefficients))] * NA
        )
      } else {
        tmp <- y$coefficients
      }

      tmp
    })
  ))
  names(res) <- names(object)
  res
}

print_confint <- function(
  object,
  parm = NULL,
  digits = max(3, getOption("digits") - 3)
) {
  methods_present <- intersect(
    names(object),
    c("RC", "ERC", "MCML", "FMA", "BMA")
  )

  for (m in methods_present) {
    if (is.null(object[[m]]$CI)) {
      next
    }

    cat(m, "confidence intervals:\n\n")

    CI <- object[[m]]$CI

    # Filter rows based on parm
    if (!is.null(parm)) {
      if (identical(parm, "dose")) {
        # Match dose-related parameters
        keep <- startsWith(rownames(CI), "dose") |
          grepl(")_dose", rownames(CI))
        CI <- CI[keep, , drop = FALSE]
      } else if (!identical(parm, "all")) {
        # Match specific parameter names
        keep <- rownames(CI) %in% parm
        CI <- CI[keep, , drop = FALSE]
      }
    }

    if (nrow(CI) == 0) {
      cat("No parameters matched parm.\n\n")
      next
    }

    print(format(CI, digits = digits), row.names = TRUE)
    cat("\n")
  }

  invisible(object)
}

summary.amerasfit <- function(object, ...) {
  object0 <- object[intersect(
    names(object),
    c("RC", "ERC", "MCML", "FMA", "BMA")
  )]

  bma <- ("BMA" %in% names(object0))

  summary_table <- do.call(
    "rbind",
    lapply(1:length(object0), function(i) {
      y <- object0[[i]]
      method <- names(object0)[i]

      coef <- y$coefficients

      se <- y$sd

      res <- data.frame(
        Method = method,
        Term = names(coef),
        Estimate = coef,
        SE = se
      )

      if (object$CI.computed) {
        CI <- y$CI
        CI.lower <- CI.upper <- coef * NA
        CI.lower[match(rownames(CI), names(coef))] <- CI$lower
        CI.upper[match(rownames(CI), names(coef))] <- CI$upper

        res <- cbind(
          res,
          data.frame(
            CI.lower = CI.lower,
            CI.upper = CI.upper
          )
        )
        # Add p-value columns for profile likelihood intervals. If not used, these will be removed after
        res <- cbind(res, data.frame(pval.lower = NA, pval.upper = NA))

        if ("pval.lower" %in% names(CI) && "pval.upper" %in% names(CI)) {
          pval.lower <- pval.upper <- coef * NA
          pval.lower[match(rownames(CI), names(coef))] <- CI$pval.lower
          pval.upper[match(rownames(CI), names(coef))] <- CI$pval.upper
          res$pval.lower <- pval.lower
          res$pval.upper <- pval.upper
          #res <- cbind(res, data.frame(pval.lower = pval.lower, pval.upper = pval.upper))
        }
      }

      if (bma) {
        if (method == "BMA") {
          res <- cbind(res, y$Rhat)
        } else {
          res <- cbind(res, data.frame(Rhat = NA, n.eff = NA))
        }
      }
      rownames(res) <- NULL
      res
    })
  )

  # Remove p-value columns if not used
  if (object$CI.computed) {
    if (
      all(is.na(summary_table$pval.lower)) &
        all(is.na(summary_table$pval.upper))
    ) {
      summary_table <- summary_table[,
        -which(names(summary_table) %in% c("pval.lower", "pval.upper"))
      ]
    }
  }
  runtime_table <- do.call(
    "rbind",
    lapply(1:length(object0), function(i) {
      y <- object0[[i]]
      method <- names(object0)[i]

      runtime <- as.numeric(strsplit(y$runtime, " seconds")[[1]])

      res <- data.frame(Method = method, Runtime = runtime)
      rownames(res) <- NULL
      res
    })
  )

  total_runtime_seconds <- sum(sapply(object0, function(x) {
    as.numeric(strsplit(x$runtime, " seconds")[[1]])
  }))

  ans <- list(
    call = object$call,
    summary_table = summary_table,
    runtime_table = runtime_table,
    total_runtime_seconds = total_runtime_seconds,
    CI.computed = object$CI.computed
  )

  class(ans) <- "summary.amerasfit"
  return(ans)
}


print.summary.amerasfit <- function(
  x,
  digits = max(3, getOption("digits") - 3),
  ...
) {
  cat("Call:\n")
  print(x$call)
  cat(paste0("\nTotal run time: ", x$total_runtime_seconds, " seconds\n\n"))

  cat("Runtime in seconds by method:\n\n")
  print(format(x$runtime_table, digits = digits, nsmall = 1), row.names = FALSE)

  cat("\nSummary of coefficients by method:\n\n")
  print(format(x$summary_table, digits = digits, nsmall = 2), row.names = FALSE)

  if (!x$CI.computed) {
    cat(
      "\nNote: confidence intervals not yet computed.",
      "Use confint() to add them.\n"
    )
  }

  cat("\n")

  invisible(x)
}
Rhat <- function(x, ...) UseMethod("Rhat")

Rhat.amerasfit <- function(x, ...) {
  if (is.null(x$BMA)) {
    stop("BMA not present in fitted object")
  }
  x$BMA$Rhat
}

included_realizations <- function(x, ...) UseMethod("included_realizations")

included_realizations.amerasfit <- function(x, methods = c("FMA", "BMA"), ...) {
  methods <- match.arg(methods, several.ok = TRUE)

  available <- intersect(methods, names(x)[names(x) %in% c("FMA", "BMA")])

  if (!length(available)) {
    stop("None of the requested methods were run")
  }

  result <- lapply(available, function(m) x[[m]]$included.realizations)
  names(result) <- available

  if (length(available) == 1) result[[1]] else result
}

traceplot <- function(object, ...) {
  UseMethod("traceplot")
}

traceplot.amerasfit <- function(
  object,
  iter = 5000,
  Rhat = TRUE,
  n.eff = TRUE,
  pdf = FALSE,
  ...
) {
  if ("BMA" %in% names(object)) {
    MCMCtrace(
      object$BMA$samples,
      iter = iter,
      Rhat = Rhat,
      n.eff = n.eff,
      pdf = pdf,
      ...
    )
  } else {
    stop(
      "ERROR: traceplot() requires BMA output in the provided 'amerasfit' object."
    )
  }
}


confint.amerasfit <- function(
  object,
  parm = "dose",
  level = 0.95,
  type = c("proflik", "percentile"),
  maxit.profCI = 20,
  tol.profCI = 1e-2,
  data = NULL,
  force = FALSE,
  digits = max(3, getOption("digits") - 3),
  ...
) {
  # Check if CIs have already been computed
  if (object$CI.computed && !force) {
    message(
      "Confidence intervals have already been computed for this object. ",
      "Returning the object unchanged. Use force=TRUE to recompute."
    )
    print_confint(object, parm = parm, digits = digits)
    return(invisible(object))
  }

  # Validate level
  if (!is.numeric(level) || length(level) != 1 || level <= 0 || level >= 1) {
    stop("level must be a single numeric value between 0 and 1")
  }

  fitobj <- object[intersect(
    names(object),
    c("RC", "ERC", "MCML", "FMA", "BMA")
  )]

  res <- NULL
  for (it in 1:length(fitobj)) {
    method <- names(fitobj)[it]

    if (method %in% c("FMA", "BMA")) {
      type.i <- match.arg(
        type,
        choices = c("percentile", "hpd"),
        several.ok = TRUE
      )
      if (length(type.i) != 1) {
        stop(
          "ERROR: Exactly one CI type (either 'percentile' or 'hpd') for FMA and/or BMA should be specified"
        )
      }

      if (method == "BMA") {
        pars <- names(fitobj[[it]]$coefficients)
        samples <- fitobj[[it]]$samples
        if (is.list(samples)) {
          samples <- do.call("rbind", samples)
        }
        samples <- samples[, pars]
      } else {
        samples <- fitobj[[it]]$samples
      }
      if (!is.null(samples)) {
        object[[method]]$CI <- compute_sample_CI(
          samples = samples,
          level = level,
          type = type.i
        )
      } else {
        if (is.null(samples) || nrow(samples) == 0) {
          warning(
            "No samples available for CI computation. ",
            "All realizations may have been excluded during fitting."
          )
        }
        object[[method]]$CI <- data.frame(
          lower = NA * object[[method]]$coefficients,
          upper = NA * object[[method]]$coefficients
        )
      }
    } else {
      # RC, ERC, MCML
      type.i <- match.arg(
        type,
        choices = c("wald.orig", "wald.transformed", "proflik"),
        several.ok = TRUE
      )
      if (length(type.i) != 1) {
        stop(
          "ERROR: Exactly one CI type (one of 'proflik', 'wald.orig', or 'wald.transformed') for RC, ERC, and/or MCML should be specified"
        )
      }
      if (!is.null(parm) && type.i != "proflik") {
        #warning("parm is only used when method='proflik' and will be ignored")
      }

      if (type.i == "wald.orig") {
        object[[method]]$CI <- compute_wald_CI(
          method_fit = fitobj[[it]],
          transform = object$transform,
          level = level,
          type = "orig",
          other.args = object$other.args
        )
      } else if (type.i == "wald.transformed") {
        object[[method]]$CI <- compute_wald_CI(
          method_fit = fitobj[[it]],
          transform = object$transform,
          level = level,
          type = "transformed",
          other.args = object$other.args
        )
      } else if (type.i == "proflik") {
        resolved_data <- resolve_data(object, data)
        object[[method]]$CI <- compute_proflik_CI(
          method_fit = fitobj[[it]],
          object = object,
          method_name = method,
          data = resolved_data,
          parm = parm,
          level = level,
          maxit.profCI = maxit.profCI,
          tol.profCI = tol.profCI,
          optim.method = object$model$optim.method,
          control = object$model$control
        )
        rm(resolved_data)
      }
    }
  }

  object$CI.computed <- TRUE

  print_confint(object, parm = parm, digits = digits)
  invisible(object)
}


CI <- function(x, ...) UseMethod("CI")

CI.amerasfit <- function(
  x,
  methods = c("RC", "ERC", "MCML", "FMA", "BMA"),
  parm = NULL,
  ...
) {
  methods <- match.arg(methods, several.ok = TRUE)

  available <- intersect(
    methods,
    names(x)[names(x) %in% c("RC", "ERC", "MCML", "FMA", "BMA")]
  )

  if (!length(available)) {
    stop("None of the requested methods were run")
  }

  result <- lapply(available, function(m) {
    ci <- x[[m]]$CI
    if (is.null(ci)) {
      return(NULL)
    }

    if (!is.null(parm)) {
      if (identical(parm, "dose")) {
        keep <- startsWith(rownames(ci), "dose") |
          grepl(")_dose", rownames(ci))
        ci <- ci[keep, , drop = FALSE]
      } else if (!identical(parm, "all")) {
        keep <- rownames(ci) %in% parm
        ci <- ci[keep, , drop = FALSE]
      }
    }
    ci
  })
  names(result) <- available

  if (length(available) == 1) result[[1]] else result
}
