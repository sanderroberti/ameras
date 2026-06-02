new_amerasfit <- function(x = list()) {
  stopifnot(is.list(x))
  structure(
    x,
    class = "amerasfit"
  )
}

vcov.amerasfit <- function(
  object,
  methods = c("RC", "ERC", "MCML", "FMA", "BMA"),
  ...
) {
  methods <- match.arg(
    methods,
    c("RC", "ERC", "MCML", "FMA", "BMA"),
    several.ok = TRUE
  )

  available <- intersect(
    methods,
    names(object)[names(object) %in% c("RC", "ERC", "MCML", "FMA", "BMA")]
  )

  if (!length(available)) {
    stop("None of the requested methods were run")
  }

  result <- lapply(available, function(m) object[[m]]$vcov)
  names(result) <- available

  if (length(result) == 1) result[[1]] else result
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

  cat(paste0("\nTotal runtime: ", total_runtime_seconds, " seconds\n\n"))

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

print_confint <- function(object, digits = max(3, getOption("digits") - 3)) {
  methods_present <- intersect(
    names(object),
    c("RC", "ERC", "MCML", "FMA", "BMA")
  )

  for (m in methods_present) {
    if (is.null(object[[m]]$CI)) {
      next()
    }

    cat(m, "confidence intervals:\n\n")
    print(format(object[[m]]$CI, digits = digits), row.names = TRUE)
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
  print = TRUE,
  digits = max(3, getOption("digits") - 3),
  ...
) {
  # Check if CIs have already been computed
  if (object$CI.computed && !force) {
    message(
      "Confidence intervals have already been computed for this object. ",
      "Returning the object unchanged. Use force=TRUE to recompute."
    )
    print_confint(object, digits = digits)
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
  if (print) {
    print_confint(object, digits = digits)
  }
  invisible(object)
}

summary_table <- function(object, ...) UseMethod("summary_table")

summary_table.amerasfit <- function(object, ...) {
  summ <- summary(object, ...)
  return(summ$summary_table)
}


residuals.amerasfit <- function(
  object,
  method = "RC",
  type = NULL,
  data = NULL,
  dose.col = NULL,
  scaled.schoenfeld = TRUE,
  ...
) {
  if (is.null(type)) {
    if (object$model$family == "prophaz") {
      type <- "schoenfeld"
    } else {
      type <- "pearson"
    }
  }
  type <- match.arg(type, c("pearson", "response", "deviance", "schoenfeld"))

  if (object$model$family == "prophaz") {
    if (type != "schoenfeld") {
      stop("Only schoenfeld residuals are supported for family 'prophaz'")
    }
  } else {
    # Families other than prophaz
    if (type == "schoenfeld") {
      stop(paste0(
        "schoenfeld residuals not supported for family='",
        object$model$family,
        "'"
      ))
    }
  }
  method <- match.arg(method, choices = c("RC", "ERC", "MCML", "FMA", "BMA"))

  if (is.null(object[[method]])) {
    stop("Method '", method, "' not present in fitted object")
  }

  resolved_data <- resolve_data(object, data = data)

  if (is.null(dose.col)) {
    dose.col <- select_dose_col(object, method, resolved_data)
  } else {
    if (!dose.col %in% colnames(resolved_data)) {
      stop("dose.col '", dose.col, "' not found in data")
    }
  }

  m <- object$model

  mus <- compute_fitted(
    object,
    method = method,
    data = resolved_data,
    dose.col = dose.col
  )
  Y <- if (m$family == "clogit") {
    resolved_data[, m$status]
  } else {
    resolved_data[, m$Y]
  }
  if (m$family == "gaussian") {
    raw <- Y - mus
    sigma <- object[[method]]$coefficients["sigma"]
    if (type == "response") {
      return(raw)
    }
    if (type == "pearson") {
      return(raw / sigma)
    }
    if (type == "deviance") return(raw)
  } else if (m$family %in% c("binomial", "clogit")) {
    raw <- Y - mus
    if (type == "response") {
      return(raw)
    }
    if (type == "pearson") {
      return(raw / sqrt(mus * (1 - mus)))
    }
    if (type == "deviance") {
      dev <- ifelse(
        Y == 1,
        sqrt(2 * log(1 / mus)),
        sqrt(2 * log(1 / (1 - mus)))
      )
      return(ifelse(Y > mus, dev, -dev))
    }
  } else if (m$family == "poisson") {
    raw <- Y - mus
    if (type == "response") {
      return(raw)
    }
    if (type == "pearson") {
      return(raw / sqrt(mus))
    }
    if (type == "deviance") {
      dev <- sqrt(2 * ifelse(Y > 0, Y * log(Y / mus) - (Y - mus), mus))
      return(ifelse(Y > mus, dev, -dev))
    }
  } else if (m$family == "multinomial") {
    # mus is N x Z matrix of fitted probabilities
    # we return an N x Z matrix of residuals rather than a vector
    Ymat <- diag(nlevels(resolved_data[, m$Y]))[
      as.integer(resolved_data[, m$Y]),
    ]
    colnames(Ymat) <- levels(resolved_data[, m$Y])

    raw <- Ymat - mus

    if (type == "response") {
      return(raw)
    }
    if (type == "pearson") {
      return(raw / sqrt(mus * (1 - mus)))
    }
    if (type == "deviance") {
      dev <- sign(raw) *
        sqrt(
          -2 *
            ifelse(Ymat > 0, log(pmax(mus, 1e-10)), log(pmax(1 - mus, 1e-10)))
        )
      return(dev)
    }
  } else if (m$family == "prophaz") {
    if (type == "schoenfeld") {
      if (m$doseRRmod == "LINEXP") {
        dose_mat <- cbind(
          resolved_data[, dose.col],
          resolved_data[, dose.col]
        )
        colnames(dose_mat) <- c("dose_linear", "dose_exponential")
      } else if (m$deg == 2) {
        dose_mat <- cbind(
          resolved_data[, dose.col],
          resolved_data[, dose.col]^2
        )
        colnames(dose_mat) <- c("dose", "dose_squared")
      } else {
        dose_mat <- matrix(resolved_data[, dose.col], ncol = 1)
        colnames(dose_mat) <- "dose"
      }

      # Modifier columns: observed M * D values
      if (!is.null(m$M)) {
        M_mat <- as.matrix(resolved_data[, m$M, drop = FALSE])
        M_names <- m$M_names

        if (m$doseRRmod == "LINEXP") {
          mod_lin <- sweep(M_mat, 1, resolved_data[, dose.col], "*")
          mod_exp <- sweep(M_mat, 1, resolved_data[, dose.col], "*")
          colnames(mod_lin) <- paste0("dose_linear:", M_names)
          colnames(mod_exp) <- paste0("dose_exponential:", M_names)
          dose_mat <- cbind(dose_mat, mod_lin, mod_exp)
        } else if (m$deg == 2) {
          mod_lin <- sweep(M_mat, 1, resolved_data[, dose.col], "*")
          mod_sq <- sweep(M_mat, 1, resolved_data[, dose.col]^2, "*")
          colnames(mod_lin) <- paste0("dose:", M_names)
          colnames(mod_sq) <- paste0("dose_squared:", M_names)
          dose_mat <- cbind(dose_mat, mod_lin, mod_sq)
        } else {
          mod_lin <- sweep(M_mat, 1, resolved_data[, dose.col], "*")
          colnames(mod_lin) <- paste0("dose:", M_names)
          dose_mat <- cbind(dose_mat, mod_lin)
        }
      }

      if (!is.null(m$X)) {
        X_mat <- resolved_data[, m$X, drop = FALSE]
      } else {
        X_mat <- NULL
      }
      covariates <- if (!is.null(X_mat)) {
        cbind(dose_mat, X_mat)
      } else {
        dose_mat
      }

      vcov.mat <- object[[method]]$vcov[
        colnames(covariates),
        colnames(covariates),
        drop = FALSE
      ]

      sf <- compute_schoenfeld_residuals(
        exit = resolved_data[, m$exit],
        status = resolved_data[, m$status],
        covariates = covariates,
        rr = mus,
        entry = if (!is.null(m$entry)) resolved_data[, m$entry] else NULL,
        scaled = scaled.schoenfeld,
        vcov.mat = vcov.mat
      )
      return(sf)
    }
  } else {
    stop("Residuals not implemented for family='", m$family, "'")
  }
}


plot.amerasfit <- function(
  x,
  methods = c("RC", "ERC", "MCML", "FMA", "BMA"),
  which = NULL,
  type = NULL,
  dose.col = NULL,
  add.smooth = getOption("add.smooth", TRUE),
  qqline = TRUE,
  id.n = 3,
  ask = NULL,
  data = NULL,
  ...
) {
  defaults <- switch(
    x$model$family,
    "prophaz" = list(which = "schoenfeld", type = "schoenfeld"),
    list(which = c("residuals-vs-fitted", "qq"), type = "pearson")
  )

  if (is.null(which)) {
    which <- defaults$which
  }
  if (is.null(type)) {
    type <- defaults$type
  }

  if (x$model$family == "prophaz") {
    which <- match.arg(which, c("schoenfeld"))
    type <- match.arg(type, c("schoenfeld"))
  } else {
    which <- match.arg(which, c("residuals-vs-fitted", "qq"), several.ok = TRUE)
    type <- match.arg(type, c("pearson", "response", "deviance"))
  }

  methods <- match.arg(
    methods,
    c("RC", "ERC", "MCML", "FMA", "BMA"),
    several.ok = TRUE
  )
  available <- intersect(
    methods,
    names(x)[names(x) %in% c("RC", "ERC", "MCML", "FMA", "BMA")]
  )

  if (!length(available)) {
    stop("None of the requested methods were run")
  }

  if (is.null(ask)) {
    n_cats <- if (x$model$family == "multinomial") {
      # Coefficient names for multinomial have format "(level)_paramname"
      # Extract unique level prefixes, excluding the reference level
      coef_names <- names(x[[available[1]]]$coefficients)
      prefixes <- unique(regmatches(
        coef_names,
        regexpr("^\\([^)]+\\)", coef_names)
      ))
      length(prefixes)
    } else if (x$model$family == "prophaz") {
      sum(!grepl("^h0\\[[0-9]+\\]$", names(x[[available[1]]]$coefficients))) # Count all variables except h0[1], h0[2], ... which is relevant when available[1] is "BMA"
    } else {
      1
    }

    n_panels <- length(which) * length(available) * n_cats
    ask <- prod(par("mfcol")) < n_panels && dev.interactive()
  }

  resolved_data <- resolve_data(x, data = data)

  id.n <- min(id.n, nrow(resolved_data))

  if (ask) {
    op <- par(ask = TRUE)
    on.exit(par(op), add = TRUE)
  }

  # Precompute for all methods before plotting
  precomputed <- lapply(available, function(method) {
    if (is.null(dose.col)) {
      dose.col <- select_dose_col(x, method, resolved_data)
      dose.col.plot <- if (dose.col == "rcdose_ameras") {
        "mean of realizations"
      } else {
        paste0("realization: ", dose.col)
      }
    } else {
      dose.col.plot <- paste0("realization: ", dose.col)
    }

    fitted_vals <- compute_fitted(
      x,
      method = method,
      data = resolved_data,
      dose.col = dose.col
    )
    resids <- residuals.amerasfit(
      x,
      method = method,
      type = type,
      data = resolved_data,
      dose.col = dose.col,
      scaled.schoenfeld = TRUE
    )
    list(
      fitted_vals = fitted_vals,
      resids = resids,
      dose.col = dose.col,
      dose.col.plot = dose.col.plot,
      residMatrix = is.matrix(resids)
    )
  })
  names(precomputed) <- available

  if (x$model$family != "prophaz") {
    for (w in which) {
      for (method in available) {
        fitted_vals <- precomputed[[method]]$fitted_vals
        resids <- precomputed[[method]]$resids
        is_matrix <- precomputed[[method]]$residMatrix

        if (is_matrix) {
          # One panel per outcome category
          cats <- colnames(resids)
          cats <- cats[-length(cats)] # exclude reference category
          for (cat in cats) {
            fv_cat <- fitted_vals[, cat]
            res_cat <- resids[, cat]

            if (w == "residuals-vs-fitted") {
              plot(
                fv_cat,
                res_cat,
                xlab = "Fitted probability",
                ylab = paste0(tools::toTitleCase(type), " residuals"),
                main = paste0(
                  "Residuals vs Fitted\n Category (",
                  cat,
                  ") \n(",
                  method,
                  ", ",
                  precomputed[[method]]$dose.col.plot,
                  ")"
                ),
                ...
              )
              abline(h = 0, lty = 2, col = "grey")
              if (add.smooth) {
                panel.smooth(fv_cat, res_cat, col.smooth = "red")
              }
              if (id.n > 0) {
                extreme <- order(abs(res_cat), decreasing = TRUE)[seq_len(id.n)]
                side <- ifelse(
                  fv_cat[extreme] < (min(fv_cat) + max(fv_cat)) / 2,
                  4,
                  2
                ) # left tail -> label right; right tail -> label left
                text(
                  fv_cat[extreme],
                  res_cat[extreme],
                  labels = extreme,
                  pos = side,
                  offset = .35,
                  cex = 0.7
                )
              }
            }

            if (w == "qq") {
              qqnorm(
                res_cat,
                main = paste0(
                  "Normal Q-Q\n Category (",
                  cat,
                  ") \n(",
                  method,
                  ", ",
                  precomputed[[method]]$dose.col.plot,
                  ")"
                ),
                ylab = paste0(tools::toTitleCase(type), " residuals"),
                ...
              )
              if (qqline) {
                stats::qqline(res_cat, lty = 2, col = "grey")
              }
              if (id.n > 0) {
                o <- order(res_cat)
                qq_x <- qnorm(ppoints(length(res_cat)))
                y <- res_cat[o]

                extreme <- order(abs(res_cat), decreasing = TRUE)[seq_len(id.n)]
                r <- match(extreme, o)
                side <- ifelse(qq_x[r] < 0, 4, 2) # left tail -> label right; right tail -> label left
                text(
                  qq_x[r],
                  y[r],
                  labels = extreme,
                  pos = side,
                  offset = .35,
                  cex = 0.7
                )
              }
            }
          }
        } else {
          if (w == "residuals-vs-fitted") {
            plot(
              fitted_vals,
              resids,
              xlab = "Fitted values",
              ylab = paste0(tools::toTitleCase(type), " residuals"),
              main = paste0(
                "Residuals vs Fitted\n(",
                method,
                ", ",
                precomputed[[method]]$dose.col.plot,
                ")"
              ),
              ...
            )
            abline(h = 0, lty = 2, col = "grey")
            if (add.smooth) {
              panel.smooth(fitted_vals, resids, col.smooth = "red")
            }
            if (id.n > 0) {
              extreme <- order(abs(resids), decreasing = TRUE)[seq_len(id.n)]
              side <- ifelse(
                fitted_vals[extreme] <
                  (min(fitted_vals) + max(fitted_vals)) / 2,
                4,
                2
              ) # left tail -> label right; right tail -> label left
              text(
                fitted_vals[extreme],
                resids[extreme],
                labels = extreme,
                pos = side,
                offset = .35,
                cex = 0.7
              )
            }
          }

          if (w == "qq") {
            qqnorm(
              resids,
              main = paste0(
                "Normal Q-Q\n(",
                method,
                ", ",
                precomputed[[method]]$dose.col.plot,
                ")"
              ),
              ylab = paste0(tools::toTitleCase(type), " residuals"),
              ...
            )
            if (qqline) {
              stats::qqline(resids, lty = 2, col = "grey")
            }
            if (id.n > 0) {
              o <- order(resids)
              qq_x <- qnorm(ppoints(length(resids)))
              y <- resids[o]

              extreme <- order(abs(resids), decreasing = TRUE)[seq_len(id.n)]
              r <- match(extreme, o)

              side <- ifelse(qq_x[r] < 0, 4, 2) # left tail -> label right; right tail -> label left

              text(
                qq_x[r],
                y[r],
                labels = extreme,
                pos = side,
                cex = 0.7,
                offset = .35
              )
            }
          }
        }
      }
    }
  } else {
    # for prophaz, resids is a data frame with columns id, time, and all covariates.
    covariate_cols <- setdiff(
      colnames(precomputed[[available[1]]]$resids),
      c("id", "time")
    )
    for (mycol in covariate_cols) {
      for (method in available) {
        resids <- precomputed[[method]]$resids

        o <- order(resids$time)
        x_sf <- resids$time[o]
        y <- resids[o, mycol]
        plot(
          x_sf,
          y,
          pch = 1,
          cex = 1,
          xlab = "Time",
          ylab = "Scaled Schoenfeld residuals",
          main = paste0(
            mycol,
            "\n(",
            method,
            ", ",
            precomputed[[method]]$dose.col.plot,
            ")"
          ),
          ...
        )

        # smooth trend, similar in spirit to cox.zph
        ok <- is.finite(x_sf) & is.finite(y)
        if (sum(ok) > 3) {
          sm <- smooth.spline(x_sf[ok], y[ok])
          lines(sm, lwd = 2, col = "black")
        }
      }
    }
  }
  invisible(x)
}
