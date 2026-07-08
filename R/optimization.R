# Internal optimizer wrapper used by FMA/BMA, MCML, and RC/ERC. It normalizes
# the one-dimensional optimize() path and the general optim() path to the same
# fit object shape, then attaches a numeric Hessian for downstream checks.
fit_objective_with_hessian <- function(
  start,
  fn,
  optim.method = "Nelder-Mead",
  control = list(reltol = 1e-10),
  lower = -20,
  upper = 5,
  use_optimize = length(start) == 1,
  gradient.check = !use_optimize,
  ...
) {
  if (use_optimize) {
    fit0 <- optimize(
      f = fn,
      lower = lower,
      upper = upper,
      ...
    )
    fit <- list(
      par = fit0$minimum,
      value = fit0$objective,
      convergence = 0
    )
  } else {
    fit <- optim(
      start,
      fn,
      method = optim.method,
      control = control,
      ...
    )
    if (optim.method == "Nelder-Mead") {
      count0 <- fit$counts
      fit <- optim(
        fit$par,
        fn,
        method = "BFGS",
        control = control,
        ...
      )
      fit$counts <- replace(count0, is.na(count0), 0) +
        replace(fit$counts, is.na(fit$counts), 0)
    }
  }

  fit$hessian <- numDeriv::hessian(
    func = fn,
    x = fit$par,
    ...
  )
  if (isTRUE(gradient.check)) {
    fit <- add_optimizer_gradient_diagnostics(fit, fn, ...)
  }
  fit
}


add_optimizer_gradient_diagnostics <- function(fit, fn, ...) {
  fit$gradient <- tryCatch(
    numDeriv::grad(
      func = fn,
      x = fit$par,
      ...
    ),
    error = function(e) {
      fit$gradient.error <- conditionMessage(e)
      NULL
    }
  )

  if (!is.null(fit$gradient) && all(is.finite(fit$gradient))) {
    fit$gradient.max.abs <- max(abs(fit$gradient))
    fit$gradient.rms <- sqrt(mean(fit$gradient^2))
    fit$gradient.scaled.rms <- fit$gradient.rms / max(1, abs(fit$value))
  }

  fit
}


optimizer_gradient_is_large <- function(
  fit,
  abs.tol = 1e-3,
  rel.tol = 1e-4
) {
  isTRUE(fit$convergence == 0) &&
    is.numeric(fit$gradient.rms) &&
    is.numeric(fit$gradient.scaled.rms) &&
    is.finite(fit$gradient.rms) &&
    is.finite(fit$gradient.scaled.rms) &&
    fit$gradient.rms > abs.tol &&
    fit$gradient.scaled.rms > rel.tol
}


warn_if_large_optimizer_gradient <- function(fit) {
  if (!optimizer_gradient_is_large(fit)) {
    return(invisible(NULL))
  }

  warning(
    paste0(
      "WARNING: optim() reported convergence, but the numerical gradient ",
      "at the solution is still large (RMS gradient = ",
      signif(fit$gradient.rms, 4),
      ", scaled RMS gradient = ",
      signif(fit$gradient.scaled.rms, 4),
      ", max absolute gradient = ",
      signif(fit$gradient.max.abs, 4),
      "). The fit may be sensitive to covariate scaling or starting values; ",
      "consider centering/scaling continuous covariates, supplying inpar, ",
      "or adjusting optimizer controls."
    ),
    call. = FALSE
  )

  invisible(NULL)
}


compute_fit_gradient_diagnostics <- function(
  object,
  method_name,
  method_fit,
  data
) {
  fit <- list(
    par = method_fit$optim$par,
    value = -1 * method_fit$loglik,
    convergence = method_fit$optim$convergence
  )
  loglik_fn <- make_loglik_fn(object, method_name, method_fit, data)
  add_optimizer_gradient_diagnostics(fit, loglik_fn)
}


stored_gradient_diagnostics <- function(method_fit) {
  optim <- method_fit$optim

  if (
    is.null(optim$gradient.rms) ||
      is.null(optim$gradient.scaled.rms) ||
      is.null(optim$gradient.max.abs)
  ) {
    return(NULL)
  }

  list(
    convergence = optim$convergence,
    gradient = optim$gradient,
    gradient.max.abs = optim$gradient.max.abs,
    gradient.rms = optim$gradient.rms,
    gradient.scaled.rms = optim$gradient.scaled.rms
  )
}


gradient_diagnostic_row <- function(method_name, fit) {
  gradient_large <- optimizer_gradient_is_large(fit)
  data.frame(
    method = method_name,
    convergence = fit$convergence %||% NA_integer_,
    gradient.rms = fit$gradient.rms %||% NA_real_,
    gradient.scaled.rms = fit$gradient.scaled.rms %||% NA_real_,
    gradient.max.abs = fit$gradient.max.abs %||% NA_real_,
    gradient.large = gradient_large,
    row.names = method_name
  )
}


is_finite_square_matrix <- function(x) {
  !is.null(x) &&
    is.numeric(x) &&
    length(dim(x)) == 2 &&
    nrow(x) == ncol(x) &&
    nrow(x) > 0 &&
    all(is.finite(x))
}

# Convergence-aware Hessian screen used when deciding whether a realization is
# admissible before model averaging.
fit_passes_hessian_check <- function(fit) {
  if (is.null(fit$hessian) || fit$convergence != 0) {
    return(FALSE)
  }

  if (!is_finite_square_matrix(fit$hessian)) {
    return(FALSE)
  }

  det(fit$hessian) != 0 &&
    rcond(fit$hessian) > .Machine$double.eps &&
    all(eigen(fit$hessian)$values > 0)
}

# Variance extraction for MCML and RC/ERC only needs to know whether the Hessian
# is usable. This deliberately does not check optimizer convergence, preserving
# the pre-refactor result assembly behavior for those methods.
hessian_supports_vcov <- function(hessian) {
  if (!is_finite_square_matrix(hessian)) {
    return(FALSE)
  }

  det(hessian) != 0 &&
    rcond(hessian) > .Machine$double.eps &&
    all(eigen(hessian)$values > 0)
}

# Shared post-fit assembly for MCML and RC/ERC. The returned list preserves the
# public result shape of those methods; method-specific fields such as ERC are
# appended through `extra`.
assemble_frequentist_fit_result <- function(
  fit,
  parnames,
  t0,
  transform = NULL,
  transform.jacobian = NULL,
  extra = list(),
  ...
) {
  # MCML and RC/ERC share this post-optimization path: optionally transform
  # coefficients, propagate variance through the Jacobian, then package results.
  if (!is.null(transform) & !is.null(transform.jacobian)) {
    if (is.function(transform) & is.function(transform.jacobian)) {
      if ("boundcheck" %in% names(formals(transform))) {
        coefs <- transform(fit$par, boundcheck = TRUE, ...)
      } else {
        coefs <- transform(fit$par, ...)
      }

      if (hessian_supports_vcov(fit$hessian)) {
        cholH <- tryCatch(chol(fit$hessian), error = function(e) NULL)
        if (!is.null(cholH)) {
          jac <- transform.jacobian(fit$par, ...)
          tmpsolve <- backsolve(cholH, t(jac), transpose = TRUE)
          vcov <- crossprod(tmpsolve)
          #vcov <- jac %*% MASS::ginv(fit$hessian) %*% t(jac)
        } else {
          warning(
            "WARNING: Hessian was not invertible or inverse was not positive definite, variance matrix could not be obtained"
          )
          vcov <- matrix(NA, ncol = length(parnames), nrow = length(parnames))
        }
      } else {
        warning(
          "WARNING: Hessian was not invertible or inverse was not positive definite, variance matrix could not be obtained"
        )
        vcov <- matrix(NA, ncol = length(parnames), nrow = length(parnames))
      }
    } else {
      stop("transform and transform.jacobian should be functions")
    }
  } else {
    coefs <- fit$par
    if (hessian_supports_vcov(fit$hessian)) {
      vcov <- chol2inv(chol(fit$hessian)) #MASS::ginv(fit$hessian)
    } else {
      warning(
        "WARNING: Hessian was not invertible or inverse was not positive definite, variance matrix could not be obtained"
      )
      vcov <- matrix(NA, ncol = length(parnames), nrow = length(parnames))
    }
  }

  warn_if_large_optimizer_gradient(fit)

  names(coefs) <- parnames
  rownames(vcov) <- colnames(vcov) <- parnames

  fit_timing <- stop_runtime_timer(t0)
  timing <- new_method_timing(fit = fit_timing)

  c(
    list(
      coefficients = coefs,
      sd = sqrt(diag(vcov)),
      vcov = vcov,
      optim = list(
        par = fit$par,
        hessian = fit$hessian,
        gradient = fit$gradient,
        gradient.max.abs = fit$gradient.max.abs,
        gradient.rms = fit$gradient.rms,
        gradient.scaled.rms = fit$gradient.scaled.rms,
        convergence = fit$convergence,
        counts = fit$counts
      ),
      loglik = -1 * fit$value,
      timing = timing,
      runtime = format_runtime(timing$total$cpu)
    ),
    extra
  )
}
