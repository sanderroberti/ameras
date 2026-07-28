# Refine a Nelder-Mead solution with BFGS when its numerical finite differences
# are usable. A failed refinement should not discard an otherwise valid
# Nelder-Mead fit.
refine_nelder_mead_with_bfgs <- function(
  fit,
  fn,
  control,
  warn = TRUE,
  ...
) {
  count0 <- fit$counts
  refined <- tryCatch(
    optim(
      fit$par,
      fn,
      method = "BFGS",
      control = control,
      ...
    ),
    error = identity
  )

  if (inherits(refined, "error")) {
    fit$bfgs.refinement.error <- conditionMessage(refined)
    if (isTRUE(warn)) {
      warning(
        "Automatic BFGS refinement failed (",
        fit$bfgs.refinement.error,
        "); retaining the Nelder-Mead solution. Review the optimizer ",
        "convergence diagnostics and consider rescaling covariates or ",
        "supplying different starting values.",
        call. = FALSE
      )
    }
    return(fit)
  }

  refined$counts <- replace(count0, is.na(count0), 0) +
    replace(refined$counts, is.na(refined$counts), 0)
  refined
}


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
  compute.hessian = TRUE,
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
      fit <- refine_nelder_mead_with_bfgs(
        fit = fit,
        fn = fn,
        control = control,
        warn = gradient.check,
        ...
      )
    }
  }

  if (isTRUE(compute.hessian)) {
    fit$hessian <- numDeriv::hessian(
      func = fn,
      x = fit$par,
      ...
    )
  } else {
    fit$hessian <- NULL
  }
  if (isTRUE(gradient.check) && isTRUE(compute.hessian)) {
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
    fit$gradient.rms.scaled <- fit$gradient.rms / max(1, abs(fit$value))

    if (
      is_finite_square_matrix(fit$hessian) &&
        length(fit$gradient) == ncol(fit$hessian)
    ) {
      cholH <- tryCatch(chol(fit$hessian), error = function(e) NULL)
      if (!is.null(cholH)) {
        # For a quadratic objective, half the squared Newton decrement is the
        # approximate objective improvement still available from one Newton
        # step. This is much less sensitive to likelihood scale and parameter
        # units than the raw gradient norm.
        z <- backsolve(cholH, fit$gradient, transpose = TRUE)
        newton_step <- backsolve(cholH, z)
        fit$newton.decrement <- sqrt(sum(z^2))
        fit$newton.improvement <- 0.5 * sum(z^2)
        fit$newton.improvement.relative <-
          fit$newton.improvement / max(1, abs(fit$value))
        fit$newton.step.max.abs <- max(abs(newton_step))
      }
    }
  }

  fit
}


optimizer_gradient_is_large <- function(
  fit,
  newton.improvement.tol = 1e-2,
  newton.relative.improvement.tol = 1e-6,
  abs.tol = 1e-3,
  rel.tol = 1e-4
) {
  if (!isTRUE(fit$convergence == 0)) {
    return(FALSE)
  }

  if (
    is.numeric(fit$newton.improvement) &&
      is.finite(fit$newton.improvement) &&
      is.numeric(fit$newton.improvement.relative) &&
      is.finite(fit$newton.improvement.relative)
  ) {
    return(
      fit$newton.improvement > newton.improvement.tol &&
        fit$newton.improvement.relative >
          newton.relative.improvement.tol
    )
  }

  is.numeric(fit$gradient.rms) &&
    is.numeric(fit$gradient.rms.scaled) &&
    is.finite(fit$gradient.rms) &&
    is.finite(fit$gradient.rms.scaled) &&
    fit$gradient.rms > abs.tol &&
    fit$gradient.rms.scaled > rel.tol
}


warn_if_large_optimizer_gradient <- function(fit) {
  if (!optimizer_gradient_is_large(fit)) {
    return(invisible(NULL))
  }

  warning(
    paste0(
      "WARNING: optim() reported convergence, but optimizer diagnostics ",
      "suggest the solution may not be fully stationary ",
      "(RMS gradient = ",
      signif(fit$gradient.rms, 4),
      ", scaled RMS gradient = ",
      signif(fit$gradient.rms.scaled, 4),
      ", max absolute gradient = ",
      signif(fit$gradient.max.abs %||% NA_real_, 4),
      if (
        is.numeric(fit$newton.decrement) &&
          is.finite(fit$newton.decrement) &&
          is.numeric(fit$newton.improvement) &&
          is.finite(fit$newton.improvement) &&
          is.numeric(fit$newton.improvement.relative) &&
          is.finite(fit$newton.improvement.relative)
      ) {
        paste0(
          ", Newton decrement = ",
          signif(fit$newton.decrement, 4),
          ", approximate objective improvement = ",
          signif(fit$newton.improvement, 4),
          ", relative improvement = ",
          signif(fit$newton.improvement.relative, 4)
        )
      } else {
        ""
      },
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
    convergence = method_fit$optim$convergence,
    hessian = method_fit$optim$hessian
  )
  loglik_fn <- make_loglik_fn(object, method_name, method_fit, data)
  add_optimizer_gradient_diagnostics(fit, loglik_fn)
}


stored_gradient_diagnostics <- function(method_fit) {
  optim <- method_fit$optim

  if (
    is.null(optim$gradient.rms) ||
      is.null(optim$gradient.rms.scaled)
  ) {
    return(NULL)
  }

  out <- list(
    convergence = optim$convergence,
    gradient = optim$gradient,
    gradient.max.abs = optim$gradient.max.abs,
    gradient.rms = optim$gradient.rms,
    gradient.rms.scaled = optim$gradient.rms.scaled,
    newton.decrement = optim$newton.decrement,
    newton.improvement = optim$newton.improvement,
    newton.improvement.relative = optim$newton.improvement.relative,
    newton.step.max.abs = optim$newton.step.max.abs
  )

  if (
    is.null(out$newton.improvement) &&
      !is.null(out$gradient) &&
      is_finite_square_matrix(optim$hessian) &&
      length(out$gradient) == ncol(optim$hessian)
  ) {
    cholH <- tryCatch(chol(optim$hessian), error = function(e) NULL)
    if (!is.null(cholH)) {
      z <- backsolve(cholH, out$gradient, transpose = TRUE)
      newton_step <- backsolve(cholH, z)
      out$newton.decrement <- sqrt(sum(z^2))
      out$newton.improvement <- 0.5 * sum(z^2)
      out$newton.improvement.relative <-
        out$newton.improvement / max(1, abs(method_fit$loglik))
      out$newton.step.max.abs <- max(abs(newton_step))
    }
  }

  if (
    is.null(out$newton.improvement.relative) &&
      !is.null(out$newton.improvement)
  ) {
    out$newton.improvement.relative <-
      out$newton.improvement / max(1, abs(method_fit$loglik))
  }

  out
}


gradient_diagnostic_row <- function(method_name, fit) {
  convergence_warning <- optimizer_gradient_is_large(fit)
  data.frame(
    method = method_name,
    optim.convergence = fit$convergence %||% NA_integer_,
    gradient.rms = fit$gradient.rms %||% NA_real_,
    gradient.rms.scaled = fit$gradient.rms.scaled %||% NA_real_,
    newton.improvement = fit$newton.improvement %||% NA_real_,
    newton.improvement.relative = fit$newton.improvement.relative %||% NA_real_,
    convergence.warning = convergence_warning,
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

  hessian_supports_vcov(fit$hessian)
}

# Variance extraction for MCML and RC/ERC only needs to know whether the Hessian
# is usable. This deliberately does not check optimizer convergence, preserving
# the pre-refactor result assembly behavior for those methods.
hessian_supports_vcov <- function(hessian) {
  if (!is_finite_square_matrix(hessian)) {
    return(FALSE)
  }

  det_hessian <- tryCatch(det(hessian), error = function(e) NA_real_)
  rcond_hessian <- tryCatch(rcond(hessian), error = function(e) NA_real_)
  eigen_values <- tryCatch(
    eigen(hessian, symmetric = TRUE, only.values = TRUE)$values,
    error = function(e) NA_real_
  )

  isTRUE(
    is.finite(det_hessian) &&
      det_hessian != 0 &&
      is.finite(rcond_hessian) &&
      rcond_hessian > .Machine$double.eps &&
      all(is.finite(eigen_values)) &&
      all(eigen_values > 0)
  )
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
        gradient.rms.scaled = fit$gradient.rms.scaled,
        newton.decrement = fit$newton.decrement,
        newton.improvement = fit$newton.improvement,
        newton.improvement.relative = fit$newton.improvement.relative,
        newton.step.max.abs = fit$newton.step.max.abs,
        convergence = fit$convergence,
        counts = fit$counts
      ),
      loglik = -1 * fit$value,
      timing = timing
    ),
    extra
  )
}
