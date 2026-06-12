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
  fit
}

# Convergence-aware Hessian screen used when deciding whether a realization is
# admissible before model averaging.
fit_passes_hessian_check <- function(fit) {
  if (is.null(fit$hessian) || fit$convergence != 0) {
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
