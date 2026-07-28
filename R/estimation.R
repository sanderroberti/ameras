compute_sample_CI <- function(samples, type = "percentile", level = .95) {
  samples <- as.matrix(samples)
  if (type == "percentile") {
    CIlower = apply(samples, 2, function(x) quantile(x, (1 - level) / 2))
    CIupper = apply(samples, 2, function(x) quantile(x, (1 + level) / 2))
    CIresult <- data.frame(lower = CIlower, upper = CIupper)
  } else if (type == "hpd") {
    CIresult <- as.data.frame(
      coda::HPDinterval(coda::as.mcmc(samples)),
      prob = level
    )
  }
  return(CIresult)
}

compute_wald_CI <- function(
  method_fit,
  level = .95,
  transform = NULL,
  other.args = NULL,
  type = "orig"
) {
  if (type == "transformed" & is.null(transform)) {
    stop(
      "No transformation specified, specify transformation or choose a different CI type"
    )
  }

  z <- qnorm(1 - (1 - level) / 2)
  hessian <- method_fit$optim$hessian

  if (type == "transformed") {
    # Check invertibility
    invertible <- hessian_supports_vcov(hessian)

    if (!invertible) {
      warning(
        "Hessian not invertible, confidence intervals could not be obtained"
      )
      na_vec <- NA * method_fit$coefficients
      return(data.frame(lower = na_vec, upper = na_vec))
    }

    hess_inv <- chol2inv(chol(hessian))

    CIlower <- do.call(
      transform,
      c(
        list(params = method_fit$optim$par - z * sqrt(diag(hess_inv))),
        other.args
      )
    )
    CIupper <- do.call(
      transform,
      c(
        list(params = method_fit$optim$par + z * sqrt(diag(hess_inv))),
        other.args
      )
    )
  } else if (type == "orig") {
    CIlower <- method_fit$coefficients - z * sqrt(diag(method_fit$vcov))
    CIupper <- method_fit$coefficients + z * sqrt(diag(method_fit$vcov))
  }

  CIresult <- data.frame(
    lower = CIlower,
    upper = CIupper,
    row.names = names(method_fit$coefficients)
  )
}

compute_proflik_CI <- function(
  method_fit,
  object,
  method_name,
  data,
  parm = "dose",
  level = 0.95,
  maxit.profCI = 20,
  tol.profCI = 1e-2,
  optim.method = "Nelder-Mead",
  control = list(reltol = 1e-10)
) {
  alpha <- 1 - level
  optval <- -method_fit$loglik # stored as positive log-likelihood
  inpar <- method_fit$optim$par
  parnames <- names(method_fit$coefficients)

  # Determine which parameters to compute CIs for
  if (identical(parm, "dose")) {
    CIindices <- which(
      startsWith(parnames, "dose") |
        grepl(")_dose", parnames)
    )
  } else if (identical(parm, "all")) {
    CIindices <- seq_along(parnames)
  } else {
    CIindices <- which(parnames %in% parm)
    if (!length(CIindices)) {
      stop(
        "No parameters matched parm. Available parameters are: ",
        paste(parnames, collapse = ", ")
      )
    }
  }

  # Reconstruct the likelihood function from stored model
  loglik_fn <- make_loglik_fn(object, method_name, method_fit, data)

  # Use stored hessian to get sensible search bounds
  hessian <- method_fit$optim$hessian
  invertible <- hessian_supports_vcov(hessian)

  if (invertible) {
    se <- sqrt(diag(chol2inv(chol(hessian))))
    lowlims <- inpar - 4 * se
    uplims <- inpar + 4 * se
  } else {
    warning(
      "Hessian not invertible, using default search bounds for profile ",
      "likelihood CI. Results may be unreliable."
    )
    lowlims <- rep(-20, length(inpar))
    uplims <- rep(10, length(inpar))
  }

  # Initialise output vectors
  CIlower <- rep(NA_real_, length(parnames))
  CIupper <- rep(NA_real_, length(parnames))
  pval_lower <- rep(NA_real_, length(parnames))
  pval_upper <- rep(NA_real_, length(parnames))
  iter_lower <- rep(NA_integer_, length(parnames))
  iter_upper <- rep(NA_integer_, length(parnames))

  for (myindex in seq_along(parnames)) {
    if (!(myindex %in% CIindices)) {
      next
    }

    message("Obtaining profile likelihood CI for ", parnames[myindex])

    profCI_one <- compute_proflik_ci_one(
      index = myindex,
      inpar = inpar,
      optval = optval,
      loglik_fn = loglik_fn,
      lowlim = lowlims[myindex],
      uplim = uplims[myindex],
      alpha = alpha,
      maxit.profCI = maxit.profCI,
      tol.profCI = tol.profCI,
      optim.method = optim.method,
      control = control,
      parname = parnames[myindex],
      transform = object$transform,
      other.args = object$other.args
    )

    CIlower[myindex] <- profCI_one$lower
    CIupper[myindex] <- profCI_one$upper
    pval_lower[myindex] <- profCI_one$pval.lower
    pval_upper[myindex] <- profCI_one$pval.upper
    iter_lower[myindex] <- profCI_one$iter.lower
    iter_upper[myindex] <- profCI_one$iter.upper
  }

  # Back-transform if needed
  if (!is.null(object$transform)) {
    CIlower <- do.call(
      object$transform,
      c(list(params = CIlower), object$other.args)
    )
    CIupper <- do.call(
      object$transform,
      c(list(params = CIupper), object$other.args)
    )
  }

  result <- data.frame(
    lower = CIlower,
    upper = CIupper,
    pval.lower = pval_lower,
    pval.upper = pval_upper,
    iter.lower = iter_lower,
    iter.upper = iter_upper,
    row.names = parnames
  )

  result[CIindices, , drop = FALSE]
}


format_profile_boundary <- function(x, digits = 4) {
  if (is.na(x) || !is.finite(x)) {
    return(as.character(x))
  }

  format(signif(x, digits), trim = TRUE, scientific = FALSE)
}


make_profile_target_fn <- function(
  index,
  inpar,
  optval,
  loglik_fn,
  alpha,
  optim.method,
  control
) {
  cache <- new.env(parent = emptyenv())

  function(par) {
    vapply(par, function(mypar) {
      key <- format(mypar, digits = 17, scientific = TRUE)
      if (exists(key, envir = cache, inherits = FALSE)) {
        return(get(key, envir = cache, inherits = FALSE))
      }

      value <- tryCatch(
        1 -
          pchisq(
            2 *
              (proflik(
                parvalue = mypar,
                index = index,
                fun = loglik_fn,
                inpar = inpar,
                optim.method = optim.method,
                control = control
              ) -
                optval),
            df = 1
          ) -
          alpha,
        error = function(e) NA_real_
      )
      assign(key, value, envir = cache)
      value
    }, numeric(1))
  }
}


profile_reported_value <- function(
  par,
  index,
  inpar,
  transform = NULL,
  other.args = NULL
) {
  if (is.null(transform)) {
    return(par)
  }

  params <- inpar
  params[index] <- par
  tryCatch(
    do.call(transform, c(list(params = params), other.args))[index],
    error = function(e) NA_real_
  )
}


warn_profile_bound_failure <- function(
  parname,
  direction,
  reason,
  par,
  target,
  alpha,
  index,
  inpar,
  transform = NULL,
  other.args = NULL
) {
  pval <- target + alpha
  reported <- profile_reported_value(
    par = par,
    index = index,
    inpar = inpar,
    transform = transform,
    other.args = other.args
  )

  warning(
    paste0(
      "Profile likelihood ",
      direction,
      " bound for ",
      parname,
      " ",
      reason,
      "; returning NA for this bound. Last attempted value: ",
      format_profile_boundary(reported),
      if (!is.null(transform)) {
        paste0(
          " (internal optimizer-scale value: ",
          format_profile_boundary(par),
          ")"
        )
      } else {
        ""
      },
      ", p-value at last attempted value: ",
      format_profile_boundary(pval),
      ", target p-value: ",
      format_profile_boundary(alpha),
      "."
    ),
    call. = FALSE
  )
}


find_profile_bracket <- function(
  target_fn,
  center,
  initial,
  direction,
  max_expand = 12L
) {
  center_value <- target_fn(center)
  if (!is.finite(center_value) || center_value <= 0) {
    return(list(
      found = FALSE,
      lower = NA_real_,
      upper = NA_real_,
      last.par = center,
      last.target = center_value
    ))
  }

  if (!is.finite(initial) || direction * (initial - center) <= 0) {
    initial <- center + direction
  }

  step <- abs(initial - center)
  if (!is.finite(step) || step <= sqrt(.Machine$double.eps)) {
    step <- 1
  }

  candidate <- initial
  last_candidate <- candidate
  last_value <- NA_real_
  for (i in 0:max_expand) {
    if (!is.finite(candidate)) {
      break
    }

    last_candidate <- candidate
    last_value <- target_fn(candidate)
    if (is.finite(last_value) && last_value <= 0) {
      if (direction < 0) {
        return(list(
          found = TRUE,
          lower = candidate,
          upper = center,
          last.par = candidate,
          last.target = last_value
        ))
      }

      return(list(
        found = TRUE,
        lower = center,
        upper = candidate,
        last.par = candidate,
        last.target = last_value
      ))
    }

    step <- step * 2
    candidate <- center + direction * step
  }

  list(
    found = FALSE,
    lower = NA_real_,
    upper = NA_real_,
    last.par = last_candidate,
    last.target = last_value
  )
}


solve_profile_bound <- function(
  target_fn,
  center,
  initial,
  direction,
  direction_name,
  parname,
  alpha,
  maxit.profCI,
  tol.profCI,
  index,
  inpar,
  transform = NULL,
  other.args = NULL,
  target.tol = 0.005
) {
  bracket <- find_profile_bracket(
    target_fn = target_fn,
    center = center,
    initial = initial,
    direction = direction
  )

  if (!isTRUE(bracket$found)) {
    warn_profile_bound_failure(
      parname = parname,
      direction = direction_name,
      reason = "could not be bracketed",
      par = bracket$last.par,
      target = bracket$last.target,
      alpha = alpha,
      index = index,
      inpar = inpar,
      transform = transform,
      other.args = other.args
    )
    return(list(
      root = NA_real_,
      f.root = bracket$last.target,
      iter = NA_integer_
    ))
  }

  root_warning <- NULL
  root_error <- NULL
  root <- tryCatch(
    withCallingHandlers(
      uniroot(
        target_fn,
        lower = bracket$lower,
        upper = bracket$upper,
        maxiter = maxit.profCI,
        tol = tol.profCI
      ),
      warning = function(w) {
        root_warning <<- conditionMessage(w)
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      root_error <<- conditionMessage(e)
      NULL
    }
  )

  if (is.null(root)) {
    warn_profile_bound_failure(
      parname = parname,
      direction = direction_name,
      reason = paste0("could not be solved", if (!is.null(root_error)) {
        paste0(" (", root_error, ")")
      } else {
        ""
      }),
      par = bracket$last.par,
      target = bracket$last.target,
      alpha = alpha,
      index = index,
      inpar = inpar,
      transform = transform,
      other.args = other.args
    )
    return(list(root = NA_real_, f.root = bracket$last.target, iter = NA_integer_))
  }

  if (
    !is.null(root_warning) ||
      is.na(root$f.root) ||
      abs(root$f.root) > target.tol
  ) {
    warn_profile_bound_failure(
      parname = parname,
      direction = direction_name,
      reason = "was not solved accurately",
      par = root$root,
      target = root$f.root,
      alpha = alpha,
      index = index,
      inpar = inpar,
      transform = transform,
      other.args = other.args
    )
    return(list(root = NA_real_, f.root = root$f.root, iter = root$iter))
  }

  root
}


compute_proflik_ci_one <- function(
  index,
  inpar,
  optval,
  loglik_fn,
  lowlim,
  uplim,
  alpha = 0.05,
  maxit.profCI = 20,
  tol.profCI = 1e-2,
  optim.method = "Nelder-Mead",
  control = list(reltol = 1e-10),
  parname = "parameter",
  transform = NULL,
  other.args = NULL
) {
  profCI_fn <- make_profile_target_fn(
    index = index,
    inpar = inpar,
    optval = optval,
    loglik_fn = loglik_fn,
    alpha = alpha,
    optim.method = optim.method,
    control = control
  )

  lowroot <- solve_profile_bound(
    target_fn = profCI_fn,
    center = inpar[index],
    initial = lowlim,
    direction = -1,
    direction_name = "lower",
    parname = parname,
    alpha = alpha,
    maxit.profCI = maxit.profCI,
    tol.profCI = tol.profCI,
    index = index,
    inpar = inpar,
    transform = transform,
    other.args = other.args
  )

  uproot <- solve_profile_bound(
    target_fn = profCI_fn,
    center = inpar[index],
    initial = uplim,
    direction = 1,
    direction_name = "upper",
    parname = parname,
    alpha = alpha,
    maxit.profCI = maxit.profCI,
    tol.profCI = tol.profCI,
    index = index,
    inpar = inpar,
    transform = transform,
    other.args = other.args
  )

  list(
    lower = lowroot$root,
    upper = uproot$root,
    pval.lower = lowroot$f.root + alpha,
    pval.upper = uproot$f.root + alpha,
    iter.lower = lowroot$iter,
    iter.upper = uproot$iter
  )
}

proflik <- function(
  parvalue,
  index,
  fun,
  inpar,
  optim.method = optim.method,
  control = list(reltol = 1e-10),
  ...
) {
  if (length(inpar) == 1) {
    # 1-parameter model -> just return likelihood itself

    return(fun(params = parvalue, ...))
  } else if (length(inpar) == 2) {
    # 2-parameter model -> 1 parameter optimization for profile likelihood -> use optimize

    innerfun <- function(params, ...) {
      if (index > 1) {
        params <- c(params[1:(index - 1)], parvalue, params[-(1:(index - 1))])
      } else {
        params <- c(parvalue, params)
      }

      fun(params = params, ...)
    }

    innerfit <- optimize(innerfun, lower = -10, upper = 10, ...)
    return(innerfit$objective)
  } else {
    # 3+ parameter model -> 2+ parameter optimization for profile likelihood -> use optim

    innerfun <- function(params, ...) {
      if (index > 1) {
        params <- c(params[1:(index - 1)], parvalue, params[-(1:(index - 1))])
      } else {
        params <- c(parvalue, params)
      }

      fun(params = params, ...)
    }
    innerfit <- optim(
      par = inpar[-index],
      fn = innerfun,
      method = optim.method,
      control = control,
      ...
    )
    if (optim.method == "Nelder-Mead") {
      innerfit <- refine_nelder_mead_with_bfgs(
        fit = innerfit,
        fn = innerfun,
        control = control,
        warn = FALSE,
        ...
      )
    }
    return(innerfit$value)
  }
}
