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

fit_passes_hessian_check <- function(fit) {
  if (is.null(fit$hessian) || fit$convergence != 0) {
    return(FALSE)
  }

  det(fit$hessian) != 0 &&
    rcond(fit$hessian) > .Machine$double.eps &&
    all(eigen(fit$hessian)$values > 0)
}
