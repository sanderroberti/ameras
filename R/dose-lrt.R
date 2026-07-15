dose_lrt <- function(object, ...) UseMethod("dose_lrt")


dose_parameter_indices <- function(parnames) {
  base_names <- sub("^\\([^)]*\\)_", "", parnames)
  dose_terms <- c("dose", "dose_squared", "dose_linear", "dose_exponential")

  which(
    base_names %in% dose_terms |
      # Subgroup-coded modifier names look like dose[M=level].
      grepl(
        paste0("^(", paste(dose_terms, collapse = "|"), ")\\["),
        base_names
      ) |
      grepl(
        paste0("^(", paste(dose_terms, collapse = "|"), "):"),
        base_names
      )
  )
}


make_reported_scale_loglik_fn <- function(object, method_name, method_fit, data) {
  reported_object <- object
  # The constrained null objective fixes tested dose coefficients after any
  # fitted transformation, so the likelihood itself must consume reported-scale
  # parameters directly.
  reported_object$transform <- NULL
  reported_object$transform.jacobian <- NULL

  make_loglik_fn(reported_object, method_name, method_fit, data)
}


apply_fitted_transform <- function(object, params) {
  if (is.null(object$transform)) {
    return(params)
  }

  other_args <- object$other.args %||% list()
  do.call(object$transform, c(list(params = params), other_args))
}


make_constrained_lrt_objective <- function(
  object,
  method_name,
  method_fit,
  data,
  fixed.idx,
  null.value = 0
) {
  inpar <- method_fit$optim$par
  parnames <- names(method_fit$coefficients)

  if (!length(fixed.idx)) {
    stop("At least one parameter must be fixed for an LRT")
  }
  if (any(!(fixed.idx %in% seq_along(inpar)))) {
    stop("Fixed parameter index out of range")
  }
  if (length(inpar) != length(parnames)) {
    stop("Optimizer parameter vector and coefficient vector have different lengths")
  }

  fixed_values <- rep(null.value, length(fixed.idx))
  free_idx <- setdiff(seq_along(inpar), fixed.idx)
  reported_loglik <- make_reported_scale_loglik_fn(
    object,
    method_name,
    method_fit,
    data
  )

  parameter_fn <- function(free_params) {
    if (length(free_params) != length(free_idx)) {
      stop("Free parameter vector has incorrect length")
    }

    # Free nuisance parameters stay on the optimizer scale until the fitted
    # transformation is applied. Tested dose parameters are then overwritten at
    # the null value on the reported scale, avoiding inverse-transform boundary
    # problems for ERR/LINEXP-style fits.
    full_internal <- inpar
    full_internal[free_idx] <- free_params
    reported_params <- apply_fitted_transform(object, full_internal)
    reported_params[fixed.idx] <- fixed_values
    names(reported_params) <- parnames
    reported_params
  }

  list(
    fn = function(free_params) reported_loglik(parameter_fn(free_params)),
    start = inpar[free_idx],
    free.idx = free_idx,
    fixed.idx = fixed.idx,
    parameter_fn = parameter_fn
  )
}


fit_constrained_lrt_objective <- function(
  objective,
  optim.method = "Nelder-Mead",
  control = list(reltol = 1e-10)
) {
  if (!length(objective$start)) {
    return(list(
      par = numeric(0),
      value = objective$fn(numeric(0)),
      convergence = 0L,
      counts = NULL,
      hessian = NULL
    ))
  }

  fit_objective_with_hessian(
    start = objective$start,
    fn = objective$fn,
    optim.method = optim.method,
    control = control,
    use_optimize = FALSE,
    gradient.check = FALSE,
    compute.hessian = FALSE
  )
}


dose_lrt_row <- function(
  method_name,
  type,
  term,
  df,
  full_loglik,
  null_fit,
  negative.tol = 1e-6
) {
  null_loglik <- -1 * null_fit$value
  statistic <- 2 * (full_loglik - null_loglik)
  null_convergence <- null_fit$convergence
  if (is.null(null_convergence) || !length(null_convergence)) {
    null_convergence <- NA_integer_
  }
  null_convergence <- null_convergence[1]

  if (!is.finite(statistic)) {
    warning("Dose LRT statistic for ", method_name, " / ", term, " is not finite")
    statistic <- NA_real_
  } else if (statistic < 0) {
    if (statistic < -negative.tol) {
      warning(
        "Dose LRT statistic for ",
        method_name,
        " / ",
        term,
        " is negative (",
        signif(statistic, 4),
        "), likely due to numerical optimization differences; reporting 0."
      )
    }
    statistic <- 0
  }
  if (!is.na(null_convergence) && null_convergence != 0) {
    warning(
      "Constrained null fit for ",
      method_name,
      " / ",
      term,
      " returned optimizer convergence code ",
      null_convergence,
      "."
    )
  }

  data.frame(
    method = method_name,
    type = type,
    term = term,
    df = df,
    logLik = full_loglik,
    logLik.null = null_loglik,
    statistic = statistic,
    p.value = pchisq(statistic, df = df, lower.tail = FALSE),
    null.optim.convergence = null_convergence,
    row.names = NULL
  )
}


compute_dose_lrt_for_method <- function(
  object,
  method_name,
  data,
  type = c("global", "individual")
) {
  method_fit <- object[[method_name]]
  if (is.null(method_fit)) {
    stop("Method not present in fitted object: ", method_name)
  }

  parnames <- names(method_fit$coefficients)
  dose_idx <- dose_parameter_indices(parnames)
  if (!length(dose_idx)) {
    stop("No dose-related parameters found for method: ", method_name)
  }

  optim_method <- object$model$optim.method %||% "Nelder-Mead"
  control <- object$model$control %||% list(reltol = 1e-10)

  fit_fixed <- function(fixed_idx) {
    objective <- make_constrained_lrt_objective(
      object = object,
      method_name = method_name,
      method_fit = method_fit,
      data = data,
      fixed.idx = fixed_idx,
      null.value = 0
    )
    fit_constrained_lrt_objective(
      objective,
      optim.method = optim_method,
      control = control
    )
  }

  rows <- list()
  if ("global" %in% type) {
    null_fit <- fit_fixed(dose_idx)
    rows <- c(rows, list(dose_lrt_row(
      method_name = method_name,
      type = "global",
      term = "all dose terms",
      df = length(dose_idx),
      full_loglik = method_fit$loglik,
      null_fit = null_fit
    )))
  }

  if ("individual" %in% type) {
    individual_rows <- lapply(dose_idx, function(idx) {
      null_fit <- fit_fixed(idx)
      dose_lrt_row(
        method_name = method_name,
        type = "individual",
        term = parnames[idx],
        df = 1L,
        full_loglik = method_fit$loglik,
        null_fit = null_fit
      )
    })
    rows <- c(rows, individual_rows)
  }

  do.call(rbind, rows)
}


dose_lrt.amerasfit <- function(
  object,
  methods = c("RC", "ERC", "MCML"),
  type = c("global", "individual"),
  data = NULL,
  ...
) {
  supported_methods <- c("RC", "ERC", "MCML")

  if (!is.character(methods) || !length(methods)) {
    stop("methods must be a non-empty character vector")
  }
  methods <- unique(methods)

  unsupported <- setdiff(methods, supported_methods)
  if (length(unsupported)) {
    stop(
      "dose_lrt() supports RC, ERC, and MCML only. Unsupported method(s): ",
      paste(unsupported, collapse = ", ")
    )
  }

  type <- match.arg(type, choices = c("global", "individual"), several.ok = TRUE)
  type <- unique(type)

  available <- methods[vapply(methods, function(m) !is.null(object[[m]]), logical(1))]
  if (!length(available)) {
    stop(
      "None of the requested supported methods are present in the fitted object"
    )
  }

  resolved_data <- resolve_data(object, data)
  out <- do.call(
    rbind,
    lapply(available, function(method_name) {
      compute_dose_lrt_for_method(
        object = object,
        method_name = method_name,
        data = resolved_data,
        type = type
      )
    })
  )
  row.names(out) <- NULL
  out
}
