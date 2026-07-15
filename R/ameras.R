ameras_main <- function(
  family = "gaussian",
  methods = "RC",
  dosevars,
  data,
  deg = 1,
  doseRRmod = "ERR",
  transform = NULL,
  transform.jacobian = NULL,
  Y = NULL,
  M = NULL,
  modifier_info = NULL,
  X = NULL,
  offset = NULL,
  inpar = NULL,
  entry = NULL,
  exit = NULL,
  status = NULL,
  setnr = NULL,
  unweightedFMA = FALSE,
  loglim = 1e-30,
  MFMA = 100000,
  future.chunk.size.FMA = NULL,
  prophaz.numints.BMA = 10,
  ERRprior.BMA = "doubleexponential",
  nburnin.BMA = 5000,
  niter.BMA = 20000,
  nchains.BMA = 2,
  thin.BMA = 10,
  included.realizations.BMA = 1:length(dosevars),
  optim.method = "Nelder-Mead",
  control = NULL,
  ...
) {
  if (is.null(control)) {
    control <- list(reltol = 1e-10)
  }
  future.chunk.size.FMA <- check_future_chunk_size_FMA(future.chunk.size.FMA)

  if (!is.null(transform) & is.null(transform.jacobian)) {
    stop("transform.jacobian is required when using a transformation")
  }

  # Family in ("gaussian", "poisson", "prophaz", "binomial")
  # gaussian requires Y and can use X and M
  # binomial requires Y and can use X and M
  # poisson requires Y and can use offset, X and M
  # prophaz requires status and exit and can use entry, X and M

  # Method in ("MCML", "RC", "ERC", "FMA", "BMA")
  # If both FMA and BMA are to be run, run FMA first so output is ordered.
  # BMA inclusions are controlled by included.realizations.BMA.
  if ("FMA" %in% methods & "BMA" %in% methods) {
    methods[methods %in% c("FMA", "BMA")] <- c("FMA", "BMA")
  }

  # If both RC and ERC are to be run, run RC first to determine initial values for ERC (not implemented yet)
  if ("RC" %in% methods & "ERC" %in% methods) {
    methods[methods %in% c("RC", "ERC")] <- c("RC", "ERC")
  }

  # CI in ("wald.orig","wald.transformed", "proflik") for method in ("RC", "ERC", "MCML")
  # CI in ("percentile", "hpd") for method in ("FMA", "BMA")

  # Input checks
  if (family != "gaussian" & is.null(doseRRmod)) {
    stop("doseRRmod is required for all families other than Gaussian")
  }

  if (family == "prophaz") {
    if (is.null(exit)) {
      stop("exit is required for family=prophaz")
    }
    if (is.null(status)) stop("status is required for family=prophaz")
  }

  # For the Gaussian family, sigma needs to be >0 and so use standard reparametrization if another is not defined
  if (family == "gaussian" & is.null(transform)) {
    transform <- make_transform(index.t = 2 + length(X) + length(M) * deg + deg)
    transform.jacobian <- make_transform.jacobian(
      index.t = 2 + length(X) + length(M) * deg + deg
    )
  }

  # If the user did not provide a transformation, add default bounds for the
  # dose-response parameters that need them. The helper returns the parameter
  # positions and lower limits used by make_transform().
  if (family != "gaussian") {
    if (doseRRmod == "ERR" & is.null(transform)) {
      if (deg == 1) {
        lwlmt <- -1 / max(data[, dosevars])
      } else {
        lwlmt <- c(0, -1 / max(data[, dosevars]^2))
      }
      transform_settings <- err_default_transform_settings(
        family = family,
        X = X,
        M = M,
        deg = deg,
        modifier_info = modifier_info,
        Y = Y,
        data = data,
        lowlimit = lwlmt
      )
      transform <- make_transform(
        index.t = transform_settings$index.t,
        lowlimit = transform_settings$lowlimit
      )
      transform.jacobian <- make_transform.jacobian(
        index.t = transform_settings$index.t,
        lowlimit = transform_settings$lowlimit
      )
    } else if (doseRRmod == "LINEXP" & is.null(transform)) {
      lwlmt <- 0 # Lower bound of 0 for beta1, no bound for beta2

      transform_settings <- linexp_default_transform_settings(
        family = family,
        X = X,
        M = M,
        deg = deg,
        modifier_info = modifier_info,
        Y = Y,
        data = data,
        lowlimit = lwlmt
      )
      transform <- make_transform(
        index.t = transform_settings$index.t,
        lowlimit = transform_settings$lowlimit
      )
      transform.jacobian <- make_transform.jacobian(
        index.t = transform_settings$index.t,
        lowlimit = transform_settings$lowlimit
      )
    }
  }

  results <- list(
    transform = transform,
    transform.jacobian = transform.jacobian
  )
  for (method in methods) {
    if (method == "MCML") {
      message("Fitting MCML")
      fit <- ameras.mcml(
        family = family,
        dosevars = dosevars,
        data = data,
        deg = deg,
        transform = transform,
        transform.jacobian = transform.jacobian,
        Y = Y,
        M = M,
        modifier_info = modifier_info,
        X = X,
        offset = offset,
        inpar = inpar,
        entry = entry,
        exit = exit,
        status = status,
        setnr = setnr,
        doseRRmod = doseRRmod,
        loglim = loglim,
        control = control,
        optim.method = optim.method,
        ...
      )
      results <- c(results, list(MCML = fit))
    } else if (method == "RC") {
      message("Fitting RC")
      fit <- ameras.rc(
        family = family,
        dosevars = dosevars,
        data = data,
        deg = deg,
        ERC = FALSE,
        transform = transform,
        transform.jacobian = transform.jacobian,
        Y = Y,
        M = M,
        modifier_info = modifier_info,
        X = X,
        offset = offset,
        inpar = inpar,
        entry = entry,
        exit = exit,
        status = status,
        setnr = setnr,
        doseRRmod = doseRRmod,
        loglim = loglim,
        control = control,
        optim.method = optim.method,
        ...
      )
      results <- c(results, list(RC = fit))
    } else if (method == "ERC") {
      message("Fitting ERC")
      fit <- ameras.rc(
        family = family,
        dosevars = dosevars,
        data = data,
        deg = deg,
        ERC = TRUE,
        transform = transform,
        transform.jacobian = transform.jacobian,
        Y = Y,
        M = M,
        modifier_info = modifier_info,
        X = X,
        offset = offset,
        inpar = inpar,
        entry = entry,
        exit = exit,
        status = status,
        setnr = setnr,
        doseRRmod = doseRRmod,
        loglim = loglim,
        control = control,
        optim.method = optim.method,
        ...
      )
      results <- c(results, list(ERC = fit))
    } else if (method == "FMA") {
      message("Fitting FMA")
      fit <- ameras.fma(
        family = family,
        dosevars = dosevars,
        data = data,
        deg = deg,
        transform = transform,
        transform.jacobian = transform.jacobian,
        Y = Y,
        M = M,
        modifier_info = modifier_info,
        X = X,
        offset = offset,
        inpar = inpar,
        entry = entry,
        exit = exit,
        status = status,
        setnr = setnr,
        doseRRmod = doseRRmod,
        unweighted = unweightedFMA,
        MFMA = MFMA,
        future.chunk.size.FMA = future.chunk.size.FMA,
        control = control,
        ...
      )
      results <- c(results, list(FMA = fit))
    } else if (method == "BMA") {
      if (modifier_is_group_coded(modifier_info)) {
        stop(
          "BMA does not yet support subgroup-coded modifiers. ",
          "Use modifier = ~ M or omit BMA."
        )
      }
      message("Fitting BMA")
      if (!is.null(included.realizations.BMA)) {
        increps <- included.realizations.BMA
      } else if (!is.null(results$FMA$included.realizations)) {
        increps <- results$FMA$included.realizations
      } else {
        increps <- NULL
      }
      fit <- ameras.bma(
        family = family,
        dosevars = dosevars,
        data = data,
        deg = deg,
        transform = transform,
        Y = Y,
        M = M,
        modifier_info = modifier_info,
        X = X,
        offset = offset,
        inpar = inpar,
        entry = entry,
        exit = exit,
        status = status,
        setnr = setnr,
        doseRRmod = doseRRmod,
        ERRprior = ERRprior.BMA,
        prophaz_numints = prophaz.numints.BMA,
        nburnin = nburnin.BMA,
        niter = niter.BMA,
        nchains = nchains.BMA,
        thin = thin.BMA,
        control = control,
        included.realizations = increps,
        optim.method = optim.method,
        ...
      )
      results <- c(results, list(BMA = fit))
    }
  }

  return(results)
}

ameras <- function(
  formula = NULL,
  data,
  family = "gaussian",
  methods = "RC",
  # Old arguments kept temporarily for deprecation
  Y = NULL,
  dosevars = NULL,
  doseRRmod = NULL,
  deg = NULL,
  M = NULL,
  X = NULL,
  offset = NULL,
  entry = NULL,
  exit = NULL,
  setnr = NULL,
  CI = NULL,
  params.profCI = NULL,
  maxit.profCI = NULL,
  tol.profCI = NULL,
  # Current arguments
  transform = NULL,
  transform.jacobian = NULL,
  inpar = NULL,
  loglim = 1e-30,
  MFMA = 100000,
  future.chunk.size.FMA = NULL,
  prophaz.numints.BMA = 10,
  ERRprior.BMA = "doubleexponential",
  nburnin.BMA = 5000,
  niter.BMA = 20000,
  nchains.BMA = 2,
  thin.BMA = 10,
  included.realizations.BMA = NULL,
  included.replicates.BMA = NULL,
  optim.method = "Nelder-Mead",
  control = NULL,
  keep.data = TRUE,
  na.action = getOption("na.action"),
  ...
) {
  call_env <- parent.frame()

  old_interface <- is.null(formula) || !inherits(formula, "formula")
  if (!is.null(included.replicates.BMA)) {
    included.realizations.BMA <- included.replicates.BMA
    lifecycle::deprecate_warn(
      when = "0.3.0",
      what = "ameras(included.replicates.BMA)",
      details = paste0(
        "Use included.realizations.BMA instead."
      )
    )
  }

  # Check for errors
  check_df(data)
  check_family(family)
  na_action_fun <- resolve_na_action(na.action, env = call_env)
  num.rows.original <- nrow(data)
  methods <- check_methods(methods)
  future.chunk.size.FMA <- check_future_chunk_size_FMA(future.chunk.size.FMA)

  if ("BMA" %in% methods) {
    message("Note: BMA may require extensive computation time")
  }

  # Detect which interface is being used and resolve parsed arguments
  if (!old_interface) {
    # New formula interface
    # Error if old arguments are also supplied
    old_supplied <- Filter(
      Negate(is.null),
      list(
        Y = Y,
        dosevars = dosevars,
        doseRRmod = doseRRmod,
        deg = deg,
        M = M,
        X = X,
        offset = offset,
        entry = entry,
        exit = exit,
        setnr = setnr
      )
    )
    if (length(old_supplied)) {
      stop(
        "Old-style arguments cannot be combined with the formula interface: ",
        paste(names(old_supplied), collapse = ", "),
        ".\n",
        "Please use the formula interface only. See ?ameras for details."
      )
    }

    # Warn for CI arguments moved to confint
    if (!is.null(CI)) {
      lifecycle::deprecate_warn(
        when = "0.2.0",
        what = "ameras(CI)",
        details = paste0(
          "CI is deprecated in ameras(). ",
          "Use confint() after fitting instead, e.g.:\n",
          "  fit <- ameras(...)\n",
          "  fit <- confint(fit, type='wald.orig')"
        )
      )
    }

    parsed <- parse_ameras_formula(formula, data, family, env = call_env)
  } else {
    # Old direct argument interface
    lifecycle::deprecate_warn(
      when = "0.2.0",
      what = "ameras(Y)",
      details = paste0(
        "The direct argument interface is deprecated and will be removed ",
        "in the next release. ",
        "Please use the formula interface instead:\n",
        "  ameras(Y ~ dose(D1:D10, model='ERR') + X1, ",
        "data=data, family='binomial')\n",
        "See ?ameras for details."
      )
    )

    if (!is.null(CI)) {
      lifecycle::deprecate_warn(
        when = "0.2.0",
        what = "ameras(CI)",
        details = paste0(
          "CI is deprecated in ameras(). ",
          "Use confint() after fitting instead, e.g.:\n",
          "  fit <- ameras(...)\n",
          "  fit <- confint(fit, method='wald.orig')"
        )
      )
    }

    # Validate old-style arguments
    if (is.null(dosevars)) {
      stop("dosevars is required")
    }
    if (is.null(Y)) {
      stop("Y is required")
    }

    # Construct parsed list from direct arguments
    # matching the structure returned by parse_ameras_formula
    if (family == "prophaz") {
      parsed_Y <- NULL
      parsed_status <- Y
    } else {
      parsed_Y <- Y
      parsed_status <- NULL
    }

    parsed <- list(
      Y = parsed_Y,
      status = parsed_status,
      entry = entry,
      exit = exit,
      dosevars = dosevars,
      doseRRmod = doseRRmod %||% if (family != "gaussian") "ERR" else NULL,
      deg = deg %||% 1,
      M = M,
      modifier_info = new_modifier_info(
        source_vars = M,
        coding = if (!is.null(M)) "legacy" else "none",
        design_vars = M,
        parameter_names = M
      ),
      X = X,
      X_formula = if (!is.null(X)) {
        as.formula(paste("~", paste(X, collapse = "+")))
      } else {
        NULL
      },
      offset = offset,
      setnr = setnr
    )
  }

  X_design <- build_X_design(parsed$X_formula, data)
  if (!is.null(X_design$matrix)) {
    X_matrix <- X_design$matrix
    X_colnames <- colnames(X_matrix)

    # Add expanded columns to data, avoiding name conflicts
    new_cols <- setdiff(X_colnames, colnames(data))

    if (length(new_cols)) {
      data[, new_cols] <- X_matrix[, new_cols, drop = FALSE]
    }

    X <- X_colnames
  } else {
    X <- NULL
  }

  check_reserved_names(data)

  if (family == "clogit") {
    Y <- NULL
    status <- parsed$Y
  } else if (family == "prophaz") {
    Y <- NULL
    status <- parsed$status
  } else {
    Y <- parsed$Y
    status <- NULL
  }

  na_vars <- model_na_vars(
    family = family,
    dosevars = parsed$dosevars,
    Y = Y,
    status = status,
    M = parsed$M,
    X = X,
    offset = parsed$offset,
    entry = parsed$entry,
    exit = parsed$exit,
    setnr = parsed$setnr
  )
  na_result <- apply_na_action_to_data(data, na_vars, na_action_fun)
  data <- na_result$data
  fit_na_action <- na_result$na.action
  if (!nrow(data)) {
    stop("ERROR: no rows remain after applying na.action")
  }

  # Formula modifiers can be factors/logicals. Convert them once to numeric
  # design columns, while keeping source-variable metadata for reconstruction.
  modifier_inputs <- prepare_modifier_inputs(data, parsed$modifier_info)
  data <- modifier_inputs$data
  parsed$modifier_info <- modifier_inputs$modifier_info
  modifier_design_vars <- modifier_inputs$design_vars

  # Checks that depend on parsed formula
  if (family != "gaussian") {
    check_doseRRmod(parsed$doseRRmod)
  }

  check_Y(Y %||% status, data, family)
  check_D(parsed$dosevars, data, methods)
  check_M(modifier_design_vars, data)
  check_X(X, data)

  if (family == "poisson") {
    check_offset(parsed$offset, data)
  }
  if (family == "prophaz") {
    check_entry_exit(parsed$entry, parsed$exit, data)
    data <- filter_prophaz_zero_followup(
      data,
      parsed$entry,
      parsed$exit,
      status
    )
  }
  if (family == "clogit") {
    check_setnr(parsed$setnr, data)
    data <- filter_clogit_sets(data, status, parsed$setnr)
  }

  if (!is.null(X)) {
    warn_if_poorly_conditioned_X(
      X_matrix = data[, X, drop = FALSE],
      family = family
    )
  }

  deg <- check_deg(parsed$deg)

  if (!is.null(parsed$doseRRmod) && parsed$doseRRmod == "LINEXP") {
    deg <- 2
  }

  X_names <- X
  M <- getVarNumbers(modifier_design_vars, data)
  X <- getVarNumbers(X, data)

  if (family != "multinomial") {
    inpar <- check_inpar(inpar, family, M, X, deg)
  } else {
    inpar <- check_inpar(
      inpar,
      family,
      M,
      X,
      deg,
      multinom_levels = length(levels(data[, parsed$Y]))
    )
  }

  # Add mean dose for RC and ERC to the data
  data$rcdose_ameras <- rowMeans(data[, parsed$dosevars, drop = FALSE])
  X_formula_to_store <- parsed$X_formula
  if (!is.null(X_formula_to_store)) {
    environment(X_formula_to_store) <- baseenv()
  }
  model_list <- list(
    data = if (keep.data) data else NULL,
    keep.data = keep.data,
    family = family,
    dosevars = parsed$dosevars,
    Y = Y,
    M = M,
    M_names = parsed$M, # names for resolve_data and required_vars
    modifier_info = parsed$modifier_info,
    X_formula = X_formula_to_store,
    X_design_info = X_design$X_design_info,
    X_names = X_names,
    X = X,
    na_action_fun = na_action_fun,
    offset = parsed$offset,
    entry = parsed$entry,
    exit = parsed$exit,
    status = status,
    setnr = parsed$setnr,
    deg = deg,
    doseRRmod = parsed$doseRRmod,
    loglim = loglim,
    optim.method = optim.method,
    control = control
  )

  result <- ameras_main(
    family = family,
    methods = methods,
    dosevars = parsed$dosevars,
    data = data,
    deg = deg,
    doseRRmod = parsed$doseRRmod,
    transform = transform,
    transform.jacobian = transform.jacobian,
    setnr = parsed$setnr,
    Y = Y,
    M = M,
    modifier_info = parsed$modifier_info,
    X = X,
    offset = parsed$offset,
    inpar = inpar,
    entry = parsed$entry,
    exit = parsed$exit,
    status = status,
    loglim = loglim,
    MFMA = MFMA,
    future.chunk.size.FMA = future.chunk.size.FMA,
    prophaz.numints.BMA = prophaz.numints.BMA,
    ERRprior.BMA = ERRprior.BMA,
    nburnin.BMA = nburnin.BMA,
    niter.BMA = niter.BMA,
    nchains.BMA = nchains.BMA,
    thin.BMA = thin.BMA,
    included.realizations.BMA = included.realizations.BMA %||%
      seq_along(parsed$dosevars),
    control = control,
    optim.method = optim.method,
    ...
  )
  formula_to_store <- formula
  if (!is.null(formula_to_store)) {
    environment(formula_to_store) <- baseenv()
  }
  ret <- c(
    list(
      call = match.call(),
      formula = formula_to_store,
      num.rows = nrow(data),
      num.rows.original = num.rows.original,
      na.action = fit_na_action,
      num.realizations = length(parsed$dosevars),
      transform = result$transform,
      transform.jacobian = result$transform.jacobian,
      other.args = list(...),
      model = model_list,
      CI.computed = FALSE
    ),
    result[setdiff(names(result), c("transform", "transform.jacobian"))]
  )

  ret <- new_amerasfit(ret)

  # Compute CIs for old interface backwards compatibility
  if (old_interface && !is.null(CI)) {
    ret <- confint(
      ret,
      type = CI,
      parm = params.profCI %||% "dose",
      maxit.profCI = maxit.profCI %||% 20,
      tol.profCI = tol.profCI %||% 1e-2
    )
  }

  ret
}
