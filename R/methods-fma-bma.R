summarize_fma_realization_fit <- function(
  fit,
  npar,
  transform = NULL,
  transform.jacobian = NULL,
  ...
) {
  out <- list(
    AIC = 2 * fit$value + 2 * npar,
    convergence = fit$convergence,
    include = FALSE,
    exclude.reason = NULL,
    sampling = NULL
  )

  if (!fit_passes_hessian_check(fit)) {
    out$exclude.reason <- "hessian"
    return(out)
  }

  sampling <- prepare_fma_sampling_inputs(
    params = fit$par,
    hessian = fit$hessian,
    transform = transform,
    transform.jacobian = transform.jacobian,
    ...
  )

  if (is.null(sampling)) {
    out$exclude.reason <- "sampling"
    return(out)
  }

  out$include <- TRUE
  out$sampling <- sampling
  out
}

prepare_fma_sampling_inputs <- function(
  params,
  hessian,
  transform = NULL,
  transform.jacobian = NULL,
  ...
) {
  # A realization can pass the optimizer/Hessian screen but still fail when
  # its Hessian is inverted or transformed for FMA sampling. Validate the
  # actual mean/covariance pair before handing it to mvtnorm::rmvnorm().
  if (!is.null(transform) && !is.null(transform.jacobian)) {
    if (!is.function(transform) || !is.function(transform.jacobian)) {
      stop("transform and transform.jacobian should be functions")
    }

    samplemeans <- transform(params, ...)
    cholH <- tryCatch(chol(hessian), error = function(e) NULL)
    if (is.null(cholH)) {
      return(NULL)
    }

    jac <- transform.jacobian(params, ...)
    tmpsolve <- tryCatch(
      backsolve(cholH, t(jac), transpose = TRUE),
      error = function(e) NULL
    )
    if (is.null(tmpsolve)) {
      return(NULL)
    }

    samplevar <- crossprod(tmpsolve)
    #samplevar <- jac %*% solve(hessian) %*% t(jac)
  } else {
    samplemeans <- params
    samplevar <- tryCatch(
      chol2inv(chol(hessian)),
      error = function(e) NULL
    )
  }

  if (
    is.null(samplevar) ||
      !is.numeric(samplemeans) ||
      !is.numeric(samplevar) ||
      length(dim(samplevar)) != 2 ||
      nrow(samplevar) != ncol(samplevar) ||
      nrow(samplevar) != length(samplemeans) ||
      any(!is.finite(samplemeans)) ||
      any(!is.finite(samplevar))
  ) {
    return(NULL)
  }

  list(mean = samplemeans, var = samplevar)
}

fma_realization_lapply <- function(
  X,
  FUN,
  ...,
  future.chunk.size.FMA = NULL
) {
  check_future_chunk_size_FMA(future.chunk.size.FMA)

  # future.apply respects the user's future::plan(). Its default plan is
  # sequential, so this preserves current behavior unless the user opts into a
  # parallel plan outside ameras().
  if (requireNamespace("future.apply", quietly = TRUE)) {
    return(
      future.apply::future_lapply(
        X,
        FUN,
        ...,
        future.seed = FALSE,
        future.chunk.size = future.chunk.size.FMA
      )
    )
  }

  # Keep future.apply optional. Without it, FMA realization fitting follows the
  # historical sequential path.
  lapply(X, FUN, ...)
}

fit_fma_realizations <- function(
  family,
  dosevars,
  data,
  inpar,
  npar,
  optim.method,
  control,
  Y = NULL,
  M = NULL,
  X = NULL,
  offset = NULL,
  entry = NULL,
  exit = NULL,
  status = NULL,
  setnr = NULL,
  doseRRmod = NULL,
  deg = 1,
  transform = NULL,
  transform.jacobian = NULL,
  modifier_info = NULL,
  designmat = NULL,
  future.chunk.size.FMA = NULL,
  ...
) {
  # Keep the per-realization FMA optimization path in one place so each
  # family branch only supplies the model-specific inputs and AIC parameter
  # count.
  fma_realization_lapply(seq_along(dosevars), function(Xi) {
    # Keep each realization fit on the reported subgroup scale when modifiers
    # use no-intercept coding, while evaluating the existing contrast likelihood.
    loglik_transform <- make_modifier_loglik_transform(
      transform = transform,
      modifier_info = modifier_info,
      family = family,
      X = X,
      M = M,
      deg = deg,
      Y = Y,
      data = data
    )
    loglik_fn <- make_single_realization_loglik(
      family = family,
      dose.col = dosevars[Xi],
      data = data,
      Y = Y,
      M = M,
      X = X,
      offset = offset,
      entry = entry,
      exit = exit,
      status = status,
      setnr = setnr,
      doseRRmod = doseRRmod,
      deg = deg,
      ERC = FALSE,
      transform = loglik_transform,
      designmat = designmat
    )
    fit <- fit_objective_with_hessian(
      inpar,
      loglik_fn,
      optim.method = optim.method,
      control = control,
      gradient.check = FALSE,
      ...
    )
    summarize_fma_realization_fit(
      fit,
      npar,
      transform = transform,
      transform.jacobian = transform.jacobian,
      ...
    )
  }, future.chunk.size.FMA = future.chunk.size.FMA)
}

assemble_fma_result <- function(
  FMAfits,
  dosevars,
  parnames,
  inpar,
  MFMA,
  unweighted = FALSE,
  t0 = proc.time()
) {
  n_total <- length(dosevars)
  exclude_reason <- vapply(
    FMAfits,
    function(x) x$exclude.reason %||% if (isTRUE(x$include)) NA_character_ else "hessian",
    character(1)
  )
  hess_excluded <- !is.na(exclude_reason) & exclude_reason == "hessian"
  sampling_excluded <- !is.na(exclude_reason) & exclude_reason == "sampling"
  included <- vapply(FMAfits, function(x) isTRUE(x$include), logical(1))

  n_hess_excluded <- sum(hess_excluded)
  n_hess_included <- n_total - n_hess_excluded

  # First drop fits that never produced a usable optimum/covariance basis.
  if (n_hess_excluded > 0) {
    warning(
      n_hess_excluded,
      " of ",
      n_total,
      " realization(s) excluded due to convergence or ",
      "Hessian issues. Using different bounds or starting values may help."
    )
  }

  # Recompute model-averaging weights after this exclusion so reported weights
  # refer only to realizations that can actually contribute draws.
  n_sampling_excluded <- sum(sampling_excluded)
  if (n_sampling_excluded > 0) {
    warning(
      n_sampling_excluded,
      " of ",
      n_hess_included,
      " Hessian-eligible realization(s) excluded because FMA sampling ",
      "means or covariance matrices contained non-finite values or could ",
      "not be computed."
    )
  }

  original_indices <- which(included)
  FMAfits <- FMAfits[included]
  sampling_inputs <- lapply(FMAfits, function(x) x$sampling)

  if (length(FMAfits) > 0) {
    minAIC <- min(sapply(FMAfits, function(x) x$AIC))

    FMAfits <- lapply(FMAfits, function(x) {
      x$AIC_cent <- x$AIC - minAIC
      x
    })

    FMAfits <- lapply(FMAfits, function(x) {
      if (unweighted) {
        x$wgt <- 1 / length(FMAfits)
      } else {
        x$wgt <- exp(-.5 * x$AIC_cent) /
          sum(sapply(FMAfits, function(y) exp(-.5 * y$AIC_cent)))
      }
      x$M <- round(x$wgt * MFMA)
      if (x$M == 0) {
        x$include <- FALSE
      }
      x
    })

    weight_included <- sapply(FMAfits, function(x) x$include)
    included.realizations <- original_indices[weight_included]
    n_weight_excluded <- sum(!weight_included)

    # Fits with tiny model-averaging weights are valid, but contribute no
    # simulated draws after integer allocation for the requested MFMA.
    if (n_weight_excluded > 0) {
      message(
        n_weight_excluded,
        " realization(s) excluded due to negligible model averaging weight ",
        "(weight < ",
        round(0.5 / MFMA, 8),
        " corresponding to 0 samples ",
        "for MFMA=",
        MFMA,
        ")."
      )
    }

    if (length(included.realizations) > 0) {
      weight_included <- sapply(FMAfits, function(x) x$include)
      FMAfits <- FMAfits[weight_included]
      sampling_inputs <- sampling_inputs[weight_included]

      FMAsamples <- Map(
        function(fit.FMAi, sampling_input) {
          rmvnorm(
            n = fit.FMAi$M,
            mean = sampling_input$mean,
            sigma = sampling_input$var
          )
        },
        FMAfits,
        sampling_inputs
      )

      FMAsamples <- do.call("rbind", FMAsamples)

      coefs <- colMeans(FMAsamples)
      sd <- apply(FMAsamples, 2, sd)
      vcov <- var(FMAsamples)

      FMAsamples <- as.data.frame(FMAsamples)
      rownames(vcov) <- colnames(vcov) <- names(coefs) <- names(sd) <- names(
        FMAsamples
      ) <- parnames

      included.samples <- nrow(FMAsamples)
      wgts <- sapply(FMAfits, function(x) x$wgt)
      names(wgts) <- dosevars[included.realizations]
    } else {
      FMAsamples <- NULL
      coefs <- sd <- NA * inpar
      vcov <- matrix(NA, ncol = length(inpar), nrow = length(inpar))
      rownames(vcov) <- colnames(vcov) <- names(coefs) <- names(sd) <- parnames

      included.samples <- 0
      wgts <- NULL
    }
  } else {
    included.realizations <- integer(0)
    FMAsamples <- NULL
    coefs <- sd <- NA * inpar
    vcov <- matrix(NA, ncol = length(inpar), nrow = length(inpar))
    rownames(vcov) <- colnames(vcov) <- names(coefs) <- names(sd) <- parnames

    included.samples <- 0
    wgts <- NULL
  }

  fit_timing <- stop_runtime_timer(t0)
  timing <- new_method_timing(fit = fit_timing)

  n_final <- length(included.realizations)

  if (n_final == 0) {
    warning(
      "No realizations contributed to FMA. ",
      "Returning NA for all estimates."
    )
  } else if (n_final == 1) {
    warning(
      "FMA is based on a single realization. ",
      "Results do not incorporate model averaging uncertainty. ",
      "Consider using RC instead."
    )
  } else if (n_final <= 5) {
    warning(
      "FMA is based on only ",
      n_final,
      " realizations. "
    )
  }

  list(
    coefficients = coefs,
    sd = sd,
    vcov = vcov,
    included.realizations = included.realizations,
    included.samples = included.samples,
    weights = wgts,
    samples = FMAsamples,
    timing = timing,
    runtime = format_runtime(timing$total$cpu)
  )
}

ameras.fma <- function(
  family,
  dosevars,
  data,
  deg,
  transform = NULL,
  transform.jacobian = NULL,
  Y = NULL,
  M = NULL,
  X = NULL,
  offset = NULL,
  inpar = NULL,
  entry = NULL,
  exit = NULL,
  status = NULL,
  setnr = setnr,
  unweighted = NULL,
  doseRRmod = NULL,
  modifier_info = NULL,
  MFMA = 100000,
  future.chunk.size.FMA = NULL,
  optim.method = "Nelder-Mead",
  control = list(reltol = 1e-10),
  ...
) {
  if (is.null(unweighted)) {
    unweighted <- FALSE
  } else if (!is.logical(unweighted)) {
    stop("unweighted should be TRUE or FALSE")
  }

  t0 <- proc.time()
  if (family == "gaussian") {
    if (is.null(Y)) {
      stop("Y is required for family=gaussian")
    }

    if (is.null(inpar)) {
      inpar <- rep(0, 2 + length(X) + length(M) * deg + deg)
    }

    FMAfits <- fit_fma_realizations(
      family = "gaussian",
      dosevars = dosevars,
      data = data,
      inpar = inpar,
      npar = 2 + length(X) + length(M) * deg + deg,
      optim.method = optim.method,
      control = control,
      Y = Y,
      M = M,
      X = X,
      deg = deg,
      transform = transform,
      transform.jacobian = transform.jacobian,
      modifier_info = modifier_info,
      future.chunk.size.FMA = future.chunk.size.FMA,
      ...
    )

  } else if (family == "binomial") {
    if (is.null(Y)) {
      stop("Y is required for family=binomial")
    }

    if (is.null(inpar)) {
      inpar <- rep(0, 1 + length(X) + length(M) * deg + deg)
    }

    FMAfits <- fit_fma_realizations(
      family = "binomial",
      dosevars = dosevars,
      data = data,
      inpar = inpar,
      npar = 1 + length(X) + length(M) * deg + deg,
      optim.method = optim.method,
      control = control,
      Y = Y,
      M = M,
      X = X,
      doseRRmod = doseRRmod,
      deg = deg,
      transform = transform,
      transform.jacobian = transform.jacobian,
      modifier_info = modifier_info,
      future.chunk.size.FMA = future.chunk.size.FMA,
      ...
    )

  } else if (family == "poisson") {
    if (is.null(Y)) {
      stop("Y is required for family=poisson")
    }

    if (is.null(inpar)) {
      inpar <- rep(0, 1 + length(X) + length(M) * deg + deg)
    }

    FMAfits <- fit_fma_realizations(
      family = "poisson",
      dosevars = dosevars,
      data = data,
      inpar = inpar,
      npar = 1 + length(X) + length(M) * deg + deg,
      optim.method = optim.method,
      control = control,
      Y = Y,
      M = M,
      X = X,
      offset = offset,
      doseRRmod = doseRRmod,
      deg = deg,
      transform = transform,
      transform.jacobian = transform.jacobian,
      modifier_info = modifier_info,
      future.chunk.size.FMA = future.chunk.size.FMA,
      ...
    )

  } else if (family == "clogit") {
    if (is.null(doseRRmod)) {
      stop("doseRRmod is required for family=clogit")
    }
    if (is.null(status)) {
      stop("status is required for family=clogit")
    }

    designmat <- t(model.matrix(~ as.factor(data[, setnr]) - 1))

    if (is.null(inpar)) {
      inpar <- rep(0, length(X) + length(M) * deg + deg)
    }

    FMAfits <- fit_fma_realizations(
      family = "clogit",
      dosevars = dosevars,
      data = data,
      inpar = inpar,
      npar = length(X) + length(M) * deg + deg,
      optim.method = optim.method,
      control = control,
      M = M,
      X = X,
      status = status,
      doseRRmod = doseRRmod,
      deg = deg,
      transform = transform,
      transform.jacobian = transform.jacobian,
      modifier_info = modifier_info,
      designmat = designmat,
      future.chunk.size.FMA = future.chunk.size.FMA,
      ...
    )

  } else if (family == "prophaz") {
    if (is.null(exit)) {
      stop("exit is required for family=prophaz")
    }
    if (is.null(doseRRmod)) {
      stop("doseRRmod is required for family=prophaz")
    }
    if (is.null(status)) {
      stop("status is required for family=prophaz")
    }

    if (is.null(inpar)) {
      inpar <- rep(0, length(X) + length(M) * deg + deg)
    }

    FMAfits <- fit_fma_realizations(
      family = "prophaz",
      dosevars = dosevars,
      data = data,
      inpar = inpar,
      npar = length(X) + length(M) * deg + deg,
      optim.method = optim.method,
      control = control,
      M = M,
      X = X,
      entry = entry,
      exit = exit,
      status = status,
      doseRRmod = doseRRmod,
      deg = deg,
      transform = transform,
      transform.jacobian = transform.jacobian,
      modifier_info = modifier_info,
      future.chunk.size.FMA = future.chunk.size.FMA,
      ...
    )

  } else if (family == "multinomial") {
    if (is.null(Y)) {
      stop("Y is required for family=multinomial")
    }

    if (is.null(inpar)) {
      inpar <- rep(
        0,
        (length(unique(data[, Y])) - 1) *
          (1 + length(X) + length(M) * deg + deg)
      )
    }

    FMAfits <- fit_fma_realizations(
      family = "multinomial",
      dosevars = dosevars,
      data = data,
      inpar = inpar,
      npar = 1 + length(X) + length(M) * deg + deg,
      optim.method = optim.method,
      control = control,
      Y = Y,
      M = M,
      X = X,
      doseRRmod = doseRRmod,
      deg = deg,
      transform = transform,
      transform.jacobian = transform.jacobian,
      modifier_info = modifier_info,
      future.chunk.size.FMA = future.chunk.size.FMA,
      ...
    )

  }
  parnames <- make_method_parnames(
    family = family,
    data = data,
    Y = Y,
    X = X,
    M = M,
    deg = deg,
    doseRRmod = doseRRmod,
    modifier_info = modifier_info
  )

  out <- assemble_fma_result(
    FMAfits = FMAfits,
    dosevars = dosevars,
    parnames = parnames,
    inpar = inpar,
    MFMA = MFMA,
    unweighted = unweighted,
    t0 = t0
  )

  return(out)
}


ameras.bma <- function(
  family,
  dosevars,
  data,
  deg,
  Y = NULL,
  M = NULL,
  X = NULL,
  offset = NULL,
  entry = NULL,
  exit = NULL,
  status = NULL,
  setnr = setnr,
  transform = NULL,
  inpar = NULL,
  doseRRmod = NULL,
  modifier_info = NULL,
  ERRprior = "doubleexponential",
  prophaz_numints = 10,
  nburnin = 1000,
  niter = 5000,
  included.realizations = 1:length(dosevars),
  nchains = 2,
  thin = 10,
  optim.method = "Nelder-Mead",
  control = list(reltol = 1e-10),
  ...
) {
  # Remove build warnings, local functions may need to be pulled out from this function
  # HPDinterval <- K <- Mlen <- Mmat <- N <- Xlen <- Xmat <- a <- as.mcmc <-
  #   b <- bm <- bm1 <- bm2 <- buildMCMC <- nsets <- setmat <-
  #   col.ind <- compileNimble <- configureMCMC <- delta <-
  #   dosemat <- equals <- h0 <- inprod <- nimbleModel <- runMCMC <- NULL

  ndoses <- length(dosevars)

  if (family != "gaussian") {
    if (doseRRmod %in% c("ERR", "LINEXP") & is.null(ERRprior)) {
      stop("Please specify prior for ERR parameters")
    } else if (doseRRmod %in% c("ERR", "LINEXP") & !is.null(ERRprior)) {
      if (
        !(ERRprior %in%
          c(
            "truncated_normal",
            "truncated_horseshoe",
            "truncated_doubleexponential",
            "normal",
            "horseshoe",
            "doubleexponential"
          ))
      ) {
        stop(
          "Incorrect prior for ERR parameters specified, should be truncated_normal, truncated_horseshoe, truncated_doubleexponential, normal, horseshoe, or doubleexponential"
        )
      }
    }
  }

  t0 <- proc.time()

  if (is.null(included.realizations)) {
    stop(
      "Automatic BMA realization screening is not currently enabled. ",
      "Supply included.realizations explicitly."
    )
  }
  check_included_realizations_BMA(included.realizations, dosevars)

  if (family == "gaussian") {
    if (is.null(Y)) {
      stop("Y is required for family=gaussian")
    }

    dosevars <- dosevars[included.realizations]

    nimbledata <- list(Y = data[, Y], dosemat = data[, dosevars, drop = FALSE])

    nimbleinits <- function() {
      L <- list(
        a0 = rnorm(1),
        b = rexp(deg),
        sigma = rexp(1),
        col.ind = sample(1:length(dosevars), 1)
      )

      if (length(X) > 0) {
        L <- c(L, list(a = rnorm(length(X))))
      }

      if (length(M) > 0) {
        if (deg == 1) {
          L <- c(L, list(bm = rexp(length(M))))
        } else if (deg > 1) {
          L <- c(L, list(bm1 = rexp(length(M)), bm2 = rexp(length(M))))
        }
      }
      L
    }
  } else if (family == "binomial") {
    if (is.null(Y)) {
      stop("Y is required for family=binomial")
    }
    if (is.null(doseRRmod)) {
      stop("doseRRmod is required for family=binomial")
    }

    dosevars <- dosevars[included.realizations]

    nimbledata <- list(Y = data[, Y], dosemat = data[, dosevars, drop = FALSE])

    nimbleinits <- function() {
      L <- list(
        a0 = rnorm(1),
        b = rexp(deg),
        col.ind = sample(1:length(dosevars), 1)
      )

      if (length(X) > 0) {
        L <- c(L, list(a = rnorm(length(X))))
      }

      if (length(M) > 0) {
        if (deg == 1) {
          L <- c(L, list(bm = rexp(length(M))))
        } else if (deg > 1) {
          L <- c(L, list(bm1 = rexp(length(M)), bm2 = rexp(length(M))))
        }
      }
      L
    }
  } else if (family == "poisson") {
    if (is.null(Y)) {
      stop("Y is required for family=poisson")
    }
    if (is.null(doseRRmod)) {
      stop("doseRRmod is required for family=poisson")
    }

    if (is.null(offset)) {
      P <- rep(1, nrow(data))
    } else {
      P <- data[, offset]
    }

    dosevars <- dosevars[included.realizations]

    nimbledata <- list(Y = data[, Y], dosemat = data[, dosevars, drop = FALSE])

    nimbleinits <- function() {
      L <- list(
        a0 = rnorm(1),
        b = rexp(deg),
        col.ind = sample(1:length(dosevars), 1)
      )

      if (length(X) > 0) {
        L <- c(L, list(a = rnorm(length(X))))
      }

      if (length(M) > 0) {
        if (deg == 1) {
          L <- c(L, list(bm = rexp(length(M))))
        } else if (deg > 1) {
          L <- c(L, list(bm1 = rexp(length(M)), bm2 = rexp(length(M))))
        }
      }
      L
    }
  } else if (family == "prophaz") {
    if (is.null(exit)) {
      stop("exit is required for family=prophaz")
    }
    if (is.null(doseRRmod)) {
      stop("doseRRmod is required for family=prophaz")
    }
    if (is.null(status)) {
      stop("status is required for family=prophaz")
    }

    dosevars <- dosevars[included.realizations]

    prophaz_timepoints <- as.numeric(quantile(
      data[data[, status] == 1, exit],
      probs = seq(0, 1, length.out = prophaz_numints + 1)
    )) # define cut points using quantiles of event time distribution
    if (is.null(entry)) {
      prophaz_timepoints[1] <- 0
    } else {
      prophaz_timepoints[1] <- min(data[, entry])
    }

    prophaz_timepoints[prophaz_numints + 1] <- max(data[, exit]) + .001

    if (!is.null(entry)) {
      int.entry <- as.numeric(cut(
        data[, entry],
        breaks = prophaz_timepoints,
        right = FALSE
      )) # indicator for which interval entry time belongs to
    } else {
      int.entry <- rep(1, nrow(data))
    }

    int.exit <- as.numeric(cut(
      data[, exit],
      breaks = prophaz_timepoints,
      right = FALSE
    )) # indicator for which interval exit time belongs to

    nimbledata <- list(
      delta = data[, status],
      exit = data[, exit],
      dosemat = data[, dosevars, drop = FALSE],
      zeros = rep(0, nrow(data))
    )
    if (!is.null(entry)) {
      nimbledata <- c(nimbledata, list(entry = data[, entry]))
    } else {
      nimbledata <- c(nimbledata, list(entry = rep(0, nrow(data))))
    }

    nimbleinits <- function() {
      L <- list(
        b = rexp(deg),
        col.ind = sample(1:length(dosevars), 1),
        h0 = runif(prophaz_numints, .1)
      )

      if (length(X) > 0) {
        L <- c(L, list(a = rnorm(length(X))))
      }

      if (length(M) > 0) {
        if (deg == 1) {
          L <- c(L, list(bm = rexp(length(M))))
        } else if (deg > 1) {
          L <- c(L, list(bm1 = rexp(length(M)), bm2 = rexp(length(M))))
        }
      }
      L
    }
  } else if (family == "clogit") {
    if (is.null(doseRRmod)) {
      stop("doseRRmod is required for family=clogit")
    }

    # Reorder and determine indexing for nimble
    data <- data[order(data[, setnr]), ]
    nsets <- length(unique(data[, setnr]))
    set_sizes <- as.numeric(table(data[, setnr]))
    set_start <- cumsum(c(1, head(set_sizes, -1)))
    set_end <- cumsum(set_sizes)

    dosevars <- dosevars[included.realizations]

    nimbledata <- list(
      Y = as.numeric(data[, status]),
      dosemat = data[, dosevars, drop = FALSE]
    )

    nimbleinits <- function() {
      L <- list(b = rexp(deg), col.ind = sample(1:length(dosevars), 1))

      if (length(X) > 0) {
        L <- c(L, list(a = rnorm(length(X))))
      }

      if (length(M) > 0) {
        if (deg == 1) {
          L <- c(L, list(bm = rexp(length(M))))
        } else if (deg > 1) {
          L <- c(L, list(bm1 = rexp(length(M)), bm2 = rexp(length(M))))
        }
      }
      L
    }
  } else if (family == "multinomial") {
    if (is.null(Y)) {
      stop("Y is required for family=multinomial")
    }
    if (is.null(doseRRmod)) {
      stop("doseRRmod is required for family=multinomial")
    }

    Z <- nlevels(data[, Y])
    dosevars <- dosevars[included.realizations]

    nimbledata <- list(
      Y = model.matrix(~ data[, Y] - 1),
      dosemat = data[, dosevars, drop = FALSE]
    )

    nimbleinits <- function() {
      L <- list(
        a0 = rnorm(Z - 1),
        b = matrix(rexp(deg * (Z - 1)), ncol = Z - 1),
        col.ind = sample(1:length(dosevars), 1)
      )
      if (deg == 1) {
        L$b <- as.vector(L$b)
      }

      if (length(X) > 0) {
        L <- c(L, list(a = matrix(rnorm((Z - 1) * length(X)), ncol = Z - 1)))
      }

      if (length(M) > 0) {
        if (deg == 1) {
          L <- c(L, list(bm = matrix(rexp((Z - 1) * length(M)), ncol = Z - 1)))
        } else if (deg > 1) {
          L <- c(
            L,
            list(
              bm1 = matrix(rexp((Z - 1) * length(M)), ncol = Z - 1),
              bm2 = matrix(rexp((Z - 1) * length(M)), ncol = Z - 1)
            )
          )
        }
      }
      L
    }
  }

  if (!(family %in% c("prophaz", "clogit"))) {
    mons <- c("a0", "b", "col.ind")
  } else {
    mons <- c("b", "col.ind")
  }

  if (family == "gaussian") {
    doseRRmod <- "EXP"
  }
  nimbleconst <- list(
    N = nrow(data),
    K = length(dosevars),
    deg = deg,
    Xlen = length(X),
    Mlen = length(M),
    w = rep(1 / length(dosevars), length(dosevars)),
    doseRRmod = doseRRmod,
    family = family
  )
  if (family == "poisson") {
    nimbleconst <- c(nimbleconst, list(P = P))
  }
  if (family == "multinomial") {
    nimbleconst <- c(nimbleconst, list(Z = Z))
  }
  if (length(X) > 0) {
    nimbleconst <- c(nimbleconst, list(Xmat = data[, X]))
    mons <- c(mons, "a")
  }
  if (length(M) > 0) {
    nimbleconst <- c(nimbleconst, list(Mmat = data[, M]))
    if (deg == 1) {
      mons <- c(mons, "bm")
    } else if (deg > 1) {
      mons <- c(mons, "bm1", "bm2")
    }
  }
  if (family == "gaussian") {
    mons <- c(mons, "sigma")
  }
  if (family == "prophaz") {
    nimbleconst <- c(
      nimbleconst,
      list(
        prophaz_timepoints = prophaz_timepoints,
        int.exit = int.exit,
        int.entry = int.entry,
        prophaz_numints = prophaz_numints
      )
    )
    mons <- c(mons, "h0")
  }
  if (family == "clogit") {
    nimbleconst <- c(
      nimbleconst,
      list(nsets = nsets, set_start = set_start, set_end = set_end)
    )
  }

  if (doseRRmod %in% c("ERR", "LINEXP")) {
    nimbleconst <- c(nimbleconst, list(ERRprior = ERRprior))
  }

  mymod <- nimbleModel(
    nimblemod,
    data = nimbledata,
    constants = nimbleconst,
    inits = list(col.ind = sample(1:length(dosevars), 1))
  )

  mymod_C <- compileNimble(mymod)

  mymod <- configureMCMC(mymod, monitors = mons, thin = thin, print = FALSE)

  mymod$removeSamplers('col.ind')
  mymod$addSampler(
    target = 'col.ind',
    type = boundedSliceSampler,
    control = list(
      K = length(dosevars),
      adaptive = FALSE,
      sliceWidth = round(length(dosevars) / 2)
    )
  )

  mymod_MCMC <- buildMCMC(mymod)
  mymod_compiled <- compileNimble(mymod_MCMC, project = mymod_C)

  nimblesamples <- runMCMC(
    mymod_compiled,
    nburnin = nburnin,
    niter = niter,
    nchains = nchains,
    inits = nimbleinits
  )

  if (family != "multinomial") {
    if (!(family %in% c("prophaz", "clogit"))) {
      pars <- "a0"
    } else {
      pars <- NULL
    }

    if (length(X) == 1) {
      pars <- c(pars, "a")
    } else if (length(X) > 1) {
      pars <- c(pars, paste0("a[", 1:length(X), "]"))
    }
    if (deg == 1) {
      pars <- c(pars, "b")
    } else if (deg > 1) {
      pars <- c(pars, "b[1]", "b[2]")
    }
    if (length(M) == 1) {
      if (deg == 1) {
        pars <- c(pars, "bm")
      } else if (deg > 1) {
        pars <- c(pars, "bm1", "bm2")
      }
    } else if (length(M) > 1) {
      if (deg == 1) {
        pars <- c(pars, paste0("bm[", 1:length(M), "]"))
      } else if (deg > 1) {
        pars <- c(
          pars,
          paste0("bm1[", 1:length(M), "]"),
          paste0("bm2[", 1:length(M), "]")
        )
      }
    }
    if (family == "gaussian") {
      pars <- c(pars, "sigma")
    }
    if (family == "prophaz") {
      pars <- c(pars, paste0("h0[", 1:prophaz_numints, "]"))
    }
  } else {
    # multinomial
    pars <- NULL
    for (z in 1:(Z - 1)) {
      pars <- c(pars, paste0("a0[", z, "]"))

      if (length(X) == 1) {
        pars <- c(pars, paste0("a[", z, "]"))
      } else if (length(X) > 1) {
        pars <- c(pars, paste0("a[", 1:length(X), ", ", z, "]"))
      }

      if (deg == 1) {
        pars <- c(pars, paste0("b[", z, "]"))
      } else {
        pars <- c(pars, paste0("b[1, ", z, "]"), paste0("b[2, ", z, "]"))
      }

      if (length(M) == 1) {
        if (deg == 1) {
          pars <- c(pars, paste0("bm[", z, "]"))
        } else {
          pars <- c(pars, paste0("bm1[", z, "]"), paste0("bm2[", z, "]"))
        }
      } else if (length(M) > 1) {
        if (deg == 1) {
          pars <- c(pars, paste0("bm[", 1:length(M), ", ", z, "]"))
        } else {
          pars <- c(
            pars,
            paste0("bm1[", 1:length(M), ", ", z, "]"),
            paste0("bm2[", 1:length(M), ", ", z, "]")
          )
        }
      }
    }
  }

  parnames <- make_method_parnames(
    family = family,
    data = data,
    Y = Y,
    X = X,
    M = M,
    deg = deg,
    doseRRmod = doseRRmod,
    prophaz_numints = if (family == "prophaz") prophaz_numints else NULL,
    modifier_info = modifier_info
  )

  format_bma_samples <- function(samples, include_col_ind = TRUE) {
    param_samples <- samples[, pars, drop = FALSE]
    param_samples <- modifier_internal_to_reported_sample_matrix(
      samples = param_samples,
      family = family,
      X = X,
      M = M,
      deg = deg,
      modifier_info = modifier_info,
      Y = Y,
      data = data
    )
    colnames(param_samples) <- parnames

    if (include_col_ind) {
      param_samples <- cbind(
        param_samples,
        col.ind = samples[, "col.ind", drop = TRUE]
      )
      colnames(param_samples)[ncol(param_samples)] <- "col.ind"
    }

    param_samples
  }

  if (nchains > 1) {
    nimblesamples.stacked <- do.call("rbind", nimblesamples)
    nimblesamples <- lapply(nimblesamples, format_bma_samples)
  } else if (nchains == 1) {
    nimblesamples.stacked <- nimblesamples
    nimblesamples <- format_bma_samples(nimblesamples)
  }
  nimblesamples.stacked <- format_bma_samples(
    nimblesamples.stacked,
    include_col_ind = FALSE
  )

  coef <- colMeans(nimblesamples.stacked)
  names(coef) <- parnames

  sd <- apply(nimblesamples.stacked, 2, sd)
  names(sd) <- parnames

  vcov <- var(nimblesamples.stacked)
  rownames(vcov) <- colnames(vcov) <- parnames

  mcmcsum <- MCMCsummary(nimblesamples)
  Rhat <- mcmcsum[-nrow(mcmcsum), c("Rhat", "n.eff")]
  if (nchains > 1) {
    if (any(Rhat$Rhat > 1.05)) {
      warning(
        "WARNING: Potential problems with MCMC convergence, consider using longer chains"
      )
    }
  } else {
    warning(
      "WARNING: MCMC convergence cannot be assessed using a single chain"
    )
  }

  fit_timing <- stop_runtime_timer(t0)
  timing <- new_method_timing(fit = fit_timing)

  prc_excluded <- round(100 * (1 - length(included.realizations) / ndoses), 1)
  if (length(included.realizations) / ndoses < .8) {
    warning(paste0(
      "WARNING: ",
      prc_excluded,
      "% of realizations excluded from model averaging. Try different bounds or starting values."
    ))
  }

  out <- list(
    coefficients = coef,
    sd = sd,
    vcov = vcov,
    #CI=CIresult,
    Rhat = Rhat,
    samples = nimblesamples,
    included.realizations = included.realizations,
    timing = timing,
    runtime = format_runtime(timing$total$cpu)
  )
  if (family == "prophaz") {
    out <- c(
      out,
      list(
        prophaz_timepoints = prophaz_timepoints
      )
    )
  }
  return(out)
}
