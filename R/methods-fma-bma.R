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
  MFMA = 100000,
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

    FMAfits <- lapply(1:length(dosevars), function(Xi) {
      fit.FMAi <- optim(
        inpar,
        loglik.gaussian,
        D = dosevars[Xi],
        X = X,
        Y = Y,
        M = M,
        data = data,
        deg = deg,
        ERC = FALSE,
        transform = transform,
        method = optim.method,
        control = control,
        ...
      )
      if (optim.method == "Nelder-Mead") {
        fit.FMAi <- optim(
          fit.FMAi$par,
          loglik.gaussian,
          D = dosevars[Xi],
          X = X,
          Y = Y,
          M = M,
          data = data,
          deg = deg,
          ERC = FALSE,
          transform = transform,
          method = "BFGS",
          control = control,
          ...
        )
      }

      fit.FMAi$hessian <- numDeriv::hessian(
        func = loglik.gaussian,
        x = fit.FMAi$par,
        D = dosevars[Xi],
        X = X,
        Y = Y,
        M = M,
        data = data,
        deg = deg,
        ERC = FALSE,
        transform = transform,
        ...
      )

      if (
        det(fit.FMAi$hessian) != 0 &
          rcond(fit.FMAi$hessian) > .Machine$double.eps &
          fit.FMAi$convergence == 0 &
          all(eigen(fit.FMAi$hessian)$values > 0)
      ) {
        include <- TRUE
      } else {
        include <- FALSE
      }

      list(
        coef = fit.FMAi$par,
        hess = fit.FMAi$hessian,
        AIC = 2 * fit.FMAi$value + 2 * (2 + length(X) + length(M) * deg + deg),
        convergence = fit.FMAi$convergence,
        include = include
      )
    })

    parnames <- c(
      "(Intercept)",
      names(data[, X, drop = FALSE]),
      c("dose", "dose_squared")[1:deg]
    )
    if (!is.null(M)) {
      parnames <- c(parnames, paste0("dose:", names(data[, M, drop = FALSE])))
      if (deg == 2) {
        parnames <- c(
          parnames,
          paste0("dose_squared:", names(data[, M, drop = FALSE]))
        )
      }
    }
    parnames <- c(parnames, "sigma")
  } else if (family == "binomial") {
    if (is.null(Y)) {
      stop("Y is required for family=binomial")
    }

    if (is.null(inpar)) {
      inpar <- rep(0, 1 + length(X) + length(M) * deg + deg)
    }

    FMAfits <- lapply(1:length(dosevars), function(Xi) {
      fit.FMAi <- optim(
        inpar,
        loglik.binomial,
        D = dosevars[Xi],
        X = X,
        Y = Y,
        M = M,
        doseRRmod = doseRRmod,
        data = data,
        deg = deg,
        ERC = FALSE,
        transform = transform,
        method = optim.method,
        control = control,
        ...
      )
      if (optim.method == "Nelder-Mead") {
        fit.FMAi <- optim(
          fit.FMAi$par,
          loglik.binomial,
          D = dosevars[Xi],
          X = X,
          Y = Y,
          M = M,
          doseRRmod = doseRRmod,
          data = data,
          deg = deg,
          ERC = FALSE,
          transform = transform,
          method = "BFGS",
          control = control,
          ...
        )
      }
      fit.FMAi$hessian <- numDeriv::hessian(
        func = loglik.binomial,
        x = fit.FMAi$par,
        D = dosevars[Xi],
        X = X,
        Y = Y,
        M = M,
        doseRRmod = doseRRmod,
        data = data,
        deg = deg,
        ERC = FALSE,
        transform = transform,
        ...
      )

      if (
        det(fit.FMAi$hessian) != 0 &
          rcond(fit.FMAi$hessian) > .Machine$double.eps &
          fit.FMAi$convergence == 0 &
          all(eigen(fit.FMAi$hessian)$values > 0)
      ) {
        include <- TRUE
      } else {
        include <- FALSE
      }

      list(
        coef = fit.FMAi$par,
        hess = fit.FMAi$hessian,
        AIC = 2 * fit.FMAi$value + 2 * (1 + length(X) + length(M) * deg + deg),
        convergence = fit.FMAi$convergence,
        include = include
      )
    })

    if (doseRRmod != "LINEXP") {
      parnames <- c(
        "(Intercept)",
        names(data[, X, drop = FALSE]),
        c("dose", "dose_squared")[1:deg]
      )
      if (!is.null(M)) {
        parnames <- c(parnames, paste0("dose:", names(data[, M, drop = FALSE])))
        if (deg == 2) {
          parnames <- c(
            parnames,
            paste0("dose_squared:", names(data[, M, drop = FALSE]))
          )
        }
      }
    } else {
      parnames <- c(
        "(Intercept)",
        names(data[, X, drop = FALSE]),
        c("dose_linear", "dose_exponential")
      )
      if (!is.null(M)) {
        parnames <- c(
          parnames,
          paste0("dose_linear:", names(data[, M, drop = FALSE]))
        )
        parnames <- c(
          parnames,
          paste0("dose_exponential:", names(data[, M, drop = FALSE]))
        )
      }
    }
  } else if (family == "poisson") {
    if (is.null(Y)) {
      stop("Y is required for family=poisson")
    }

    if (is.null(inpar)) {
      inpar <- rep(0, 1 + length(X) + length(M) * deg + deg)
    }

    FMAfits <- lapply(1:length(dosevars), function(Xi) {
      fit.FMAi <- optim(
        inpar,
        loglik.poisson,
        D = dosevars[Xi],
        X = X,
        Y = Y,
        M = M,
        offset = offset,
        doseRRmod = doseRRmod,
        data = data,
        deg = deg,
        transform = transform,
        method = optim.method,
        control = control,
        ...
      )
      if (optim.method == "Nelder-Mead") {
        fit.FMAi <- optim(
          fit.FMAi$par,
          loglik.poisson,
          D = dosevars[Xi],
          X = X,
          Y = Y,
          M = M,
          offset = offset,
          doseRRmod = doseRRmod,
          data = data,
          deg = deg,
          transform = transform,
          method = "BFGS",
          control = control,
          ...
        )
      }
      fit.FMAi$hessian <- numDeriv::hessian(
        func = loglik.poisson,
        x = fit.FMAi$par,
        D = dosevars[Xi],
        X = X,
        Y = Y,
        M = M,
        offset = offset,
        doseRRmod = doseRRmod,
        data = data,
        deg = deg,
        transform = transform,
        ...
      )

      if (
        det(fit.FMAi$hessian) != 0 &
          rcond(fit.FMAi$hessian) > .Machine$double.eps &
          fit.FMAi$convergence == 0 &
          all(eigen(fit.FMAi$hessian)$values > 0)
      ) {
        include <- TRUE
      } else {
        include <- FALSE
      }

      list(
        coef = fit.FMAi$par,
        hess = fit.FMAi$hessian,
        AIC = 2 * fit.FMAi$value + 2 * (1 + length(X) + length(M) * deg + deg),
        convergence = fit.FMAi$convergence,
        include = include
      )
    })

    if (doseRRmod != "LINEXP") {
      parnames <- c(
        "(Intercept)",
        names(data[, X, drop = FALSE]),
        c("dose", "dose_squared")[1:deg]
      )
      if (!is.null(M)) {
        parnames <- c(parnames, paste0("dose:", names(data[, M, drop = FALSE])))
        if (deg == 2) {
          parnames <- c(
            parnames,
            paste0("dose_squared:", names(data[, M, drop = FALSE]))
          )
        }
      }
    } else {
      parnames <- c(
        "(Intercept)",
        names(data[, X, drop = FALSE]),
        c("dose_linear", "dose_exponential")
      )
      if (!is.null(M)) {
        parnames <- c(
          parnames,
          paste0("dose_linear:", names(data[, M, drop = FALSE]))
        )
        parnames <- c(
          parnames,
          paste0("dose_exponential:", names(data[, M, drop = FALSE]))
        )
      }
    }
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

    FMAfits <- lapply(1:length(dosevars), function(Xi) {
      if (length(X) + length(M) * deg + deg == 1) {
        # Optimize 1-dimensional model: use optimize instead of optim
        fit0 <- optimize(
          f = loglik.clogit,
          lower = -20,
          upper = 5,
          D = dosevars[Xi],
          status = status,
          X = X,
          M = M,
          designmat = designmat,
          doseRRmod = doseRRmod,
          data = data,
          deg = deg,
          ERC = FALSE,
          transform = transform,
          ...
        )
        fit.FMAi <- list(
          par = fit0$minimum,
          value = fit0$objective,
          convergence = 0,
          hessian = numDeriv::hessian(
            func = loglik.clogit,
            x = fit0$minimum,
            D = dosevars[Xi],
            status = status,
            X = X,
            M = M,
            designmat = designmat,
            doseRRmod = doseRRmod,
            data = data,
            deg = deg,
            ERC = FALSE,
            transform = transform,
            ...
          )
        )
      } else {
        fit.FMAi <- optim(
          inpar,
          loglik.clogit,
          D = dosevars[Xi],
          status = status,
          X = X,
          M = M,
          designmat = designmat,
          doseRRmod = doseRRmod,
          data = data,
          deg = deg,
          ERC = FALSE,
          transform = transform,
          method = optim.method,
          control = control,
          ...
        )
        if (optim.method == "Nelder-Mead") {
          fit.FMAi <- optim(
            fit.FMAi$par,
            loglik.clogit,
            D = dosevars[Xi],
            status = status,
            X = X,
            M = M,
            designmat = designmat,
            doseRRmod = doseRRmod,
            data = data,
            deg = deg,
            ERC = FALSE,
            transform = transform,
            method = "BFGS",
            control = control,
            ...
          )
        }
        fit.FMAi$hessian <- numDeriv::hessian(
          func = loglik.clogit,
          x = fit.FMAi$par,
          D = dosevars[Xi],
          status = status,
          X = X,
          M = M,
          designmat = designmat,
          doseRRmod = doseRRmod,
          data = data,
          deg = deg,
          ERC = FALSE,
          transform = transform,
          ...
        )
      }

      if (
        det(fit.FMAi$hessian) != 0 &
          rcond(fit.FMAi$hessian) > .Machine$double.eps &
          fit.FMAi$convergence == 0 &
          all(eigen(fit.FMAi$hessian)$values > 0)
      ) {
        include <- TRUE
      } else {
        include <- FALSE
      }

      list(
        coef = fit.FMAi$par,
        hess = fit.FMAi$hessian,
        AIC = 2 * fit.FMAi$value + 2 * (length(X) + length(M) * deg + deg),
        convergence = fit.FMAi$convergence,
        include = include
      )
    })

    if (doseRRmod != "LINEXP") {
      parnames <- c(
        names(data[, X, drop = FALSE]),
        c("dose", "dose_squared")[1:deg]
      )
      if (!is.null(M)) {
        parnames <- c(parnames, paste0("dose:", names(data[, M, drop = FALSE])))
        if (deg == 2) {
          parnames <- c(
            parnames,
            paste0("dose_squared:", names(data[, M, drop = FALSE]))
          )
        }
      }
    } else {
      parnames <- c(
        names(data[, X, drop = FALSE]),
        c("dose_linear", "dose_exponential")
      )
      if (!is.null(M)) {
        parnames <- c(
          parnames,
          paste0("dose_linear:", names(data[, M, drop = FALSE]))
        )
        parnames <- c(
          parnames,
          paste0("dose_exponential:", names(data[, M, drop = FALSE]))
        )
      }
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

    if (is.null(inpar)) {
      inpar <- rep(0, length(X) + length(M) * deg + deg)
    }

    FMAfits <- lapply(1:length(dosevars), function(Xi) {
      if (length(X) + length(M) * deg + deg == 1) {
        # Optimize 1-dimensional model: use optimize instead of optim
        fit0 <- optimize(
          f = loglik.prophaz,
          lower = -20,
          upper = 5,
          D = dosevars[Xi],
          status = status,
          X = X,
          M = M,
          entry = entry,
          exit = exit,
          doseRRmod = doseRRmod,
          data = data,
          deg = deg,
          ERC = FALSE,
          transform = transform,
          ...
        )
        fit.FMAi <- list(
          par = fit0$minimum,
          value = fit0$objective,
          convergence = 0,
          hessian = numDeriv::hessian(
            func = loglik.prophaz,
            x = fit0$minimum,
            D = dosevars[Xi],
            status = status,
            X = X,
            M = M,
            entry = entry,
            exit = exit,
            doseRRmod = doseRRmod,
            data = data,
            deg = deg,
            ERC = FALSE,
            transform = transform,
            ...
          )
        )
      } else {
        fit.FMAi <- optim(
          inpar,
          loglik.prophaz,
          D = dosevars[Xi],
          status = status,
          X = X,
          M = M,
          entry = entry,
          exit = exit,
          doseRRmod = doseRRmod,
          data = data,
          deg = deg,
          ERC = FALSE,
          transform = transform,
          method = optim.method,
          control = control,
          ...
        )
        if (optim.method == "Nelder-Mead") {
          fit.FMAi <- optim(
            fit.FMAi$par,
            loglik.prophaz,
            D = dosevars[Xi],
            status = status,
            X = X,
            M = M,
            entry = entry,
            exit = exit,
            doseRRmod = doseRRmod,
            data = data,
            deg = deg,
            ERC = FALSE,
            transform = transform,
            method = "BFGS",
            control = control,
            ...
          )
        }
        fit.FMAi$hessian <- numDeriv::hessian(
          func = loglik.prophaz,
          x = fit.FMAi$par,
          D = dosevars[Xi],
          status = status,
          X = X,
          M = M,
          entry = entry,
          exit = exit,
          doseRRmod = doseRRmod,
          data = data,
          deg = deg,
          ERC = FALSE,
          transform = transform,
          ...
        )
      }

      if (
        det(fit.FMAi$hessian) != 0 &
          rcond(fit.FMAi$hessian) > .Machine$double.eps &
          fit.FMAi$convergence == 0 &
          all(eigen(fit.FMAi$hessian)$values > 0)
      ) {
        include <- TRUE
      } else {
        include <- FALSE
      }

      list(
        coef = fit.FMAi$par,
        hess = fit.FMAi$hessian,
        AIC = 2 * fit.FMAi$value + 2 * (length(X) + length(M) * deg + deg),
        convergence = fit.FMAi$convergence,
        include = include
      )
    })

    if (doseRRmod != "LINEXP") {
      parnames <- c(
        names(data[, X, drop = FALSE]),
        c("dose", "dose_squared")[1:deg]
      )
      if (!is.null(M)) {
        parnames <- c(parnames, paste0("dose:", names(data[, M, drop = FALSE])))
        if (deg == 2) {
          parnames <- c(
            parnames,
            paste0("dose_squared:", names(data[, M, drop = FALSE]))
          )
        }
      }
    } else {
      parnames <- c(
        names(data[, X, drop = FALSE]),
        c("dose_linear", "dose_exponential")
      )
      if (!is.null(M)) {
        parnames <- c(
          parnames,
          paste0("dose_linear:", names(data[, M, drop = FALSE]))
        )
        parnames <- c(
          parnames,
          paste0("dose_exponential:", names(data[, M, drop = FALSE]))
        )
      }
    }
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

    FMAfits <- lapply(1:length(dosevars), function(Xi) {
      fit.FMAi <- optim(
        inpar,
        loglik.multinomial,
        D = dosevars[Xi],
        X = X,
        Y = Y,
        M = M,
        doseRRmod = doseRRmod,
        data = data,
        deg = deg,
        ERC = FALSE,
        transform = transform,
        method = optim.method,
        control = control,
        ...
      )
      if (optim.method == "Nelder-Mead") {
        fit.FMAi <- optim(
          fit.FMAi$par,
          loglik.multinomial,
          D = dosevars[Xi],
          X = X,
          Y = Y,
          M = M,
          doseRRmod = doseRRmod,
          data = data,
          deg = deg,
          ERC = FALSE,
          transform = transform,
          method = "BFGS",
          control = control,
          ...
        )
      }
      fit.FMAi$hessian <- numDeriv::hessian(
        func = loglik.multinomial,
        x = fit.FMAi$par,
        D = dosevars[Xi],
        X = X,
        Y = Y,
        M = M,
        doseRRmod = doseRRmod,
        data = data,
        deg = deg,
        ERC = FALSE,
        transform = transform,
        ...
      )

      if (
        det(fit.FMAi$hessian) != 0 &
          rcond(fit.FMAi$hessian) > .Machine$double.eps &
          fit.FMAi$convergence == 0 &
          all(eigen(fit.FMAi$hessian)$values > 0)
      ) {
        include <- TRUE
      } else {
        include <- FALSE
      }

      list(
        coef = fit.FMAi$par,
        hess = fit.FMAi$hessian,
        AIC = 2 * fit.FMAi$value + 2 * (1 + length(X) + length(M) * deg + deg),
        convergence = fit.FMAi$convergence,
        include = include
      )
    })

    if (doseRRmod != "LINEXP") {
      parnames <- c(
        "(Intercept)",
        names(data[, X, drop = FALSE]),
        c("dose", "dose_squared")[1:deg]
      )
      if (!is.null(M)) {
        parnames <- c(parnames, paste0("dose:", names(data[, M, drop = FALSE])))
        if (deg == 2) {
          parnames <- c(
            parnames,
            paste0("dose_squared:", names(data[, M, drop = FALSE]))
          )
        }
      }
    } else {
      parnames <- c(
        "(Intercept)",
        names(data[, X, drop = FALSE]),
        c("dose_linear", "dose_exponential")
      )
      if (!is.null(M)) {
        parnames <- c(
          parnames,
          paste0("dose_linear:", names(data[, M, drop = FALSE]))
        )
        parnames <- c(
          parnames,
          paste0("dose_exponential:", names(data[, M, drop = FALSE]))
        )
      }
    }

    mylv <- levels(data[, Y])

    mylv <- mylv[-length(mylv)]

    parnames <- do.call(
      "c",
      lapply(mylv, function(lv) paste0("(", lv, ")_", parnames))
    )
  }

  #Extra exclusion check for asymmetric variance matrix
  for (iii in 1:length(FMAfits)) {
    fit.FMAi <- FMAfits[[iii]]
    if (fit.FMAi$include) {
      if (!is.null(transform) & !is.null(transform.jacobian)) {
        if (is.function(transform) & is.function(transform.jacobian)) {
          jac <- transform.jacobian(fit.FMAi$coef, ...)
          samplevar <- jac %*% solve(fit.FMAi$hess) %*% t(jac)
        } else {
          stop("transform and transform.jacobian should be functions")
        }
      } else {
        samplevar <- solve(fit.FMAi$hess)
      }
      if (!isSymmetric(samplevar)) {
        FMAfits[[iii]]$include <- FALSE
      }
    }
  }

  #allfits <- FMAfits
  included.replicates <- which(sapply(FMAfits, function(x) x$include))

  if (length(included.replicates) > 0) {
    FMAfits <- FMAfits[sapply(FMAfits, function(x) x$include)]

    meanAIC <- mean(sapply(FMAfits, function(x) x$AIC))

    FMAfits <- lapply(FMAfits, function(x) {
      x$AIC_cent <- x$AIC - meanAIC
      x
    })

    FMAfits <- lapply(FMAfits, function(x) {
      if (unweighted) {
        x$wgt <- 1 / length(FMAfits)
      } else {
        x$wgt <- exp(-.5 * x$AIC_cent) /
          sum(sapply(FMAfits, function(y) exp(-.5 * y$AIC_cent)))
      }
      x$M <- max(round(x$wgt * MFMA), 1)
      x
    })

    FMAsamples <- lapply(
      FMAfits,
      function(fit.FMAi, ...) {
        if (!is.null(transform) & !is.null(transform.jacobian)) {
          if (is.function(transform) & is.function(transform.jacobian)) {
            samplemeans <- transform(fit.FMAi$coef, ...)
            jac <- transform.jacobian(fit.FMAi$coef, ...)
            samplevar <- jac %*% solve(fit.FMAi$hess) %*% t(jac)
          } else {
            stop("transform and transform.jacobian should be functions")
          }
        } else {
          samplemeans <- fit.FMAi$coef
          samplevar <- solve(fit.FMAi$hess)
        }

        return(rmvnorm(n = fit.FMAi$M, mean = samplemeans, sigma = samplevar))

        # if(isSymmetric(samplevar)){
        #   return(rmvnorm(n=fit.FMAi$M, mean=samplemeans, sigma=samplevar))
        # } else{
        #   return(NULL)
        # }
      },
      ...
    )

    FMAsamples <- do.call("rbind", FMAsamples)

    coefs <- colMeans(FMAsamples)
    sd <- apply(FMAsamples, 2, sd)

    FMAsamples <- as.data.frame(FMAsamples)
    names(coefs) <- names(sd) <- names(FMAsamples) <- parnames

    included.samples <- nrow(FMAsamples)
  } else {
    coefs <- sd <- NA * FMAfits[[1]]$coef
    names(coefs) <- names(sd) <- parnames

    included.samples <- 0
  }

  t1 <- proc.time()
  timedif <- t1 - t0
  runtime <- paste(
    round(as.numeric(as.difftime(timedif["elapsed"], units = "secs")), 1),
    "seconds"
  )

  prc_excluded <- round(
    100 * (1 - length(included.replicates) / length(dosevars)),
    1
  )
  if (length(included.replicates) / length(dosevars) < .8) {
    warning(paste0(
      "WARNING: ",
      prc_excluded,
      "% of replicates excluded from model averaging. Try different bounds or starting values."
    ))
  }

  out <- list(
    coefficients = coefs,
    sd = sd,
    included.replicates = included.replicates,
    included.samples = included.samples,
    samples = FMAsamples,
    #includedfits=FMAfits,
    #allfits=allfits,
    runtime = runtime
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
  ERRprior = "doubleexponential",
  prophaz_numints = 10,
  nburnin = 1000,
  niter = 5000,
  included.replicates = 1:length(dosevars),
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

  if (family == "gaussian") {
    if (is.null(Y)) {
      stop("Y is required for family=gaussian")
    }

    if (is.null(included.replicates)) {
      if (is.null(inpar)) {
        inpar <- rep(0, 2 + length(X) + length(M) * deg + deg)
      }

      toInclude <- sapply(1:length(dosevars), function(Xi) {
        fit.FMAi <- optim(
          inpar,
          loglik.gaussian,
          D = dosevars[Xi],
          X = X,
          Y = Y,
          M = M,
          data = data,
          deg = deg,
          ERC = FALSE,
          transform = transform,
          method = optim.method,
          control = control,
          ...
        )
        if (optim.method == "Nelder-Mead") {
          fit.FMAi <- optim(
            fit.FMAi$par,
            loglik.gaussian,
            D = dosevars[Xi],
            X = X,
            Y = Y,
            M = M,
            data = data,
            deg = deg,
            ERC = FALSE,
            transform = transform,
            method = "BFGS",
            control = control,
            ...
          )
        }
        fit.FMAi$hessian <- numDeriv::hessian(
          func = loglik.gaussian,
          x = fit.FMAi$par,
          D = dosevars[Xi],
          X = X,
          Y = Y,
          M = M,
          data = data,
          deg = deg,
          ERC = FALSE,
          transform = transform,
          ...
        )

        if (
          det(fit.FMAi$hessian) != 0 &
            rcond(fit.FMAi$hessian) > .Machine$double.eps &
            fit.FMAi$convergence == 0 &
            all(eigen(fit.FMAi$hessian)$values > 0)
        ) {
          include <- TRUE
        } else {
          include <- FALSE
        }

        return(include)
      })
      included.replicates <- which(toInclude == TRUE)
    }
    dosevars <- dosevars[included.replicates]

    nimbledata <- list(Y = data[, Y], dosemat = data[, dosevars])

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

    if (is.null(included.replicates)) {
      if (is.null(inpar)) {
        inpar <- rep(0, 1 + length(X) + length(M) * deg + deg)
      }

      toInclude <- sapply(1:length(dosevars), function(Xi) {
        fit.FMAi <- optim(
          inpar,
          loglik.binomial,
          D = dosevars[Xi],
          X = X,
          Y = Y,
          M = M,
          doseRRmod = doseRRmod,
          data = data,
          deg = deg,
          ERC = FALSE,
          transform = transform,
          method = optim.method,
          control = control,
          ...
        )
        if (optim.method == "Nelder-Mead") {
          fit.FMAi <- optim(
            fit.FMAi$par,
            loglik.binomial,
            D = dosevars[Xi],
            X = X,
            Y = Y,
            M = M,
            doseRRmod = doseRRmod,
            data = data,
            deg = deg,
            ERC = FALSE,
            transform = transform,
            method = "BFGS",
            control = control,
            ...
          )
        }
        fit.FMAi$hessian <- numDeriv::hessian(
          func = loglik.binomial,
          x = fit.FMAi$par,
          D = dosevars[Xi],
          X = X,
          Y = Y,
          M = M,
          doseRRmod = doseRRmod,
          data = data,
          deg = deg,
          ERC = FALSE,
          transform = transform,
          ...
        )

        if (
          det(fit.FMAi$hessian) != 0 &
            rcond(fit.FMAi$hessian) > .Machine$double.eps &
            fit.FMAi$convergence == 0 &
            all(eigen(fit.FMAi$hessian)$values > 0)
        ) {
          include <- TRUE
        } else {
          include <- FALSE
        }

        return(include)
      })
      included.replicates <- which(toInclude == TRUE)
    }

    dosevars <- dosevars[included.replicates]

    nimbledata <- list(Y = data[, Y], dosemat = data[, dosevars])

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

    if (is.null(included.replicates)) {
      if (is.null(inpar)) {
        inpar <- rep(0, 1 + length(X) + length(M) * deg + deg)
      }

      toInclude <- sapply(1:length(dosevars), function(Xi) {
        fit.FMAi <- optim(
          inpar,
          loglik.poisson,
          D = dosevars[Xi],
          X = X,
          Y = Y,
          M = M,
          offset = offset,
          doseRRmod = doseRRmod,
          data = data,
          deg = deg,
          transform = transform,
          method = optim.method,
          control = control,
          ...
        )
        if (optim.method == "Nelder-Mead") {
          fit.FMAi <- optim(
            fit.FMAi$par,
            loglik.poisson,
            D = dosevars[Xi],
            X = X,
            Y = Y,
            M = M,
            offset = offset,
            doseRRmod = doseRRmod,
            data = data,
            deg = deg,
            transform = transform,
            method = "BFGS",
            control = control,
            ...
          )
        }
        fit.FMAi$hessian <- numDeriv::hessian(
          func = loglik.poisson,
          x = fit.FMAi$par,
          D = dosevars[Xi],
          X = X,
          Y = Y,
          M = M,
          offset = offset,
          doseRRmod = doseRRmod,
          data = data,
          deg = deg,
          transform = transform,
          ...
        )

        if (
          det(fit.FMAi$hessian) != 0 &
            rcond(fit.FMAi$hessian) > .Machine$double.eps &
            fit.FMAi$convergence == 0 &
            all(eigen(fit.FMAi$hessian)$values > 0)
        ) {
          include <- TRUE
        } else {
          include <- FALSE
        }

        return(include)
      })
      included.replicates <- which(toInclude == TRUE)
    }
    dosevars <- dosevars[included.replicates]

    nimbledata <- list(Y = data[, Y], dosemat = data[, dosevars])

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

    if (is.null(included.replicates)) {
      if (is.null(inpar)) {
        inpar <- rep(0, length(X) + length(M) * deg + deg)
      }

      toInclude <- sapply(1:length(dosevars), function(Xi) {
        if (length(X) + length(M) * deg + deg == 1) {
          # Optimize 1-dimensional model: use optimize instead of optim
          fit0 <- optimize(
            f = loglik.prophaz,
            lower = -20,
            upper = 5,
            D = dosevars[Xi],
            status = status,
            X = X,
            M = M,
            entry = entry,
            exit = exit,
            doseRRmod = doseRRmod,
            data = data,
            deg = deg,
            ERC = FALSE,
            transform = transform,
            ...
          )
          fit.FMAi <- list(
            par = fit0$minimum,
            value = fit0$objective,
            convergence = 0,
            hessian = numDeriv::hessian(
              func = loglik.prophaz,
              x = fit0$minimum,
              D = dosevars[Xi],
              status = status,
              X = X,
              M = M,
              entry = entry,
              exit = exit,
              doseRRmod = doseRRmod,
              data = data,
              deg = deg,
              ERC = FALSE,
              transform = transform,
              ...
            )
          )
        } else {
          fit.FMAi <- optim(
            inpar,
            loglik.prophaz,
            D = dosevars[Xi],
            status = status,
            X = X,
            M = M,
            entry = entry,
            exit = exit,
            doseRRmod = doseRRmod,
            data = data,
            deg = deg,
            ERC = FALSE,
            transform = transform,
            method = optim.method,
            control = control,
            ...
          )
          if (optim.method == "Nelder-Mead") {
            fit.FMAi <- optim(
              fit.FMAi$par,
              loglik.prophaz,
              D = dosevars[Xi],
              status = status,
              X = X,
              M = M,
              entry = entry,
              exit = exit,
              doseRRmod = doseRRmod,
              data = data,
              deg = deg,
              ERC = FALSE,
              transform = transform,
              method = "BFGS",
              control = control,
              ...
            )
          }
          fit.FMAi$hessian <- numDeriv::hessian(
            func = loglik.prophaz,
            x = fit.FMAi$par,
            D = dosevars[Xi],
            status = status,
            X = X,
            M = M,
            entry = entry,
            exit = exit,
            doseRRmod = doseRRmod,
            data = data,
            deg = deg,
            ERC = FALSE,
            transform = transform,
            ...
          )
        }

        if (
          det(fit.FMAi$hessian) != 0 &
            rcond(fit.FMAi$hessian) > .Machine$double.eps &
            fit.FMAi$convergence == 0 &
            all(eigen(fit.FMAi$hessian)$values > 0)
        ) {
          include <- TRUE
        } else {
          include <- FALSE
        }

        return(include)
      })
      included.replicates <- which(toInclude == TRUE)
    }

    dosevars <- dosevars[included.replicates]

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
      dosemat = data[, dosevars],
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

    # Remove sets of size 1
    set_counts <- table(data[, setnr])
    valid_sets <- as.numeric(names(set_counts[set_counts > 1]))
    data <- data[data[, setnr] %in% valid_sets, ]

    # Reorder and determine indexing for nimble
    data <- data[order(data[, setnr]), ]
    nsets <- length(unique(data[, setnr]))
    set_sizes <- as.numeric(table(data[, setnr]))
    set_start <- cumsum(c(1, head(set_sizes, -1)))
    set_end <- cumsum(set_sizes)

    if (is.null(included.replicates)) {
      designmat <- t(model.matrix(~ as.factor(data[, setnr]) - 1))

      if (is.null(inpar)) {
        inpar <- rep(0, length(X) + length(M) * deg + deg)
      }

      toInclude <- sapply(1:length(dosevars), function(Xi) {
        if (length(X) + length(M) * deg + deg == 1) {
          # Optimize 1-dimensional model: use optimize instead of optim
          fit0 <- optimize(
            f = loglik.clogit,
            lower = -20,
            upper = 5,
            D = dosevars[Xi],
            status = status,
            X = X,
            M = M,
            designmat = designmat,
            doseRRmod = doseRRmod,
            data = data,
            deg = deg,
            ERC = FALSE,
            transform = transform,
            ...
          )
          fit.FMAi <- list(
            par = fit0$minimum,
            value = fit0$objective,
            convergence = 0,
            hessian = numDeriv::hessian(
              func = loglik.clogit,
              x = fit0$minimum,
              D = dosevars[Xi],
              status = status,
              X = X,
              M = M,
              designmat = designmat,
              doseRRmod = doseRRmod,
              data = data,
              deg = deg,
              ERC = FALSE,
              transform = transform,
              ...
            )
          )
        } else {
          fit.FMAi <- optim(
            inpar,
            loglik.clogit,
            D = dosevars[Xi],
            status = status,
            X = X,
            M = M,
            designmat = designmat,
            doseRRmod = doseRRmod,
            data = data,
            deg = deg,
            ERC = FALSE,
            transform = transform,
            method = optim.method,
            control = control,
            ...
          )
          if (optim.method == "Nelder-Mead") {
            fit.FMAi <- optim(
              fit.FMAi$par,
              loglik.clogit,
              D = dosevars[Xi],
              status = status,
              X = X,
              M = M,
              designmat = designmat,
              doseRRmod = doseRRmod,
              data = data,
              deg = deg,
              ERC = FALSE,
              transform = transform,
              method = "BFGS",
              control = control,
              ...
            )
          }
          fit.FMAi$hessian <- numDeriv::hessian(
            func = loglik.clogit,
            x = fit.FMAi$par,
            D = dosevars[Xi],
            status = status,
            X = X,
            M = M,
            designmat = designmat,
            doseRRmod = doseRRmod,
            data = data,
            deg = deg,
            ERC = FALSE,
            transform = transform,
            ...
          )
        }

        if (
          det(fit.FMAi$hessian) != 0 &
            rcond(fit.FMAi$hessian) > .Machine$double.eps &
            fit.FMAi$convergence == 0 &
            all(eigen(fit.FMAi$hessian)$values > 0)
        ) {
          include <- TRUE
        } else {
          include <- FALSE
        }

        return(include)
      })
      included.replicates <- which(toInclude == TRUE)
    }

    dosevars <- dosevars[included.replicates]

    nimbledata <- list(
      Y = as.numeric(data[, status]),
      dosemat = data[, dosevars]
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
    if (is.null(included.replicates)) {
      if (is.null(inpar)) {
        inpar <- rep(0, (Z - 1) * (1 + length(X) + length(M) * deg + deg))
      }

      toInclude <- sapply(1:length(dosevars), function(Xi) {
        fit.FMAi <- optim(
          inpar,
          loglik.multinomial,
          D = dosevars[Xi],
          X = X,
          Y = Y,
          M = M,
          doseRRmod = doseRRmod,
          data = data,
          deg = deg,
          ERC = FALSE,
          transform = transform,
          method = optim.method,
          control = control,
          ...
        )
        if (optim.method == "Nelder-Mead") {
          fit.FMAi <- optim(
            fit.FMAi$par,
            loglik.multinomial,
            D = dosevars[Xi],
            X = X,
            Y = Y,
            M = M,
            doseRRmod = doseRRmod,
            data = data,
            deg = deg,
            ERC = FALSE,
            transform = transform,
            method = "BFGS",
            control = control,
            ...
          )
        }
        fit.FMAi$hessian <- numDeriv::hessian(
          func = loglik.multinomial,
          x = fit.FMAi$par,
          D = dosevars[Xi],
          X = X,
          Y = Y,
          M = M,
          doseRRmod = doseRRmod,
          data = data,
          deg = deg,
          ERC = FALSE,
          transform = transform,
          ...
        )

        if (
          det(fit.FMAi$hessian) != 0 &
            rcond(fit.FMAi$hessian) > .Machine$double.eps &
            fit.FMAi$convergence == 0 &
            all(eigen(fit.FMAi$hessian)$values > 0)
        ) {
          include <- TRUE
        } else {
          include <- FALSE
        }

        return(include)
      })
      included.replicates <- which(toInclude == TRUE)
    }

    dosevars <- dosevars[included.replicates]

    nimbledata <- list(
      Y = model.matrix(~ data[, Y] - 1),
      dosemat = data[, dosevars]
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

  if (!(family %in% c("prophaz", "clogit"))) {
    parnames <- "(Intercept)"
  } else {
    parnames <- NULL
  }
  if (!is.null(doseRRmod)) {
    if (doseRRmod == "LINEXP") {
      parnames <- c(
        parnames,
        names(data[, X, drop = FALSE]),
        c("dose_linear", "dose_exponential")
      )
      if (!is.null(M)) {
        parnames <- c(
          parnames,
          paste0("dose_linear:", names(data[, M, drop = FALSE]))
        )
        parnames <- c(
          parnames,
          paste0("dose_exponential:", names(data[, M, drop = FALSE]))
        )
      }
    } else {
      parnames <- c(
        parnames,
        names(data[, X, drop = FALSE]),
        c("dose", "dose_squared")[1:deg]
      )
      if (!is.null(M)) {
        parnames <- c(parnames, paste0("dose:", names(data[, M, drop = FALSE])))
        if (deg == 2) {
          parnames <- c(
            parnames,
            paste0("dose_squared:", names(data[, M, drop = FALSE]))
          )
        }
      }
    }
  } else {
    parnames <- c(
      parnames,
      names(data[, X, drop = FALSE]),
      c("dose", "dose_squared")[1:deg]
    )
    if (!is.null(M)) {
      parnames <- c(parnames, paste0("dose:", names(data[, M, drop = FALSE])))
      if (deg == 2) {
        parnames <- c(
          parnames,
          paste0("dose_squared:", names(data[, M, drop = FALSE]))
        )
      }
    }
  }
  if (family == "gaussian") {
    parnames <- c(parnames, "sigma")
  }
  if (family == "prophaz") {
    parnames <- c(parnames, paste0("h0[", 1:prophaz_numints, "]"))
  }

  if (doseRRmod == "LINEXP") {
    parnames <- sub("\\bdose\\b", "dose_linear", parnames)
    parnames <- sub("\\bdose_squared\\b", "dose_exponential", parnames)
  }

  if (family == "multinomial") {
    mylv <- levels(data[, Y])

    mylv <- mylv[-length(mylv)]

    parnames <- do.call(
      "c",
      lapply(mylv, function(lv) paste0("(", lv, ")_", parnames))
    )
  }

  if (nchains > 1) {
    nimblesamples.stacked <- do.call("rbind", nimblesamples)
    for (ichain in 1:length(nimblesamples)) {
      nimblesamples[[ichain]] <- nimblesamples[[ichain]][, c(pars, "col.ind")]
      colnames(nimblesamples[[ichain]]) <- c(parnames, "col.ind")
    }
  } else if (nchains == 1) {
    nimblesamples.stacked <- nimblesamples
    nimblesamples <- nimblesamples[, c(pars, "col.ind")]
    colnames(nimblesamples) <- c(parnames, "col.ind")
  }
  nimblesamples.stacked <- nimblesamples.stacked[, pars]

  coef <- colMeans(nimblesamples.stacked)
  names(coef) <- parnames

  sd <- apply(nimblesamples.stacked, 2, sd)
  names(sd) <- parnames

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
      "WARNING: MCMC convergence cannot be assessed using a single chains"
    )
  }

  # if(CI=="percentile"){
  #   CIlower=apply(nimblesamples.stacked, 2, function(x) quantile(x, .025))
  #   CIupper=apply(nimblesamples.stacked, 2, function(x) quantile(x, .975))
  #   CIresult <- data.frame(lower=CIlower, upper=CIupper)
  # } else if(CI=="hpd"){
  #   CIresult <- as.data.frame(HPDinterval(as.mcmc(nimblesamples.stacked)))
  # }

  # rownames(CIresult) <- parnames

  t1 <- proc.time()
  timedif <- t1 - t0
  runtime <- paste(
    round(as.numeric(as.difftime(timedif["elapsed"], units = "secs")), 1),
    "seconds"
  )

  prc_excluded <- round(100 * (1 - length(included.replicates) / ndoses), 1)
  if (length(included.replicates) / ndoses < .8) {
    warning(paste0(
      "WARNING: ",
      prc_excluded,
      "% of replicates excluded from model averaging. Try different bounds or starting values."
    ))
  }

  out <- list(
    coefficients = coef,
    sd = sd,
    #CI=CIresult,
    Rhat = Rhat,
    samples = nimblesamples,
    included.replicates = included.replicates,
    runtime = runtime
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
