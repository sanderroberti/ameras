make_base_loglik_fn <- function(object, method_fit, data) {
  m <- object$model
  ERC <- method_fit$ERC
  if (is.null(ERC)) {
    ERC <- FALSE
  }
  if (ERC && m$family == "prophaz") {
    ord_exit <- order(data[[m$exit]])
    dosemat_ord <- as.matrix(data[ord_exit, m$dosevars, drop = FALSE])
    storage.mode(dosemat_ord) <- "double"
    rcdose_ord <- data[ord_exit, "rcdose_ameras"]
    Xc_ord <- dosemat_ord - rcdose_ord
    Kmat_diag_ord <- rowSums(Xc_ord^2) / (ncol(dosemat_ord) - 1)
    rm(dosemat_ord)
    gc()

    Kmat <- NULL
  } else if (ERC && m$family == "poisson") {
    dosemat_poisson <- as.matrix(data[, m$dosevars, drop = FALSE])
    storage.mode(dosemat_poisson) <- "double"
    Xc_poisson <- dosemat_poisson - data[, "rcdose_ameras"]
    Kmat_diag_poisson <- rowSums(Xc_poisson^2) / (ncol(dosemat_poisson) - 1)
    rm(dosemat_poisson)
    gc()

    Kmat <- NULL
  } else if (ERC) {
    Kmat <- cov(t(data[, m$dosevars, drop = FALSE]))
  } else {
    Kmat <- NULL
  }

  if (m$family == "prophaz") {
    if (ERC) {
      # Xc_ord and Kmat_diag_ord captured in closure
      function(params, D) {
        do.call(
          loglik.prophaz.erc,
          c(
            list(
              params = params,
              D = D,
              X = m$X,
              status = m$status,
              entry = m$entry,
              exit = m$exit,
              M = m$M,
              doseRRmod = m$doseRRmod,
              data = data,
              deg = m$deg,
              loglim = m$loglim,
              transform = object$transform,
              Xc_ord = Xc_ord,
              Kmat_diag_ord = Kmat_diag_ord
            ),
            object$transform.args
          )
        )
      }
    } else {
      function(params, D) {
        do.call(
          loglik.prophaz,
          c(
            list(
              params = params,
              D = D,
              X = m$X,
              status = m$status,
              entry = m$entry,
              exit = m$exit,
              M = m$M,
              doseRRmod = m$doseRRmod,
              data = data,
              deg = m$deg,
              loglim = m$loglim,
              transform = object$transform
            ),
            object$transform.args
          )
        )
      }
    }
  } else if (m$family == "poisson") {
    if (ERC) {
      # Xc_poisson and Kmat_diag_poisson captured in closure
      function(params, D) {
        do.call(
          loglik.poisson.erc,
          c(
            list(
              params = params,
              D = D,
              X = m$X,
              Y = m$Y,
              M = m$M,
              offset = m$offset,
              doseRRmod = m$doseRRmod,
              data = data,
              deg = m$deg,
              loglim = m$loglim,
              transform = object$transform,
              Xc = Xc_poisson,
              Kmat_diag = Kmat_diag_poisson
            ),
            object$transform.args
          )
        )
      }
    } else {
      function(params, D) {
        do.call(
          loglik.poisson,
          c(
            list(
              params = params,
              D = D,
              X = m$X,
              Y = m$Y,
              M = m$M,
              offset = m$offset,
              doseRRmod = m$doseRRmod,
              data = data,
              deg = m$deg,
              loglim = m$loglim,
              transform = object$transform
            ),
            object$transform.args
          )
        )
      }
    }
  } else if (m$family == "gaussian") {
    function(params, D) {
      do.call(
        loglik.gaussian,
        c(
          list(
            params = params,
            D = D,
            X = m$X,
            Y = m$Y,
            M = m$M,
            data = data,
            deg = m$deg,
            ERC = ERC,
            Kmat = Kmat,
            loglim = m$loglim,
            transform = object$transform
          ),
          object$other.args
        )
      )
    }
  } else if (m$family == "binomial") {
    function(params, D) {
      do.call(
        loglik.binomial,
        c(
          list(
            params = params,
            D = D,
            X = m$X,
            Y = m$Y,
            M = m$M,
            doseRRmod = m$doseRRmod,
            data = data,
            deg = m$deg,
            ERC = ERC,
            Kmat = Kmat,
            loglim = m$loglim,
            transform = object$transform
          ),
          object$other.args
        )
      )
    }
  } else if (m$family == "clogit") {
    designmat <- t(model.matrix(~ as.factor(m$data[, m$setnr]) - 1))
    function(params, D) {
      do.call(
        loglik.clogit,
        c(
          list(
            params = params,
            D = D,
            X = m$X,
            status = m$status,
            designmat = designmat,
            M = m$M,
            doseRRmod = m$doseRRmod,
            data = data,
            deg = m$deg,
            ERC = ERC,
            Kmat = Kmat,
            loglim = m$loglim,
            transform = object$transform
          ),
          object$other.args
        )
      )
    }
  } else if (m$family == "multinomial") {
    function(params, D) {
      do.call(
        loglik.multinomial,
        c(
          list(
            params = params,
            D = D,
            X = m$X,
            Y = m$Y,
            M = m$M,
            doseRRmod = m$doseRRmod,
            data = data,
            deg = m$deg,
            ERC = ERC,
            Kmat = Kmat,
            loglim = m$loglim,
            transform = object$transform
          ),
          object$other.args
        )
      )
    }
  } else {
    stop("Unknown family: ", m$family)
  }
}

make_base_loglik_fn_single <- function(object, method_fit, data, dose.col) {
  m <- object$model
  ERC <- FALSE # always non-ERC since we are evaluating at a single column

  if (m$family == "gaussian") {
    function(params) {
      do.call(
        loglik.gaussian,
        c(
          list(
            params = params,
            D = dose.col,
            X = m$X,
            Y = m$Y,
            M = m$M,
            data = data,
            deg = m$deg,
            ERC = ERC,
            Kmat = NULL,
            loglim = m$loglim,
            transform = object$transform
          ),
          object$transform.args
        )
      )
    }
  } else if (m$family == "binomial") {
    function(params) {
      do.call(
        loglik.binomial,
        c(
          list(
            params = params,
            D = dose.col,
            X = m$X,
            Y = m$Y,
            M = m$M,
            doseRRmod = m$doseRRmod,
            data = data,
            deg = m$deg,
            ERC = ERC,
            Kmat = NULL,
            loglim = m$loglim,
            transform = object$transform
          ),
          object$transform.args
        )
      )
    }
  } else if (m$family == "poisson") {
    function(params) {
      do.call(
        loglik.poisson,
        c(
          list(
            params = params,
            D = dose.col,
            X = m$X,
            Y = m$Y,
            M = m$M,
            offset = m$offset,
            doseRRmod = m$doseRRmod,
            data = data,
            deg = m$deg,
            loglim = m$loglim,
            transform = object$transform
          ),
          object$transform.args
        )
      )
    }
  } else if (m$family == "multinomial") {
    function(params) {
      do.call(
        loglik.multinomial,
        c(
          list(
            params = params,
            D = dose.col,
            X = m$X,
            Y = m$Y,
            M = m$M,
            doseRRmod = m$doseRRmod,
            data = data,
            deg = m$deg,
            ERC = FALSE,
            Kmat = NULL,
            loglim = m$loglim,
            transform = object$transform
          ),
          object$transform.args
        )
      )
    }
  } else if (m$family == "clogit") {
    designmat <- t(model.matrix(~ as.factor(data[, m$setnr]) - 1))
    set_members <- lapply(sort(unique(data[, m$setnr])), function(s) {
      which(data[, m$setnr] == s) - 1L
    })
    function(params) {
      do.call(
        loglik.clogit,
        c(
          list(
            params = params,
            D = dose.col,
            X = m$X,
            status = m$status,
            M = m$M,
            doseRRmod = m$doseRRmod,
            designmat = designmat,
            set_members = set_members,
            data = data,
            deg = m$deg,
            ERC = FALSE,
            Kmat = NULL,
            loglim = m$loglim,
            transform = object$transform
          ),
          object$transform.args
        )
      )
    }
  } else if (m$family == "prophaz") {
    function(params) {
      do.call(
        loglik.prophaz,
        c(
          list(
            params = params,
            D = dose.col,
            X = m$X,
            status = m$status,
            entry = m$entry,
            exit = m$exit,
            M = m$M,
            doseRRmod = m$doseRRmod,
            data = data,
            deg = m$deg,
            loglim = m$loglim,
            transform = object$transform
          ),
          object$transform.args
        )
      )
    }
  } else {
    stop("Unknown family: ", m$family)
  }
}

make_loglik_fn <- function(object, method_name, method_fit, data) {
  m <- object$model
  base_fn <- make_base_loglik_fn(object, method_fit, data)

  if (method_name == "MCML") {
    function(params) {
      logliks <- base_fn(params, D = m$dosevars)
      shifted <- pmax(pmin(-logliks - max(-logliks), 7e1), -7e1)
      -log(mean(exp(shifted))) - max(-logliks)
    }
  } else if (method_name %in% c("RC", "ERC")) {
    # Poisson ERC uses dosevars directly inside loglik.poisson.erc
    # so D argument is ignored in that case, but we pass it for consistency
    function(params) {
      base_fn(params, D = "rcdose_ameras")
    }
  } else {
    stop("make_loglik_fn not applicable for method: ", method_name)
  }
}


resolve_data <- function(object, data = NULL) {
  if (!is.null(object$model$data)) {
    return(object$model$data)
  }

  if (is.null(data)) {
    stop(
      "Data not stored on object (keep.data=FALSE). ",
      "Please supply data argument to confint()."
    )
  }

  check_df(data)

  m <- object$model

  # Check dose columns separately for a clearer error message
  missing_dose <- setdiff(m$dosevars, colnames(data))
  if (length(missing_dose)) {
    stop(
      "Supplied data is missing dose columns present during fitting: ",
      paste(missing_dose, collapse = ", ")
    )
  }

  # Check remaining required variables excluding X and dose
  # X is checked separately via X_formula, dose already checked above
  needed <- setdiff(required_vars(m), c(m$dosevars, m$X))
  missing_vars <- setdiff(needed, colnames(data))
  if (length(missing_vars)) {
    stop(
      "Supplied data is missing columns present during fitting: ",
      paste(missing_vars, collapse = ", ")
    )
  }

  # Check raw X variables referenced in the X formula
  if (!is.null(m$X_formula)) {
    X_raw_vars <- all.vars(m$X_formula)
    missing_X <- setdiff(X_raw_vars, colnames(data))
    if (length(missing_X)) {
      stop(
        "Supplied data is missing X variables present during fitting: ",
        paste(missing_X, collapse = ", ")
      )
    }
  }

  # Check dose columns are numeric and finite
  for (v in m$dosevars) {
    check_num_vec(data[, v, drop = TRUE], nm = paste0("dosevars:", v))
  }

  if (nrow(data) != object$num.rows) {
    warning(
      "Supplied data has ",
      nrow(data),
      " rows but original data had ",
      object$num.rows,
      " rows."
    )
  }

  # Re-expand X formula on supplied data
  if (!is.null(m$X_formula)) {
    X_matrix <- model.matrix(m$X_formula, data = data)[, -1, drop = FALSE]
    X_colnames <- colnames(X_matrix)
    new_cols <- setdiff(X_colnames, colnames(data))
    if (length(new_cols)) {
      data[, new_cols] <- X_matrix[, new_cols, drop = FALSE]
    }
  }

  data$rcdose_ameras <- rowMeans(data[, m$dosevars, drop = FALSE])

  data
}


check_pkgs <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Please install required packages: ", paste(missing, collapse = ", "))
  }
}


`%||%` <- function(x, y) if (is.null(x)) y else x


compute_fitted <- function(
  object,
  method = "RC",
  data,
  dose.col = "rcdose_ameras"
) {
  m <- object$model
  coefs <- object[[method]]$coefficients

  if (m$family %in% c("prophaz", "clogit")) {
    # For prophaz we return the vector relative risks
    # For clogit, the conditional probabilities
    if (!is.null(m$X)) {
      a <- coefs[1:length(m$X)]
      Xlinpred <- c(as.matrix(data[, m$X]) %*% a)
    } else {
      Xlinpred <- 0
    }

    b1 <- coefs[length(m$X) + 1]

    if (m$deg == 2) {
      b2 <- coefs[length(m$X) + 2]
    } else {
      b2 <- NULL
    }

    if (!is.null(m$M)) {
      bm1 <- coefs[
        (1 + m$deg + length(m$X)):(length(m$M) + m$deg + length(m$X))
      ]
      if (m$deg == 2) {
        bm2 <- coefs[
          (length(m$M) + 1 + m$deg + length(m$X)):(2 *
            length(m$M) +
            m$deg +
            length(m$X))
        ]
      } else {
        bm2 <- NULL
      }
    } else {
      bm1 <- bm2 <- NULL
    }

    betavec <- c(b1, b2, bm1, bm2)

    RR <- exposureRR(
      params = betavec,
      D = dose.col,
      M = m$M,
      data = data,
      doseRRmod = m$doseRRmod,
      deg = m$deg
    )
    fullRR <- drop(exp(pmin(Xlinpred, 7e1)) * RR)
    if (m$family == "prophaz") {
      return(fullRR)
    } else {
      # clogit
      probvec <- rep(NA, nrow(data))
      setnrs <- data[, m$setnr]
      for (sn in unique(setnrs)) {
        probvec[setnrs == sn] <- fullRR[setnrs == sn] /
          sum(fullRR[setnrs == sn])
      }
      return(probvec)
    }
  } else if (m$family == "multinomial") {
    # For multinomial we return the probability matrix
    Z <- nlevels(data[, m$Y])
    params <- matrix(coefs, ncol = Z - 1)
    params <- cbind(params, rep(0, nrow(params)))

    a0 <- params[1, ]

    if (!is.null(m$X)) {
      a <- params[2:(length(m$X) + 1), , drop = FALSE]
      Xlinpred <- sweep(as.matrix(data[, m$X]) %*% a, 2, a0, "+")
    } else {
      Xlinpred <- matrix(a0, nrow = nrow(data), ncol = Z, byrow = TRUE)
    }

    betamat <- params[
      (length(m$X) + 2):(m$deg * length(m$M) + m$deg + length(m$X) + 1),
      ,
      drop = FALSE
    ]

    RRmat <- matrix(1, nrow = nrow(data), ncol = Z)
    for (ii in 1:(Z - 1)) {
      RRmat[, ii] <- exposureRR(
        params = betamat[, ii],
        D = dose.col,
        M = m$M,
        data = data,
        doseRRmod = m$doseRRmod,
        deg = m$deg
      )
    }

    RRmat <- RRmat * exp(pmin(Xlinpred, 7e1))
    probmat <- RRmat / rowSums(RRmat)
    colnames(probmat) <- levels(data[, m$Y])
    return(probmat)
  } else {
    a0 <- coefs[1]

    if (!is.null(m$X)) {
      a <- coefs[2:(length(m$X) + 1)]
      Xlinpred <- c(as.matrix(data[, m$X]) %*% a)
    } else {
      Xlinpred <- 0
    }

    b1 <- coefs[length(m$X) + 2]

    if (m$deg == 2) {
      b2 <- coefs[length(m$X) + 3]
    } else {
      b2 <- NULL
    }

    if (!is.null(m$M)) {
      bm1 <- coefs[
        (2 + m$deg + length(m$X)):(length(m$M) + m$deg + length(m$X) + 1)
      ]
      if (m$deg == 2) {
        bm2 <- coefs[
          (length(m$M) + 2 + m$deg + length(m$X)):(2 *
            length(m$M) +
            m$deg +
            length(m$X) +
            1)
        ]
      } else {
        bm2 <- NULL
      }
    } else {
      bm1 <- bm2 <- NULL
    }

    betavec <- c(b1, b2, bm1, bm2)

    RR <- exposureRR(
      params = betavec,
      D = dose.col,
      M = m$M,
      data = data,
      doseRRmod = m$doseRRmod,
      deg = m$deg
    )

    if (m$family == "gaussian") {
      return(drop(a0 + Xlinpred + RR - 1))
    } else if (m$family == "binomial") {
      A <- exp(pmin(a0 + Xlinpred, 7e1)) * RR
      return(drop(A / (1 + A)))
    } else if (m$family == "poisson") {
      offs <- if (!is.null(m$offset)) data[, m$offset] else 1
      return(drop(exp(pmin(a0 + Xlinpred, 7e1)) * RR * offs))
    }
  }
}


select_dose_col <- function(object, method, data) {
  m <- object$model
  dosevars <- m$dosevars

  if (method %in% c("RC", "ERC")) {
    return("rcdose_ameras")
  } else if (method == "MCML") {
    # Realization yielding highest likelihood at final parameter estimates
    params <- object$MCML$optim$par
    logliks <- sapply(dosevars, function(dv) {
      fn <- make_base_loglik_fn_single(
        object,
        method_fit = object$MCML,
        data = data,
        dose.col = dv
      )
      fn(params)
    })
    return(dosevars[which.max(-logliks)])
  } else if (method == "FMA") {
    weights <- object$FMA$weights
    inc <- object$FMA$included.realizations

    if (is.null(weights) || is.null(inc)) {
      return("rcdose_ameras")
    }

    # Highest weight = highest likelihood since p is same for all realizations
    best_idx <- which.max(weights)
    return(dosevars[inc[best_idx]])
  } else if (method == "BMA") {
    # Realization with highest posterior probability
    # identified from col.ind column in MCMC samples
    samples <- object$BMA$samples

    if (is.list(samples)) {
      # Multiple chains: rbind all chains
      col_ind <- do.call(
        c,
        lapply(samples, function(chain) {
          chain[, "col.ind"]
        })
      )
    } else {
      # Single chain or already stacked matrix
      col_ind <- samples[, "col.ind"]
    }
    best <- as.integer(names(which.max(table(col_ind))))
    inc <- object$BMA$included.realizations
    return(dosevars[inc[best]])
  }
}


compute_schoenfeld_residuals <- function(
  exit,
  status,
  covariates,
  rr,
  entry = NULL,
  scaled = TRUE,
  vcov.mat = NULL
) {
  exit <- as.numeric(exit)
  status <- as.integer(status)
  X <- as.matrix(covariates)
  rr <- as.numeric(rr)

  n <- length(exit)
  stopifnot(length(status) == n, nrow(X) == n, length(rr) == n)

  if (is.null(entry)) {
    entry <- rep(-Inf, n)
  } else {
    entry <- as.numeric(entry)
    stopifnot(length(entry) == n)
  }

  if (any(entry > exit, na.rm = TRUE)) {
    stop("Each subject must satisfy entry <= exit.")
  }

  event_rows <- which(status == 1)
  event_times <- sort(unique(exit[event_rows]))

  p <- ncol(X)
  raw_out <- matrix(NA_real_, nrow = length(event_rows), ncol = p)
  colnames(raw_out) <- colnames(covariates)
  in_risk_set <- function(t) {
    which(entry < t & exit >= t)
  }

  for (t in event_times) {
    fail <- which(status == 1 & exit == t)
    risk <- in_risk_set(t)
    d <- length(fail)

    if (length(risk) == 0L) {
      stop(sprintf("Empty risk set at time %g.", t))
    }

    S0 <- sum(rr[risk])
    S1 <- colSums(X[risk, , drop = FALSE] * rr[risk])

    # Ties handled using Breslow approximation
    xbar <- S1 / S0
    raw_out[match(fail, event_rows), ] <-
      X[fail, , drop = FALSE] - matrix(xbar, nrow = d, ncol = p, byrow = TRUE)
  }

  scaled_out <- NULL
  if (scaled) {
    if (is.null(vcov.mat)) {
      stop("vcov.mat is required when scaled = TRUE.")
    }
    vcov.mat <- as.matrix(vcov.mat)
    if (!all(dim(vcov.mat) == c(p, p))) {
      stop(
        "vcov.mat must be a p x p covariance matrix matching ncol(covariates)."
      )
    }

    scaled_out <- raw_out %*% vcov.mat
    if (!is.null(colnames(raw_out))) {
      colnames(scaled_out) <- colnames(covariates)
    }
  }

  if (scaled) {
    data.frame(
      id = event_rows,
      time = exit[event_rows],
      scaled_out,
      row.names = NULL,
      check.names = FALSE
    )
  } else {
    data.frame(
      id = event_rows,
      time = exit[event_rows],
      raw_out,
      row.names = NULL,
      check.names = FALSE
    )
  }
}
