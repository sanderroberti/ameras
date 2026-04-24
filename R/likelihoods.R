loglik.binomial <- function(
  params,
  Y,
  D,
  M = NULL,
  X = NULL,
  data,
  doseRRmod,
  ERC = FALSE,
  Kmat = NULL,
  deg = 1,
  loglim = 1e-30,
  transform = NULL,
  ...
) {
  # params = c(a0, a1, ..., ap, b1, b2, (bm1), (bm2))

  if (length(params) != deg * length(M) + deg + length(X) + 1) {
    stop("Parameter vector length mismatch")
  }
  if (ERC & length(D) > 1) {
    stop(
      "ERC only works with one supplied dose (i.e., the mean across replicates)"
    )
  }
  # doseRRmod = ERR or EXP
  if (!is.null(transform)) {
    if (is.function(transform)) {
      params <- transform(params = params, ...)
    } else {
      stop("transform should be a function")
    }
  }

  a0 <- params[1]

  if (!is.null(X)) {
    a <- params[2:(length(X) + 1)]
    Xlinpred <- c(as.matrix(data[, X]) %*% a)
  } else {
    Xlinpred <- 0
  }

  b1 <- params[length(X) + 2]

  if (deg == 2) {
    b2 <- params[length(X) + 3]
  } else if (deg == 1) {
    b2 <- NULL
  }

  if (!is.null(M)) {
    bm1 <- params[(2 + deg + length(X)):(length(M) + deg + length(X) + 1)]
    if (deg == 2) {
      bm2 <- params[
        (length(M) + 2 + deg + length(X)):(2 * length(M) + deg + length(X) + 1)
      ]
    } else {
      bm2 <- NULL
    }
  } else {
    bm1 <- bm2 <- NULL
  }

  if (deg == 2) {
    betavec <- c(b1, b2, bm1, bm2)
  } else {
    betavec <- c(b1, bm1)
  }

  A <- exp(pmin(a0 + Xlinpred, 7e1)) *
    exposureRR(
      params = betavec,
      D = D,
      M = M,
      data = data,
      doseRRmod = doseRRmod,
      deg = deg
    )

  if (any(A < 0)) {
    stop("RR <0, please check supplied transformation")
  }
  p <- A / (1 + A)
  p <- pmin(pmax(p, 1e-10), 1 - 1e-10)

  if (length(D) > 1) {
    ls <- colSums(log(p^data[, Y] * (1 - p)^(1 - data[, Y])))
    ls <- as.numeric(ls)
  } else {
    ls <- sum(log(p^data[, Y] * (1 - p)^(1 - data[, Y])))
  }
  if (ERC) {
    derivs <- dRRdD(
      params = betavec,
      D = D,
      M = M,
      data = data,
      doseRRmod = doseRRmod,
      deg = deg
    )
    dpdk <- exp(pmin(a0 + Xlinpred, 7e1)) * derivs$first / (1 + A)^2
    dpdk2 <- exp(pmin(a0 + Xlinpred, 7e1)) *
      derivs$second /
      (1 + A)^2 -
      2 * exp(pmin(a0 + Xlinpred, 7e1))^2 * derivs$first^2 / (1 + A)^3

    mymat <- tcrossprod((-1)^(1 - data[, Y])) *
      tcrossprod(dpdk) /
      tcrossprod(p^data[, Y] * (1 - p)^(1 - data[, Y]))

    diag(mymat) <- (-1)^(1 - data[, Y]) *
      dpdk2 /
      (p^data[, Y] * (1 - p)^(1 - data[, Y]))

    return(-1 * (ls + log(max(1 + .5 * sum(mymat * Kmat), loglim))))
  } else {
    return(-1 * ls)
  }
}


loglik.poisson <- function(
  params,
  Y,
  D,
  X = NULL,
  offset = NULL,
  M = NULL,
  doseRRmod,
  data,
  deg = 1,
  loglim = 1e-30,
  transform = NULL,
  ...
) {
  # C++ version
  # params = c(a0, a1, ..., ap, b1, b2, (bm1), (bm2))

  if (length(params) != deg * length(M) + deg + length(X) + 1) {
    stop("Parameter vector length mismatch")
  }

  if (!is.null(transform)) {
    if (is.function(transform)) {
      params <- transform(params = params, ...)
    } else {
      stop("transform should be a function")
    }
  }

  if (is.null(offset)) {
    offset <- 1
  } else {
    offset <- data[, offset]
  }

  a0 <- params[1]

  if (!is.null(X)) {
    a <- params[2:(length(X) + 1)]
    Xlinpred <- c(as.matrix(data[, X]) %*% a)
  } else {
    Xlinpred <- 0
  }

  b1 <- params[length(X) + 2]

  if (deg == 2) {
    b2 <- params[length(X) + 3]
  } else if (deg == 1) {
    b2 <- 0
  }

  if (!is.null(M)) {
    bm1 <- params[(2 + deg + length(X)):(length(M) + deg + length(X) + 1)]
    if (deg == 2) {
      bm2 <- params[
        (length(M) + 2 + deg + length(X)):(2 * length(M) + deg + length(X) + 1)
      ]
    } else {
      bm2 <- NULL
    }
  } else {
    bm1 <- bm2 <- NULL
  }

  if (deg == 2) {
    betavec <- c(b1, b2, bm1, bm2)
  } else {
    betavec <- c(b1, bm1)
  }

  mus <- pmax(
    exp(pmin(a0 + Xlinpred, 7e1)) *
      pmax(
        exposureRR(
          params = betavec,
          D = D,
          M = M,
          data = data,
          doseRRmod = doseRRmod,
          deg = deg
        ),
        loglim
      ) *
      offset,
    loglim
  )
  if (any(mus < 0)) {
    stop("RR <0, please check supplied transformation")
  }

  if (length(D) > 1) {
    #ls <- colSums(data[,Y]*log(mus)-mus-lfactorial(data[,Y]))
    ls <- colSums(
      data[, Y] *
        log(mus) -
        mus -
        ifelse(data[, Y] > 0, data[, Y] * log(data[, Y]) - data[, Y], 0)
    ) # Stirling's approximation like Mark
    ls <- as.numeric(ls)
  } else {
    #ls <- sum(data[,Y]*log(mus)-mus-lfactorial(data[,Y]))
    ls <- sum(
      data[, Y] *
        log(mus) -
        mus -
        ifelse(data[, Y] > 0, data[, Y] * log(data[, Y]) - data[, Y], 0)
    ) # Stirling's approximation like Mark
  }

  return(-1 * ls)
}


loglik.poisson.erc <- function(
  params,
  Y,
  D,
  X = NULL,
  offset = NULL,
  M = NULL,
  doseRRmod,
  data,
  deg = 1,
  loglim = 1e-30,
  transform = NULL,
  ...
) {
  # C++ version
  # params = c(a0, a1, ..., ap, b1, b2, (bm1), (bm2))

  if (length(params) != deg * length(M) + deg + length(X) + 1) {
    stop("Parameter vector length mismatch")
  }
  #if(ERC & is.null(Kmat)) stop("Kmat necessary for ERC")
  #if(ERC & length(D)>1) stop("ERC only works with one supplied dose (i.e., the mean across replicates)")

  dosemat <- as.matrix(data[, D])

  if (!is.null(transform)) {
    if (is.function(transform)) {
      params <- transform(params = params, ...)
    } else {
      stop("transform should be a function")
    }
  }

  if (is.null(offset)) {
    offset <- 1
  } else {
    offset <- data[, offset]
  }

  a0 <- params[1]

  if (!is.null(X)) {
    a <- params[2:(length(X) + 1)]
    Xlinpred <- c(as.matrix(data[, X]) %*% a)
  } else {
    Xlinpred <- 0
  }

  b1 <- params[length(X) + 2]

  if (deg == 2) {
    b2 <- params[length(X) + 3]
  } else if (deg == 1) {
    b2 <- 0
  }

  if (!is.null(M)) {
    bm1 <- params[(2 + deg + length(X)):(length(M) + deg + length(X) + 1)]
    if (deg == 2) {
      bm2 <- params[
        (length(M) + 2 + deg + length(X)):(2 * length(M) + deg + length(X) + 1)
      ]
    } else {
      bm2 <- NULL
    }
  } else {
    bm1 <- bm2 <- NULL
  }

  if (deg == 2) {
    betavec <- c(b1, b2, bm1, bm2)
  } else {
    betavec <- c(b1, bm1)
  }

  #data$rcdose_ameras <- rowMeans(dosemat)

  mus <- pmax(
    exp(pmin(a0 + Xlinpred, 7e1)) *
      pmax(
        exposureRR(
          params = betavec,
          D = "rcdose_ameras",
          M = M,
          data = data,
          doseRRmod = doseRRmod,
          deg = deg
        ),
        loglim
      ) *
      offset,
    loglim
  )
  if (any(mus < 0)) {
    stop("RR <0, please check supplied transformation")
  }

  #ls <- sum(data[,Y]*log(mus)-mus-lfactorial(data[,Y]))
  ls <- sum(
    data[, Y] *
      log(mus) -
      mus -
      ifelse(data[, Y] > 0, data[, Y] * log(data[, Y]) - data[, Y], 0)
  ) # Stirling's approximation like Mark

  derivs <- dRRdD(
    params = betavec,
    D = "rcdose_ameras",
    M = M,
    data = data,
    doseRRmod = doseRRmod,
    deg = deg
  )
  dmdd <- exp(pmin(a0 + Xlinpred, 7e1)) * derivs$first * offset
  dmdd2 <- exp(pmin(a0 + Xlinpred, 7e1)) * derivs$second * offset

  v <- (data[, Y] / mus - 1) * dmdd

  # center dose for covariance
  Xc <- dosemat - data[, "rcdose_ameras"]
  Xt_v <- crossprod(Xc, v)

  term1 <- sum(Xt_v^2) / (ncol(dosemat) - 1)

  # diagonal correction
  corrterm <- (data[, Y] / mus - 1) * dmdd2 - data[, Y] / mus^2 * dmdd^2

  # diag(Kmat) = row variances of dosemat
  row_var <- rowSums(Xc^2) / (ncol(dosemat) - 1)

  term2 <- sum(row_var * corrterm)

  val <- term1 + term2

  return(-1 * (ls + log(max(1 + .5 * val, loglim))))
}


loglik.gaussian <- function(
  params,
  D,
  Y,
  X = NULL,
  M = NULL,
  data,
  ERC = FALSE,
  Kmat = NULL,
  deg = 1,
  loglim = 1e-30,
  transform = NULL,
  ...
) {
  # params = c(a0, a1, ..., ap, b1, b2, (bm1), (bm2), sigma)

  if (length(params) != deg * length(M) + deg + length(X) + 2) {
    stop("Parameter vector length mismatch")
  }
  if (ERC & is.null(Kmat)) {
    stop("Kmat necessary for ERC")
  }
  if (ERC & length(D) > 1) {
    stop(
      "ERC only works with one supplied dose (i.e., the mean across replicates)"
    )
  }

  if (!is.null(transform)) {
    if (is.function(transform)) {
      params <- transform(params = params, ...)
    } else {
      stop("transform should be a function")
    }
  }

  a0 <- params[1]

  if (!is.null(X)) {
    a <- params[2:(length(X) + 1)]
    Xlinpred <- c(as.matrix(data[, X]) %*% a)
  } else {
    Xlinpred <- 0
  }

  b1 <- params[length(X) + 2]

  if (deg == 2) {
    b2 <- params[length(X) + 3]
  } else if (deg == 1) {
    b2 <- 0
  }

  if (!is.null(M)) {
    bm1 <- params[(2 + deg + length(X)):(length(M) + deg + length(X) + 1)]
    if (deg == 2) {
      bm2 <- params[
        (length(M) + 2 + deg + length(X)):(2 * length(M) + deg + length(X) + 1)
      ]
    } else {
      bm2 <- NULL
    }
  } else {
    bm1 <- bm2 <- NULL
  }

  if (deg == 2) {
    betavec <- c(b1, b2, bm1, bm2)
  } else {
    betavec <- c(b1, bm1)
  }

  sigma <- params[deg * length(M) + deg + length(X) + 2]

  #L <- (1/(2*pi*sigma^2))^(length(Y)/2)*exp(-sum((Y-alpha-beta1*D-beta2*D^2)^2)/(2*sigma^2))

  # Here I use mu = a0 + X^Ta + RR - 1 with RR using the linear ERR
  mus <- a0 +
    Xlinpred +
    exposureRR(
      params = betavec,
      D = D,
      M = M,
      data = data,
      doseRRmod = "ERR",
      deg = deg
    ) -
    1

  ls <- -log(2 * pi * sigma^2) *
    nrow(data) /
    2 -
    colSums((data[, Y] - as.matrix(mus, ncol = length(D)))^2) / (2 * sigma^2)

  ls <- as.numeric(ls)

  if (ERC) {
    # See above, dmu/dD = dRR/dD

    derivs <- dRRdD(
      params = betavec,
      D = D,
      M = M,
      data = data,
      doseRRmod = "ERR",
      deg = deg
    )
    dmdd <- derivs$first
    dmdd2 <- derivs$second

    mymat <- (tcrossprod(dmdd) / sigma^2) *
      (tcrossprod(data[, Y] - mus) /
        sigma^2 -
        diag(1, nrow = nrow(data), ncol = nrow(data)))
    diag(mymat) <- diag(mymat) + dmdd2 / sigma^2 * (data[, Y] - mus)

    #res <- compute_weighted_sum_kahan(dmdd, dmdd2, data[,Y]-mus, Kmat, sigma)

    return(-1 * (ls + log(max(1 + .5 * sum(mymat * Kmat), loglim))))
    #return(-1*(l+log(max(1+.5*res,loglim))))
  } else {
    return(-1 * ls)
  }
}


loglik.clogit <- function(
  params,
  D,
  status,
  X = NULL,
  M = NULL,
  doseRRmod,
  designmat,
  data,
  deg = 1,
  ERC = FALSE,
  Kmat = NULL,
  loglim = 1e-30,
  transform = NULL,
  ...
) {
  # params = c(a1, ..., ap, b1, b2, (bm1), (bm2))
  #print(params)
  if (length(params) != deg * length(M) + deg + length(X)) {
    stop("Parameter vector length mismatch")
  }
  if (ERC & is.null(Kmat)) {
    stop("Kmat necessary for ERC")
  }
  if (ERC & length(D) > 1) {
    stop(
      "ERC only works with one supplied dose (i.e., the mean across replicates)"
    )
  }

  if (!is.null(transform)) {
    if (is.function(transform)) {
      params <- transform(params = params, ...)
    } else {
      stop("transform should be a function")
    }
  }

  if (!is.null(X)) {
    a <- params[1:length(X)]
    Xlinpred <- c(as.matrix(data[, X]) %*% a)
  } else {
    Xlinpred <- 0
  }

  b1 <- params[length(X) + 1]

  if (deg == 2) {
    b2 <- params[length(X) + 2]
  } else if (deg == 1) {
    b2 <- 0
  }

  if (!is.null(M)) {
    bm1 <- params[(1 + deg + length(X)):(length(M) + deg + length(X))]
    if (deg == 2) {
      bm2 <- params[
        (length(M) + 1 + deg + length(X)):(2 * length(M) + deg + length(X))
      ]
    } else {
      bm2 <- NULL
    }
  } else {
    bm1 <- bm2 <- NULL
  }

  if (deg == 2) {
    betavec <- c(b1, b2, bm1, bm2)
  } else {
    betavec <- c(b1, bm1)
  }

  RRs <- exp(pmin(Xlinpred, 7e1)) *
    exposureRR(
      params = betavec,
      D = D,
      M = M,
      data = data,
      doseRRmod = doseRRmod,
      deg = deg
    )

  if (any(RRs < 0)) {
    stop("RR <0, please check supplied transformation")
  }
  ls <- colSums(
    log(RRs[data[, status] == 1]) -
      log(designmat %*% as.matrix(RRs, ncol = length(D)))
  )
  ls <- as.numeric(ls)

  if (ERC) {
    derivs <- dRRdD(
      params = betavec,
      D = D,
      M = M,
      data = data,
      doseRRmod = doseRRmod,
      deg = deg
    )
    drdd <- exp(pmin(Xlinpred, 7e1)) * derivs$first
    drdd2 <- exp(pmin(Xlinpred, 7e1)) * derivs$second

    tmp <- c(dldd_clogit(designmat, RRs))
    dldd <- drdd / RRs * (data[, status] == 1) - tmp * drdd

    mymat <- compute_ERCmatrix_clogit(
      designmat,
      RRs,
      drdd,
      drdd2,
      as.integer(data[, status])
    )

    return(
      -1 *
        (ls + log(max(1 + .5 * sum((tcrossprod(dldd) + mymat) * Kmat), loglim)))
    )
  } else {
    return(-1 * ls)
  }
}

loglik.prophaz <- function(
  params,
  D,
  status,
  X = NULL,
  M = NULL,
  doseRRmod,
  data,
  deg = 1,
  entry = NULL,
  exit,
  ERC = FALSE,
  Kmat = NULL,
  loglim = 1e-30,
  transform = NULL,
  ...
) {
  # params = c(a1, ..., ap, b1, b2, (bm1), (bm2))
  #print(params)
  if (length(params) != deg * length(M) + deg + length(X)) {
    stop("Parameter vector length mismatch")
  }
  if (ERC & is.null(Kmat)) {
    stop("Kmat necessary for ERC")
  }
  if (ERC & length(D) > 1) {
    stop(
      "ERC only works with one supplied dose (i.e., the mean across replicates)"
    )
  }

  if (!is.null(transform)) {
    if (is.function(transform)) {
      params <- transform(params = params, ...)
    } else {
      stop("transform should be a function")
    }
  }

  if (!is.null(X)) {
    a <- params[1:length(X)]
    Xlinpred <- c(as.matrix(data[, X]) %*% a)
  } else {
    Xlinpred <- 0
  }

  b1 <- params[length(X) + 1]

  if (deg == 2) {
    b2 <- params[length(X) + 2]
  } else if (deg == 1) {
    b2 <- 0
  }

  if (!is.null(M)) {
    bm1 <- params[(1 + deg + length(X)):(length(M) + deg + length(X))]
    if (deg == 2) {
      bm2 <- params[
        (length(M) + 1 + deg + length(X)):(2 * length(M) + deg + length(X))
      ]
    } else {
      bm2 <- NULL
    }
  } else {
    bm1 <- bm2 <- NULL
  }

  if (deg == 2) {
    betavec <- c(b1, b2, bm1, bm2)
  } else {
    betavec <- c(b1, bm1)
  }

  #RRs <- exp(pmin(Xlinpred, 7e1))*exposureRR(params=betavec, D=D, M=M, data=data, doseRRmod=doseRRmod, deg=deg)
  e1 <- exposureRR(
    params = betavec,
    D = D,
    M = M,
    data = data,
    doseRRmod = doseRRmod,
    deg = deg
  )
  e1 <- as.matrix(e1, ncol = length(D))
  if (any(e1 < 0)) {
    stop("RR <0, please check supplied transformation")
  }

  logRRs <- Xlinpred + log(e1)
  #logRRs <- logRRs - max(logRRs) # Is not compatible with ERC calculation!!

  RRs <- exp(pmin(logRRs, 7e1))

  ord_exit <- order(data[[exit]])
  exit_t <- data[[exit]][ord_exit]
  status_ord <- data[[status]][ord_exit]
  RR_exit <- RRs[ord_exit, , drop = FALSE]

  if (!is.null(entry)) {
    ord_entry <- order(data[[entry]])
    entry_t <- data[[entry]][ord_entry]
    RR_entry <- RRs[ord_entry, , drop = FALSE]
  } else {
    entry_t <- rep(min(exit_t) - 1, nrow(data))
    RR_entry <- RRs
  }

  ls <- loglik_prophaz_rcpp(
    exit_t,
    entry_t,
    RR_entry,
    RR_exit,
    status_ord,
    loglim
  )

  if (ERC) {
    derivs <- dRRdD(
      params = betavec,
      D = D,
      M = M,
      data = data,
      doseRRmod = doseRRmod,
      deg = deg
    )
    drdd <- exp(pmin(Xlinpred, 7e1)) * derivs$first
    drdd2 <- exp(pmin(Xlinpred, 7e1)) * derivs$second

    drdd <- drdd[ord_exit]
    drdd2 <- drdd2[ord_exit]

    Kmat <- Kmat[ord_exit, ord_exit]

    if (!is.null(entry)) {
      entry_t2 <- data[[entry]][ord_exit]
    } else {
      entry_t2 <- entry_t
    }
    tmp <- c(dldd_prophaz(entry_t2, exit_t, status_ord, RR_exit))
    dldd <- drdd / RR_exit * (status_ord == 1) - tmp * drdd

    mymat <- compute_ERCmatrix_prophaz(
      entry_t2,
      exit_t,
      status_ord,
      RR_exit,
      drdd,
      drdd2
    )

    return(
      -1 *
        (ls + log(max(1 + .5 * sum((tcrossprod(dldd) + mymat) * Kmat), loglim)))
    )
  } else {
    return(-1 * ls)
  }
}


loglik.multinomial <- function(
  params,
  Y,
  D,
  M = NULL,
  X = NULL,
  data,
  doseRRmod,
  ERC = FALSE,
  Kmat = NULL,
  deg = 1,
  loglim = 1e-30,
  transform = NULL,
  ...
) {
  Z <- length(unique(data[, Y])) # number of outcome categories
  # params = c(a0^(1), (a)^(1), b1^(1), b2^(1), (bm1)^(1), (bm2)^(1), ..., a0^(Z-1), (a)^(Z-1), b1^(Z-1), b2^(Z-1), (bm1)^(Z-1), (bm2)^(Z-1) )
  if (length(params) != (Z - 1) * (deg * length(M) + deg + length(X) + 1)) {
    stop("Parameter vector length mismatch")
  }
  if (ERC & length(D) > 1) {
    stop(
      "ERC only works with one supplied dose (i.e., the mean across replicates)"
    )
  }

  if (!is.null(transform)) {
    if (is.function(transform)) {
      params <- transform(params = params, ...)
    } else {
      stop("transform should be a function")
    }
  }

  params <- matrix(params, ncol = Z - 1)
  params <- cbind(params, rep(0, nrow(params)))

  a0 <- params[1, ]

  if (!is.null(X)) {
    a <- params[2:(length(X) + 1), ]
    Xlinpred <- sweep(as.matrix(data[, X]) %*% a, 2, a0, "+")
  } else {
    Xlinpred <- matrix(a0, nrow = nrow(data), ncol = Z, byrow = TRUE)
  }

  betamat <- params[
    (length(X) + 2):(deg * length(M) + deg + length(X) + 1),
    ,
    drop = FALSE
  ]

  if (length(D) > 1) {
    myarray <- array(NA_real_, dim = c(nrow(data), length(D), Z))

    for (ii in 1:(Z - 1)) {
      myarray[,, ii] <- exposureRR(
        params = betamat[, ii],
        D = D,
        M = M,
        data = data,
        doseRRmod = doseRRmod,
        deg = deg
      )
    }
    myarray[,, Z] <- 1

    if (any(myarray < 0)) {
      stop("RR <0, please check supplied transformation")
    }

    Xlinpred_array <- array(exp(pmin(Xlinpred, 7e1)), dim = c(nrow(data), 1, Z))
    Xlinpred_array <- Xlinpred_array[, rep(1, length(D)), , drop = FALSE]

    myarray <- myarray * Xlinpred_array # RR(d) * exp(X^t a)

    Asums <- colSums(aperm(myarray, c(3, 1, 2)))
    Asums_array <- array(Asums, dim = c(nrow(data), length(D), 1))
    Asums_array <- Asums_array[,, rep(1, Z), drop = FALSE]

    prob_array <- myarray / Asums_array # array with probabilities, indexed by individual x dose replicate x outcome category

    Ymat <- as.matrix(model.matrix(~ data[, Y] - 1))
    Y_array <- array(Ymat, dim = c(nrow(data), 1, Z))
    Y_array <- Y_array[, rep(1, length(D)), , drop = FALSE]

    myarray2 <- log(prob_array) * Y_array
    ls <- colSums(aperm(myarray2, c(1, 3, 2)), dims = 2)
  } else {
    RRmat <- matrix(1, nrow = nrow(data), ncol = Z)

    for (ii in 1:(Z - 1)) {
      RRmat[, ii] <- exposureRR(
        params = betamat[, ii],
        D = D,
        M = M,
        data = data,
        doseRRmod = doseRRmod,
        deg = deg
      )
    }

    if (any(RRmat < 0)) {
      stop("RR <0, please check supplied transformation")
    }

    RRmat <- RRmat * exp(pmin(Xlinpred, 7e1))
    probmat <- RRmat / rowSums(RRmat)[, drop = FALSE]

    Ymat <- diag(Z)[as.integer(data[, Y]), ] #as.matrix(model.matrix(~data[,Y]-1))
    ls <- sum(log(probmat) * Ymat)
  }

  if (ERC) {
    # Only one dose supplied for ERC

    A <- RRmat

    drdd <- apply(betamat, 2, function(betavec) {
      dRRdD(
        params = betavec,
        D = D,
        M = M,
        data = data,
        doseRRmod = doseRRmod,
        deg = deg
      )$first
    })
    drdd2 <- apply(betamat, 2, function(betavec) {
      dRRdD(
        params = betavec,
        D = D,
        M = M,
        data = data,
        doseRRmod = doseRRmod,
        deg = deg
      )$second
    })

    exp_X <- exp(pmin(Xlinpred, 7e1))
    rowSums_A <- rowSums(A)

    probvec <- rowSums(Ymat * probmat) #diag(tcrossprod(Ymat, probmat))

    drdd_expX <- drdd * exp_X
    Ymat_drdd_expX <- rowSums(Ymat * drdd_expX) #diag(tcrossprod(Ymat, drdd_expX))
    cross_drdd_exp_X <- rowSums(drdd * exp_X) #diag(tcrossprod(drdd, exp_X))

    dpdd <- (Ymat_drdd_expX - probvec * cross_drdd_exp_X) / rowSums_A
    dpdd2 <- (rowSums(Ymat * exp_X * drdd2) -
      dpdd * cross_drdd_exp_X -
      probvec * rowSums(drdd2 * exp_X)) /
      rowSums_A +
      (cross_drdd_exp_X^2 * probvec - cross_drdd_exp_X * Ymat_drdd_expX) /
        rowSums_A^2

    mymat <- tcrossprod(dpdd / probvec)
    diag(mymat) <- dpdd2 / probvec

    return(-1 * (ls + log(max(1 + .5 * sum(mymat * Kmat), loglim))))
  } else {
    return(-1 * ls)
  }
}


exposureRR <- function(params, D, M, data, doseRRmod, deg) {
  # params = c(beta1, beta2, (beta_m1), (beta_m2))
  params[is.infinite(params) & params > 0] <- 7e1
  dosemat <- as.matrix(data[, D, drop = FALSE])
  b1 <- params[1]

  if (deg == 2) {
    b2 <- params[2]
  } else {
    b2 <- 0
  }

  if (!is.null(M)) {
    bm1 <- params[(1 + deg):(length(M) + deg)]
    Mlinpred1 <- c(as.matrix(data[, M]) %*% bm1)
    if (deg == 2) {
      bm2 <- params[(length(M) + 1 + deg):(2 * length(M) + deg)]
      Mlinpred2 <- c(as.matrix(data[, M]) %*% bm2)
    } else {
      Mlinpred2 <- 0
    }
  } else {
    Mlinpred1 <- Mlinpred2 <- 0
  }

  if (doseRRmod == "ERR") {
    return(
      1 +
        b1 * dosemat +
        b2 * dosemat^2 +
        Mlinpred1 * dosemat +
        Mlinpred2 * dosemat^2
    )
  } else if (doseRRmod == "EXP") {
    val <- b1 *
      dosemat +
      b2 * dosemat^2 +
      Mlinpred1 * dosemat +
      Mlinpred2 * dosemat^2
    val <- pmin(val, 7e1) # cap large finite values
    val[!is.finite(val)] <- 7e1 # replace NaN, Inf, -Inf
    return(exp(val))
  } else if (doseRRmod == "LINEXP") {
    val <- 1 +
      (b1 + Mlinpred1) * dosemat * exp(pmin((b2 + Mlinpred2) * dosemat, 7e1))
    val <- pmin(val, exp(7e1))
    val[!is.finite(val)] <- exp(7e1) # replace Inf/NaN with capped value
    return(val)
    #return(pmin(1+(b1+Mlinpred1)*dosemat*exp(pmin((b2+Mlinpred2)*dosemat, 8e1)), exp(8e1)))
  }
}


dRRdD <- function(params, D, M, data, doseRRmod, deg) {
  # params = c(beta1, beta2, (beta_m1), (beta_m2))
  params[is.infinite(params) & params > 0] <- 7e1
  b1 <- params[1]

  if (deg == 2) {
    b2 <- params[2]
  } else {
    b2 <- 0
  }

  if (!is.null(M)) {
    bm1 <- params[(1 + deg):(length(M) + deg)]
    Mlinpred1 <- c(as.matrix(data[, M]) %*% bm1)
    if (deg == 2) {
      bm2 <- params[(length(M) + 1 + deg):(2 * length(M) + deg)]
      Mlinpred2 <- c(as.matrix(data[, M]) %*% bm2)
    } else {
      Mlinpred2 <- 0
    }
  } else {
    Mlinpred1 <- Mlinpred2 <- 0
  }

  # First derivative
  if (doseRRmod == "ERR") {
    first <- b1 + 2 * b2 * data[, D] + Mlinpred1 + 2 * Mlinpred2 * data[, D]
  } else if (doseRRmod == "EXP") {
    first <- (b1 + 2 * b2 * data[, D] + Mlinpred1 + 2 * Mlinpred2 * data[, D]) *
      exp(pmin(
        b1 *
          data[, D] +
          b2 * data[, D]^2 +
          Mlinpred1 * data[, D] +
          Mlinpred2 * data[, D]^2,
        8e1
      ))
  } else if (doseRRmod == "LINEXP") {
    first <- (b1 + Mlinpred1) *
      exp((b2 + Mlinpred2) * data[, D]) +
      (b2 + Mlinpred2) *
        (b1 + Mlinpred1) *
        data[, D] *
        exp((b2 + Mlinpred2) * data[, D])
  }

  # Second derivative
  if (doseRRmod == "ERR") {
    second <- 2 * (b2 + Mlinpred2)

    if (length(second) == 1) second <- rep(second, length(first))
  } else if (doseRRmod == "EXP") {
    second <- 2 *
      (b2 + Mlinpred2) *
      exp(pmin(
        b1 *
          data[, D] +
          b2 * data[, D]^2 +
          Mlinpred1 * data[, D] +
          Mlinpred2 * data[, D]^2,
        8e1
      )) +
      (b1 + 2 * b2 * data[, D] + Mlinpred1 + 2 * Mlinpred2 * data[, D])^2 *
        exp(pmin(
          b1 *
            data[, D] +
            b2 * data[, D]^2 +
            Mlinpred1 * data[, D] +
            Mlinpred2 * data[, D]^2,
          8e1
        ))

    if (length(second) == 1) second <- rep(second, length(first))
  } else if (doseRRmod == "LINEXP") {
    second <- (b1 + Mlinpred1) *
      (b2 + Mlinpred2) *
      exp((b2 + Mlinpred2) * data[, D]) +
      (b2 + Mlinpred2) * first

    if (length(second) == 1) second <- rep(second, length(first))
  }

  return(list(first = first, second = second))
}
