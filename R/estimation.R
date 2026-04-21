compute_sample_CI <- function(samples, type="percentile", level=.95){
  if(type=="percentile"){
    CIlower=apply(samples, 2, function(x) quantile(x, (1-level)/2))
    CIupper=apply(samples, 2, function(x) quantile(x, (1+level)/2))
    CIresult <- data.frame(lower=CIlower, upper=CIupper)
  } else if(type=="hpd"){
    CIresult <- as.data.frame(coda::HPDinterval(coda::as.mcmc(samples)), prob=level)
  }
  return(CIresult)
}

compute_wald_CI <- function(method_fit, level=.95, transform=NULL, other.args=NULL, type="orig"){
  
  if(type=="transformed" & is.null(transform)) stop("No transformation specified, specify transformation or choose a different CI type")
  
  z     <- qnorm(1 - (1-level)/2)
  hessian <- method_fit$optim$hessian
  
  if(type=="transformed"){
    # Check invertibility
    invertible <- det(hessian) != 0 &&
      rcond(hessian) > .Machine$double.eps &&
      all(eigen(hessian)$values > 0)
    
    if (!invertible) {
      warning("Hessian not invertible, confidence intervals could not be obtained")
      na_vec <- NA * method_fit$coefficients
      return(data.frame(lower=na_vec, upper=na_vec))
    }
    
    hess_inv <- solve(hessian)
    
    CIlower <- do.call(transform, c(list(params=method_fit$optim$par-z*sqrt(diag(hess_inv))), other.args))
    CIupper <- do.call(transform, c(list(params=method_fit$optim$par+z*sqrt(diag(hess_inv))), other.args))
  } else if(type=="orig"){
    CIlower <- method_fit$coefficients-z*sqrt(diag(method_fit$vcov))
    CIupper <- method_fit$coefficients+z*sqrt(diag(method_fit$vcov))
  }
  
  CIresult <- data.frame(lower=CIlower, upper=CIupper, row.names=names(method_fit$coefficients))
  
}

compute_proflik_CI <- function(method_fit, object, method_name, data,
                               parm="dose", level=0.95,
                               maxit.profCI=20, tol.profCI=1e-2,
                               optim.method="Nelder-Mead",
                               control=list(reltol=1e-10)) {
  
  alpha    <- 1 - level
  optval   <- -method_fit$loglik  # stored as positive log-likelihood
  inpar    <- method_fit$optim$par
  parnames <- names(method_fit$coefficients)
  
  # Determine which parameters to compute CIs for
  if (identical(parm, "dose")) {
    CIindices <- which(startsWith(parnames, "dose") |
                         grepl(")_dose", parnames))
  } else if (identical(parm, "all")) {
    CIindices <- seq_along(parnames)
  } else {
    CIindices <- which(parnames %in% parm)
    if (!length(CIindices)) {
      stop(
        "No parameters matched parm. Available parameters are: ",
        paste(parnames, collapse=", ")
      )
    }
  }
  
  # Reconstruct the likelihood function from stored model
  loglik_fn <- make_loglik_fn(object, method_name, method_fit, data)
  
  # Use stored hessian to get sensible search bounds
  hessian    <- method_fit$optim$hessian
  invertible <- det(hessian) != 0 &&
    rcond(hessian) > .Machine$double.eps &&
    all(eigen(hessian)$values > 0)
  
  if (invertible) {
    se      <- sqrt(diag(solve(hessian)))
    lowlims <- inpar - 4 * se
    uplims  <- inpar + 4 * se
  } else {
    warning(
      "Hessian not invertible, using default search bounds for profile ",
      "likelihood CI. Results may be unreliable."
    )
    lowlims <- rep(-20, length(inpar))
    uplims  <- rep(10,  length(inpar))
  }
  
  # Initialise output vectors
  CIlower    <- rep(NA_real_,    length(parnames))
  CIupper    <- rep(NA_real_,    length(parnames))
  pval_lower <- rep(NA_real_,    length(parnames))
  pval_upper <- rep(NA_real_,    length(parnames))
  iter_lower <- rep(NA_integer_, length(parnames))
  iter_upper <- rep(NA_integer_, length(parnames))
  
  for (myindex in seq_along(parnames)) {
    
    if (!(myindex %in% CIindices)) next
    
    message("Obtaining profile likelihood CI for ", parnames[myindex])
    
    profCI_one <- compute_proflik_ci_one(
      index        = myindex,
      inpar        = inpar,
      optval       = optval,
      loglik_fn    = loglik_fn,
      lowlim       = lowlims[myindex],
      uplim        = uplims[myindex],
      alpha        = alpha,
      maxit.profCI = maxit.profCI,
      tol.profCI   = tol.profCI,
      optim.method = optim.method,
      control      = control,
      parname      = parnames[myindex],
      transform    = object$transform,
      other.args   = object$other.args
    )
    
    CIlower[myindex]    <- profCI_one$lower
    CIupper[myindex]    <- profCI_one$upper
    pval_lower[myindex] <- profCI_one$pval.lower
    pval_upper[myindex] <- profCI_one$pval.upper
    iter_lower[myindex] <- profCI_one$iter.lower
    iter_upper[myindex] <- profCI_one$iter.upper
  }
  
  # Back-transform if needed
  if (!is.null(object$transform)) {
    CIlower <- do.call(object$transform, 
                       c(list(params=CIlower), object$other.args))
    CIupper <- do.call(object$transform,
                       c(list(params=CIupper), object$other.args))
  }
  
  result <- data.frame(
    lower      = CIlower,
    upper      = CIupper,
    pval.lower = pval_lower,
    pval.upper = pval_upper,
    iter.lower = iter_lower,
    iter.upper = iter_upper,
    row.names  = parnames
  )
  
  result[CIindices, , drop=FALSE]
}


compute_proflik_ci_one <- function(index, inpar, optval, loglik_fn,
                                   lowlim, uplim, alpha=0.05,
                                   maxit.profCI=20, tol.profCI=1e-2,
                                   optim.method="Nelder-Mead",
                                   control=list(reltol=1e-10),
                                   parname="parameter", transform=NULL,
                                   other.args=NULL) {
  
  profCI_fn <- function(par) {
    sapply(par, function(mypar) {
      1 - pchisq(
        2 * (proflik(parvalue      = mypar,
                     index         = index,
                     fun           = loglik_fn,
                     inpar         = inpar,
                     optim.method  = optim.method,
                     control       = control) - optval),
        df = 1
      ) - alpha
    })
  }
  
  # Upper bound
  val_at_15 <- profCI_fn(15)
  if (val_at_15 > 0) {
    if(is.null(transform)){
      warning(
        "Upper bound for ", parname, " is > 15 and may not exist. ",
        "Consider rescaling the variable."
      )
    } else{
      warning(
        "Lower bound for ", parname, " is > ",round(do.call(transform, c(list(params=rep(15, length(inpar))), other.args))[index], 1), " and may not exist. ",
        "Consider rescaling the variable."
      )
    }

    uproot <- list(root=Inf, f.root=val_at_15, iter=NA_integer_)
  } else {
    uproot <- tryCatch(
      uniroot(profCI_fn, 
              lower       = inpar[index], 
              upper       = uplim,
              extendInt   = "downX", 
              maxiter     = maxit.profCI, 
              tol         = tol.profCI),
      error = function(e) list(root=NA_real_, f.root=NA_real_, iter=NA_integer_)
    )
    if (!is.na(uproot$f.root) && abs(uproot$f.root) > 0.005) {
      warning(
        "P-value for ", parname, " upper bound is more than 0.005 from ",
        alpha, ". Consider reducing tol.profCI or increasing maxit.profCI."
      )
    }
  }
  
  # Lower bound
  val_at_minus10 <- profCI_fn(-10)
  if (val_at_minus10 > 0) {
    if(is.null(transform)){
      warning(
        "Lower bound for ", parname, " is < -10 and may not exist. ",
        "Consider rescaling the variable."
      )
    } else{
      warning(
        "Lower bound for ", parname, " is < ",round(do.call(transform, c(list(params=rep(-10, length(inpar))), other.args))[index], 1), " and may not exist. ",
        "Consider rescaling the variable or using a different transformation."
      )
    }

    lowroot <- list(root=-Inf, f.root=val_at_minus10, iter=NA_integer_)
  } else {
    lowroot <- tryCatch(
      uniroot(profCI_fn, 
              lower       = lowlim, 
              upper       = inpar[index],
              extendInt   = "upX", 
              maxiter     = maxit.profCI, 
              tol         = tol.profCI),
      error = function(e) list(root=NA_real_, f.root=NA_real_, iter=NA_integer_)
    )
    if (!is.na(lowroot$f.root) && abs(lowroot$f.root) > 0.005) {
      warning(
        "P-value for ", parname, " lower bound is more than 0.005 from ",
        alpha, ". Consider reducing tol.profCI or increasing maxit.profCI."
      )
    }
  }
  
  list(
    lower      = lowroot$root,
    upper      = uproot$root,
    pval.lower = lowroot$f.root + alpha,
    pval.upper = uproot$f.root + alpha,
    iter.lower = lowroot$iter,
    iter.upper = uproot$iter
  )
}

proflik <- memoise(function(parvalue, index, fun, inpar, optim.method=optim.method, control=list(reltol=1e-10), ...){
  
  if(length(inpar)==1){ # 1-parameter model -> just return likelihood itself
    
    return(fun(params=parvalue, ...))
    
  } else if(length(inpar)==2){ # 2-parameter model -> 1 parameter optimization for profile likelihood -> use optimize
    
    innerfun <- function(params, ...){
      if(index>1){
        params <- c(params[1:(index-1)], parvalue, params[-(1:(index-1))])
      } else{
        params <- c(parvalue, params)
      }
      
      fun(params=params, ... )
    }
    
    innerfit <- optimize(innerfun, lower=-10, upper=10, ...) 
    return(innerfit$objective)
    
  } else{ # 3+ parameter model -> 2+ parameter optimization for profile likelihood -> use optim
    
    innerfun <- function(params, ...){
      if(index>1){
        params <- c(params[1:(index-1)], parvalue, params[-(1:(index-1))])
      } else{
        params <- c(parvalue, params)
      }
      
      fun(params=params, ... )
    }
    innerfit <- optim(par=inpar[-index],fn=innerfun, method=optim.method, control=control, ...) 
    if(optim.method=="Nelder-Mead"){
      count0 <- innerfit$counts
      innerfit <- optim(par=innerfit$par,fn=innerfun, method="BFGS", control=control, ...)
      innerfit$counts <- replace(count0, is.na(count0), 0) + replace(innerfit$counts, is.na(innerfit$counts), 0)
    }
    return(innerfit$value)
  }
  
})
