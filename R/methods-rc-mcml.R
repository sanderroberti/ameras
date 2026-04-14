
ameras.mcml <- function(family, dosevars, data, deg, transform=NULL,transform.jacobian=NULL, Y=NULL, M=NULL, X=NULL, offset=NULL, inpar=NULL, entry=NULL, exit=NULL, status=NULL,setnr=setnr, CI=NULL,params.profCI=NULL, maxit.profCI=NULL,tol.profCI=NULL, doseRRmod=NULL, loglim=1e-30, optim.method="Nelder-Mead", ...){
  if(length(CI)==0) stop("No CI method specified")
  CI <- CI[CI %in% c("wald.orig","wald.transformed", "proflik")]
  if(length(CI)>1) stop("Provide one type of CI for MCML: one of wald.orig, wald.transformed, and proflik")
  if(length(CI)==0) stop("Incorrect CI method specified, should be one of wald.orig, wald.transformed, and proflik")
  if(CI=="proflik" & !(params.profCI %in% c("dose","all"))) stop("Incorrect choice of parameters for profile likelihood CI supplied, should be either all or dose")
  if(CI=="wald.transformed" & is.null(transform)) stop("No transformation specified, specify transformation or choose a different CI type")
  
  if(CI=="proflik") message("Note: computation times for profile likelihood intervals for MCML may be extensive with large datasets or complex models")
  
  if(family=="gaussian"){
    if(is.null(Y)) stop("Y is required for family=gaussian")
    
    if(is.null(inpar)){
      inpar <- rep(0, 2+length(X)+length(M)*deg+deg)
    }
    
    loglik.mcml <- function(params, ...){
      logliks <- loglik.gaussian(params, D=dosevars,X=X, Y=Y, M=M, data=data, deg=deg, ERC=FALSE, loglim=loglim, transform=transform, ...)
      return(-1*log(mean(exp(pmax(pmin(-1*logliks-max(-1*logliks), 7e1), -7e1))))-max(-1*logliks)) # with explim=7e1
      #return(mean(logliks)) # mean of log likelihoods
    }
    
    parnames <- c("(Intercept)",names(data[,X,drop=FALSE]),c("dose","dose_squared")[1:deg])
    if(!is.null(M)){
      parnames <- c(parnames, paste0("dose:",names(data[, M,drop=FALSE])))
      if(deg==2){
        parnames <- c(parnames, paste0("dose_squared:",names(data[, M,drop=FALSE])))
      }
    }
    parnames <- c(parnames, "sigma")
    
  } else if(family=="binomial"){
    
    if(is.null(Y)) stop("Y is required for family=binomial")
    
    if(is.null(inpar)){
      inpar <- rep(0, 1+length(X)+length(M)*deg+deg)
    }
    
    loglik.mcml <- function(params, ...){
      logliks <- loglik.binomial(params, D=dosevars,X=X, Y=Y, M=M,doseRRmod=doseRRmod, data=data, deg=deg, ERC=FALSE, loglim=loglim, transform=transform, ...)
      #return(log(mean(exp(logliks-mean(logliks))))+mean(logliks)) # log of mean of likelihoods
      return(-1*log(mean(exp(pmax(pmin(-1*logliks-max(-1*logliks), 7e1), -7e1))))-max(-1*logliks)) # with explim=7e1
      #return(mean(logliks)) # mean of log likelihoods
    }
    
    
    if(doseRRmod!="LINEXP"){
      parnames <- c("(Intercept)",names(data[,X,drop=FALSE]),c("dose","dose_squared")[1:deg])
      if(!is.null(M)){
        parnames <- c(parnames, paste0("dose:",names(data[, M,drop=FALSE])))
        if(deg==2){
          parnames <- c(parnames, paste0("dose_squared:",names(data[, M,drop=FALSE])))
        }
      }
    } else{
      parnames <- c("(Intercept)",names(data[,X,drop=FALSE]),c("dose_linear","dose_exponential"))
      if(!is.null(M)){
        parnames <- c(parnames, paste0("dose_linear:",names(data[, M,drop=FALSE])))
        parnames <- c(parnames, paste0("dose_exponential:",names(data[, M,drop=FALSE])))
      }
    }
    
    
    
  } else if(family=="poisson"){
    if(is.null(Y)) stop("Y is required for family=poisson")
    
    if(is.null(inpar)){
      inpar <- rep(0, 1+length(X)+length(M)*deg+deg)
    }
    
    loglik.mcml <- function(params, ...){
      logliks <- loglik.poisson(params, D=dosevars,X=X, Y=Y,offset=offset, M=M,doseRRmod=doseRRmod, data=data, deg=deg, loglim=loglim, transform=transform, ...)
      #return(log(mean(exp(logliks-mean(logliks))))+mean(logliks)) # log of mean of likelihoods
      return(-1*log(mean(exp(pmax(pmin(-1*logliks-max(-1*logliks), 7e1), -7e1))))-max(-1*logliks)) # with explim=7e1
      #return(mean(logliks)) # mean of log likelihoods
    }
    
    if(doseRRmod!="LINEXP"){
      parnames <- c("(Intercept)",names(data[,X,drop=FALSE]),c("dose","dose_squared")[1:deg])
      if(!is.null(M)){
        parnames <- c(parnames, paste0("dose:",names(data[, M,drop=FALSE])))
        if(deg==2){
          parnames <- c(parnames, paste0("dose_squared:",names(data[, M,drop=FALSE])))
        }
      }
    } else{
      parnames <- c("(Intercept)",names(data[,X,drop=FALSE]),c("dose_linear","dose_exponential"))
      if(!is.null(M)){
        parnames <- c(parnames, paste0("dose_linear:",names(data[, M,drop=FALSE])))
        parnames <- c(parnames, paste0("dose_exponential:",names(data[, M,drop=FALSE])))
      }
    }
    
    
  } else if(family=="clogit"){
    
    
    if(is.null(doseRRmod)) stop("doseRRmod is required for family=prophaz")
    if(is.null(status)) stop("status is required for family=prophaz")
    
    designmat <- t(model.matrix(~as.factor(data[,setnr])-1))
    
    
    
    if(is.null(inpar)){
      inpar <- rep(0, length(X)+length(M)*deg+deg)
    }
    
    loglik.mcml <- function(params, ...){
      logliks <- loglik.clogit(params, D=dosevars,status=status, X=X, M=M,doseRRmod=doseRRmod, designmat=designmat, data=data, deg=deg, ERC=FALSE, loglim=loglim, transform=transform,...)
      return(-1*log(mean(exp(pmax(pmin(-1*logliks-max(-1*logliks), 7e1), -7e1))))-max(-1*logliks)) # with explim=7e1
      #return(mean(logliks)) # mean of log likelihoods
    }
    
    if(doseRRmod!="LINEXP"){
      parnames <- c(names(data[,X,drop=FALSE]),c("dose","dose_squared")[1:deg])
      if(!is.null(M)){
        parnames <- c(parnames, paste0("dose:",names(data[, M,drop=FALSE])))
        if(deg==2){
          parnames <- c(parnames, paste0("dose_squared:",names(data[, M,drop=FALSE])))
        }
      }
    } else{
      parnames <- c(names(data[,X,drop=FALSE]),c("dose_linear","dose_exponential"))
      if(!is.null(M)){
        parnames <- c(parnames, paste0("dose_linear:",names(data[, M,drop=FALSE])))
        parnames <- c(parnames, paste0("dose_exponential:",names(data[, M,drop=FALSE])))
      }
    }
    
    
  } else if(family=="prophaz"){
    
    if(is.null(exit)) stop("exit is required for family=prophaz")
    if(is.null(doseRRmod)) stop("doseRRmod is required for family=prophaz")
    if(is.null(status)) stop("status is required for family=prophaz")
    
    
    if(is.null(inpar)){
      inpar <- rep(0, length(X)+length(M)*deg+deg)
    }
    
    loglik.mcml <- function(params, ...){
      logliks <- loglik.prophaz(params, D=dosevars,status=status, X=X, M=M,doseRRmod=doseRRmod, entry=entry, exit=exit, data=data, deg=deg, ERC=FALSE, loglim=loglim, transform=transform,...)
      return(-1*log(mean(exp(pmax(pmin(-1*logliks-max(-1*logliks), 7e1), -7e1))))-max(-1*logliks)) # with explim=7e1
      #return(mean(logliks)) # mean of log likelihoods
    }
    
    if(doseRRmod!="LINEXP"){
      parnames <- c(names(data[,X,drop=FALSE]),c("dose","dose_squared")[1:deg])
      if(!is.null(M)){
        parnames <- c(parnames, paste0("dose:",names(data[, M,drop=FALSE])))
        if(deg==2){
          parnames <- c(parnames, paste0("dose_squared:",names(data[, M,drop=FALSE])))
        }
      }
    } else{
      parnames <- c(names(data[,X,drop=FALSE]),c("dose_linear","dose_exponential"))
      if(!is.null(M)){
        parnames <- c(parnames, paste0("dose_linear:",names(data[, M,drop=FALSE])))
        parnames <- c(parnames, paste0("dose_exponential:",names(data[, M,drop=FALSE])))
      }
    }
    
    
  } else if(family=="multinomial"){
    
    if(is.null(Y)) stop("Y is required for family=multinomial")
    
    if(is.null(inpar)){
      inpar <- rep(0, (length(unique(data[,Y]))-1)*(1+length(X)+length(M)*deg+deg))
    }
    
    loglik.mcml <- function(params, ...){
      logliks <- loglik.multinomial(params, D=dosevars,X=X, Y=Y, M=M,doseRRmod=doseRRmod, data=data, deg=deg, ERC=FALSE, loglim=loglim, transform=transform, ...)
      #return(log(mean(exp(logliks-mean(logliks))))+mean(logliks)) # log of mean of likelihoods
      return(-1*log(mean(exp(pmax(pmin(-1*logliks-max(-1*logliks), 7e1), -7e1))))-max(-1*logliks)) # with explim=7e1
      #return(mean(logliks)) # mean of log likelihoods
    }
    
    
    if(doseRRmod!="LINEXP"){
      parnames <- c("(Intercept)",names(data[,X,drop=FALSE]),c("dose","dose_squared")[1:deg])
      if(!is.null(M)){
        parnames <- c(parnames, paste0("dose:",names(data[, M,drop=FALSE])))
        if(deg==2){
          parnames <- c(parnames, paste0("dose_squared:",names(data[, M,drop=FALSE])))
        }
      }
    } else{
      parnames <- c("(Intercept)",names(data[,X,drop=FALSE]),c("dose_linear","dose_exponential"))
      if(!is.null(M)){
        parnames <- c(parnames, paste0("dose_linear:",names(data[, M,drop=FALSE])))
        parnames <- c(parnames, paste0("dose_exponential:",names(data[, M,drop=FALSE])))
      }
    }
    
    mylv <- levels(data[,Y])
    
    mylv <- mylv[-length(mylv)]
    
    parnames <- do.call("c",lapply(mylv, function(lv) paste0("(",lv, ")_", parnames)))
    
  }
  t0 <- proc.time()
  if(length(parnames)==1){ # Optimize 1-dimensional model: use optimize instead of optim
    fit0 <- optimize(f=loglik.mcml, lower=-20, upper=5, ...)
    fit <- list(par=fit0$minimum, value=fit0$objective, convergence=0, hessian=numDeriv::hessian(func=loglik.mcml, x=fit0$minimum, ...))
  } else{
    fit <- optim(inpar, loglik.mcml, method=optim.method, ...)
    if(optim.method=="Nelder-Mead"){
      count0 <- fit$counts
      fit <- optim(fit$par, loglik.mcml, method="BFGS", ...)
      fit$counts <- replace(count0, is.na(count0), 0) + replace(fit$counts, is.na(fit$counts), 0)
    }
    fit$hessian <- numDeriv::hessian(func=loglik.mcml, x=fit$par, ...)
  }
  
  if(!is.null(transform) & !is.null(transform.jacobian)){
    if(is.function(transform) & is.function(transform.jacobian)){
      if("boundcheck" %in% names(formals(transform))){
        coefs <- transform(fit$par, boundcheck=TRUE, ...)
      } else {
        coefs <- transform(fit$par, ...)
      }
      if(det(fit$hessian)!=0 & rcond(fit$hessian)>.Machine$double.eps &  all(eigen(fit$hessian)$values > 0)){
        jac <- transform.jacobian(fit$par, ...)
        vcov <- jac %*% solve(fit$hessian) %*% t(jac)
      } else{
        warning("WARNING: Hessian was not invertible or inverse was not positive definite, variance matrix could not be obtained")
        vcov <- matrix(NA, ncol=length(parnames), nrow=length(parnames))
      }
    } else{
      stop("transform and transform.jacobian should be functions")
    } 
  } else{
    coefs <- fit$par
    if(det(fit$hessian)!=0 & rcond(fit$hessian)>.Machine$double.eps &  all(eigen(fit$hessian)$values > 0)){
      vcov <- solve(fit$hessian)
    } else{
      warning("WARNING: Hessian was not invertible or inverse was not positive definite, variance matrix could not be obtained")
      vcov <- matrix(NA, ncol=length(parnames), nrow=length(parnames))
    }
  }
  
  names(coefs) <- parnames
  rownames(vcov) <- colnames(vcov) <- parnames
  
  if(CI=="wald.transformed"){
    if(det(fit$hessian)!=0 & rcond(fit$hessian)>.Machine$double.eps &  all(eigen(fit$hessian)$values > 0)){
      CIlower <- transform(fit$par-1.96*sqrt(diag(solve(fit$hessian))), ...)
      CIupper <- transform(fit$par+1.96*sqrt(diag(solve(fit$hessian))), ...)
    } else{
      warning("WARNING: Hessian was not invertible or inverse was not positive definite, confidence intervals could not be obtained")
      CIlower <- CIupper <- NA * fit$par
    }
  } else if(CI=="wald.orig"){
    CIlower <- coefs-1.96*sqrt(diag(vcov))
    CIupper <- coefs+1.96*sqrt(diag(vcov))
  } else if(CI=="proflik"){
    if(det(fit$hessian)!=0 & rcond(fit$hessian)>.Machine$double.eps &  all(eigen(fit$hessian)$values > 0)){
      lowlims <- fit$par-4*sqrt(diag(solve(fit$hessian)))
      uplims <- fit$par+4*sqrt(diag(solve(fit$hessian)))
    } else {
      lowlims <- rep(-20, length(coefs))
      uplims <- rep(10, length(coefs))
    }
    
    CIlower <- rep(NA, length(coefs))
    CIupper <- rep(NA, length(coefs))
    
    pval_lower <- pval_upper <- rep(NA, length(coefs))
    iter_lower <- iter_upper <- rep(NA, length(coefs))
    
    if(params.profCI=="all"){
      CIindices <- 1:length(parnames)
    } else if(params.profCI=="dose"){
      CIindices <- (1:length(parnames))[startsWith(parnames, "dose") | grepl(")_dose", parnames)]
    }
    
    if(is.null(transform)){
      p15 <- rep(15, length(parnames))
      pminus10 <- rep(-10, length(parnames))
    } else{
      p15 <- transform(rep(15, length(parnames)), ...)
      pminus10 <- transform(rep(-10, length(parnames)), ...)
    }
    
    for(myindex in 1:length(parnames)){
      if(myindex %in% CIindices){  # Implemented this way rather than a loop over a subset of parameters so that the full covariate vector length is maintained, keeping transformation functionality intact
        message(paste0("Obtaining profile likelihood CI for ", parnames[myindex]))
        
        val15 <- profCIfun(par=15,optval=fit$value, index=myindex, fun=loglik.mcml, inpar=fit$par, optim.method=optim.method, ...)
        if(val15 > 0){
          if(is.null(transform)){
            warning(paste0("WARNING: Upper bound for ", parnames[myindex], " is >15 and may not exist if rescaling the variable does not help"))
          } else{
            warning(paste0("WARNING: Upper bound for ", parnames[myindex], " is >",round(p15[myindex],1)," and may not exist if rescaling the variable does not help"))
          }
          uproot <- list(root=Inf, f.root=val15, iter=NA)
        } else{
          uproot <- tryCatch(uniroot(profCIfun, lower=fit$par[myindex], upper=uplims[myindex], optval=fit$value, index=myindex, fun=loglik.mcml, inpar=fit$par, optim.method=optim.method,extendInt="downX", maxiter=maxit.profCI, tol=tol.profCI, ...), error=function(e) list(root=NA, f.root=NA, iter=NA))
          if(!is.na(uproot$f.root)){
            if(abs(uproot$f.root)>.005) warning(paste0("P-value for ", parnames[myindex]," upper bound more than 0.005 away from 0.05, reducing tol.profCI and/or increasing maxit.profCI is recommended"))
          }
        }
        
        valminus10 <- profCIfun(par=-10,optval=fit$value, index=myindex, fun=loglik.mcml, inpar=fit$par, optim.method=optim.method, ...)
        if(valminus10 > 0){
          if(is.null(transform)){
            warning(paste0("WARNING: Lower bound for ", parnames[myindex], " is < -10 and may not exist if rescaling the variable does not help"))
          } else{
            warning(paste0("WARNING: Lower bound for ", parnames[myindex], " is < ",round(pminus10[myindex],1)," and may not exist if rescaling the variable does not help"))
          }
          lowroot <- list(root=-Inf, f.root=valminus10, iter=NA)
        } else{
          lowroot <- tryCatch(uniroot(profCIfun, lower=lowlims[myindex], upper=fit$par[myindex], optval=fit$value, index=myindex, fun=loglik.mcml, optim.method=optim.method, inpar=fit$par, extendInt="upX", maxiter=maxit.profCI,tol=tol.profCI, ...), error=function(e) list(root=NA, f.root=NA, iter=NA))
          if(!is.na(lowroot$f.root)){
            if(abs(lowroot$f.root)>.005) warning(paste0("P-value for ", parnames[myindex]," lower bound more than 0.005 away from 0.05, reducing tol.profCI and/or increasing maxit.profCI is recommended"))
          }
        }
        
        CIlower[myindex] <- lowroot$root
        CIupper[myindex] <- uproot$root
        
        pval_lower[myindex] <- lowroot$f.root+.05
        pval_upper[myindex] <- uproot$f.root+.05
        
        iter_lower[myindex] <- lowroot$iter
        iter_upper[myindex] <- uproot$iter
      } else{
        CIlower[myindex] <- CIupper[myindex] <- pval_lower[myindex] <- pval_upper[myindex] <- iter_lower[myindex] <- iter_upper[myindex] <- NA
      }
    }
    
    if(!is.null(transform)){
      CIlower <- transform(CIlower, ...)
      CIupper <- transform(CIupper, ...)
    }
    
    profCI_pval <- data.frame(pval.lower=pval_lower, pval.upper=pval_upper)
    profCI_iter <- data.frame(iter.lower=iter_lower, iter.upper=iter_upper)
  }
  
  CIresult <- data.frame(lower=CIlower, upper=CIupper)
  if(CI=="proflik"){
    CIresult <- cbind(CIresult, profCI_pval, profCI_iter)
    rownames(CIresult) <- parnames
    CIresult <- CIresult[CIindices,]
  } else{
    rownames(CIresult) <- parnames
  }
  
  
  t1 <- proc.time()
  timedif <- t1-t0
  runtime <- paste(round(as.numeric(as.difftime(timedif["elapsed"], units="secs")),1), "seconds")
  
  
  out <- list(coefficients=coefs,
              sd=sqrt(diag(vcov)),
              CI=CIresult, 
              vcov=vcov, 
              convergence.optim=fit$convergence, 
              counts.optim=fit$counts,
              loglik=-1*fit$value,
              runtime=runtime)
  return(out)
}

ameras.rc <- function(family, dosevars, data, deg, ERC=FALSE, transform=NULL, transform.jacobian=NULL, Y=NULL, M=NULL, X=NULL, offset=NULL, inpar=NULL, entry=NULL, exit=NULL, status=NULL, setnr=NULL, CI=NULL, params.profCI=NULL, maxit.profCI=NULL, tol.profCI=NULL, doseRRmod=NULL, loglim=1e-30, optim.method="Nelder-Mead", ...){
  if(length(CI)==0) stop("No CI method specified")
  CI <- CI[CI %in% c("wald.orig","wald.transformed", "proflik")]
  if(length(CI)>1) stop("Provide one type of CI for (E)RC: one of wald.orig, wald.transformed, and proflik")
  if(length(CI)==0) stop("Incorrect CI method specified, should be one of wald.orig, wald.transformed, and proflik") 
  if(CI=="proflik" & !(params.profCI %in% c("dose","all"))) stop("Incorrect choice of parameters for profile likelihood CI supplied, should be either all or dose")
  if(CI=="wald.transformed" & is.null(transform)) stop("No transformation specified, specify transformation or choose a different CI type")
  
  if(CI=="proflik" & ERC==TRUE) message("Note: computation times for profile likelihood intervals for ERC may be extensive with large datasets or complex models")
  
  #data.rc <- data#[,-dosevars]
  data$rcdose_ameras <- rowMeans(data[,dosevars, drop=FALSE])
  
  if(ERC & family!="poisson"){
    Kmat <- cov(t(data[,dosevars]))
  } else{
    Kmat <- NULL
  }
  
  t0 <- proc.time()
  if(family=="gaussian"){
    if(is.null(Y)) stop("Y is required for family=gaussian")
    
    if(is.null(inpar)){
      inpar <- rep(0, 2+length(X)+length(M)*deg+deg)
    }
    
    fit <- optim(inpar, loglik.gaussian, D="rcdose_ameras",X=X, Y=Y, M=M, data=data, deg=deg, ERC=ERC, Kmat=Kmat, loglim=loglim, transform=transform, method=optim.method, ...)
    if(optim.method=="Nelder-Mead"){
      count0 <- fit$counts
      fit <- optim(fit$par, loglik.gaussian, D="rcdose_ameras",X=X, Y=Y, M=M, data=data, deg=deg, ERC=ERC, Kmat=Kmat, loglim=loglim, transform=transform, method="BFGS", ...)
      fit$counts <- replace(count0, is.na(count0), 0) + replace(fit$counts, is.na(fit$counts), 0)
    }
    fit$hessian <- numDeriv::hessian(func=loglik.gaussian, x=fit$par, D="rcdose_ameras",X=X, Y=Y, M=M, data=data, deg=deg, ERC=ERC, Kmat=Kmat, loglim=loglim, transform=transform, ...)
    
    parnames <- c("(Intercept)",names(data[,X,drop=FALSE]),c("dose","dose_squared")[1:deg])
    if(!is.null(M)){
      parnames <- c(parnames, paste0("dose:",names(data[, M,drop=FALSE])))
      
      if(deg==2){
        parnames <- c(parnames, paste0("dose_squared:",names(data[, M,drop=FALSE])))
      }
    }
    parnames <- c(parnames, "sigma")
    
    
    
  } else if(family=="binomial"){
    
    if(is.null(Y)) stop("Y is required for family=binomial")
    
    if(is.null(inpar)){
      inpar <- rep(0, 1+length(X)+length(M)*deg+deg)
    }
    
    fit <- optim(inpar, loglik.binomial, D="rcdose_ameras",X=X, Y=Y, M=M, doseRRmod=doseRRmod, data=data, deg=deg, ERC=ERC, Kmat=Kmat, loglim=loglim, transform=transform, method=optim.method, ...)
    if(optim.method=="Nelder-Mead"){
      count0 <- fit$counts
      fit <- optim(fit$par, loglik.binomial, D="rcdose_ameras",X=X, Y=Y, M=M, doseRRmod=doseRRmod, data=data, deg=deg, ERC=ERC, Kmat=Kmat, loglim=loglim, transform=transform, method="BFGS", ...)
      fit$counts <- replace(count0, is.na(count0), 0) + replace(fit$counts, is.na(fit$counts), 0)
    }
    
    fit$hessian <- numDeriv::hessian(func=loglik.binomial, x=fit$par, D="rcdose_ameras",X=X, Y=Y, M=M, doseRRmod=doseRRmod, data=data, deg=deg, ERC=ERC, Kmat=Kmat, loglim=loglim, transform=transform, ...)
    
    
    if(doseRRmod!="LINEXP"){
      parnames <- c("(Intercept)",names(data[,X,drop=FALSE]),c("dose","dose_squared")[1:deg])
      if(!is.null(M)){
        parnames <- c(parnames, paste0("dose:",names(data[, M,drop=FALSE])))
        if(deg==2){
          parnames <- c(parnames, paste0("dose_squared:",names(data[, M,drop=FALSE])))
        }
      }
    } else{
      parnames <- c("(Intercept)",names(data[,X,drop=FALSE]),c("dose_linear","dose_exponential"))
      if(!is.null(M)){
        parnames <- c(parnames, paste0("dose_linear:",names(data[, M,drop=FALSE])))
        parnames <- c(parnames, paste0("dose_exponential:",names(data[, M,drop=FALSE])))
      }
    }
    
    
    
  } else if(family=="poisson"){
    
    
    if(is.null(Y)) stop("Y is required for family=poisson")
    
    if(is.null(inpar)){
      inpar <- rep(0, 1+length(X)+length(M)*deg+deg)
    }
    if(ERC){
      fit <- optim(inpar, loglik.poisson.erc, D=dosevars,X=X, Y=Y, offset=offset, M=M, doseRRmod=doseRRmod, data=data, deg=deg, loglim=loglim, transform=transform, method=optim.method, ...)
    } else{
      fit <- optim(inpar, loglik.poisson, D="rcdose_ameras",X=X, Y=Y, offset=offset, M=M, doseRRmod=doseRRmod, data=data, deg=deg, loglim=loglim, transform=transform, method=optim.method, ...)
    }
    if(optim.method=="Nelder-Mead"){
      count0 <- fit$counts
      if(ERC){
        fit <- optim(fit$par, loglik.poisson.erc, D=dosevars,X=X, Y=Y, offset=offset, M=M, doseRRmod=doseRRmod, data=data, deg=deg, loglim=loglim, transform=transform, method="BFGS", ...)
      } else {
        fit <- optim(fit$par, loglik.poisson, D="rcdose_ameras",X=X, Y=Y, offset=offset, M=M, doseRRmod=doseRRmod, data=data, deg=deg, loglim=loglim, transform=transform, method="BFGS", ...)
      }
      fit$counts <- replace(count0, is.na(count0), 0) + replace(fit$counts, is.na(fit$counts), 0)
    }
    
    if(ERC){
      fit$hessian <- numDeriv::hessian(func=loglik.poisson.erc, x=fit$par, D=dosevars,X=X, Y=Y, offset=offset, M=M, doseRRmod=doseRRmod, data=data, deg=deg, loglim=loglim, transform=transform, ...)
    } else {
      fit$hessian <- numDeriv::hessian(func=loglik.poisson, x=fit$par, D="rcdose_ameras",X=X, Y=Y, offset=offset, M=M, doseRRmod=doseRRmod, data=data, deg=deg, loglim=loglim, transform=transform, ...)
    }
    if(doseRRmod!="LINEXP"){
      parnames <- c("(Intercept)",names(data[,X,drop=FALSE]),c("dose","dose_squared")[1:deg])
      if(!is.null(M)){
        parnames <- c(parnames, paste0("dose:",names(data[, M,drop=FALSE])))
        if(deg==2){
          parnames <- c(parnames, paste0("dose_squared:",names(data[, M,drop=FALSE])))
        }
      }
    } else{
      parnames <- c("(Intercept)",names(data[,X,drop=FALSE]),c("dose_linear","dose_exponential"))
      if(!is.null(M)){
        parnames <- c(parnames, paste0("dose_linear:",names(data[, M,drop=FALSE])))
        parnames <- c(parnames, paste0("dose_exponential:",names(data[, M,drop=FALSE])))
      }
    }
    
  } else if(family %in% c("clogit")){
    if(is.null(doseRRmod)) stop("doseRRmod is required for family=prophaz")
    if(is.null(status)) stop("status is required for family=prophaz")
    
    designmat <- t(model.matrix(~as.factor(data[,setnr])-1))
    
    if(is.null(inpar)){
      inpar <- rep(0, length(X)+length(M)*deg+deg)
    }
    
    if(length(X)+length(M)*deg+deg == 1){ # Optimize 1-dimensional model: use optimize instead of optim
      
      fit0 <- optimize(f=loglik.clogit, lower=-20, upper=5, D="rcdose_ameras", status=status,X=X, M=M, doseRRmod=doseRRmod, designmat=designmat,entry=entry,exit=exit, data=data, deg=deg, ERC=ERC, Kmat=Kmat, loglim=loglim, transform=transform, ...)
      
      
      fit <- list(par=fit0$minimum, value=fit0$objective, convergence=0, hessian=numDeriv::hessian(func=loglik.clogit, x=fit0$minimum, D="rcdose_ameras", status=status,X=X, M=M, doseRRmod=doseRRmod, designmat=designmat,entry=entry,exit=exit, data=data, deg=deg, ERC=ERC, Kmat=Kmat, loglim=loglim, transform=transform, ...))
    } else {
      
      fit <- optim(inpar, loglik.clogit, D="rcdose_ameras", status=status,X=X, M=M, doseRRmod=doseRRmod, designmat=designmat,entry=entry,exit=exit, data=data, deg=deg, ERC=ERC, Kmat=Kmat, loglim=loglim, transform=transform, method=optim.method, ...)
      
      if(optim.method=="Nelder-Mead"){
        count0 <- fit$counts
        fit <- optim(fit$par, loglik.clogit, D="rcdose_ameras", status=status,X=X, M=M, doseRRmod=doseRRmod, designmat=designmat,entry=entry,exit=exit, data=data, deg=deg, ERC=ERC, Kmat=Kmat, loglim=loglim, transform=transform, method="BFGS", ...)
        fit$counts <- replace(count0, is.na(count0), 0) + replace(fit$counts, is.na(fit$counts), 0)
        
      }
      fit$hessian <- numDeriv::hessian(func=loglik.clogit, x=fit$par, D="rcdose_ameras", status=status,X=X, M=M, doseRRmod=doseRRmod, designmat=designmat,entry=entry,exit=exit, data=data, deg=deg, ERC=ERC, Kmat=Kmat, loglim=loglim, transform=transform, ...)
    }
    
    if(doseRRmod!="LINEXP"){
      parnames <- c(names(data[,X,drop=FALSE]),c("dose","dose_squared")[1:deg])
      if(!is.null(M)){
        parnames <- c(parnames, paste0("dose:",names(data[, M,drop=FALSE])))
        if(deg==2){
          parnames <- c(parnames, paste0("dose_squared:",names(data[, M,drop=FALSE])))
        }
      }
    } else{
      parnames <- c(names(data[,X,drop=FALSE]),c("dose_linear","dose_exponential"))
      if(!is.null(M)){
        parnames <- c(parnames, paste0("dose_linear:",names(data[, M,drop=FALSE])))
        parnames <- c(parnames, paste0("dose_exponential:",names(data[, M,drop=FALSE])))
      }
    }
    
  } else if(family == "prophaz"){
    if(is.null(exit)) stop("exit is required for family=prophaz")
    if(is.null(doseRRmod)) stop("doseRRmod is required for family=prophaz")
    if(is.null(status)) stop("status is required for family=prophaz")
    
    
    if(is.null(inpar)){
      inpar <- rep(0, length(X)+length(M)*deg+deg)
    }
    
    if(length(X)+length(M)*deg+deg == 1){ # Optimize 1-dimensional model: use optimize instead of optim
      fit0 <- optimize(f=loglik.prophaz, lower=-20, upper=5, D="rcdose_ameras", status=status,X=X, M=M, doseRRmod=doseRRmod,entry=entry,exit=exit, data=data, deg=deg, ERC=ERC, Kmat=Kmat, loglim=loglim, transform=transform, ...)
      fit <- list(par=fit0$minimum, value=fit0$objective, convergence=0, hessian=numDeriv::hessian(func=loglik.prophaz, x=fit0$minimum, D="rcdose_ameras", status=status,X=X, M=M, doseRRmod=doseRRmod,entry=entry,exit=exit, data=data, deg=deg, ERC=ERC, Kmat=Kmat, loglim=loglim, transform=transform, ...))
    } else {
      
      fit <- optim(inpar, loglik.prophaz, D="rcdose_ameras", status=status,X=X, M=M, doseRRmod=doseRRmod,entry=entry,exit=exit, data=data, deg=deg, ERC=ERC, Kmat=Kmat, loglim=loglim, transform=transform, method=optim.method, ...)
      
      if(optim.method=="Nelder-Mead"){
        count0 <- fit$counts
        fit <- optim(fit$par, loglik.prophaz, D="rcdose_ameras", status=status,X=X, M=M, doseRRmod=doseRRmod, entry=entry,exit=exit, data=data, deg=deg, ERC=ERC, Kmat=Kmat, loglim=loglim, transform=transform, method="BFGS", ...)
        fit$counts <- replace(count0, is.na(count0), 0) + replace(fit$counts, is.na(fit$counts), 0)
        
      }
      fit$hessian <- numDeriv::hessian(func=loglik.prophaz, x=fit$par, D="rcdose_ameras", status=status,X=X, M=M, doseRRmod=doseRRmod, entry=entry,exit=exit, data=data, deg=deg, ERC=ERC, Kmat=Kmat, loglim=loglim, transform=transform, ...)
    }
    
    if(doseRRmod!="LINEXP"){
      parnames <- c(names(data[,X,drop=FALSE]),c("dose","dose_squared")[1:deg])
      if(!is.null(M)){
        parnames <- c(parnames, paste0("dose:",names(data[, M,drop=FALSE])))
        if(deg==2){
          parnames <- c(parnames, paste0("dose_squared:",names(data[, M,drop=FALSE])))
        }
      }
    } else{
      parnames <- c(names(data[,X,drop=FALSE]),c("dose_linear","dose_exponential"))
      if(!is.null(M)){
        parnames <- c(parnames, paste0("dose_linear:",names(data[, M,drop=FALSE])))
        parnames <- c(parnames, paste0("dose_exponential:",names(data[, M,drop=FALSE])))
      }
    }
    
  } else if(family=="multinomial"){
    
    if(is.null(Y)) stop("Y is required for family=multinomial")
    
    if(is.null(inpar)){
      inpar <- rep(0, (length(unique(data[,Y]))-1)*(1+length(X)+length(M)*deg+deg))
    }
    
    fit <- optim(inpar, loglik.multinomial, D="rcdose_ameras", X=X, Y=Y, M=M, doseRRmod=doseRRmod, data=data, deg=deg, ERC=ERC, Kmat=Kmat, loglim=loglim, transform=transform, method=optim.method, ...)
    if(optim.method=="Nelder-Mead"){
      count0 <- fit$counts
      fit <- optim(fit$par, loglik.multinomial, D="rcdose_ameras", X=X, Y=Y, M=M, doseRRmod=doseRRmod, data=data, deg=deg, ERC=ERC, Kmat=Kmat, loglim=loglim, transform=transform, method="BFGS", ...)
      fit$counts <- replace(count0, is.na(count0), 0) + replace(fit$counts, is.na(fit$counts), 0)
    }
    fit$hessian <- numDeriv::hessian(func=loglik.multinomial, x=fit$par, D="rcdose_ameras",X=X, Y=Y, M=M, doseRRmod=doseRRmod, data=data, deg=deg, ERC=ERC, Kmat=Kmat, loglim=loglim, transform=transform, ...)
    
    
    if(doseRRmod!="LINEXP"){
      parnames <- c("(Intercept)",names(data[,X,drop=FALSE]),c("dose","dose_squared")[1:deg])
      if(!is.null(M)){
        parnames <- c(parnames, paste0("dose:",names(data[, M,drop=FALSE])))
        if(deg==2){
          parnames <- c(parnames, paste0("dose_squared:",names(data[, M,drop=FALSE])))
        }
      }
    } else{
      parnames <- c("(Intercept)",names(data[,X,drop=FALSE]),c("dose_linear","dose_exponential"))
      if(!is.null(M)){
        parnames <- c(parnames, paste0("dose_linear:",names(data[, M,drop=FALSE])))
        parnames <- c(parnames, paste0("dose_exponential:",names(data[, M,drop=FALSE])))
      }
    }
    
    mylv <- levels(data[,Y])
    
    mylv <- mylv[-length(mylv)]
    
    parnames <- do.call("c",lapply(mylv, function(lv) paste0("(",lv, ")_", parnames)))
    
    
  }
  
  if(!is.null(transform) & !is.null(transform.jacobian)){
    if(is.function(transform) & is.function(transform.jacobian)){
      if("boundcheck" %in% names(formals(transform))){
        coefs <- transform(fit$par, boundcheck=TRUE, ...)
      } else {
        coefs <- transform(fit$par, ...)
      }
      
      if(det(fit$hessian)!=0 & rcond(fit$hessian)>.Machine$double.eps &  all(eigen(fit$hessian)$values > 0)){
        jac <- transform.jacobian(fit$par, ...)
        vcov <- jac %*% solve(fit$hessian) %*% t(jac)
      } else{
        warning("WARNING: Hessian was not invertible or inverse was not positive definite, variance matrix could not be obtained")
        vcov <- matrix(NA, ncol=length(parnames), nrow=length(parnames))
      }
    } else{
      stop("transform and transform.jacobian should be functions")
    } 
  } else{
    coefs <- fit$par
    if(det(fit$hessian)!=0 & rcond(fit$hessian)>.Machine$double.eps &  all(eigen(fit$hessian)$values > 0)){
      vcov <- solve(fit$hessian)
    } else{
      warning("WARNING: Hessian was not invertible or inverse was not positive definite, variance matrix could not be obtained")
      vcov <- matrix(NA, ncol=length(parnames), nrow=length(parnames))
    }
  }
  
  
  
  names(coefs) <- parnames
  rownames(vcov) <- colnames(vcov) <- parnames
  
  if(CI=="wald.transformed"){
    if(det(fit$hessian)!=0 & rcond(fit$hessian)>.Machine$double.eps &  all(eigen(fit$hessian)$values > 0)){
      CIlower <- transform(fit$par-1.96*sqrt(diag(solve(fit$hessian))), ...)
      CIupper <- transform(fit$par+1.96*sqrt(diag(solve(fit$hessian))), ...)
    } else {
      warning("WARNING: Hessian was not invertible or inverse was not positive definite, confidence intervals could not be obtained")
      CIlower <- CIupper <- NA * fit$par
    }
  } else if(CI=="wald.orig"){
    CIlower <- coefs-1.96*sqrt(diag(vcov))
    CIupper <- coefs+1.96*sqrt(diag(vcov))
  } else if(CI=="proflik"){
    if(det(fit$hessian)!=0 & rcond(fit$hessian)>.Machine$double.eps &  all(eigen(fit$hessian)$values > 0)){
      lowlims <- fit$par-4*sqrt(diag(solve(fit$hessian)))
      uplims <- fit$par+4*sqrt(diag(solve(fit$hessian)))
    } else{
      lowlims <- rep(-20, length(fit$par))
      uplims <- rep(10, length(fit$par))
    }
    
    if(family=="gaussian"){
      funforprofCI <- function(params, ...){
        loglik.gaussian(params=params, D="rcdose_ameras",X=X, Y=Y, M=M, data=data, deg=deg, ERC=ERC, Kmat=Kmat, loglim=loglim, transform=transform, optim.method=optim.method, ...)
      }
    } else if(family=="binomial"){
      funforprofCI <- function(params, ...){
        loglik.binomial(params=params, D="rcdose_ameras",X=X, Y=Y, M=M, data=data, deg=deg,doseRRmod=doseRRmod, ERC=ERC, Kmat=Kmat, loglim=loglim, transform=transform, optim.method=optim.method,...)
      }
    } else if(family=="poisson"){
      if(ERC){
        funforprofCI <- function(params, ...){
          loglik.poisson.erc(params=params, D=dosevars,X=X, Y=Y, M=M, offset=offset, data=data, doseRRmod=doseRRmod, deg=deg, loglim=loglim, transform=transform, optim.method=optim.method,...)
        }
      } else{
        funforprofCI <- function(params, ...){
          loglik.poisson(params=params, D="rcdose_ameras",X=X, Y=Y, M=M, offset=offset, data=data, doseRRmod=doseRRmod, deg=deg, loglim=loglim, transform=transform, optim.method=optim.method,...)
        }
      }
      
    } else if(family=="prophaz"){
      funforprofCI <- function(params, ...){
        loglik.prophaz(params=params, D="rcdose_ameras",X=X, status=status,entry=entry, exit=exit, M=M, data=data, doseRRmod=doseRRmod, deg=deg, ERC=ERC, Kmat=Kmat, loglim=loglim, transform=transform, optim.method=optim.method,...)
      } 
    } else if(family=="clogit"){
      funforprofCI <- function(params, ...){
        loglik.clogit(params=params, D="rcdose_ameras",X=X, status=status, designmat=designmat, M=M, data=data, doseRRmod=doseRRmod, deg=deg, ERC=ERC, Kmat=Kmat, loglim=loglim, transform=transform, optim.method=optim.method,...)
      } 
    } else if(family=="multinomial"){
      funforprofCI <- function(params, ...){
        loglik.multinomial(params=params, D="rcdose_ameras",X=X, Y=Y, M=M, data=data, deg=deg,doseRRmod=doseRRmod, ERC=ERC, Kmat=Kmat, loglim=loglim, transform=transform, optim.method=optim.method,...)
      }
    }
    
    CIlower <- rep(NA, length(coefs))
    CIupper <- rep(NA, length(coefs))
    
    pval_lower <- pval_upper <- rep(NA, length(coefs))
    iter_lower <- iter_upper <- rep(NA, length(coefs))
    
    if(params.profCI=="all"){
      CIindices <- 1:length(parnames)
    } else if(params.profCI=="dose"){
      CIindices <- (1:length(parnames))[startsWith(parnames, "dose") | grepl(")_dose", parnames)]
    }
    
    if(is.null(transform)){
      p15 <- rep(15, length(parnames))
      pminus10 <- rep(-10, length(parnames))
    } else{
      p15 <- transform(rep(15, length(parnames)), ...)
      pminus10 <- transform(rep(-10, length(parnames)), ...)
    }
    
    for(myindex in 1:length(parnames)){
      if(myindex %in% CIindices){  # Implemented this way rather than a loop over a subset of parameters so that the full covariate vector length is maintained, keeping transformation functionality intact
        
        message(paste0("Obtaining profile likelihood CI for ", parnames[myindex]))
        val15 <- profCIfun(par=15,optval=fit$value, index=myindex, fun=funforprofCI, inpar=fit$par, optim.method=optim.method, ...)
        if(val15 > 0){
          if(is.null(transform)){
            warning(paste0("WARNING: Upper bound for ", parnames[myindex], " is >15 and may not exist if rescaling the variable does not help"))
          } else{
            warning(paste0("WARNING: Upper bound for ", parnames[myindex], " is >",round(p15[myindex],1)," and may not exist if rescaling the variable does not help"))
          }
          uproot <- list(root=Inf, f.root=val15, iter=NA)
        } else{
          uproot <- tryCatch(uniroot(profCIfun, lower=fit$par[myindex], upper=uplims[myindex], optval=fit$value, index=myindex, fun=funforprofCI, inpar=fit$par, optim.method=optim.method,extendInt="downX", maxiter=maxit.profCI, tol=tol.profCI, ...), error=function(e) list(root=NA, f.root=NA, iter=NA))
          if(!is.na(uproot$f.root)){
            if(abs(uproot$f.root)>.005) warning(paste0("P-value for ", parnames[myindex]," upper bound more than 0.005 away from 0.05, reducing tol.profCI and/or increasing maxit.profCI is recommended"))
          }
        }
        
        valminus10 <- profCIfun(par=-10,optval=fit$value, index=myindex, fun=funforprofCI, inpar=fit$par, optim.method=optim.method, ...)
        if(valminus10 > 0){
          if(is.null(transform)){
            warning(paste0("WARNING: Lower bound for ", parnames[myindex], " is < -10 and may not exist if rescaling the variable does not help"))
          } else{
            warning(paste0("WARNING: Lower bound for ", parnames[myindex], " is < ",round(pminus10[myindex],1)," and may not exist if rescaling the variable does not help"))
          }
          lowroot <- list(root=-Inf, f.root=valminus10, iter=NA)
        } else{
          lowroot <- tryCatch(uniroot(profCIfun, lower=lowlims[myindex], upper=fit$par[myindex], optval=fit$value, index=myindex, fun=funforprofCI, inpar=fit$par, optim.method=optim.method,extendInt="upX", maxiter=maxit.profCI, tol=tol.profCI, ...), error=function(e) list(root=NA, f.root=NA, iter=NA))
          if(!is.na(lowroot$f.root)){
            if(abs(lowroot$f.root)>.005) warning(paste0("P-value for ", parnames[myindex]," lower bound more than 0.005 away from 0.05, reducing tol.profCI and/or increasing maxit.profCI is recommended"))
          }
        }
        CIlower[myindex] <- lowroot$root
        CIupper[myindex] <- uproot$root
        
        pval_lower[myindex] <- lowroot$f.root+.05
        pval_upper[myindex] <- uproot$f.root+.05
        
        
        iter_lower[myindex] <- lowroot$iter
        iter_upper[myindex] <- uproot$iter
      } else{
        CIlower[myindex] <- CIupper[myindex] <- pval_lower[myindex] <- pval_upper[myindex] <- iter_lower[myindex] <- iter_upper[myindex] <- NA
      }
    }
    
    if(!is.null(transform)){
      CIlower <- transform(CIlower, ...)
      CIupper <- transform(CIupper, ...)
    }
    
    profCI_pval <- data.frame(pval.lower=pval_lower, pval.upper=pval_upper)
    profCI_iter <- data.frame(iter.lower=iter_lower, iter.upper=iter_upper)
    
  }
  
  CIresult <- data.frame(lower=CIlower, upper=CIupper)
  if(CI=="proflik"){
    CIresult <- cbind(CIresult, profCI_pval, profCI_iter)
    rownames(CIresult) <- parnames
    CIresult <- CIresult[CIindices,]
  } else{
    rownames(CIresult) <- parnames
  }
  
  t1 <- proc.time()
  timedif <- t1-t0
  runtime <- paste(round(as.numeric(as.difftime(timedif["elapsed"], units="secs")),1), "seconds")
  
  
  
  
  
  
  out <- list(coefficients=coefs, 
              sd=sqrt(diag(vcov)),
              vcov=vcov,
              CI=CIresult, 
              convergence.optim=fit$convergence, 
              counts.optim=fit$counts,
              loglik=-1*fit$value,
              runtime=runtime)
  
  return(out)
}
