
ameras.mcml <- function(family, dosevars, data, deg, transform=NULL,transform.jacobian=NULL, Y=NULL, M=NULL, X=NULL, offset=NULL, inpar=NULL, entry=NULL, exit=NULL, status=NULL,setnr=setnr, doseRRmod=NULL, loglim=1e-30, optim.method="Nelder-Mead", ...){
  # if(length(CI)==0) stop("No CI method specified")
  # CI <- CI[CI %in% c("wald.orig","wald.transformed", "proflik")]
  # if(length(CI)>1) stop("Provide one type of CI for MCML: one of wald.orig, wald.transformed, and proflik")
  # if(length(CI)==0) stop("Incorrect CI method specified, should be one of wald.orig, wald.transformed, and proflik")
  # if(CI=="proflik" & !(params.profCI %in% c("dose","all"))) stop("Incorrect choice of parameters for profile likelihood CI supplied, should be either all or dose")
  # if(CI=="wald.transformed" & is.null(transform)) stop("No transformation specified, specify transformation or choose a different CI type")
  # 
  # if(CI=="proflik") message("Note: computation times for profile likelihood intervals for MCML may be extensive with large datasets or complex models")
  
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
  
  t1 <- proc.time()
  timedif <- t1-t0
  runtime <- paste(round(as.numeric(as.difftime(timedif["elapsed"], units="secs")),1), "seconds")
  
  
  out <- list(coefficients=coefs,
              sd=sqrt(diag(vcov)),
              vcov=vcov, 
              optim = list(
                par = fit$par,
                hessian = fit$hessian,
                convergence = fit$convergence,
                counts = fit$counts
              ),
              loglik=-1*fit$value,
              runtime=runtime)
  return(out)
}

ameras.rc <- function(family, dosevars, data, deg, ERC=FALSE, transform=NULL, transform.jacobian=NULL, Y=NULL, M=NULL, X=NULL, offset=NULL, inpar=NULL, entry=NULL, exit=NULL, status=NULL, setnr=NULL, doseRRmod=NULL, loglim=1e-30, optim.method="Nelder-Mead", ...){
 
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
  
  
  t1 <- proc.time()
  timedif <- t1-t0
  runtime <- paste(round(as.numeric(as.difftime(timedif["elapsed"], units="secs")),1), "seconds")
  
  
  out <- list(coefficients=coefs, 
              sd=sqrt(diag(vcov)),
              vcov=vcov,
              optim = list(
                par = fit$par,
                hessian = fit$hessian,
                convergence = fit$convergence,
                counts = fit$counts
              ),
              loglik=-1*fit$value,
              runtime=runtime,
              ERC = ERC)
  
  return(out)
}
