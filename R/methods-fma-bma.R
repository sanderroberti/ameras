

ameras.fma <- function(family, dosevars, data, deg, transform=NULL,transform.jacobian=NULL, Y=NULL, M=NULL, X=NULL, offset=NULL, inpar=NULL, entry=NULL, exit=NULL, status=NULL, setnr=setnr, unweighted=NULL, doseRRmod=NULL, MFMA=100000, optim.method="Nelder-Mead", ...){
  

  if(is.null(unweighted)){
    unweighted <- FALSE
  } else if(!is.logical(unweighted)){
    stop("unweighted should be TRUE or FALSE")
  }
  
  t0 <- proc.time()
  if(family=="gaussian"){
    if(is.null(Y)) stop("Y is required for family=gaussian")
    
    if(is.null(inpar)){
      inpar <- rep(0, 2+length(X)+length(M)*deg+deg)
    }
    
    FMAfits <-  lapply(1:length(dosevars), function(Xi){
      
      fit.FMAi <- optim(inpar, loglik.gaussian, D=dosevars[Xi], X=X, Y=Y, M=M, data=data, deg=deg, ERC=FALSE, transform=transform, method=optim.method, ...)
      if(optim.method=="Nelder-Mead"){
        fit.FMAi <- optim(fit.FMAi$par, loglik.gaussian, D=dosevars[Xi], X=X, Y=Y, M=M, data=data, deg=deg, ERC=FALSE, transform=transform, method="BFGS", ...)
      }
      
      fit.FMAi$hessian <- numDeriv::hessian(func=loglik.gaussian, x=fit.FMAi$par, D=dosevars[Xi], X=X, Y=Y, M=M, data=data, deg=deg, ERC=FALSE, transform=transform, ...)
      
      if(det(fit.FMAi$hessian)!=0 & rcond(fit.FMAi$hessian)>.Machine$double.eps & fit.FMAi$convergence==0 &  all(eigen(fit.FMAi$hessian)$values > 0)){
        include <- TRUE
      } else{
        include <- FALSE
      }
      
      list(coef=fit.FMAi$par, hess=fit.FMAi$hessian, AIC=2*fit.FMAi$value+2*(2+length(X)+length(M)*deg+deg), convergence=fit.FMAi$convergence, include=include)
      
    })
    
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
    
    FMAfits <-  lapply(1:length(dosevars), function(Xi){
      
      fit.FMAi <- optim(inpar, loglik.binomial, D=dosevars[Xi], X=X, Y=Y, M=M, doseRRmod=doseRRmod, data=data, deg=deg, ERC=FALSE, transform=transform, method=optim.method, ...)
      if(optim.method=="Nelder-Mead"){
        fit.FMAi <- optim(fit.FMAi$par, loglik.binomial, D=dosevars[Xi], X=X, Y=Y, M=M, doseRRmod=doseRRmod, data=data, deg=deg, ERC=FALSE, transform=transform, method="BFGS", ...)
      }
      fit.FMAi$hessian <- numDeriv::hessian(func=loglik.binomial, x=fit.FMAi$par, D=dosevars[Xi], X=X, Y=Y, M=M, doseRRmod=doseRRmod, data=data, deg=deg, ERC=FALSE, transform=transform, ...)
      
      if(det(fit.FMAi$hessian)!=0 & rcond(fit.FMAi$hessian)>.Machine$double.eps & fit.FMAi$convergence==0 &  all(eigen(fit.FMAi$hessian)$values > 0)){
        include <- TRUE
      } else{
        include <- FALSE
      }
      
      list(coef=fit.FMAi$par, hess=fit.FMAi$hessian, AIC=2*fit.FMAi$value+2*(1+length(X)+length(M)*deg+deg), convergence=fit.FMAi$convergence, include=include)
      
    })
    
    
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
    
    FMAfits <-  lapply(1:length(dosevars), function(Xi){
      
      fit.FMAi <- optim(inpar, loglik.poisson, D=dosevars[Xi], X=X, Y=Y, M=M, offset=offset, doseRRmod=doseRRmod, data=data, deg=deg, transform=transform, method=optim.method, ...)
      if(optim.method=="Nelder-Mead"){
        fit.FMAi <- optim(fit.FMAi$par, loglik.poisson, D=dosevars[Xi], X=X, Y=Y, M=M, offset=offset, doseRRmod=doseRRmod, data=data, deg=deg, transform=transform, method="BFGS", ...)
      }
      fit.FMAi$hessian <- numDeriv::hessian(func=loglik.poisson, x=fit.FMAi$par, D=dosevars[Xi], X=X, Y=Y, M=M, offset=offset, doseRRmod=doseRRmod, data=data, deg=deg, transform=transform, ...)
      
      if(det(fit.FMAi$hessian)!=0 & rcond(fit.FMAi$hessian)>.Machine$double.eps & fit.FMAi$convergence==0 &  all(eigen(fit.FMAi$hessian)$values > 0)){
        include <- TRUE
      } else{
        include <- FALSE
      }
      
      list(coef=fit.FMAi$par, hess=fit.FMAi$hessian, AIC=2*fit.FMAi$value+2*(1+length(X)+length(M)*deg+deg), convergence=fit.FMAi$convergence, include=include)
      
    })
    
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
    if(family=="prophaz" & is.null(exit)) stop("exit is required for family=prophaz")
    if(is.null(doseRRmod)) stop("doseRRmod is required for family=prophaz")
    if(is.null(status)) stop("status is required for family=prophaz")
    
    
    designmat <- t(model.matrix(~as.factor(data[,setnr])-1))
    
    
    
    if(is.null(inpar)){
      inpar <- rep(0, length(X)+length(M)*deg+deg)
    }
    
    FMAfits <-  lapply(1:length(dosevars), function(Xi){
      
      if(length(X)+length(M)*deg+deg == 1){ # Optimize 1-dimensional model: use optimize instead of optim
        fit0 <- optimize(f=loglik.clogit, lower=-20, upper=5, D=dosevars[Xi], status=status, X=X, M=M, designmat=designmat, doseRRmod=doseRRmod, data=data, deg=deg, ERC=FALSE, transform=transform, ...)
        fit.FMAi <- list(par=fit0$minimum, value=fit0$objective, convergence=0, hessian=numDeriv::hessian(func=loglik.clogit, x=fit0$minimum, D=dosevars[Xi], status=status, X=X, M=M, designmat=designmat, doseRRmod=doseRRmod, data=data, deg=deg, ERC=FALSE, transform=transform, ...))
      } else {
        fit.FMAi <- optim(inpar, loglik.clogit, D=dosevars[Xi], status=status, X=X, M=M, designmat=designmat, doseRRmod=doseRRmod, data=data, deg=deg, ERC=FALSE, transform=transform, method=optim.method, ...)
        if(optim.method=="Nelder-Mead"){
          fit.FMAi <- optim(fit.FMAi$par, loglik.clogit, D=dosevars[Xi], status=status, X=X, M=M, designmat=designmat, doseRRmod=doseRRmod, data=data, deg=deg, ERC=FALSE, transform=transform, method="BFGS", ...)
        }
        fit.FMAi$hessian <- numDeriv::hessian(func=loglik.clogit, x=fit.FMAi$par, D=dosevars[Xi], status=status, X=X, M=M, designmat=designmat, doseRRmod=doseRRmod, data=data, deg=deg, ERC=FALSE, transform=transform, ...)
      }
      
      
      if(det(fit.FMAi$hessian)!=0 & rcond(fit.FMAi$hessian)>.Machine$double.eps & fit.FMAi$convergence==0 &  all(eigen(fit.FMAi$hessian)$values > 0)){
        include <- TRUE
      } else{
        include <- FALSE
      }
      
      list(coef=fit.FMAi$par, hess=fit.FMAi$hessian, AIC=2*fit.FMAi$value+2*(length(X)+length(M)*deg+deg), convergence=fit.FMAi$convergence, include=include)
      
    })
    
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
    
    FMAfits <-  lapply(1:length(dosevars), function(Xi){
      
      if(length(X)+length(M)*deg+deg == 1){ # Optimize 1-dimensional model: use optimize instead of optim
        fit0 <- optimize(f=loglik.prophaz, lower=-20, upper=5, D=dosevars[Xi], status=status, X=X, M=M, entry=entry, exit=exit, doseRRmod=doseRRmod, data=data, deg=deg, ERC=FALSE, transform=transform, ...)
        fit.FMAi <- list(par=fit0$minimum, value=fit0$objective, convergence=0, hessian=numDeriv::hessian(func=loglik.prophaz, x=fit0$minimum, D=dosevars[Xi], status=status, X=X, M=M, entry=entry, exit=exit, doseRRmod=doseRRmod, data=data, deg=deg, ERC=FALSE, transform=transform, ...))
      } else {
        fit.FMAi <- optim(inpar, loglik.prophaz, D=dosevars[Xi], status=status, X=X, M=M, entry=entry, exit=exit, doseRRmod=doseRRmod, data=data, deg=deg, ERC=FALSE, transform=transform, method=optim.method, ...)
        if(optim.method=="Nelder-Mead"){
          fit.FMAi <- optim(fit.FMAi$par, loglik.prophaz, D=dosevars[Xi], status=status, X=X, M=M, entry=entry, exit=exit, doseRRmod=doseRRmod, data=data, deg=deg, ERC=FALSE, transform=transform, method="BFGS", ...)
        }
        fit.FMAi$hessian <- numDeriv::hessian(func=loglik.prophaz, x=fit.FMAi$par, D=dosevars[Xi], status=status, X=X, M=M, entry=entry, exit=exit, doseRRmod=doseRRmod, data=data, deg=deg, ERC=FALSE, transform=transform, ...)
      }
      
      
      
      if(det(fit.FMAi$hessian)!=0 & rcond(fit.FMAi$hessian)>.Machine$double.eps & fit.FMAi$convergence==0 &  all(eigen(fit.FMAi$hessian)$values > 0)){
        include <- TRUE
      } else{
        include <- FALSE
      }
      
      list(coef=fit.FMAi$par, hess=fit.FMAi$hessian, AIC=2*fit.FMAi$value+2*(length(X)+length(M)*deg+deg), convergence=fit.FMAi$convergence, include=include)
      
    })
    
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
    
    FMAfits <-  lapply(1:length(dosevars), function(Xi){
      
      fit.FMAi <- optim(inpar, loglik.multinomial, D=dosevars[Xi], X=X, Y=Y, M=M, doseRRmod=doseRRmod, data=data, deg=deg, ERC=FALSE, transform=transform, method=optim.method, ...)
      if(optim.method=="Nelder-Mead"){
        fit.FMAi <- optim(fit.FMAi$par, loglik.multinomial, D=dosevars[Xi], X=X, Y=Y, M=M, doseRRmod=doseRRmod, data=data, deg=deg, ERC=FALSE, transform=transform, method="BFGS", ...)
      }
      fit.FMAi$hessian <- numDeriv::hessian(func=loglik.multinomial, x=fit.FMAi$par, D=dosevars[Xi], X=X, Y=Y, M=M, doseRRmod=doseRRmod, data=data, deg=deg, ERC=FALSE, transform=transform, ...)
      
      if(det(fit.FMAi$hessian)!=0 & rcond(fit.FMAi$hessian)>.Machine$double.eps & fit.FMAi$convergence==0 &  all(eigen(fit.FMAi$hessian)$values > 0)){
        include <- TRUE
      } else{
        include <- FALSE
      }
      
      list(coef=fit.FMAi$par, hess=fit.FMAi$hessian, AIC=2*fit.FMAi$value+2*(1+length(X)+length(M)*deg+deg), convergence=fit.FMAi$convergence, include=include)
      
    })
    
    
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
  
  #Extra exclusion check for asymmetric variance matrix
  for(iii in 1:length(FMAfits)){
    fit.FMAi <- FMAfits[[iii]]
    if(fit.FMAi$include){
      if(!is.null(transform) & !is.null(transform.jacobian)){
        if(is.function(transform) & is.function(transform.jacobian)){
          jac <- transform.jacobian(fit.FMAi$coef, ...)
          samplevar <- jac %*% solve(fit.FMAi$hess) %*% t(jac)
        } else{
          stop("transform and transform.jacobian should be functions")
        } 
      } else{
        samplevar <- solve(fit.FMAi$hess)
      }
      if(!isSymmetric(samplevar)){
        FMAfits[[iii]]$include <- FALSE
      }
    }
  }
  
  #allfits <- FMAfits
  included.replicates <- which(sapply(FMAfits, function(x) x$include))
  
  if(length(included.replicates)>0){
    FMAfits <- FMAfits[sapply(FMAfits, function(x) x$include)]
    
    meanAIC <- mean(sapply(FMAfits, function(x) x$AIC))
    
    FMAfits <- lapply(FMAfits, function(x){
      x$AIC_cent <- x$AIC-meanAIC
      x
    })
    
    FMAfits <- lapply(FMAfits, function(x){
      if(unweighted){
        x$wgt <- 1/length(FMAfits)
      } else{
        x$wgt <- exp(-.5*x$AIC_cent)/sum(sapply(FMAfits, function(y) exp(-.5*y$AIC_cent)))
      }
      x$M <- max(round(x$wgt*MFMA),1)
      x
    })
    
    
    
    FMAsamples <- lapply(FMAfits, function(fit.FMAi, ...){
      
      if(!is.null(transform) & !is.null(transform.jacobian)){
        if(is.function(transform) & is.function(transform.jacobian)){
          samplemeans <- transform(fit.FMAi$coef, ...)
          jac <- transform.jacobian(fit.FMAi$coef, ...)
          samplevar <- jac %*% solve(fit.FMAi$hess) %*% t(jac)
        } else{
          stop("transform and transform.jacobian should be functions")
        } 
      } else{
        samplemeans <- fit.FMAi$coef
        samplevar <- solve(fit.FMAi$hess)
      }
      
      return(rmvnorm(n=fit.FMAi$M, mean=samplemeans, sigma=samplevar))
      
      # if(isSymmetric(samplevar)){
      #   return(rmvnorm(n=fit.FMAi$M, mean=samplemeans, sigma=samplevar))
      # } else{
      #   return(NULL) 
      # }
      
      
      
    }, ...)
    
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
  timedif <- t1-t0
  runtime <- paste(round(as.numeric(as.difftime(timedif["elapsed"], units="secs")),1), "seconds")
  
  prc_excluded <- round(100*(1-length(included.replicates)/length(dosevars)), 1)
  if(length(included.replicates)/length(dosevars) < .8) warning(paste0("WARNING: ",prc_excluded, "% of replicates excluded from model averaging. Try different bounds or starting values."))
  
  
  out <- list(
    coefficients=coefs,
    sd=sd,
    included.replicates=included.replicates,
    included.samples=included.samples,
    samples=FMAsamples,
    #includedfits=FMAfits,
    #allfits=allfits,
    runtime=runtime
  )
  
  return(out)
  
}

boundedSliceSampler <- nimbleFunction(
  contains = nimble::sampler_BASE,
  setup = function(model, mvSaved, target, control) {
    # Create internal slice sampler safely
    sliceSampler <- nimble::sampler_slice(
      model = model,
      mvSaved = mvSaved,
      target = target,
      control = control
    )
    
    # Optionally accept user-supplied K from control
    if (!is.null(control$K)) {
      K <- control$K
    } else {
      # Default: try to detect length of p if it exists
      if ("p" %in% model$getNodeNames()) {
        K <- length(model[['p']])
      } else {
        stop("Please provide K (number of categories) in 'control' for boundedSliceSampler.")
      }
    }
  },
  run = function() {
    sliceSampler$run()
    
    value <- model[[target]]
    
    # Apply integer bounds
    if (value < 1) value <- 1
    if (value > K) value <- K
    
    model[[target]] <<- value
  },
  methods = list(
    reset = function() {
      sliceSampler$reset()
    }
  )
)

nimblemod <- nimbleCode({
  
  col.ind~dcat(w[1:K])
  for (k in 1:K) {
    vec[k] <- equals(col.ind, k)
  }
  
  
  if(family != "multinomial"){
    
    if(!(family%in%c("prophaz","clogit"))){
      a0~dnorm(0,0.001)
    } else{
      a0 <- 0
    }
    
    if(Xlen>1){
      for(k in 1:Xlen){
        a[k]~dnorm(0, 0.001)
      }
    } else if(Xlen==1){
      a~dnorm(0, .001)
    }
    
    if(doseRRmod=="EXP" | family == "gaussian"){ # Normal priors with large variance
      if(deg>1){
        for(k in 1:deg){
          b[k]~dnorm(0, 0.001)
        }
      } else{
        b~dnorm(0, .001)
      }
      
      if(Mlen>1){
        if(deg==1){
          for(k in 1:Mlen){
            bm[k]~dnorm(0, 0.001)
          }
        } else if(deg>1){
          for(k in 1:Mlen){
            bm1[k]~dnorm(0, 0.001)
            bm2[k]~dnorm(0, 0.001)
          }
        }
      } else if(Mlen==1){
        if(deg==1){
          bm~dnorm(0, 0.001)
        } else if(deg>1){
          bm1~dnorm(0, 0.001)
          bm2~dnorm(0, 0.001)
        }
      }
    } else if(doseRRmod=="ERR"){
      if(ERRprior=="truncated_horseshoe"){
        tau~T(dt(0,1, df=1), 0,)
        if(deg>1){
          for(kk in 1:2){ # horseshoe prior
            lambda[kk]~T(dt(0,1, df=1), 0,)
            b[kk]~T(dnorm(0,sd=lambda[kk]*tau), 0,)
          }
        } else{
          lambda~T(dt(0,1, df=1), 0,)
          b~T(dnorm(0,sd=lambda*tau), 0,)
        }
        
        if(Mlen>1){
          if(deg==1){
            for(zz in 1:Mlen){
              lambda_m[zz]~T(dt(0,1, df=1), 0,)
              bm[zz]~T(dnorm(0,sd=lambda_m[zz]*tau), 0,)
            }
          } else if(deg>1){
            for(zz in 1:Mlen){
              # horseshoe prior
              lambda_m1[zz]~T(dt(0,1, df=1), 0,)
              bm1[zz]~T(dnorm(0,sd=lambda_m1[zz]*tau), 0,)
              lambda_m2[zz]~T(dt(0,1, df=1), 0,)
              bm2[zz]~T(dnorm(0,sd=lambda_m2[zz]*tau), 0,)
            }
          }
        } else if(Mlen==1){
          if(deg==1){
            lambda_m~T(dt(0,1, df=1), 0,)
            bm~T(dnorm(0, sd=lambda_m*tau), 0,)
          } else if(deg>1){
            lambda_m1~T(dt(0,1, df=1), 0,)
            bm1~T(dnorm(0, sd=lambda_m1*tau), 0,)
            lambda_m2~T(dt(0,1, df=1), 0,)
            bm2~T(dnorm(0, sd=lambda_m2*tau), 0,)
          }
        }
      } else if(ERRprior=="truncated_normal"){
        if(deg>1){
          for(k in 1:deg){
            b[k]~T(dnorm(0,0.001), 0, )
          }
        } else{
          b~T(dnorm(0, .001), 0, )
        }
        
        if(Mlen>1){
          if(deg==1){
            for(k in 1:Mlen){
              bm[k]~T(dnorm(0, 0.001), 0, )
            }
          } else if(deg>1){
            for(k in 1:Mlen){
              bm1[k]~T(dnorm(0, 0.001), 0, )
              bm2[k]~T(dnorm(0, 0.001), 0, )
            }
          }
        } else if(Mlen==1){
          if(deg==1){
            bm~T(dnorm(0, 0.001), 0, )
          } else if(deg>1){
            bm1~T(dnorm(0, 0.001), 0, )
            bm2~T(dnorm(0, 0.001), 0, )
          }
        }
      } else if(ERRprior=="truncated_doubleexponential"){
        if(deg>1){
          for(kk in 1:2){
            lambda[kk]~T(dt(0,1, df=1), 0,)
            b[kk]~T(ddexp(0.0, lambda[kk]) , 0,)
          }
        } else{
          lambda~T(dt(0,1, df=1), 0,)
          b~T(ddexp(0.0, lambda) , 0,)
        }
        
        if(Mlen>1){
          if(deg==1){
            for(zz in 1:Mlen){
              lambda_m[zz]~T(dt(0,1, df=1), 0,)
              bm[zz]~T(ddexp(0.0, lambda_m[zz]), 0,)
            }
          } else if(deg>1){
            for(zz in 1:Mlen){
              lambda_m1[zz]~T(dt(0,1, df=1), 0,)
              bm1[zz]~T(ddexp(0.0, lambda_m1[zz]), 0,)
              lambda_m2[zz]~T(dt(0,1, df=1), 0,)
              bm2[zz]~T(ddexp(0.0, lambda_m2[zz]), 0,)
            }
          }
        } else if(Mlen==1){
          if(deg==1){
            lambda_m~T(dt(0,1, df=1), 0,)
            bm~T(ddexp(0.0, lambda_m), 0,)
          } else if(deg>1){
            lambda_m1~T(dt(0,1, df=1), 0,)
            bm1~T(ddexp(0.0, lambda_m1), 0,)
            lambda_m2~T(dt(0,1, df=1), 0,)
            bm2~T(ddexp(0.0, lambda_m2), 0,)
          }
        }
      } else if(ERRprior=="horseshoe"){
        tau~T(dt(0,1, df=1), 0,)
        if(deg>1){
          for(kk in 1:2){ # horseshoe prior
            lambda[kk]~T(dt(0,1, df=1), 0,)
            b[kk]~dnorm(0,sd=lambda[kk]*tau)
          }
        } else{
          lambda~T(dt(0,1, df=1), 0,)
          b~dnorm(0,sd=lambda*tau)
        }
        
        if(Mlen>1){
          if(deg==1){
            for(zz in 1:Mlen){
              lambda_m[zz]~T(dt(0,1, df=1), 0,)
              bm[zz]~dnorm(0,sd=lambda_m[zz]*tau)
            }
          } else if(deg>1){
            for(zz in 1:Mlen){
              # horseshoe prior
              lambda_m1[zz]~T(dt(0,1, df=1), 0,)
              bm1[zz]~dnorm(0,sd=lambda_m1[zz]*tau)
              lambda_m2[zz]~T(dt(0,1, df=1), 0,)
              bm2[zz]~dnorm(0,sd=lambda_m2[zz]*tau)
            }
          }
        } else if(Mlen==1){
          if(deg==1){
            lambda_m~T(dt(0,1, df=1), 0,)
            bm~dnorm(0, sd=lambda_m*tau)
          } else if(deg>1){
            lambda_m1~T(dt(0,1, df=1), 0,)
            bm1~dnorm(0, sd=lambda_m1*tau)
            lambda_m2~T(dt(0,1, df=1), 0,)
            bm2~dnorm(0, sd=lambda_m2*tau)
          }
        }
      } else if(ERRprior=="normal"){
        if(deg>1){
          for(k in 1:deg){
            b[k]~dnorm(0,0.001)
          }
        } else{
          b~dnorm(0, .001)
        }
        
        if(Mlen>1){
          if(deg==1){
            for(k in 1:Mlen){
              bm[k]~dnorm(0, 0.001)
            }
          } else if(deg>1){
            for(k in 1:Mlen){
              bm1[k]~dnorm(0, 0.001)
              bm2[k]~dnorm(0, 0.001)
            }
          }
        } else if(Mlen==1){
          if(deg==1){
            bm~dnorm(0, 0.001)
          } else if(deg>1){
            bm1~dnorm(0, 0.001)
            bm2~dnorm(0, 0.001)
          }
        }
      } else if(ERRprior=="doubleexponential"){
        if(deg>1){
          for(kk in 1:2){
            lambda[kk]~T(dt(0,1, df=1), 0,)
            b[kk]~ddexp(0.0, lambda[kk])
          }
        } else{
          lambda~T(dt(0,1, df=1), 0,)
          b~ddexp(0.0, lambda)
        }
        
        if(Mlen>1){
          if(deg==1){
            for(zz in 1:Mlen){
              lambda_m[zz]~T(dt(0,1, df=1), 0,)
              bm[zz]~ddexp(0.0, lambda_m[zz])
            }
          } else if(deg>1){
            for(zz in 1:Mlen){
              lambda_m1[zz]~T(dt(0,1, df=1), 0,)
              bm1[zz]~ddexp(0.0, lambda_m1[zz])
              lambda_m2[zz]~T(dt(0,1, df=1), 0,)
              bm2[zz]~ddexp(0.0, lambda_m2[zz])
            }
          }
        } else if(Mlen==1){
          if(deg==1){
            lambda_m~T(dt(0,1, df=1), 0,)
            bm~ddexp(0.0, lambda_m)
          } else if(deg>1){
            lambda_m1~T(dt(0,1, df=1), 0,)
            bm1~ddexp(0.0, lambda_m1)
            lambda_m2~T(dt(0,1, df=1), 0,)
            bm2~ddexp(0.0, lambda_m2)
          }
        }
      }
    } else if(doseRRmod=="LINEXP"){ 
      if(ERRprior=="truncated_horseshoe"){
        tau~T(dt(0,1, df=1), 0,)
        
        # horseshoe prior
        lambda[1]~T(dt(0,1, df=1), 0,)
        b[1]~T(dnorm(0,sd=lambda[1]*tau), 0,)
        lambda[2]~T(dt(0,1, df=1), 0,)
        b[2]~dnorm(0,sd=lambda[2]*tau)
        
        
        
        if(Mlen>1){
          
          for(zz in 1:Mlen){
            # horseshoe prior
            lambda_m1[zz]~T(dt(0,1, df=1), 0,)
            bm1[zz]~T(dnorm(0,sd=lambda_m1[zz]*tau), 0,)
            lambda_m2[zz]~T(dt(0,1, df=1), 0,)
            bm2[zz]~dnorm(0,sd=lambda_m2[zz]*tau)
          }
          
        } else if(Mlen==1){
          
          lambda_m1~T(dt(0,1, df=1), 0,)
          bm1~T(dnorm(0, sd=lambda_m1*tau), 0,)
          lambda_m2~T(dt(0,1, df=1), 0,)
          bm2~dnorm(0, sd=lambda_m2*tau)
          
        }
      } else if(ERRprior=="truncated_normal"){
        
        
        b[1]~T(dnorm(0,0.001), 0, )
        b[2]~dnorm(0,0.001)
        
        if(Mlen>1){
          
          for(k in 1:Mlen){
            bm1[k]~T(dnorm(0, 0.001), 0, )
            bm2[k]~dnorm(0, 0.001)
          }
          
        } else if(Mlen==1){
          
          bm1~T(dnorm(0, 0.001), 0, )
          bm2~dnorm(0, 0.001)
          
        }
      } else if(ERRprior=="truncated_doubleexponential"){
        
        lambda[1]~T(dt(0,1, df=1), 0,)
        b[1]~T(ddexp(0.0, lambda[1]) , 0,)
        lambda[2]~T(dt(0,1, df=1), 0,)
        b[2]~ddexp(0.0, lambda[2])
        
        if(Mlen>1){
          
          for(zz in 1:Mlen){
            lambda_m1[zz]~T(dt(0,1, df=1), 0,)
            bm1[zz]~T(ddexp(0.0, lambda_m1[zz]), 0,)
            lambda_m2[zz]~T(dt(0,1, df=1), 0,)
            bm2[zz]~ddexp(0.0, lambda_m2[zz])
          }
          
        } else if(Mlen==1){
          
          lambda_m1~T(dt(0,1, df=1), 0,)
          bm1~T(ddexp(0.0, lambda_m1), 0,)
          lambda_m2~T(dt(0,1, df=1), 0,)
          bm2~ddexp(0.0, lambda_m2)
          
        }
      } else if(ERRprior=="horseshoe"){
        tau~T(dt(0,1, df=1), 0,)
        
        for(kk in 1:2){ # horseshoe prior
          lambda[kk]~T(dt(0,1, df=1), 0,)
          b[kk]~dnorm(0,sd=lambda[kk]*tau)
        }
        
        
        if(Mlen>1){
          
          for(zz in 1:Mlen){
            # horseshoe prior
            lambda_m1[zz]~T(dt(0,1, df=1), 0,)
            bm1[zz]~dnorm(0,sd=lambda_m1[zz]*tau)
            lambda_m2[zz]~T(dt(0,1, df=1), 0,)
            bm2[zz]~dnorm(0,sd=lambda_m2[zz]*tau)
          }
          
        } else if(Mlen==1){
          
          lambda_m1~T(dt(0,1, df=1), 0,)
          bm1~dnorm(0, sd=lambda_m1*tau)
          lambda_m2~T(dt(0,1, df=1), 0,)
          bm2~dnorm(0, sd=lambda_m2*tau)
          
        }
      } else if(ERRprior=="normal"){
        
        for(k in 1:2){
          b[k]~dnorm(0,0.001)
        }
        
        
        if(Mlen>1){
          
          for(k in 1:Mlen){
            bm1[k]~dnorm(0, 0.001)
            bm2[k]~dnorm(0, 0.001)
          }
          
        } else if(Mlen==1){
          
          bm1~dnorm(0, 0.001)
          bm2~dnorm(0, 0.001)
          
        }
      } else if(ERRprior=="doubleexponential"){
        
        for(kk in 1:2){
          lambda[kk]~T(dt(0,1, df=1), 0,)
          b[kk]~ddexp(0.0, lambda[kk])
        }
        
        
        if(Mlen>1){
          
          for(zz in 1:Mlen){
            lambda_m1[zz]~T(dt(0,1, df=1), 0,)
            bm1[zz]~ddexp(0.0, lambda_m1[zz])
            lambda_m2[zz]~T(dt(0,1, df=1), 0,)
            bm2[zz]~ddexp(0.0, lambda_m2[zz])
          }
          
        } else if(Mlen==1){
          
          lambda_m1~T(dt(0,1, df=1), 0,)
          bm1~ddexp(0.0, lambda_m1)
          lambda_m2~T(dt(0,1, df=1), 0,)
          bm2~ddexp(0.0, lambda_m2)
          
        }
      }
    }
  } else { # Multinomial
    
    for(z in 1:(Z-1)){
      a0[z]~dnorm(0,0.001)
      
      if(Xlen>1){
        for(k in 1:Xlen){
          a[k,z]~dnorm(0, 0.001)
        }
      } else if(Xlen==1){
        a[z]~dnorm(0, .001)
      }
    }
    
    if(doseRRmod=="EXP"){ # Normal priors with large variance
      if(deg>1){
        for(z in 1:(Z-1)){
          for(k in 1:deg){
            b[k, z]~dnorm(0, 0.001)
          }
        }
      } else{
        for(z in 1:(Z-1)){
          b[z]~dnorm(0, .001)
        }
      }
      
      if(Mlen>1){
        if(deg==1){
          for(z in 1:(Z-1)){
            for(k in 1:Mlen){
              bm[k,z]~dnorm(0, 0.001)
            }
          }
        } else if(deg>1){
          for(z in 1:(Z-1)){
            for(k in 1:Mlen){
              bm1[k,z]~dnorm(0, 0.001)
              bm2[k,z]~dnorm(0, 0.001)
            }
          }
        }
      } else if(Mlen==1){
        if(deg==1){
          for(z in 1:(Z-1)){
            bm[z]~dnorm(0, 0.001)
          }
        } else if(deg>1){
          for(z in 1:(Z-1)){
            bm1[z]~dnorm(0, 0.001)
            bm2[z]~dnorm(0, 0.001)
          }
        }
      }
    } else if(doseRRmod=="ERR"){
      if(ERRprior=="truncated_horseshoe"){
        tau~T(dt(0,1, df=1), 0,)
        if(deg>1){
          for(z in 1:(Z-1)){
            for(kk in 1:2){ # horseshoe prior
              lambda[kk,z]~T(dt(0,1, df=1), 0,)
              b[kk,z]~T(dnorm(0,sd=lambda[kk,z]*tau), 0,)
            }
          }
        } else{
          for(z in 1:(Z-1)){
            lambda[z]~T(dt(0,1, df=1), 0,)
            b[z]~T(dnorm(0,sd=lambda[z]*tau), 0,)
          }
        }
        
        if(Mlen>1){
          if(deg==1){
            for(z in 1:(Z-1)){
              for(zz in 1:Mlen){
                lambda_m[zz,z]~T(dt(0,1, df=1), 0,)
                bm[zz,z]~T(dnorm(0,sd=lambda_m[zz,z]*tau), 0,)
              }
            }
          } else if(deg>1){
            for(z in 1:(Z-1)){
              for(zz in 1:Mlen){
                # horseshoe prior
                lambda_m1[zz,z]~T(dt(0,1, df=1), 0,)
                bm1[zz,z]~T(dnorm(0,sd=lambda_m1[zz,z]*tau), 0,)
                lambda_m2[zz,z]~T(dt(0,1, df=1), 0,)
                bm2[zz,z]~T(dnorm(0,sd=lambda_m2[zz,z]*tau), 0,)
              }
            }
          }
        } else if(Mlen==1){
          if(deg==1){
            for(z in 1:(Z-1)){
              lambda_m[z]~T(dt(0,1, df=1), 0,)
              bm[z]~T(dnorm(0, sd=lambda_m[z]*tau), 0,)
            }
          } else if(deg>1){
            for(z in 1:(Z-1)){
              lambda_m1[z]~T(dt(0,1, df=1), 0,)
              bm1[z]~T(dnorm(0, sd=lambda_m1[z]*tau), 0,)
              lambda_m2[z]~T(dt(0,1, df=1), 0,)
              bm2[z]~T(dnorm(0, sd=lambda_m2[z]*tau), 0,)
            }
          }
        }
      } else if(ERRprior=="truncated_normal"){
        if(deg>1){
          for(z in 1:(Z-1)){
            for(k in 1:deg){
              b[k,z]~T(dnorm(0,0.001), 0, )
            }
          }
        } else{
          for(z in 1:(Z-1)){
            b[z]~T(dnorm(0, .001), 0, )
          }
        }
        
        if(Mlen>1){
          if(deg==1){
            for(z in 1:(Z-1)){
              for(k in 1:Mlen){
                bm[k,z]~T(dnorm(0, 0.001), 0, )
              }
            }
          } else if(deg>1){
            for(z in 1:(Z-1)){
              for(k in 1:Mlen){
                bm1[k,z]~T(dnorm(0, 0.001), 0, )
                bm2[k,z]~T(dnorm(0, 0.001), 0, )
              }
            }
          }
        } else if(Mlen==1){
          if(deg==1){
            for(z in 1:(Z-1)){
              bm[z]~T(dnorm(0, 0.001), 0, )
            }
          } else if(deg>1){
            for(z in 1:(Z-1)){
              bm1[z]~T(dnorm(0, 0.001), 0, )
              bm2[z]~T(dnorm(0, 0.001), 0, )
            }
          }
        }
      } else if(ERRprior=="truncated_doubleexponential"){
        if(deg>1){
          for(z in 1:(Z-1)){
            for(kk in 1:2){
              lambda[kk,z]~T(dt(0,1, df=1), 0,)
              b[kk,z]~T(ddexp(0.0, lambda[kk,z]) , 0,)
            }
          }
        } else{
          for(z in 1:(Z-1)){
            lambda[z]~T(dt(0,1, df=1), 0,)
            b[z]~T(ddexp(0.0, lambda) , 0,)
          }
        }
        
        if(Mlen>1){
          if(deg==1){
            for(z in 1:(Z-1)){
              for(zz in 1:Mlen){
                lambda_m[zz,z]~T(dt(0,1, df=1), 0,)
                bm[zz,z]~T(ddexp(0.0, lambda_m[zz,z]), 0,)
              }
            }
          } else if(deg>1){
            for(z in 1:(Z-1)){
              for(zz in 1:Mlen){
                lambda_m1[zz,z]~T(dt(0,1, df=1), 0,)
                bm1[zz,z]~T(ddexp(0.0, lambda_m1[zz,z]), 0,)
                lambda_m2[zz,z]~T(dt(0,1, df=1), 0,)
                bm2[zz,z]~T(ddexp(0.0, lambda_m2[zz,z]), 0,)
              }
            }
          }
        } else if(Mlen==1){
          if(deg==1){
            for(z in 1:(Z-1)){
              lambda_m[z]~T(dt(0,1, df=1), 0,)
              bm[z]~T(ddexp(0.0, lambda_m[z]), 0,)
            }
          } else if(deg>1){
            for(z in 1:(Z-1)){
              lambda_m1[z]~T(dt(0,1, df=1), 0,)
              bm1[z]~T(ddexp(0.0, lambda_m1[z]), 0,)
              lambda_m2[z]~T(dt(0,1, df=1), 0,)
              bm2[z]~T(ddexp(0.0, lambda_m2[z]), 0,)
            }
          }
        }
      } else if(ERRprior=="horseshoe"){
        tau~T(dt(0,1, df=1), 0,)
        if(deg>1){
          for(z in 1:(Z-1)){
            for(kk in 1:2){ # horseshoe prior
              lambda[kk,z]~T(dt(0,1, df=1), 0,)
              b[kk,z]~dnorm(0,sd=lambda[kk,z]*tau)
            }
          }
        } else{
          for(z in 1:(Z-1)){
            lambda[z]~T(dt(0,1, df=1), 0,)
            b[z]~dnorm(0,sd=lambda[z]*tau)
          }
        }
        
        if(Mlen>1){
          if(deg==1){
            for(z in 1:(Z-1)){
              for(zz in 1:Mlen){
                lambda_m[zz,z]~T(dt(0,1, df=1), 0,)
                bm[zz,z]~dnorm(0,sd=lambda_m[zz,z]*tau)
              }
            }
          } else if(deg>1){
            for(z in 1:(Z-1)){
              for(zz in 1:Mlen){
                # horseshoe prior
                lambda_m1[zz,z]~T(dt(0,1, df=1), 0,)
                bm1[zz,z]~dnorm(0,sd=lambda_m1[zz,z]*tau)
                lambda_m2[zz,z]~T(dt(0,1, df=1), 0,)
                bm2[zz,z]~dnorm(0,sd=lambda_m2[zz,z]*tau)
              }
            }
          }
        } else if(Mlen==1){
          if(deg==1){
            for(z in 1:(Z-1)){
              lambda_m[z]~T(dt(0,1, df=1), 0,)
              bm[z]~dnorm(0, sd=lambda_m[z]*tau)
            }
          } else if(deg>1){
            for(z in 1:(Z-1)){
              lambda_m1[z]~T(dt(0,1, df=1), 0,)
              bm1[z]~dnorm(0, sd=lambda_m1[z]*tau)
              lambda_m2[z]~T(dt(0,1, df=1), 0,)
              bm2[z]~dnorm(0, sd=lambda_m2[z]*tau)
            }
          }
        }
      } else if(ERRprior=="normal"){
        if(deg>1){
          for(z in 1:(Z-1)){
            for(k in 1:deg){
              b[k,z]~dnorm(0,0.001)
            }
          }
        } else{
          for(z in 1:(Z-1)){
            b[z]~dnorm(0, .001)
          }
        }
        
        if(Mlen>1){
          if(deg==1){
            for(z in 1:(Z-1)){
              for(k in 1:Mlen){
                bm[k,z]~dnorm(0, 0.001)
              }
            }
          } else if(deg>1){
            for(z in 1:(Z-1)){
              for(k in 1:Mlen){
                bm1[k,z]~dnorm(0, 0.001)
                bm2[k,z]~dnorm(0, 0.001)
              }
            }
          }
        } else if(Mlen==1){
          if(deg==1){
            for(z in 1:(Z-1)){
              bm[z]~dnorm(0, 0.001)
            }
          } else if(deg>1){
            for(z in 1:(Z-1)){
              bm1[z]~dnorm(0, 0.001)
              bm2[z]~dnorm(0, 0.001)
            }
          }
        }
      } else if(ERRprior=="doubleexponential"){
        if(deg>1){
          for(z in 1:(Z-1)){
            for(kk in 1:2){
              lambda[kk,z]~T(dt(0,1, df=1), 0,)
              b[kk,z]~ddexp(0.0, lambda[kk,z])
            }
          }
        } else{
          for(z in 1:(Z-1)){
            lambda[z]~T(dt(0,1, df=1), 0,)
            b[z]~ddexp(0.0, lambda[z])
          }
        }
        
        if(Mlen>1){
          if(deg==1){
            for(z in 1:(Z-1)){
              for(zz in 1:Mlen){
                lambda_m[zz,z]~T(dt(0,1, df=1), 0,)
                bm[zz,z]~ddexp(0.0, lambda_m[zz,z])
              }
            }
          } else if(deg>1){
            for(z in 1:(Z-1)){
              for(zz in 1:Mlen){
                lambda_m1[zz,z]~T(dt(0,1, df=1), 0,)
                bm1[zz,z]~ddexp(0.0, lambda_m1[zz,z])
                lambda_m2[zz,z]~T(dt(0,1, df=1), 0,)
                bm2[zz,z]~ddexp(0.0, lambda_m2[zz,z])
              }
            }
          }
        } else if(Mlen==1){
          if(deg==1){
            for(z in 1:(Z-1)){
              lambda_m[z]~T(dt(0,1, df=1), 0,)
              bm[z]~ddexp(0.0, lambda_m[z])
            }
          } else if(deg>1){
            for(z in 1:(Z-1)){
              lambda_m1[z]~T(dt(0,1, df=1), 0,)
              bm1[z]~ddexp(0.0, lambda_m1[z])
              lambda_m2[z]~T(dt(0,1, df=1), 0,)
              bm2[z]~ddexp(0.0, lambda_m2[z])
            }
          }
        }
      }
    } else if(doseRRmod=="LINEXP"){ 
      if(ERRprior=="truncated_horseshoe"){
        tau~T(dt(0,1, df=1), 0,)
        for(z in 1:(Z-1)){
          # horseshoe prior
          lambda[1,z]~T(dt(0,1, df=1), 0,)
          b[1,z]~T(dnorm(0,sd=lambda[1,z]*tau), 0,)
          lambda[2,z]~T(dt(0,1, df=1), 0,)
          b[2,z]~dnorm(0,sd=lambda[2,z]*tau)
          
        }
        
        if(Mlen>1){
          for(z in 1:(Z-1)){
            for(zz in 1:Mlen){
              # horseshoe prior
              lambda_m1[zz,z]~T(dt(0,1, df=1), 0,)
              bm1[zz,z]~T(dnorm(0,sd=lambda_m1[zz,z]*tau), 0,)
              lambda_m2[zz,z]~T(dt(0,1, df=1), 0,)
              bm2[zz,z]~dnorm(0,sd=lambda_m2[zz,z]*tau)
            }
          }
        } else if(Mlen==1){
          for(z in 1:(Z-1)){
            lambda_m1[z]~T(dt(0,1, df=1), 0,)
            bm1[z]~T(dnorm(0, sd=lambda_m1[z]*tau), 0,)
            lambda_m2[z]~T(dt(0,1, df=1), 0,)
            bm2[z]~dnorm(0, sd=lambda_m2[z]*tau)
          }
        }
      } else if(ERRprior=="truncated_normal"){
        for(z in 1:(Z-1)){
          
          b[1,z]~T(dnorm(0,0.001), 0, )
          b[2,z]~dnorm(0,0.001)
          
        }
        if(Mlen>1){
          for(z in 1:(Z-1)){
            for(k in 1:Mlen){
              bm1[k,z]~T(dnorm(0, 0.001), 0, )
              bm2[k,z]~dnorm(0, 0.001)
            }
          }
        } else if(Mlen==1){
          for(z in 1:(Z-1)){
            bm1[z]~T(dnorm(0, 0.001), 0, )
            bm2[z]~dnorm(0, 0.001)
          }
        }
      } else if(ERRprior=="truncated_doubleexponential"){
        for(z in 1:(Z-1)){
          
          lambda[1,z]~T(dt(0,1, df=1), 0,)
          b[1,z]~T(ddexp(0.0, lambda[1,z]) , 0,)
          lambda[2,z]~T(dt(0,1, df=1), 0,)
          b[2,z]~ddexp(0.0, lambda[2,z])
          
        }
        if(Mlen>1){
          for(z in 1:(Z-1)){
            for(zz in 1:Mlen){
              lambda_m1[zz,z]~T(dt(0,1, df=1), 0,)
              bm1[zz,z]~T(ddexp(0.0, lambda_m1[zz,z]), 0,)
              lambda_m2[zz,z]~T(dt(0,1, df=1), 0,)
              bm2[zz,z]~ddexp(0.0, lambda_m2[zz,z])
            }
          }
        } else if(Mlen==1){
          for(z in 1:(Z-1)){
            lambda_m1[z]~T(dt(0,1, df=1), 0,)
            bm1[z]~T(ddexp(0.0, lambda_m1[z]), 0,)
            lambda_m2[z]~T(dt(0,1, df=1), 0,)
            bm2[z]~ddexp(0.0, lambda_m2[z])
          }
        }
      } else if(ERRprior=="horseshoe"){
        tau~T(dt(0,1, df=1), 0,)
        for(z in 1:(Z-1)){
          for(kk in 1:2){ # horseshoe prior
            lambda[kk,z]~T(dt(0,1, df=1), 0,)
            b[kk,z]~dnorm(0,sd=lambda[kk,z]*tau)
          }
        }
        
        if(Mlen>1){
          for(z in 1:(Z-1)){
            for(zz in 1:Mlen){
              # horseshoe prior
              lambda_m1[zz,z]~T(dt(0,1, df=1), 0,)
              bm1[zz,z]~dnorm(0,sd=lambda_m1[zz,z]*tau)
              lambda_m2[zz,z]~T(dt(0,1, df=1), 0,)
              bm2[zz,z]~dnorm(0,sd=lambda_m2[zz,z]*tau)
            }
          }
        } else if(Mlen==1){
          for(z in 1:(Z-1)){
            lambda_m1[z]~T(dt(0,1, df=1), 0,)
            bm1[z]~dnorm(0, sd=lambda_m1[z]*tau)
            lambda_m2[z]~T(dt(0,1, df=1), 0,)
            bm2[z]~dnorm(0, sd=lambda_m2[z]*tau)
          }
        }
      } else if(ERRprior=="normal"){
        for(z in 1:(Z-1)){
          for(k in 1:2){
            b[k,z]~dnorm(0,0.001)
          }
        }
        
        if(Mlen>1){
          for(z in 1:(Z-1)){
            for(k in 1:Mlen){
              bm1[k,z]~dnorm(0, 0.001)
              bm2[k,z]~dnorm(0, 0.001)
            }
          }
        } else if(Mlen==1){
          for(z in 1:(Z-1)){
            bm1[z]~dnorm(0, 0.001)
            bm2[z]~dnorm(0, 0.001)
          }
        }
      } else if(ERRprior=="doubleexponential"){
        for(z in 1:(Z-1)){
          for(kk in 1:2){
            lambda[kk,z]~T(dt(0,1, df=1), 0,)
            b[kk,z]~ddexp(0.0, lambda[kk,z])
          }
        }
        
        if(Mlen>1){
          for(z in 1:(Z-1)){
            for(zz in 1:Mlen){
              lambda_m1[zz,z]~T(dt(0,1, df=1), 0,)
              bm1[zz,z]~ddexp(0.0, lambda_m1[zz,z])
              lambda_m2[zz,z]~T(dt(0,1, df=1), 0,)
              bm2[zz,z]~ddexp(0.0, lambda_m2[zz,z])
            }
          }
        } else if(Mlen==1){
          for(z in 1:(Z-1)){
            lambda_m1[z]~T(dt(0,1, df=1), 0,)
            bm1[z]~ddexp(0.0, lambda_m1[z])
            lambda_m2[z]~T(dt(0,1, df=1), 0,)
            bm2[z]~ddexp(0.0, lambda_m2[z])
          }
        }
      }
    }
  }
  
  if(family=="gaussian"){
    sigma~T(dnorm(0,0.001),0,)
  }
  if(family=="prophaz"){
    for(j in 1:prophaz_numints) {
      h0[j] ~ dgamma(0.01, 0.01)
    }
  }
  
  # Likelihood
  
  if(family=="prophaz"){
    for (i in 1:N){
      # Pieces of the cumulative baseline hazard function
      for(k in int.entry[i]:int.exit[i]) {
        start_t[i,k] <- max(prophaz_timepoints[k], entry[i])
        end_t[i,k]   <- min(prophaz_timepoints[k+1], exit[i])
        HH[i,k] <- max(end_t[i,k] - start_t[i,k], 0) * h0[k]
      }
      # Cumulative baseline hazard
      H[i] <- sum(HH[i, int.entry[i]:int.exit[i]])
    }
  }
  
  for (i in 1:N){
    dose[i] <- inprod(dosemat[i,1:K], vec[1:K])
    if(family!="multinomial"){
      if(doseRRmod != "LINEXP"){
        if(Mlen==0){
          if(deg==2){
            dosepart[i] <-b[1]*dose[i]+b[2]*dose[i]*dose[i]
          } else if(deg==1){
            dosepart[i] <-b*dose[i]
          }
        } else if (Mlen==1){
          if(deg==2){
            dosepart[i] <-(b[1]+bm1*Mmat[i])*dose[i]+(b[2]+bm2*Mmat[i])*dose[i]*dose[i]
          } else if(deg==1){
            dosepart[i] <-(b+bm*Mmat[i])*dose[i]
          }
        } else{
          if(deg==2){
            dosepart[i] <- (b[1]+inprod(Mmat[i,], bm1[1:Mlen]))*dose[i]+(b[2]+inprod(Mmat[i,], bm2[1:Mlen]))*dose[i]*dose[i]
          } else if(deg==1){
            dosepart[i] <- (b+inprod(Mmat[i,], bm[1:Mlen]))*dose[i]
          }
        }
      } else { # LINEXP
        if(Mlen==0){
          dosepart[i] <-b[1]*dose[i] * exp(b[2]*dose[i])
        } else if (Mlen==1){
          dosepart[i] <-(b[1]+bm1*Mmat[i])*dose[i]*exp((b[2]+bm2*Mmat[i])*dose[i])
        } else{
          dosepart[i] <- (b[1]+inprod(Mmat[i,], bm1[1:Mlen]))*dose[i]*exp((b[2]+inprod(Mmat[i,], bm2[1:Mlen]))*dose[i])
        }
      }
      
      if(Xlen>1){
        Xlinpred[i] <- a0+inprod(Xmat[i,], a[1:Xlen])
      } else if(Xlen==1){
        Xlinpred[i] <- a0+Xmat[i]*a
      } else if(Xlen==0){
        Xlinpred[i] <- a0
      }
      if(family=="gaussian"){
        mu[i] <- dosepart[i] + Xlinpred[i]
        Y[i]~dnorm(mu[i], sd=sigma)
      } else if(family=="binomial"){
        if(doseRRmod%in%c("ERR", "LINEXP")){
          mu[i] <- (1+dosepart[i])*exp(Xlinpred[i])
        } else if(doseRRmod=="EXP"){
          mu[i] <- exp(dosepart[i]+Xlinpred[i])
        }
        prob[i] <- mu[i]/(1+mu[i])
        Y[i] ~ dbinom(prob[i], size = 1)
      } else if(family=="poisson"){
        if(doseRRmod%in%c("ERR", "LINEXP")){
          mu[i] <- (1+dosepart[i])*exp(Xlinpred[i])*P[i]
        } else if(doseRRmod=="EXP"){
          mu[i] <- exp(dosepart[i]+Xlinpred[i])*P[i]
        }
        Y[i]~dpois(mu[i])
      } else if(family=="prophaz"){
        # Relative risk
        if(doseRRmod%in%c("ERR", "LINEXP")){
          R[i] <- (1+dosepart[i])*exp(Xlinpred[i])
        } else if(doseRRmod=="EXP"){
          R[i] <- exp(dosepart[i]+Xlinpred[i])
        }
        
        # Log-hazard
        logHaz[i] <- log(h0[int.exit[i]] * R[i])
        # Log-survival
        logSurv[i] <- -H[i]*R[i]
        
        # Using the zeros trick
        phi[i] <- 100000-delta[i]*logHaz[i]-logSurv[i]
        zeros[i] ~ dpois(phi[i])
      } else if(family=="clogit"){
        # Relative risk
        if(doseRRmod%in%c("ERR", "LINEXP")){
          R[i] <- (1+dosepart[i])*exp(Xlinpred[i])
        } else if(doseRRmod=="EXP"){
          R[i] <- exp(dosepart[i]+Xlinpred[i])
        }
      }
    } else{ # multinomial
      for(z in 1:(Z-1)){
        if(doseRRmod != "LINEXP"){
          if(Mlen==0){
            if(deg==2){
              dosepart[i,z] <-b[1,z]*dose[i]+b[2,z]*dose[i]*dose[i]
            } else if(deg==1){
              dosepart[i,z] <-b[z]*dose[i]
            }
          } else if (Mlen==1){
            if(deg==2){
              dosepart[i,z] <-(b[1,z]+bm1[z]*Mmat[i])*dose[i]+(b[2,z]+bm2[z]*Mmat[i])*dose[i]*dose[i]
            } else if(deg==1){
              dosepart[i,z] <-(b[z]+bm[z]*Mmat[i])*dose[i]
            }
          } else{
            if(deg==2){
              dosepart[i,z] <- (b[1,z]+inprod(Mmat[i,], bm1[1:Mlen,z]))*dose[i]+(b[2,z]+inprod(Mmat[i,], bm2[1:Mlen,z]))*dose[i]*dose[i]
            } else if(deg==1){
              dosepart[i,z] <- (b[z]+inprod(Mmat[i,], bm[1:Mlen,z]))*dose[i]
            }
          }
        } else { # LINEXP
          if(Mlen==0){
            dosepart[i,z] <-b[1,z]*dose[i] * exp(b[2,z]*dose[i])
          } else if (Mlen==1){
            dosepart[i,z] <-(b[1,z]+bm1[z]*Mmat[i])*dose[i]*exp((b[2,z]+bm2[z]*Mmat[i])*dose[i])
          } else{
            dosepart[i,z] <- (b[1,z]+inprod(Mmat[i,], bm1[1:Mlen,z]))*dose[i]*exp((b[2,z]+inprod(Mmat[i,], bm2[1:Mlen,z]))*dose[i])
          }
        }
        
        if(Xlen>1){
          Xlinpred[i,z] <- a0[z]+inprod(Xmat[i,], a[1:Xlen,z])
        } else if(Xlen==1){
          Xlinpred[i,z] <- a0[z]+Xmat[i]*a[z]
        } else if(Xlen==0){
          Xlinpred[i,z] <- a0[z]
        }
        
        
        if(doseRRmod%in%c("ERR", "LINEXP")){
          mu[i,z] <- (1+dosepart[i,z])*exp(Xlinpred[i,z])
        } else if(doseRRmod=="EXP"){
          mu[i,z] <- exp(dosepart[i,z]+Xlinpred[i,z])
        }
        
      }
      
      mu[i, Z] <- 1
      probs[i,1:Z] <- mu[i,1:Z]/sum(mu[i, 1:Z])
      Y[i,1:Z] ~ dmulti(prob=probs[i,1:Z], size = 1)
      
    }
  }
  
  if(family=="clogit"){
    
    for (setnr in 1:nsets){
      #probs[setnr, 1:N] <- R[1:N]*setmat[setnr,1:N]/inprod(setmat[setnr, 1:N], R[1:N])
      #Ymat[setnr, 1:N] ~ dmulti(prob=probs[setnr, 1:N], size=1)
      probs[set_start[setnr]:set_end[setnr]] <- R[set_start[setnr]:set_end[setnr]]/sum(R[set_start[setnr]:set_end[setnr]])
      Y[set_start[setnr]:set_end[setnr]] ~ dmulti(probs[set_start[setnr]:set_end[setnr]], size = 1)
    }
  }
})

ameras.bma <- function(family, dosevars, data, deg, Y=NULL, M=NULL, X=NULL, offset=NULL, entry=NULL, exit=NULL, status=NULL, setnr=setnr, transform=NULL, inpar=NULL, doseRRmod=NULL, ERRprior="doubleexponential",prophaz_numints=10, nburnin=1000, niter=5000, included.replicates=1:length(dosevars), nchains=2, thin=10, optim.method="Nelder-Mead", ...){
  
  # Remove build warnings, local functions may need to be pulled out from this function
  # HPDinterval <- K <- Mlen <- Mmat <- N <- Xlen <- Xmat <- a <- as.mcmc <- 
  #   b <- bm <- bm1 <- bm2 <- buildMCMC <- nsets <- setmat <-
  #   col.ind <- compileNimble <- configureMCMC <- delta <- 
  #   dosemat <- equals <- h0 <- inprod <- nimbleModel <- runMCMC <- NULL
  
  
  
  ndoses <- length(dosevars)
  
  if(family!="gaussian"){
    if(doseRRmod%in%c("ERR","LINEXP") & is.null(ERRprior)) {
      stop("Please specify prior for ERR parameters")
    } else if(doseRRmod%in%c("ERR","LINEXP") & !is.null(ERRprior)){
      if(!(ERRprior %in% c("truncated_normal", "truncated_horseshoe", "truncated_doubleexponential", "normal", "horseshoe", "doubleexponential"))) stop("Incorrect prior for ERR parameters specified, should be truncated_normal, truncated_horseshoe, truncated_doubleexponential, normal, horseshoe, or doubleexponential")
    }
  }   
  
  
  
  
  
  
  t0 <- proc.time()
  
  if(family=="gaussian"){
    
    if(is.null(Y)) stop("Y is required for family=gaussian")
    
    
    if(is.null(included.replicates)){
      if(is.null(inpar)){
        inpar <- rep(0, 2+length(X)+length(M)*deg+deg)
      }
      
      toInclude <-  sapply(1:length(dosevars), function(Xi){
        
        fit.FMAi <- optim(inpar, loglik.gaussian, D=dosevars[Xi], X=X, Y=Y, M=M, data=data, deg=deg, ERC=FALSE, transform=transform, method=optim.method, ...)
        if(optim.method=="Nelder-Mead"){
          fit.FMAi <- optim(fit.FMAi$par, loglik.gaussian, D=dosevars[Xi], X=X, Y=Y, M=M, data=data, deg=deg, ERC=FALSE, transform=transform, method="BFGS", ...)
        }
        fit.FMAi$hessian <- numDeriv::hessian(func=loglik.gaussian, x=fit.FMAi$par, D=dosevars[Xi], X=X, Y=Y, M=M, data=data, deg=deg, ERC=FALSE, transform=transform, ...)
        
        
        if(det(fit.FMAi$hessian)!=0 & rcond(fit.FMAi$hessian)>.Machine$double.eps & fit.FMAi$convergence==0 &  all(eigen(fit.FMAi$hessian)$values > 0)){
          include <- TRUE
        } else{
          include <- FALSE
        }
        
        return(include)
        
      })
      included.replicates <- which(toInclude==TRUE)
    } 
    dosevars <- dosevars[included.replicates]
    
    
    nimbledata <- list(Y=data[,Y], dosemat=data[,dosevars])
    
    
    nimbleinits <- function(){
      L <- list(a0=rnorm(1), b=rexp(deg), sigma=rexp(1), col.ind=sample(1:length(dosevars),1))
      
      if(length(X)>0){
        L <- c(L, list(a=rnorm(length(X))))
      }
      
      if(length(M)>0){
        if(deg==1){
          L <- c(L, list(bm=rexp(length(M))))
        } else if(deg>1){
          L <- c(L, list(bm1=rexp(length(M)), bm2=rexp(length(M))))
        }
      }
      L
    }
    
    
  } else if(family=="binomial"){
    
    if(is.null(Y)) stop("Y is required for family=binomial")
    if(is.null(doseRRmod)) stop("doseRRmod is required for family=binomial")
    
    if(is.null(included.replicates)){
      if(is.null(inpar)){
        inpar <- rep(0, 1+length(X)+length(M)*deg+deg)
      }
      
      toInclude <-  sapply(1:length(dosevars), function(Xi){
        
        fit.FMAi <- optim(inpar, loglik.binomial, D=dosevars[Xi], X=X, Y=Y, M=M, doseRRmod=doseRRmod, data=data, deg=deg, ERC=FALSE, transform=transform, method=optim.method, ...)
        if(optim.method=="Nelder-Mead"){
          fit.FMAi <- optim(fit.FMAi$par, loglik.binomial, D=dosevars[Xi], X=X, Y=Y, M=M, doseRRmod=doseRRmod, data=data, deg=deg, ERC=FALSE, transform=transform, method="BFGS", ...)
        }
        fit.FMAi$hessian <- numDeriv::hessian(func=loglik.binomial, x=fit.FMAi$par, D=dosevars[Xi], X=X, Y=Y, M=M, doseRRmod=doseRRmod, data=data, deg=deg, ERC=FALSE, transform=transform, ...)
        
        if(det(fit.FMAi$hessian)!=0 & rcond(fit.FMAi$hessian)>.Machine$double.eps & fit.FMAi$convergence==0 &  all(eigen(fit.FMAi$hessian)$values > 0)){
          include <- TRUE
        } else{
          include <- FALSE
        }
        
        return(include)
      })
      included.replicates <- which(toInclude==TRUE)
    }
    
    dosevars <- dosevars[included.replicates]
    
    
    nimbledata <- list(Y=data[,Y], dosemat=data[,dosevars])
    
    
    
    nimbleinits <- function(){
      L <- list(a0=rnorm(1), b=rexp(deg), col.ind=sample(1:length(dosevars),1))
      
      if(length(X)>0){
        L <- c(L, list(a=rnorm(length(X))))
      }
      
      if(length(M)>0){
        if(deg==1){
          L <- c(L, list(bm=rexp(length(M))))
        } else if(deg>1){
          L <- c(L, list(bm1=rexp(length(M)), bm2=rexp(length(M))))
        }
      }
      L
    }
    
  } else if(family=="poisson"){
    
    if(is.null(Y)) stop("Y is required for family=poisson")
    if(is.null(doseRRmod)) stop("doseRRmod is required for family=poisson")
    
    if(is.null(offset)){
      P <- rep(1, nrow(data))
    } else{
      P <- data[,offset]
    }
    
    
    if(is.null(included.replicates)){
      if(is.null(inpar)){
        inpar <- rep(0, 1+length(X)+length(M)*deg+deg)
      }
      
      toInclude <-  sapply(1:length(dosevars), function(Xi){
        
        fit.FMAi <- optim(inpar, loglik.poisson, D=dosevars[Xi], X=X, Y=Y, M=M, offset=offset, doseRRmod=doseRRmod, data=data, deg=deg, transform=transform, method=optim.method, ...)
        if(optim.method=="Nelder-Mead"){
          fit.FMAi <- optim(fit.FMAi$par, loglik.poisson, D=dosevars[Xi], X=X, Y=Y, M=M, offset=offset, doseRRmod=doseRRmod, data=data, deg=deg, transform=transform, method="BFGS", ...)
        }
        fit.FMAi$hessian <- numDeriv::hessian(func=loglik.poisson, x=fit.FMAi$par, D=dosevars[Xi], X=X, Y=Y, M=M, offset=offset, doseRRmod=doseRRmod, data=data, deg=deg, transform=transform, ...)
        
        if(det(fit.FMAi$hessian)!=0 & rcond(fit.FMAi$hessian)>.Machine$double.eps & fit.FMAi$convergence==0 &  all(eigen(fit.FMAi$hessian)$values > 0)){
          include <- TRUE
        } else{
          include <- FALSE
        }
        
        return(include)
      })
      included.replicates <- which(toInclude==TRUE)
    }
    dosevars <- dosevars[included.replicates]
    
    
    nimbledata <- list(Y=data[,Y], dosemat=data[,dosevars])
    
    
    
    nimbleinits <- function(){
      L <- list(a0=rnorm(1), b=rexp(deg), col.ind=sample(1:length(dosevars),1))
      
      if(length(X)>0){
        L <- c(L, list(a=rnorm(length(X))))
      }
      
      if(length(M)>0){
        if(deg==1){
          L <- c(L, list(bm=rexp(length(M))))
        } else if(deg>1){
          L <- c(L, list(bm1=rexp(length(M)), bm2=rexp(length(M))))
        }
      }
      L
    }
    
    
  } else if(family=="prophaz"){
    
    if(is.null(exit)) stop("exit is required for family=prophaz")
    if(is.null(doseRRmod)) stop("doseRRmod is required for family=prophaz")
    if(is.null(status)) stop("status is required for family=prophaz")
    
    
    if(is.null(included.replicates)){
      
      if(is.null(inpar)){
        inpar <- rep(0, length(X)+length(M)*deg+deg)
      }
      
      toInclude <- sapply(1:length(dosevars), function(Xi){
        
        if(length(X)+length(M)*deg+deg == 1){ # Optimize 1-dimensional model: use optimize instead of optim
          fit0 <- optimize(f=loglik.prophaz, lower=-20, upper=5, D=dosevars[Xi], status=status, X=X, M=M, entry=entry,exit=exit, doseRRmod=doseRRmod, data=data, deg=deg, ERC=FALSE, transform=transform, ...)
          fit.FMAi <- list(par=fit0$minimum, value=fit0$objective, convergence=0, hessian=numDeriv::hessian(func=loglik.prophaz, x=fit0$minimum, D=dosevars[Xi], status=status, X=X, M=M, entry=entry,exit=exit, doseRRmod=doseRRmod, data=data, deg=deg, ERC=FALSE, transform=transform, ...))
        } else {
          fit.FMAi <- optim(inpar, loglik.prophaz, D=dosevars[Xi], status=status, X=X, M=M, entry=entry,exit=exit, doseRRmod=doseRRmod, data=data, deg=deg, ERC=FALSE, transform=transform, method=optim.method, ...)
          if(optim.method=="Nelder-Mead"){
            fit.FMAi <- optim(fit.FMAi$par, loglik.prophaz, D=dosevars[Xi], status=status, X=X, M=M, entry=entry,exit=exit, doseRRmod=doseRRmod, data=data, deg=deg, ERC=FALSE, transform=transform, method="BFGS", ...)
          }
          fit.FMAi$hessian <- numDeriv::hessian(func=loglik.prophaz, x=fit.FMAi$par, D=dosevars[Xi], status=status, X=X, M=M, entry=entry,exit=exit, doseRRmod=doseRRmod, data=data, deg=deg, ERC=FALSE, transform=transform, ...)
        }
        
        
        if(det(fit.FMAi$hessian)!=0 & rcond(fit.FMAi$hessian)>.Machine$double.eps & fit.FMAi$convergence==0 &  all(eigen(fit.FMAi$hessian)$values > 0)){
          include <- TRUE
        } else{
          include <- FALSE
        }
        
        return(include)
        
      })
      included.replicates <- which(toInclude==TRUE)
    }
    
    dosevars <- dosevars[included.replicates]
    
    prophaz_timepoints <- as.numeric(quantile(data[data[,status]==1,exit], probs=seq(0,1, length.out=prophaz_numints+1))) # define cut points using quantiles of event time distribution
    if(is.null(entry)){
      prophaz_timepoints[1] <- 0
    } else{
      prophaz_timepoints[1] <- min(data[,entry])
    }
    
    prophaz_timepoints[prophaz_numints+1] <- max(data[,exit])+.001
    
    if(!is.null(entry)){
      int.entry <- as.numeric(cut(data[,entry], breaks=prophaz_timepoints, right=FALSE)) # indicator for which interval entry time belongs to
    } else{
      int.entry <- rep(1, nrow(data))
    }
    
    int.exit <- as.numeric(cut(data[,exit], breaks=prophaz_timepoints, right=FALSE)) # indicator for which interval exit time belongs to
    
    
    nimbledata <- list(delta=data[,status], exit=data[,exit], dosemat=data[,dosevars], zeros=rep(0, nrow(data)))
    if(!is.null(entry)){
      nimbledata <- c(nimbledata, list(entry=data[,entry]))
    } else{
      nimbledata <- c(nimbledata, list(entry=rep(0, nrow(data))))
    }
    
    nimbleinits <- function(){
      L <- list(b=rexp(deg), col.ind=sample(1:length(dosevars),1), h0=runif(prophaz_numints, .1))
      
      if(length(X)>0){
        L <- c(L, list(a=rnorm(length(X))))
      }
      
      if(length(M)>0){
        if(deg==1){
          L <- c(L, list(bm=rexp(length(M))))
        } else if(deg>1){
          L <- c(L, list(bm1=rexp(length(M)), bm2=rexp(length(M))))
        }
      }
      L
    }
    
  } else if(family=="clogit"){
    
    
    if(is.null(doseRRmod)) stop("doseRRmod is required for family=clogit")
    
    # Remove sets of size 1
    set_counts <- table(data[,setnr])
    valid_sets <- as.numeric(names(set_counts[set_counts > 1]))
    data <- data[data[,setnr] %in% valid_sets, ]
    
    # Reorder and determine indexing for nimble
    data <- data[order(data[,setnr]),]
    nsets <- length(unique(data[,setnr]))
    set_sizes <- as.numeric(table(data[,setnr]))
    set_start <- cumsum(c(1, head(set_sizes, -1)))
    set_end <- cumsum(set_sizes)
    
    
    if(is.null(included.replicates)){
      
      designmat <- t(model.matrix(~as.factor(data[,setnr])-1))
      
      if(is.null(inpar)){
        inpar <- rep(0, length(X)+length(M)*deg+deg)
      }
      
      toInclude <- sapply(1:length(dosevars), function(Xi){
        
        if(length(X)+length(M)*deg+deg == 1){ # Optimize 1-dimensional model: use optimize instead of optim
          fit0 <- optimize(f=loglik.clogit, lower=-20, upper=5, D=dosevars[Xi], status=status, X=X, M=M, designmat=designmat, doseRRmod=doseRRmod, data=data, deg=deg, ERC=FALSE, transform=transform, ...)
          fit.FMAi <- list(par=fit0$minimum, value=fit0$objective, convergence=0, hessian=numDeriv::hessian(func=loglik.clogit, x=fit0$minimum, D=dosevars[Xi], status=status, X=X, M=M, designmat=designmat, doseRRmod=doseRRmod, data=data, deg=deg, ERC=FALSE, transform=transform, ...))
        } else {
          fit.FMAi <- optim(inpar, loglik.clogit, D=dosevars[Xi], status=status, X=X, M=M, designmat=designmat, doseRRmod=doseRRmod, data=data, deg=deg, ERC=FALSE, transform=transform, method=optim.method, ...)
          if(optim.method=="Nelder-Mead"){
            fit.FMAi <- optim(fit.FMAi$par, loglik.clogit, D=dosevars[Xi], status=status, X=X, M=M, designmat=designmat, doseRRmod=doseRRmod, data=data, deg=deg, ERC=FALSE, transform=transform, method="BFGS", ...)
          }
          fit.FMAi$hessian <- numDeriv::hessian(func=loglik.clogit, x=fit.FMAi$par, D=dosevars[Xi], status=status, X=X, M=M, designmat=designmat, doseRRmod=doseRRmod, data=data, deg=deg, ERC=FALSE, transform=transform, ...)
        }
        
        
        if(det(fit.FMAi$hessian)!=0 & rcond(fit.FMAi$hessian)>.Machine$double.eps & fit.FMAi$convergence==0 &  all(eigen(fit.FMAi$hessian)$values > 0)){
          include <- TRUE
        } else{
          include <- FALSE
        }
        
        return(include)
        
      })
      included.replicates <- which(toInclude==TRUE)
    }
    
    dosevars <- dosevars[included.replicates]
    
    
    
    nimbledata <- list(Y=as.numeric(data[,status]), dosemat=data[,dosevars])
    
    
    nimbleinits <- function(){
      L <- list(b=rexp(deg), col.ind=sample(1:length(dosevars),1))
      
      if(length(X)>0){
        L <- c(L, list(a=rnorm(length(X))))
      }
      
      if(length(M)>0){
        if(deg==1){
          L <- c(L, list(bm=rexp(length(M))))
        } else if(deg>1){
          L <- c(L, list(bm1=rexp(length(M)), bm2=rexp(length(M))))
        }
      }
      L
    }
    
  }  else if(family=="multinomial"){
    if(is.null(Y)) stop("Y is required for family=multinomial")
    if(is.null(doseRRmod)) stop("doseRRmod is required for family=multinomial")
    
    Z <- nlevels(data[,Y])
    if(is.null(included.replicates)){
      if(is.null(inpar)){
        inpar <- rep(0, (Z-1)*(1+length(X)+length(M)*deg+deg))
      }
      
      toInclude <-  sapply(1:length(dosevars), function(Xi){
        
        fit.FMAi <- optim(inpar, loglik.multinomial, D=dosevars[Xi], X=X, Y=Y, M=M, doseRRmod=doseRRmod, data=data, deg=deg, ERC=FALSE, transform=transform, method=optim.method, ...)
        if(optim.method=="Nelder-Mead"){
          fit.FMAi <- optim(fit.FMAi$par, loglik.multinomial, D=dosevars[Xi], X=X, Y=Y, M=M, doseRRmod=doseRRmod, data=data, deg=deg, ERC=FALSE, transform=transform, method="BFGS", ...)
        }
        fit.FMAi$hessian <- numDeriv::hessian(func=loglik.multinomial, x=fit.FMAi$par, D=dosevars[Xi], X=X, Y=Y, M=M, doseRRmod=doseRRmod, data=data, deg=deg, ERC=FALSE, transform=transform, ...)
        
        if(det(fit.FMAi$hessian)!=0 & rcond(fit.FMAi$hessian)>.Machine$double.eps & fit.FMAi$convergence==0 &  all(eigen(fit.FMAi$hessian)$values > 0)){
          include <- TRUE
        } else{
          include <- FALSE
        }
        
        return(include)
      })
      included.replicates <- which(toInclude==TRUE)
    }
    
    dosevars <- dosevars[included.replicates]
    
    
    nimbledata <- list(Y=model.matrix(~data[,Y]-1), dosemat=data[,dosevars])
    
    
    
    nimbleinits <- function(){
      L <- list(a0=rnorm(Z-1), b=matrix(rexp(deg*(Z-1)), ncol=Z-1), col.ind=sample(1:length(dosevars),1))
      if(deg==1){
        L$b <- as.vector(L$b)
      }
      
      if(length(X)>0){
        L <- c(L, list(a=matrix(rnorm((Z-1)*length(X)), ncol=Z-1)))
      }
      
      if(length(M)>0){
        if(deg==1){
          L <- c(L, list(bm=matrix(rexp((Z-1)*length(M)), ncol=Z-1)))
        } else if(deg>1){
          L <- c(L, list(bm1=matrix(rexp((Z-1)*length(M)), ncol=Z-1), bm2=matrix(rexp((Z-1)*length(M)), ncol=Z-1)))
        }
      }
      L
    }
  }
  
  
  if(!(family%in%c("prophaz", "clogit"))){
    mons <- c("a0", "b", "col.ind")
  } else{
    mons <- c("b", "col.ind")
  }
  
  if(family=="gaussian") doseRRmod <- "EXP"
  nimbleconst <- list(N=nrow(data),K=length(dosevars), deg=deg, Xlen=length(X), Mlen=length(M), w=rep(1/length(dosevars), length(dosevars)), doseRRmod=doseRRmod, family=family)
  if(family=="poisson") nimbleconst <- c(nimbleconst, list(P=P))
  if(family=="multinomial") nimbleconst <- c(nimbleconst, list(Z=Z))
  if(length(X)>0){
    nimbleconst <- c(nimbleconst, list(Xmat=data[,X]))
    mons <- c(mons, "a")
  }
  if(length(M)>0){
    nimbleconst <- c(nimbleconst, list(Mmat=data[,M]))
    if(deg==1){
      mons <- c(mons, "bm")
    } else if(deg>1){
      mons <- c(mons, "bm1", "bm2")
    }
  }
  if(family=="gaussian"){
    mons <- c(mons, "sigma")
  }
  if(family=="prophaz"){
    nimbleconst <- c(nimbleconst, list(prophaz_timepoints=prophaz_timepoints, int.exit=int.exit, int.entry=int.entry, prophaz_numints=prophaz_numints))
    mons <- c(mons, "h0")
  }
  if(family=="clogit"){
    nimbleconst <- c(nimbleconst, list(nsets = nsets, set_start = set_start, set_end = set_end))
  }
  
  if(doseRRmod%in%c("ERR", "LINEXP")) nimbleconst <- c(nimbleconst, list(ERRprior=ERRprior))
  
  mymod <- nimbleModel(nimblemod, data=nimbledata, constants=nimbleconst, inits=list(col.ind=sample(1:length(dosevars),1)))
  
  mymod_C <- compileNimble(mymod)
  
  mymod <- configureMCMC(mymod, monitors=mons, thin=thin, print=FALSE)
  
  mymod$removeSamplers('col.ind') 
  mymod$addSampler(target='col.ind', type=boundedSliceSampler, control=list(K=length(dosevars), adaptive=FALSE, sliceWidth=round(length(dosevars)/2)))
  
  mymod_MCMC <- buildMCMC(mymod)
  mymod_compiled <- compileNimble(mymod_MCMC, project=mymod_C)
  
  
  
  nimblesamples <- runMCMC(mymod_compiled, nburnin=nburnin, niter=niter, nchains=nchains, inits=nimbleinits)
  
  if(family!="multinomial"){
    if(!(family%in%c("prophaz", "clogit"))){
      pars <- "a0"
    } else{
      pars <- NULL
    }
    
    if(length(X)==1) {
      pars <- c(pars,"a")
    } else if(length(X)>1){
      pars <- c(pars, paste0("a[",1:length(X),"]"))
    }
    if(deg==1){
      pars <- c(pars, "b")
    } else if(deg>1){
      pars <- c(pars, "b[1]","b[2]")
    }
    if(length(M)==1){
      if(deg==1){
        pars <- c(pars, "bm")
      } else if(deg>1){
        pars <- c(pars, "bm1","bm2")
      }
    } else if(length(M)>1){
      if(deg==1){
        pars <- c(pars, paste0("bm[",1:length(M),"]"))
      } else if(deg>1){
        pars <- c(pars, paste0("bm1[",1:length(M),"]"),paste0("bm2[",1:length(M),"]"))
      }
    }
    if(family=="gaussian") pars <- c(pars, "sigma")
    if(family=="prophaz") pars <- c(pars, paste0("h0[", 1:prophaz_numints,"]"))
  } else{ # multinomial
    pars <- NULL
    for(z in 1:(Z-1)){
      pars <- c(pars,paste0("a0[",z,"]"))
      
      if(length(X)==1) {
        pars <- c(pars,paste0("a[",z,"]"))
      } else if(length(X)>1){
        pars <- c(pars, paste0("a[",1:length(X),", ",z,"]"))
      }
      
      if(deg==1){
        pars <- c(pars, paste0("b[",z,"]"))
      } else{
        pars <- c(pars, paste0("b[1, ",z,"]"),paste0("b[2, ",z,"]"))
      }
      
      
      if(length(M)==1){
        if(deg==1){
          pars <- c(pars, paste0("bm[",z,"]"))
        } else{
          pars <- c(pars, paste0("bm1[",z,"]"),paste0("bm2[",z,"]"))
        }
        
        
      } else if(length(M)>1){
        if(deg==1){
          pars <- c(pars, paste0("bm[",1:length(M),", ",z,"]"))
        } else{
          pars <- c(pars, paste0("bm1[",1:length(M),", ",z,"]"),paste0("bm2[",1:length(M),", ",z,"]"))
        }
        
        
      }
    }
  }
  
  if(!(family%in%c("prophaz", "clogit"))){
    parnames <- "(Intercept)"
  } else{
    parnames <- NULL
  }
  if(!is.null(doseRRmod)){
    if(doseRRmod=="LINEXP"){
      parnames <- c(parnames,names(data[,X,drop=FALSE]),c("dose_linear","dose_exponential"))
      if(!is.null(M)){
        parnames <- c(parnames, paste0("dose_linear:",names(data[, M,drop=FALSE])))
        parnames <- c(parnames, paste0("dose_exponential:",names(data[, M,drop=FALSE])))
      }
    } else{
      parnames <- c(parnames,names(data[,X,drop=FALSE]),c("dose","dose_squared")[1:deg])
      if(!is.null(M)){
        parnames <- c(parnames, paste0("dose:",names(data[, M,drop=FALSE])))
        if(deg==2){
          parnames <- c(parnames, paste0("dose_squared:",names(data[, M,drop=FALSE])))
        }
      }
    }
  } else{
    parnames <- c(parnames,names(data[,X,drop=FALSE]),c("dose","dose_squared")[1:deg])
    if(!is.null(M)){
      parnames <- c(parnames, paste0("dose:",names(data[, M,drop=FALSE])))
      if(deg==2){
        parnames <- c(parnames, paste0("dose_squared:",names(data[, M,drop=FALSE])))
      }
    }
  }
  if(family=="gaussian") parnames <- c(parnames, "sigma")
  if(family=="prophaz") parnames <- c(parnames, paste0("h0[", 1:prophaz_numints,"]"))
  
  if(doseRRmod=="LINEXP"){
    parnames <- sub("\\bdose\\b", "dose_linear", parnames)
    parnames <- sub("\\bdose_squared\\b", "dose_exponential", parnames)
  }
  
  if(family=="multinomial"){
    mylv <- levels(data[,Y])
    
    mylv <- mylv[-length(mylv)]
    
    parnames <- do.call("c",lapply(mylv, function(lv) paste0("(",lv, ")_", parnames)))
    
  }
  
  if(nchains>1) {
    nimblesamples.stacked <- do.call("rbind", nimblesamples)
    for(ichain in 1:length(nimblesamples)){
      nimblesamples[[ichain]] <- nimblesamples[[ichain]][,c(pars, "col.ind")]
      colnames(nimblesamples[[ichain]]) <- c(parnames, "col.ind")
    }
    
  } else if(nchains==1){
    nimblesamples.stacked <- nimblesamples
    nimblesamples <- nimblesamples[,c(pars,"col.ind")]
    colnames(nimblesamples) <- c(parnames, "col.ind")
  }
  nimblesamples.stacked <- nimblesamples.stacked[,pars]
  
  coef <- colMeans(nimblesamples.stacked)
  names(coef) <- parnames
  
  sd <- apply(nimblesamples.stacked, 2, sd)
  names(sd) <- parnames
  
  mcmcsum <- MCMCsummary(nimblesamples)
  Rhat <- mcmcsum[-nrow(mcmcsum),c("Rhat", "n.eff")]
  if(nchains>1){
    if(any(Rhat$Rhat > 1.05)) warning("WARNING: Potential problems with MCMC convergence, consider using longer chains")
  } else {
    warning("WARNING: MCMC convergence cannot be assessed using a single chains")
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
  timedif <- t1-t0
  runtime <- paste(round(as.numeric(as.difftime(timedif["elapsed"], units="secs")),1), "seconds")
  
  prc_excluded <- round(100*(1-length(included.replicates)/ndoses), 1)
  if(length(included.replicates)/ndoses < .8) warning(paste0("WARNING: ",prc_excluded, "% of replicates excluded from model averaging. Try different bounds or starting values."))
  
  out <- list(coefficients=coef,
              sd=sd,
              #CI=CIresult,
              Rhat=Rhat,
              samples=nimblesamples,
              included.replicates=included.replicates,
              runtime=runtime)
  if(family=="prophaz"){
    out <- c(out, list(
      prophaz_timepoints=prophaz_timepoints
    ))
  }
  return(out)
}
