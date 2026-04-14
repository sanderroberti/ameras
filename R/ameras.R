
ameras_main <- function(family="gaussian", methods="RC", dosevars, data, deg=1, doseRRmod="ERR", transform=NULL,transform.jacobian=NULL, Y=NULL, M=NULL, X=NULL, offset=NULL, inpar=NULL, entry=NULL, exit=NULL, status=NULL, setnr=NULL, CI=c("proflik","percentile"), params.profCI="dose",maxit.profCI=20, tol.profCI=1e-2, unweightedFMA=FALSE, loglim=1e-30, MFMA=100000, prophaz.numints.BMA=10, ERRprior.BMA="doubleexponential", nburnin.BMA=5000, niter.BMA=20000, nchains.BMA=2, thin.BMA=10, included.replicates.BMA=1:length(dosevars), optim.method="Nelder-Mead", control=NULL, ... ){
  if(is.null(control)) control <- list(reltol=1e-10)
  
  if(!is.null(transform) & is.null(transform.jacobian)) stop("transform.jacobian is required when using a transformation")
  
  # Family in ("gaussian", "poisson", "prophaz", "binomial")
  # gaussian requires Y and can use X and M
  # binomial requires Y and can use X and M
  # poisson requires Y and can use offset, X and M
  # prophaz requires status and exit and can use entry, X and M
  
  # Method in ("MCML", "RC", "ERC", "FMA", "BMA")
  # If both FMA and BMA are to be run, run FMA first to determine the included replicates and skip that step in BMA
  if("FMA" %in% methods & "BMA" %in% methods) methods[methods %in% c("FMA","BMA")] <- c("FMA","BMA")
  
  # If both RC and ERC are to be run, run RC first to determine initial values for ERC (not implemented yet)
  if("RC" %in% methods & "ERC" %in% methods) methods[methods %in% c("RC","ERC")] <- c("RC","ERC")
  
  # CI in ("wald.orig","wald.transformed", "proflik") for method in ("RC", "ERC", "MCML")
  # CI in ("percentile", "hpd") for method in ("FMA", "BMA")
  
  # Input checks
  if(family!="gaussian" & is.null(doseRRmod)) stop("doseRRmod is required for all families other than Gaussian")
  
  if(family=="prophaz"){
    if(is.null(exit)) stop("exit is required for family=prophaz")
    if(is.null(status)) stop("status is required for family=prophaz")
  }
  
  # For the Gaussian family, sigma needs to be >0 and so use standard reparametrization if another is not defined
  if(family=="gaussian" & is.null(transform)){
    transform <- function(params,boundcheck=FALSE,...) {
      return(transform1(params=params, index.t=2+length(X)+length(M)*deg+deg, boundcheck=boundcheck, ...))
    }
    transform.jacobian <- function(params, ...){
      return(transform1.jacobian(params=params, index.t=2+length(X)+length(M)*deg+deg,...))
    }
  }
  
  # For linear ERR models, if no transformation is specified, use transform1 with lower limits -1/max(dose) for linear dose-response and (0,-1/max(dose^2)) for linear and linear-quadratic parameters, respectively. If effect modifiers M are specified, no transformation is used for those parameters. If negative RRs are evaluated, the program will produce an error and the user should specify a different transformation or bounds
  if(family!="gaussian"){
    if(doseRRmod=="ERR" & is.null(transform)){ 
      if(deg==1){
        lwlmt <- -1/max(data[,dosevars])
      } else{
        lwlmt <- c(0, -1/max(data[,dosevars]^2))
      }
      if(family != "multinomial"){
        transform <- function(params, boundcheck=FALSE, ...) {
          return(transform1(params=params, index.t=(1*(!(family%in%c("prophaz", "clogit")))+length(X)+1):(1*(!(family%in%c("prophaz", "clogit")))+length(X)+deg),lowlimit=lwlmt, boundcheck=boundcheck, ...))
        }
        transform.jacobian <- function(params, ...){
          return(transform1.jacobian(params=params, index.t=(1*(!(family%in%c("prophaz", "clogit")))+length(X)+1):(1*(!(family%in%c("prophaz", "clogit")))+length(X)+deg), lowlimit=lwlmt,...))
        }
      } else{
        lwlmt <- rep(lwlmt, length(levels(data[,Y]))-1)
        indx <- do.call("c",lapply(0:(length(levels(data[,Y]))-2), function(xx) xx*(1+length(X)+length(M)*deg+deg)+((1+length(X)+1):(1+length(X)+deg))))
        
        transform <- function(params, boundcheck=FALSE, ...) {
          return(transform1(params=params, index.t=indx,lowlimit=lwlmt, boundcheck=boundcheck, ...))
        }
        transform.jacobian <- function(params, ...){
          return(transform1.jacobian(params=params, index.t=indx, lowlimit=lwlmt,...))
        }
      }
    } else if(doseRRmod=="LINEXP" & is.null(transform)){ 
      lwlmt <- 0 # Lower bound of 0 for beta1, no bound for beta2
      
      if(family != "multinomial"){
        transform <- function(params, boundcheck=FALSE, ...) {
          return(transform1(params=params, index.t=(1*(!(family%in%c("prophaz", "clogit")))+length(X)+1),lowlimit=lwlmt, boundcheck=boundcheck, ...))
        }
        transform.jacobian <- function(params, ...){
          return(transform1.jacobian(params=params, index.t=(1*(!(family%in%c("prophaz", "clogit")))+length(X)+1), lowlimit=lwlmt,...))
        }
      } else{
        lwlmt <- rep(lwlmt, length(levels(data[,Y]))-1)
        indx <- do.call("c",lapply(0:(length(levels(data[,Y]))-2), function(xx) xx*(1+length(X)+length(M)*deg+deg)+(1+length(X)+1)))
        
        transform <- function(params, boundcheck=FALSE, ...) {
          return(transform1(params=params, index.t=indx,lowlimit=lwlmt, boundcheck=boundcheck, ...))
        }
        transform.jacobian <- function(params, ...){
          return(transform1.jacobian(params=params, index.t=indx, lowlimit=lwlmt,...))
        }
      }
    }
  }
  
  
  results <- NULL
  for(method in methods){
    if(method=="MCML"){
      message("Fitting MCML")
      fit <- ameras.mcml(family=family, dosevars=dosevars, data=data, deg=deg, transform=transform, transform.jacobian=transform.jacobian, Y=Y, M=M, X=X, offset=offset, inpar=inpar, entry=entry, exit=exit, status=status, setnr=setnr, CI=CI, params.profCI=params.profCI, maxit.profCI=maxit.profCI, tol.profCI=tol.profCI, doseRRmod=doseRRmod, loglim=loglim, control=control, optim.method=optim.method, ...)
      results <- c(results, list(MCML=fit))
    } else if(method=="RC"){
      message("Fitting RC")
      fit <- ameras.rc(family=family, dosevars=dosevars, data=data, deg=deg, ERC=FALSE, transform=transform, transform.jacobian=transform.jacobian, Y=Y, M=M, X=X, offset=offset, inpar=inpar, entry=entry, exit=exit, status=status, setnr=setnr, CI=CI, params.profCI=params.profCI, maxit.profCI=maxit.profCI, tol.profCI=tol.profCI, doseRRmod=doseRRmod, loglim=loglim,control=control, optim.method=optim.method, ...)
      results <- c(results, list(RC=fit))
    } else if(method=="ERC"){
      message("Fitting ERC")
      fit <- ameras.rc(family=family, dosevars=dosevars, data=data, deg=deg, ERC=TRUE, transform=transform, transform.jacobian=transform.jacobian, Y=Y, M=M, X=X, offset=offset, inpar=inpar, entry=entry, exit=exit, status=status,setnr=setnr,CI=CI, params.profCI=params.profCI, maxit.profCI=maxit.profCI, tol.profCI=tol.profCI, doseRRmod=doseRRmod, loglim=loglim,control=control, optim.method=optim.method, ...)
      results <- c(results, list(ERC=fit))
    } else if(method=="FMA"){
      message("Fitting FMA")
      fit <- ameras.fma(family=family, dosevars=dosevars, data=data, deg=deg, transform=transform,transform.jacobian=transform.jacobian, Y=Y, M=M, X=X, offset=offset, inpar=inpar, entry=entry, exit=exit, status=status, setnr=setnr,doseRRmod=doseRRmod,CI=CI, unweighted=unweightedFMA, MFMA=MFMA, control=control, ...)
      results <- c(results, list(FMA=fit))
    } else if(method=="BMA"){
      message("Fitting BMA")
      if(!is.null(included.replicates.BMA)){
        increps <- included.replicates.BMA
      } else if(!is.null(results$FMA$included.replicates)){
        increps <- results$FMA$included.replicates
      } else{
        increps <- NULL
      }
      fit <- ameras.bma(family=family, dosevars=dosevars, data=data, deg=deg, transform=transform, Y=Y, M=M, X=X, offset=offset, inpar=inpar, entry=entry, exit=exit, status=status,setnr=setnr, doseRRmod=doseRRmod,ERRprior=ERRprior.BMA, prophaz_numints=prophaz.numints.BMA, nburnin=nburnin.BMA, niter=niter.BMA, nchains=nchains.BMA, thin=thin.BMA, CI=CI, control=control, included.replicates=increps, optim.method=optim.method, ...)
      results <- c(results, list(BMA=fit))
    }
    
  }
  
  return(results)
}

ameras <- function(data, family="gaussian", Y, dosevars, M=NULL, X=NULL, offset=NULL, entry=NULL, exit=NULL, setnr=NULL,
                   methods="RC", deg=1, doseRRmod="ERR", transform=NULL, transform.jacobian=NULL, inpar=NULL, CI=c("proflik","percentile"),
                   params.profCI="dose", maxit.profCI=20, tol.profCI=1e-2, loglim=1e-30, MFMA=100000, prophaz.numints.BMA=10, 
                   ERRprior.BMA="doubleexponential", nburnin.BMA=5000, niter.BMA=20000, nchains.BMA=2, thin.BMA=10, included.replicates.BMA=1:length(dosevars), 
                   optim.method="Nelder-Mead", control=NULL, ... ){
  
  # Check for errors
  check_df(data) 
  check_family(family)
  if(family!="gaussian") check_doseRRmod(doseRRmod)
  check_Y(Y, data, family)
  
  methods <- check_methods(methods)
  
  if("BMA" %in% methods) message("Note: BMA may require extensive computation time in the order of multiple hours")
  
  check_D(dosevars, data, methods)
  check_M(M, data)
  check_X(X, data)
  if (family == "poisson") {
    check_offset(offset, data)
  }
  if (family == "prophaz") {
    check_entry_exit(entry, exit, data)
  }
  if (family == "clogit") {
    check_setnr(setnr, data)
  }
  
  
  deg     <- check_deg(deg)
  if(family != "multinomial"){
    inpar   <- check_inpar(inpar, family, M, X, deg)
  } else {
    inpar   <- check_inpar(inpar, family, M, X, deg, multinom_levels = length(levels(data[,Y])))
  }
  
  
  status <- NULL
  
  if(family=="clogit"){
    #data$exit <- 1
    status <- Y
    Y <- NULL
  }
  
  
  if (family == "prophaz") {
    status <- Y
    Y      <- NULL
  }
  
  if(!is.null(doseRRmod)){
    if(doseRRmod=="LINEXP"){
      deg    <- 2
    }
  }
  
  
  
  
  # Need variable numbers for M
  M <- getVarNumbers(M, data)
  ret <- c(list(call=match.call(),
                num.individuals=nrow(data),
                num.replicates=length(dosevars)),
           ameras_main(family, methods=methods, dosevars, data, deg, doseRRmod=doseRRmod, 
                       transform=transform, transform.jacobian=transform.jacobian, setnr=setnr,
                       Y=Y, M=M, X=X, offset=offset, inpar=inpar, entry=entry, exit=exit, status=status, 
                       CI=CI, params.profCI=params.profCI, maxit.profCI=maxit.profCI, tol.profCI=tol.profCI, 
                       loglim=loglim, MFMA=MFMA, prophaz.numints.BMA=prophaz.numints.BMA, 
                       ERRprior.BMA=ERRprior.BMA, nburnin.BMA=nburnin.BMA, niter.BMA=niter.BMA, 
                       nchains.BMA=nchains.BMA, thin.BMA=thin.BMA, 
                       included.replicates.BMA=included.replicates.BMA, control=control, 
                       optim.method=optim.method, ... )
  )
  
  ret <- new_amerasfit(ret)
  
  ret
}

