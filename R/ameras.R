
ameras_main <- function(family="gaussian", methods="RC", dosevars, data, deg=1, doseRRmod="ERR", transform=NULL,transform.jacobian=NULL, Y=NULL, M=NULL, X=NULL, offset=NULL, inpar=NULL, entry=NULL, exit=NULL, status=NULL, setnr=NULL, unweightedFMA=FALSE, loglim=1e-30, MFMA=100000, prophaz.numints.BMA=10, ERRprior.BMA="doubleexponential", nburnin.BMA=5000, niter.BMA=20000, nchains.BMA=2, thin.BMA=10, included.replicates.BMA=1:length(dosevars), optim.method="Nelder-Mead", control=NULL, ... ){
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
  
  
  results <- list(transform=transform,
                  transform.jacobian=transform.jacobian)
  for(method in methods){
    if(method=="MCML"){
      message("Fitting MCML")
      fit <- ameras.mcml(family=family, dosevars=dosevars, data=data, deg=deg, transform=transform, transform.jacobian=transform.jacobian, Y=Y, M=M, X=X, offset=offset, inpar=inpar, entry=entry, exit=exit, status=status, setnr=setnr, doseRRmod=doseRRmod, loglim=loglim, control=control, optim.method=optim.method, ...)
      results <- c(results, list(MCML=fit))
    } else if(method=="RC"){
      message("Fitting RC")
      fit <- ameras.rc(family=family, dosevars=dosevars, data=data, deg=deg, ERC=FALSE, transform=transform, transform.jacobian=transform.jacobian, Y=Y, M=M, X=X, offset=offset, inpar=inpar, entry=entry, exit=exit, status=status, setnr=setnr, doseRRmod=doseRRmod, loglim=loglim,control=control, optim.method=optim.method, ...)
      results <- c(results, list(RC=fit))
    } else if(method=="ERC"){
      message("Fitting ERC")
      fit <- ameras.rc(family=family, dosevars=dosevars, data=data, deg=deg, ERC=TRUE, transform=transform, transform.jacobian=transform.jacobian, Y=Y, M=M, X=X, offset=offset, inpar=inpar, entry=entry, exit=exit, status=status,setnr=setnr, doseRRmod=doseRRmod, loglim=loglim,control=control, optim.method=optim.method, ...)
      results <- c(results, list(ERC=fit))
    } else if(method=="FMA"){
      message("Fitting FMA")
      fit <- ameras.fma(family=family, dosevars=dosevars, data=data, deg=deg, transform=transform,transform.jacobian=transform.jacobian, Y=Y, M=M, X=X, offset=offset, inpar=inpar, entry=entry, exit=exit, status=status, setnr=setnr,doseRRmod=doseRRmod, unweighted=unweightedFMA, MFMA=MFMA, control=control, ...)
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
      fit <- ameras.bma(family=family, dosevars=dosevars, data=data, deg=deg, transform=transform, Y=Y, M=M, X=X, offset=offset, inpar=inpar, entry=entry, exit=exit, status=status,setnr=setnr, doseRRmod=doseRRmod,ERRprior=ERRprior.BMA, prophaz_numints=prophaz.numints.BMA, nburnin=nburnin.BMA, niter=niter.BMA, nchains=nchains.BMA, thin=thin.BMA, control=control, included.replicates=increps, optim.method=optim.method, ...)
      results <- c(results, list(BMA=fit))
    }
    
  }
  
  return(results)
}

ameras <- function(formula, data, family="gaussian",
                   methods="RC", transform=NULL, transform.jacobian=NULL, inpar=NULL, loglim=1e-30, 
                   MFMA=100000, prophaz.numints.BMA=10, ERRprior.BMA="doubleexponential", nburnin.BMA=5000, niter.BMA=20000, 
                   nchains.BMA=2, thin.BMA=10, included.replicates.BMA=NULL, 
                   optim.method="Nelder-Mead", control=NULL, keep.data=TRUE, ... ){
  
  # Check for errors
  check_df(data) 
  check_family(family)
  methods <- check_methods(methods)
  
  if("BMA" %in% methods) message("Note: BMA may require extensive computation time")
  
  parsed <- parse_ameras_formula(formula, data, family)

  if (!is.null(parsed$X_formula)) {
    X_matrix  <- model.matrix(parsed$X_formula, data=data)[, -1, drop=FALSE]
    X_colnames <- colnames(X_matrix)
    
    # Add expanded columns to data, avoiding name conflicts
    existing  <- intersect(X_colnames, colnames(data))
    new_cols  <- setdiff(X_colnames, colnames(data))
    
    if (length(new_cols)) {
      data[, new_cols] <- X_matrix[, new_cols, drop=FALSE]
    }
    
    X <- X_colnames
  } else {
    X <- NULL
  }
  
  if (family == "clogit") {
    Y      <- NULL
    status <- parsed$Y
  } else if (family == "prophaz") {
    Y      <- NULL
    status <- parsed$status
  } else {
    Y      <- parsed$Y
    status <- NULL
  }

  # Checks that depend on parsed formula
  if (family != "gaussian") check_doseRRmod(parsed$doseRRmod)
  check_Y(Y %||% status, data, family)
  check_D(parsed$dosevars, data, methods)
  check_M(parsed$M, data)
  check_X(X, data)
  
  if (family == "poisson")  check_offset(parsed$offset, data)
  if (family == "prophaz")  check_entry_exit(parsed$entry, parsed$exit, data)
  if (family == "clogit")   check_setnr(parsed$setnr, data)
  
  
  deg   <- check_deg(parsed$deg)
  
  if (family != "multinomial") {
    inpar <- check_inpar(inpar, family, parsed$M, parsed$X, deg)
  } else {
    inpar <- check_inpar(inpar, family, parsed$M, parsed$X, deg,
                         multinom_levels=length(levels(data[, parsed$Y])))
  }
  
  
  
  if (!is.null(parsed$doseRRmod) && parsed$doseRRmod == "LINEXP") {
    deg <- 2
  }
  


  M <- getVarNumbers(parsed$M, data)
  X <- getVarNumbers(X, data)
  
  # Add mean dose for RC and ERC to the data
  data$rcdose_ameras <- rowMeans(data[, parsed$dosevars, drop=FALSE])
  
  model_list <- list(
    data      = if (keep.data) data else NULL,
    keep.data = keep.data,
    family    = family,
    dosevars  = parsed$dosevars,
    Y         = Y,
    M         = M,
    X_formula = parsed$X_formula,
    X         = X,
    offset    = parsed$offset,
    entry     = parsed$entry,
    exit      = parsed$exit,
    status    = status,
    setnr     = parsed$setnr,
    deg       = deg,
    doseRRmod = parsed$doseRRmod,
    loglim    = loglim,
    optim.method = optim.method,
    control   = control
  )
  
  memoise::forget(proflik) # To avoid conflicts with existing cache when determining profile likelihood CI's
  
  result <- ameras_main(
    family     = family,
    methods    = methods,
    dosevars   = parsed$dosevars,
    data       = data,
    deg        = deg,
    doseRRmod  = parsed$doseRRmod,
    transform  = transform,
    transform.jacobian = transform.jacobian,
    setnr      = parsed$setnr,
    Y          = Y,
    M          = M,
    X          = X,
    offset     = parsed$offset,
    inpar      = inpar,
    entry      = parsed$entry,
    exit       = parsed$exit,
    status     = status,
    loglim     = loglim,
    MFMA       = MFMA,
    prophaz.numints.BMA  = prophaz.numints.BMA,
    ERRprior.BMA         = ERRprior.BMA,
    nburnin.BMA          = nburnin.BMA,
    niter.BMA            = niter.BMA,
    nchains.BMA          = nchains.BMA,
    thin.BMA             = thin.BMA,
    included.replicates.BMA = included.replicates.BMA %||% seq_along(parsed$dosevars),
    control      = control,
    optim.method = optim.method,
    ...
  )
  ret <- c(list(
    call               = match.call(),
    formula            = formula,
    num.rows           = nrow(data),
    num.replicates     = length(parsed$dosevars),
    transform          = result$transform,
    transform.jacobian = result$transform.jacobian,
    other.args         = list(...),
    model              = model_list,
    CI.computed        = FALSE
  ),
  result[setdiff(names(result), c("transform", "transform.jacobian"))]
  )
  
  ret <- new_amerasfit(ret)
  
  ret
}


parse_ameras_formula <- function(formula, data, family) {
  
  specials <- c("dose", "strata", "offset")
  X_formula        <- collect_X(formula)

  if (family == "prophaz") {
    surv         <- parse_surv_term(formula)
    formula[[2]] <- quote(.response)
    tt           <- terms(formula, specials=specials)
    dose         <- parse_dose_term(tt, data)
    
    return(list(
      Y         = NULL,
      status    = surv$status,
      entry     = surv$entry,
      exit      = surv$exit,
      dosevars  = dose$dosevars,
      doseRRmod = dose$doseRRmod,
      deg       = dose$deg,
      M         = dose$M,
      X_formula = X_formula,
      offset    = NULL,
      setnr     = NULL
    ))
    
  } else if (family == "clogit") {
    
    tt     <- terms(formula, specials=specials)
    dose   <- parse_dose_term(tt, data)
    strata <- parse_strata_term(tt)
    
    if (is.null(strata$setnr)) {
      stop("Formula for family='clogit' must contain a strata() term")
    }
    
    return(list(
      Y         = as.character(formula[[2]]),
      status    = NULL,
      entry     = NULL,
      exit      = NULL,
      dosevars  = dose$dosevars,
      doseRRmod = dose$doseRRmod,
      deg       = dose$deg,
      M         = dose$M,
      X_formula = X_formula,
      offset    = NULL,
      setnr     = strata$setnr
    ))
    
  } else if (family == "poisson") {
    
    tt   <- terms(formula, specials=specials)
    dose <- parse_dose_term(tt, data)
    off  <- parse_offset_term(tt)
    
    return(list(
      Y         = as.character(formula[[2]]),
      status    = NULL,
      entry     = NULL,
      exit      = NULL,
      dosevars  = dose$dosevars,
      doseRRmod = dose$doseRRmod,
      deg       = dose$deg,
      M         = dose$M,
      X_formula = X_formula,
      offset    = off$offset,
      setnr     = NULL
    ))
    
  } else {
    
    tt   <- terms(formula, specials=specials)
    dose <- parse_dose_term(tt, data)
    
    return(list(
      Y         = as.character(formula[[2]]),
      status    = NULL,
      entry     = NULL,
      exit      = NULL,
      dosevars  = dose$dosevars,
      doseRRmod = dose$doseRRmod,
      deg       = dose$deg,
      M         = dose$M,
      X_formula = X_formula,
      offset    = NULL,
      setnr     = NULL
    ))
  }
}

