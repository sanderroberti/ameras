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

getVarNumbers <- function(vars, data) {
  
  if (!length(vars)) return(NULL)
  if (is.numeric(vars)) return(vars)
  cx  <- colnames(data)
  ret <- match(vars, cx)
  ret
}


dldd_clogit <- function(designmat, RRs) {
  
  nr  <- nrow(designmat)
  nc  <- ncol(designmat)
  ret <- as.numeric(rep(-9999, nc))
  
  # designmat must be passed in as a vector with columns stacked
  tmp <- .C("C_dldd_clogit", as.numeric(designmat), as.integer(nr), as.integer(nc), 
            as.numeric(RRs), ret=ret, PACKAGE="ameras") 
  ret <- tmp$ret
  ret
}

dldd_prophaz <- function(entry, exit, status, RRs) {
  
  n <- length(entry)  # Number of individuals
  ret <- as.numeric(rep(-9999, n))  # Initialize return vector with placeholders
  
  # Ensure all inputs are numeric vectors (entry, exit, status, RRs)
  entry <- as.numeric(entry)
  exit <- as.numeric(exit)
  status <- as.integer(status)
  RRs <- as.numeric(RRs)
  
  # Call the C function using .C
  tmp <- .C("C_dldd_prophaz", 
            as.numeric(entry), as.numeric(exit), as.integer(status), 
            as.numeric(RRs), as.integer(n), ret=ret, 
            PACKAGE="ameras")
  
  # Return the result as a numeric vector
  ret <- tmp$ret
  return(ret)
}


compute_ERCmatrix_clogit <- function(designmat, RRs, drdd, drdd2, status) {
  
  nr  <- nrow(designmat)
  nc  <- ncol(designmat)
  ret <- as.numeric(rep(-9999, nc*nc))
  
  tmp <- .C("C_compute_ERCmatrix_clogit", as.numeric(designmat), as.integer(nr), as.integer(nc), 
            as.numeric(RRs), as.numeric(drdd), as.numeric(drdd2), as.integer(status), 
            ret=ret, PACKAGE="ameras") 
  ret <- matrix(tmp$ret, nrow=nc, ncol=nc, byrow=FALSE)
  ret
}


compute_ERCmatrix_prophaz <- function(entry, exit, status, RRs, drdd, drdd2) {
  
  n <- length(entry)  # Number of individuals
  ret <- as.numeric(rep(-9999, n*n))  # Initialize return vector with placeholders
  
  # Ensure all inputs are numeric vectors (entry, exit, status, RRs)
  entry <- as.numeric(entry)
  exit <- as.numeric(exit)
  status <- as.integer(status)
  RRs <- as.numeric(RRs)
  drdd <- as.numeric(drdd)
  drdd2 <- as.numeric(drdd2)
  
  # Call the C function using .C
  tmp <- .C("C_compute_ERCmatrix_prophaz", 
            as.numeric(entry), as.numeric(exit), as.integer(status), 
            as.numeric(RRs), as.numeric(drdd), as.numeric(drdd2),
            as.integer(n), ret=ret, 
            PACKAGE="ameras")
  
  # Return the result as a numeric vector
  ret <- matrix(tmp$ret, nrow=n, ncol=n, byrow=FALSE)
  ret
  return(ret)
}

loglik_prophaz_rcpp <- function(exit_t,
                                entry_t,
                                RR_entry,
                                RR_exit,
                                status_ord,
                                loglim = 1e-30) {
  
  n <- length(exit_t)
  K <- ncol(RR_exit)
  
  # sanity checks (minimal)
  stopifnot(
    length(entry_t)    == n,
    length(status_ord) == n,
    nrow(RR_entry)     == n,
    nrow(RR_exit)      == n,
    ncol(RR_entry)     == K
  )
  
  ret <- double(K)
  
  tmp <- .C(
    "C_loglik_prophaz_rcpp",
    exit_t     = as.double(exit_t),
    entry_t    = as.double(entry_t),
    RR_entry   = as.double(RR_entry), 
    RR_exit    = as.double(RR_exit),
    status_ord = as.integer(status_ord),
    n          = as.integer(n),
    K          = as.integer(K),
    loglim     = as.double(loglim),
    ret        = ret,
    PACKAGE    = "ameras"
  )
  
  tmp$ret
}





