## ----include = FALSE----------------------------------------------------------
  knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
  )

## ----setup--------------------------------------------------------------------
  library(ameras)

## -----------------------------------------------------------------------------
transform1 <- function(params, index.t=1:length(params), lowlimit=rep(0,length(index.t)), 
                       boundcheck=FALSE, boundtol=1e-3, ...){
  
  if(length(index.t)!=length(lowlimit)) 
    stop("Length mismatch between index.t and lowlimit")
  if(any(!(index.t %in% 1:length(params)))) 
    stop("Incorrect indices for transformation specified")
  params[index.t] <- exp(params[index.t]) + lowlimit
  if(boundcheck){
    if(any(params[index.t]-lowlimit < boundtol)) 
      warning(paste0("WARNING: one or multiple parameter estimates within ", boundtol, " of 
        lower bounds. Try different bounds or starting values."))
  }
  return(params)
}

transform1.jacobian <- function(params, index.t=1:length(params), ...){ 
  if(any(!(index.t %in% 1:length(params)))) 
    stop("Incorrect indices for transformation specified")
  grad <- rep(1, length(params))
  grad[index.t] <- exp(params[index.t])
  if(length(params)>1){
    return(diag(grad))
  } else{
    return(matrix(grad))
  }
}

## -----------------------------------------------------------------------------
transform.sigmoid <- function(params, index.t=1:length(params), a=rep(0,length(index.t)), 
                              b=rep(1,length(index.t)), boundcheck=FALSE, boundtol=1e-3, ...){
  
  if(length(index.t)!=length(a) | length(index.t) != length(b)) 
    stop("Length mismatch between index.t, a, and b")
  if(any(!(index.t %in% 1:length(params)))) 
    stop("Incorrect indices for transformation specified")
  
  params[index.t] <- a + (b-a) * 1/(1+exp(-1*params[index.t]))
  if(boundcheck){
    if(any( (params[index.t]-a < boundtol) | (b-params[index.t] < boundtol))) 
      warning(paste0("WARNING: one or multiple parameter estimates within ", boundtol, 
      " of bounds. Try different bounds or starting values."))
  }
  return(params)
}

## -----------------------------------------------------------------------------
transform.sigmoid.jacobian <- function(params, index.t=1:length(params), 
                                       a=rep(0,length(index.t)), b=rep(1,length(index.t)), ...){ 
  if(length(index.t)!=length(a) | length(index.t) != length(b)) 
    stop("Length mismatch between index.t, a, and b")
  if(any(!(index.t %in% 1:length(params)))) 
    stop("Incorrect indices for transformation specified")
  grad <- rep(1, length(params))
  grad[index.t] <- (b-a)*exp(-1*params[index.t])/(1+exp(-1*params[index.t]))^2
  if(length(params)>1){
    return(diag(grad))
  } else{
    return(matrix(grad))
  }
}

## ----modelfit.sigmoid---------------------------------------------------------
data(data, package="ameras")
dosevars <- paste0("V", 1:10)
fit.ameras.sigmoid <- ameras(Y="Y.binomial", dosevars=dosevars, X=c("X1","X2"), data=data, 
                            family="binomial", deg=2, doseRRmod = "ERR", methods="RC",
                         transform=transform.sigmoid, transform.jacobian=transform.sigmoid.jacobian,
                         index.t=4:5)
summary(fit.ameras.sigmoid)

## ----modelfit.transform1------------------------------------------------------
fit.ameras.transform1 <- ameras(Y="Y.binomial", dosevars=dosevars, X=c("X1","X2"), data=data, 
                            family="binomial", deg=2, doseRRmod = "ERR", methods="RC")
summary(fit.ameras.transform1)

