transform1 <- function(params, index.t=1:length(params), lowlimit=rep(0,length(index.t)), boundcheck=FALSE, boundtol=1e-3, ...){ # Transforms from reparametrized to original scale, i.e., xi to beta
  
  if(length(index.t)!=length(lowlimit)) stop("Length mismatch between index.t and lowlimit")
  if(any(!(index.t %in% 1:length(params)))) stop("Incorrect indices for transformation specified")
  params[index.t] <- exp(params[index.t]) + lowlimit
  if(boundcheck){
    if(any(params[index.t]-lowlimit < boundtol)) warning(paste0("WARNING: one or multiple parameter estimates within ", boundtol, " of lower bounds. Try different bounds or starting values."))
  }
  return(params)
}
transform1.inv <- function(params, index.t=1:length(params), lowlimit=rep(0,length(index.t)), ...){ # Transforms from original scale to reparametrized, i.e., beta to xi
  if(length(index.t)!=length(lowlimit)) stop("Length mismatch between index.t and lowlimit")
  if(any(!(index.t %in% 1:length(params)))) stop("Incorrect indices for transformation specified")
  params[index.t] <- log(params[index.t]- lowlimit)
  return(params)
}
transform1.jacobian <- function(params, index.t=1:length(params), ...){ # Transforms from reparametrized to original scale, i.e., xi to beta
  if(any(!(index.t %in% 1:length(params)))) stop("Incorrect indices for transformation specified")
  grad <- rep(1, length(params))
  grad[index.t] <- exp(params[index.t])
  if(length(params)>1){
    return(diag(grad))
  } else{
    return(matrix(grad))
  }
}



exposureRR <- function(params, D, M, data, doseRRmod, deg){
  # params = c(beta1, beta2, (beta_m1), (beta_m2))
  params[is.infinite(params) & params>0] <- 7e1
  dosemat <- as.matrix(data[, D, drop=FALSE])
  b1 <- params[1]
  
  if(deg==2){
    b2 <- params[2]
  } else{
    b2 <- 0
  }
  
  if(!is.null(M)){
    bm1 <- params[(1+deg):(length(M)+deg)]
    Mlinpred1 <- c(as.matrix(data[,M])%*%bm1)
    if(deg==2){
      bm2 <- params[(length(M)+1+deg):(2*length(M)+deg)]
      Mlinpred2 <- c(as.matrix(data[,M])%*%bm2)
    } else{
      Mlinpred2 <- 0
    }
  } else{
    Mlinpred1 <- Mlinpred2 <- 0
  }
  
  if(doseRRmod=="ERR"){
    return(1+b1*dosemat+b2*dosemat^2+Mlinpred1*dosemat+Mlinpred2*dosemat^2)
  } else if(doseRRmod=="EXP"){
    val <- b1*dosemat + b2*dosemat^2 + Mlinpred1*dosemat + Mlinpred2*dosemat^2
    val <- pmin(val, 7e1)       # cap large finite values
    val[!is.finite(val)] <- 7e1 # replace NaN, Inf, -Inf
    return(exp(val))
  } else if(doseRRmod=="LINEXP"){
    val <- 1 + (b1 + Mlinpred1) * dosemat * exp(pmin((b2 + Mlinpred2) * dosemat, 7e1))
    val <- pmin(val, exp(7e1))
    val[!is.finite(val)] <- exp(7e1)  # replace Inf/NaN with capped value
    return(val)
    #return(pmin(1+(b1+Mlinpred1)*dosemat*exp(pmin((b2+Mlinpred2)*dosemat, 8e1)), exp(8e1)))
  }
}


dRRdD <- function(params, D, M, data, doseRRmod, deg){
  # params = c(beta1, beta2, (beta_m1), (beta_m2))
  params[is.infinite(params) & params>0] <- 7e1
  b1 <- params[1]
  
  if(deg==2){
    b2 <- params[2]
  } else{
    b2 <- 0
  }
  
  if(!is.null(M)){
    bm1 <- params[(1+deg):(length(M)+deg)]
    Mlinpred1 <- c(as.matrix(data[,M])%*%bm1)
    if(deg==2){
      bm2 <- params[(length(M)+1+deg):(2*length(M)+deg)]
      Mlinpred2 <- c(as.matrix(data[,M])%*%bm2)
    } else{
      Mlinpred2 <- 0
    }
  } else{
    Mlinpred1 <- Mlinpred2 <- 0
  }
  
  # First derivative
  if(doseRRmod=="ERR"){
    first <- b1+2*b2*data[,D]+Mlinpred1+2*Mlinpred2*data[,D]
  } else if(doseRRmod=="EXP"){
    first <- (b1+2*b2*data[,D]+Mlinpred1+2*Mlinpred2*data[,D])*exp(pmin(b1*data[,D]+b2*data[,D]^2+Mlinpred1*data[,D]+Mlinpred2*data[,D]^2, 8e1))
  } else if(doseRRmod=="LINEXP"){
    first <- (b1+Mlinpred1)*exp((b2+Mlinpred2)*data[,D])+(b2+Mlinpred2)*(b1+Mlinpred1)*data[,D]*exp((b2+Mlinpred2)*data[,D])
  }
  
  # Second derivative
  if(doseRRmod=="ERR"){
    second <- 2*(b2+Mlinpred2)
    
    if(length(second)==1) second <- rep(second, length(first))
  } else if(doseRRmod=="EXP"){
    second <- 2*(b2+Mlinpred2)*exp(pmin(b1*data[,D]+b2*data[,D]^2+Mlinpred1*data[,D]+Mlinpred2*data[,D]^2,8e1))+(b1+2*b2*data[,D]+Mlinpred1+2*Mlinpred2*data[,D])^2*exp(pmin(b1*data[,D]+b2*data[,D]^2+Mlinpred1*data[,D]+Mlinpred2*data[,D]^2, 8e1))
    
    if(length(second)==1) second <- rep(second, length(first))
  } else if(doseRRmod=="LINEXP"){
    second <- (b1+Mlinpred1)*(b2+Mlinpred2)*exp((b2+Mlinpred2)*data[,D])+(b2+Mlinpred2)*first
    
    if(length(second)==1) second <- rep(second, length(first))
  }
  
  return(list(first=first, second=second))
}



loglik.binomial <- function(params, Y, D, M=NULL, X=NULL, data, doseRRmod, ERC=FALSE, Kmat=NULL, deg=1, loglim=1e-30, transform=NULL, ...){
  # params = c(a0, a1, ..., ap, b1, b2, (bm1), (bm2))
  
  if(length(params) != deg*length(M)+deg+length(X)+1) stop("Parameter vector length mismatch")
  if(ERC & length(D)>1) stop("ERC only works with one supplied dose (i.e., the mean across replicates)")
  # doseRRmod = ERR or EXP
  if(!is.null(transform)){
    if(is.function(transform)){
      params <- transform(params=params, ...)
    } else{
      stop("transform should be a function")
    }
  }
  
  a0 <- params[1]
  
  
  
  if(!is.null(X)){
    a <- params[2:(length(X)+1)]
    Xlinpred <- c(as.matrix(data[,X])%*%a)
  } else{
    Xlinpred <- 0
  }
  
  b1 <- params[length(X)+2]
  
  
  if(deg==2){
    b2 <- params[length(X)+3]
  } else if(deg==1){
    b2 <- NULL
  }
  
  if(!is.null(M)){
    bm1 <- params[(2+deg+length(X)):(length(M)+deg+length(X)+1)]
    if(deg==2){
      bm2 <- params[(length(M)+2+deg+length(X)):(2*length(M)+deg+length(X)+1)]
    } else{
      bm2 <- NULL
    }
  } else{
    bm1 <- bm2 <- NULL
  }
  
  if(deg==2){
    betavec <- c(b1, b2, bm1, bm2)
  } else{
    betavec <- c(b1, bm1)
  }
  
  A <- exp(pmin(a0+Xlinpred,7e1))*exposureRR(params=betavec, D=D, M=M, data=data, doseRRmod=doseRRmod, deg=deg)
  
  if(any(A<0)) stop("RR <0, please check supplied transformation")
  p <- A/(1+A)
  p <- pmin(pmax(p, 1e-10), 1-1e-10)
  
  if(length(D)>1){
    ls <- colSums(log(p^data[,Y]*(1-p)^(1-data[,Y])))
    ls <- as.numeric(ls)
  } else{
    ls <- sum(log(p^data[,Y]*(1-p)^(1-data[,Y])))
  }
  if(ERC){
    derivs <- dRRdD(params=betavec, D=D, M=M, data=data, doseRRmod=doseRRmod, deg=deg)
    dpdk <- exp(pmin(a0+Xlinpred,7e1))*derivs$first/(1+A)^2
    dpdk2 <- exp(pmin(a0+Xlinpred,7e1))*derivs$second/(1+A)^2-2*exp(pmin(a0+Xlinpred,7e1))^2*derivs$first^2/(1+A)^3
    
    mymat <- tcrossprod((-1)^(1-data[,Y]))*tcrossprod(dpdk)/tcrossprod(p^data[,Y]*(1-p)^(1-data[,Y]))
    
    diag(mymat) <- (-1)^(1-data[,Y])*dpdk2/(p^data[,Y]*(1-p)^(1-data[,Y]))
    
    
    return(-1*(ls+log(max(1+.5*sum(mymat*Kmat), loglim))))
    
  } else{
    return(-1*ls)
  }
  
}





loglik.poisson <- function(params, Y, D, X=NULL, offset=NULL,M=NULL, doseRRmod, data, deg=1, loglim=1e-30, transform=NULL, ...){ # C++ version
  # params = c(a0, a1, ..., ap, b1, b2, (bm1), (bm2))
  
  if(length(params) != deg*length(M)+deg+length(X)+1) stop("Parameter vector length mismatch")
  
  if(!is.null(transform)){
    if(is.function(transform)){
      params <- transform(params=params, ...)
    } else{
      stop("transform should be a function")
    }
  }
  
  if(is.null(offset)){
    offset <- 1
  } else{
    offset <- data[,offset]
  }
  
  a0 <- params[1]
  
  if(!is.null(X)){
    a <- params[2:(length(X)+1)]
    Xlinpred <- c(as.matrix(data[,X])%*%a)
  } else{
    Xlinpred <- 0
  }
  
  b1 <- params[length(X)+2]
  
  
  if(deg==2){
    b2 <- params[length(X)+3]
  } else if(deg==1){
    b2 <- 0
  }
  
  if(!is.null(M)){
    bm1 <- params[(2+deg+length(X)):(length(M)+deg+length(X)+1)]
    if(deg==2){
      bm2 <- params[(length(M)+2+deg+length(X)):(2*length(M)+deg+length(X)+1)]
    } else{
      bm2 <- NULL
    }
  } else{
    bm1 <- bm2 <- NULL
  }
  
  if(deg==2){
    betavec <- c(b1, b2, bm1, bm2)
  } else{
    betavec <- c(b1, bm1)
  }
  
  mus <- pmax(exp(pmin(a0+Xlinpred,7e1))*pmax(exposureRR(params=betavec, D=D, M=M, data=data, doseRRmod=doseRRmod, deg=deg), loglim)*offset, loglim)
  if(any(mus<0)) stop("RR <0, please check supplied transformation")
  
  if(length(D)>1){
    #ls <- colSums(data[,Y]*log(mus)-mus-lfactorial(data[,Y]))
    ls <- colSums(data[,Y]*log(mus)-mus-ifelse(data[,Y]>0,data[,Y]*log(data[,Y])-data[,Y],0)) # Stirling's approximation like Mark
    ls <- as.numeric(ls)
  } else{
    #ls <- sum(data[,Y]*log(mus)-mus-lfactorial(data[,Y]))
    ls <- sum(data[,Y]*log(mus)-mus-ifelse(data[,Y]>0,data[,Y]*log(data[,Y])-data[,Y],0))  # Stirling's approximation like Mark
  }
  
  return(-1*ls)
  
}


loglik.poisson.erc <- function(params, Y, D, X=NULL, offset=NULL,M=NULL, doseRRmod, data,deg=1, loglim=1e-30, transform=NULL, ...){ # C++ version
  # params = c(a0, a1, ..., ap, b1, b2, (bm1), (bm2))
  
  if(length(params) != deg*length(M)+deg+length(X)+1) stop("Parameter vector length mismatch")
  #if(ERC & is.null(Kmat)) stop("Kmat necessary for ERC")
  #if(ERC & length(D)>1) stop("ERC only works with one supplied dose (i.e., the mean across replicates)")
  
  dosemat <- as.matrix(data[,D])
  
  if(!is.null(transform)){
    if(is.function(transform)){
      params <- transform(params=params, ...)
    } else{
      stop("transform should be a function")
    }
  }
  
  if(is.null(offset)){
    offset <- 1
  } else{
    offset <- data[,offset]
  }
  
  a0 <- params[1]
  
  if(!is.null(X)){
    a <- params[2:(length(X)+1)]
    Xlinpred <- c(as.matrix(data[,X])%*%a)
  } else{
    Xlinpred <- 0
  }
  
  b1 <- params[length(X)+2]
  
  
  if(deg==2){
    b2 <- params[length(X)+3]
  } else if(deg==1){
    b2 <- 0
  }
  
  if(!is.null(M)){
    bm1 <- params[(2+deg+length(X)):(length(M)+deg+length(X)+1)]
    if(deg==2){
      bm2 <- params[(length(M)+2+deg+length(X)):(2*length(M)+deg+length(X)+1)]
    } else{
      bm2 <- NULL
    }
  } else{
    bm1 <- bm2 <- NULL
  }
  
  if(deg==2){
    betavec <- c(b1, b2, bm1, bm2)
  } else{
    betavec <- c(b1, bm1)
  }
  
  #data$rcdose_ameras <- rowMeans(dosemat)
  
  mus <- pmax(exp(pmin(a0+Xlinpred,7e1))*pmax(exposureRR(params=betavec, D="rcdose_ameras", M=M, data=data, doseRRmod=doseRRmod, deg=deg), loglim)*offset, loglim)
  if(any(mus<0)) stop("RR <0, please check supplied transformation")
  
  
  #ls <- sum(data[,Y]*log(mus)-mus-lfactorial(data[,Y]))
  ls <- sum(data[,Y]*log(mus)-mus-ifelse(data[,Y]>0,data[,Y]*log(data[,Y])-data[,Y],0))  # Stirling's approximation like Mark
  
  
  derivs <- dRRdD(params=betavec, D="rcdose_ameras", M=M, data=data, doseRRmod=doseRRmod, deg=deg)
  dmdd <- exp(pmin(a0+Xlinpred,7e1))*derivs$first*offset
  dmdd2 <- exp(pmin(a0+Xlinpred,7e1))*derivs$second*offset
  
  v <- (data[,Y]/mus - 1) * dmdd
  
  
  
  
  
  
  
  # center dose for covariance
  Xc <- dosemat - rowMeans(dosemat)
  Xt_v <- crossprod(Xc, v)
  
  term1 <- sum(Xt_v^2) / (ncol(dosemat) - 1)
  
  # diagonal correction
  corrterm <- (data[,Y]/mus - 1) * dmdd2 - data[,Y]/mus^2 * dmdd^2
  
  # diag(Kmat) = row variances of dosemat
  row_var <- rowSums(Xc^2) / (ncol(dosemat) - 1)
  
  term2 <- sum(row_var * corrterm)
  
  val <- term1 + term2
  
  
  return(-1*(ls+log(max(1+.5*val, loglim))))
  
  
}


loglik.gaussian <- function(params, D, Y,X=NULL,M=NULL, data, ERC=FALSE, Kmat=NULL,deg=1, loglim=1e-30, transform=NULL, ...){
  # params = c(a0, a1, ..., ap, b1, b2, (bm1), (bm2), sigma)
  
  if(length(params) != deg*length(M)+deg+length(X)+2) stop("Parameter vector length mismatch")
  if(ERC & is.null(Kmat)) stop("Kmat necessary for ERC")
  if(ERC & length(D)>1) stop("ERC only works with one supplied dose (i.e., the mean across replicates)")
  
  if(!is.null(transform)){
    if(is.function(transform)){
      params <- transform(params=params, ...)
    } else{
      stop("transform should be a function")
    }
  }
  
  a0 <- params[1]
  
  if(!is.null(X)){
    a <- params[2:(length(X)+1)]
    Xlinpred <- c(as.matrix(data[,X])%*%a)
  } else{
    Xlinpred <- 0
  }
  
  b1 <- params[length(X)+2]
  
  
  if(deg==2){
    b2 <- params[length(X)+3]
  } else if(deg==1){
    b2 <- 0
  }
  
  if(!is.null(M)){
    bm1 <- params[(2+deg+length(X)):(length(M)+deg+length(X)+1)]
    if(deg==2){
      bm2 <- params[(length(M)+2+deg+length(X)):(2*length(M)+deg+length(X)+1)]
    } else{
      bm2 <- NULL
    }
  } else{
    bm1 <- bm2 <- NULL
  }
  
  if(deg==2){
    betavec <- c(b1, b2, bm1, bm2)
  } else{
    betavec <- c(b1, bm1)
  }
  
  sigma <- params[deg*length(M)+deg+length(X)+2]
  
  #L <- (1/(2*pi*sigma^2))^(length(Y)/2)*exp(-sum((Y-alpha-beta1*D-beta2*D^2)^2)/(2*sigma^2))
  
  # Here I use mu = a0 + X^Ta + RR - 1 with RR using the linear ERR 
  mus <- a0+Xlinpred+exposureRR(params=betavec, D=D, M=M, data=data, doseRRmod="ERR", deg=deg)-1
  
  ls <- -log(2*pi*sigma^2)*nrow(data)/2-colSums((data[,Y]-as.matrix(mus, ncol=length(D)))^2)/(2*sigma^2)
  
  ls <- as.numeric(ls)
  
  if(ERC){
    # See above, dmu/dD = dRR/dD
    
    derivs <- dRRdD(params=betavec, D=D, M=M, data=data, doseRRmod="ERR", deg=deg)
    dmdd <- derivs$first
    dmdd2 <- derivs$second
    
    mymat <- (tcrossprod(dmdd)/sigma^2)*(tcrossprod(data[,Y]-mus)/sigma^2-diag(1, nrow=nrow(data), ncol=nrow(data)))
    diag(mymat) <- diag(mymat) + dmdd2/sigma^2*(data[,Y]-mus)
    
    #res <- compute_weighted_sum_kahan(dmdd, dmdd2, data[,Y]-mus, Kmat, sigma)
    
    return(-1*(ls+log(max(1+.5*sum(mymat*Kmat),loglim))))
    #return(-1*(l+log(max(1+.5*res,loglim))))
  } else{
    return(-1*ls)
  }
}



loglik.clogit <- function(params, D, status,X=NULL, M=NULL, doseRRmod, designmat, data, deg=1, ERC=FALSE, Kmat=NULL, loglim=1e-30, transform=NULL, ...){
  # params = c(a1, ..., ap, b1, b2, (bm1), (bm2))
  #print(params)
  if(length(params) != deg*length(M)+deg+length(X)) stop("Parameter vector length mismatch")
  if(ERC & is.null(Kmat)) stop("Kmat necessary for ERC")
  if(ERC & length(D)>1) stop("ERC only works with one supplied dose (i.e., the mean across replicates)")
  
  if(!is.null(transform)){
    if(is.function(transform)){
      params <- transform(params=params, ...)
    } else{
      stop("transform should be a function")
    }
  }
  
  
  if(!is.null(X)){
    a <- params[1:length(X)]
    Xlinpred <- c(as.matrix(data[,X])%*%a)
  } else{
    Xlinpred <- 0
  }
  
  b1 <- params[length(X)+1]
  
  
  if(deg==2){
    b2 <- params[length(X)+2]
  } else if(deg==1){
    b2 <- 0
  }
  
  if(!is.null(M)){
    bm1 <- params[(1+deg+length(X)):(length(M)+deg+length(X))]
    if(deg==2){
      bm2 <- params[(length(M)+1+deg+length(X)):(2*length(M)+deg+length(X))]
    } else{
      bm2 <- NULL
    }
  } else{
    bm1 <- bm2 <- NULL
  }
  
  if(deg==2){
    betavec <- c(b1, b2, bm1, bm2)
  } else{
    betavec <- c(b1, bm1)
  }
  
  RRs <- exp(pmin(Xlinpred, 7e1))*exposureRR(params=betavec, D=D, M=M, data=data, doseRRmod=doseRRmod, deg=deg)
  
  if(any(RRs<0)) stop("RR <0, please check supplied transformation")
  ls <- colSums(log(RRs[data[,status]==1])-log(designmat %*% as.matrix(RRs, ncol=length(D))))
  ls <- as.numeric(ls)
  
  if(ERC){
    derivs <- dRRdD(params=betavec, D=D, M=M, data=data, doseRRmod=doseRRmod, deg=deg)
    drdd <- exp(pmin(Xlinpred, 7e1))*derivs$first
    drdd2 <- exp(pmin(Xlinpred, 7e1))*derivs$second
    
    tmp <- c(dldd_clogit(designmat, RRs))
    dldd <- drdd/RRs*(data[,status]==1) - tmp * drdd
    
    mymat <- compute_ERCmatrix_clogit(designmat, RRs, drdd, drdd2, as.integer(data[,status]))
    
    return(-1*(ls+log(max(1+.5*sum((tcrossprod(dldd)+mymat)*Kmat), loglim))))
  } else{
    return(-1*ls)
  }
  
}

loglik.prophaz <- function(params, D, status,X=NULL, M=NULL, doseRRmod, data, deg=1,entry=NULL, exit, ERC=FALSE, Kmat=NULL, loglim=1e-30, transform=NULL, ...){
  # params = c(a1, ..., ap, b1, b2, (bm1), (bm2))
  #print(params)
  if(length(params) != deg*length(M)+deg+length(X)) stop("Parameter vector length mismatch")
  if(ERC & is.null(Kmat)) stop("Kmat necessary for ERC")
  if(ERC & length(D)>1) stop("ERC only works with one supplied dose (i.e., the mean across replicates)")
  
  if(!is.null(transform)){
    if(is.function(transform)){
      params <- transform(params=params, ...)
    } else{
      stop("transform should be a function")
    }
  }
  
  
  if(!is.null(X)){
    a <- params[1:length(X)]
    Xlinpred <- c(as.matrix(data[,X])%*%a)
  } else{
    Xlinpred <- 0
  }
  
  b1 <- params[length(X)+1]
  
  
  if(deg==2){
    b2 <- params[length(X)+2]
  } else if(deg==1){
    b2 <- 0
  }
  
  if(!is.null(M)){
    bm1 <- params[(1+deg+length(X)):(length(M)+deg+length(X))]
    if(deg==2){
      bm2 <- params[(length(M)+1+deg+length(X)):(2*length(M)+deg+length(X))]
    } else{
      bm2 <- NULL
    }
  } else{
    bm1 <- bm2 <- NULL
  }
  
  if(deg==2){
    betavec <- c(b1, b2, bm1, bm2)
  } else{
    betavec <- c(b1, bm1)
  }
  
  #RRs <- exp(pmin(Xlinpred, 7e1))*exposureRR(params=betavec, D=D, M=M, data=data, doseRRmod=doseRRmod, deg=deg)
  e1 <- exposureRR(params=betavec, D=D, M=M, data=data, doseRRmod=doseRRmod, deg=deg)
  e1 <- as.matrix(e1, ncol=length(D))
  if(any(e1<0)) stop("RR <0, please check supplied transformation")
  
  logRRs <- Xlinpred + log(e1)
  #logRRs <- logRRs - max(logRRs) # Is not compatible with ERC calculation!!
  
  RRs <- exp(pmin(logRRs,7e1))
  
  ord_exit <- order(data[[exit]])
  exit_t <- data[[exit]][ord_exit]
  status_ord <- data[[status]][ord_exit]
  RR_exit <- RRs[ord_exit,, drop=FALSE]
  
  if(!is.null(entry)){
    ord_entry <- order(data[[entry]])
    entry_t <- data[[entry]][ord_entry]
    RR_entry <- RRs[ord_entry,, drop=FALSE]
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
  
  if(ERC){
    derivs <- dRRdD(params=betavec, D=D, M=M, data=data, doseRRmod=doseRRmod, deg=deg)
    drdd <- exp(pmin(Xlinpred, 7e1))*derivs$first
    drdd2 <- exp(pmin(Xlinpred, 7e1))*derivs$second
    
    drdd <- drdd[ord_exit]
    drdd2 <- drdd2[ord_exit]
    
    Kmat <- Kmat[ord_exit, ord_exit]
    
    if(!is.null(entry)){
      entry_t2 <- data[[entry]][ord_exit]
    } else{
      entry_t2 <- entry_t
    }
    tmp <- c(dldd_prophaz(entry_t2,exit_t, status_ord, RR_exit))
    dldd <- drdd/RR_exit*(status_ord==1) - tmp * drdd
    
    
    
    mymat <- compute_ERCmatrix_prophaz(entry_t2, exit_t, status_ord, RR_exit, drdd, drdd2)
    
    return(-1*(ls+log(max(1+.5*sum((tcrossprod(dldd)+mymat)*Kmat), loglim))))
  } else{
    return(-1*ls)
  }
  
}


loglik.multinomial <- function(params, Y, D, M=NULL, X=NULL, data, doseRRmod, ERC=FALSE, Kmat=NULL, deg=1, loglim=1e-30, transform=NULL, ...){
  
  Z <- length(unique(data[,Y])) # number of outcome categories
  # params = c(a0^(1), (a)^(1), b1^(1), b2^(1), (bm1)^(1), (bm2)^(1), ..., a0^(Z-1), (a)^(Z-1), b1^(Z-1), b2^(Z-1), (bm1)^(Z-1), (bm2)^(Z-1) )
  if(length(params) != (Z-1)*(deg*length(M)+deg+length(X)+1)) stop("Parameter vector length mismatch")
  if(ERC & length(D)>1) stop("ERC only works with one supplied dose (i.e., the mean across replicates)")
  
  
  
  if(!is.null(transform)){
    if(is.function(transform)){
      params <- transform(params=params, ...)
    } else{
      stop("transform should be a function")
    }
  }
  
  params <- matrix(params, ncol = Z-1)
  params <- cbind(params, rep(0, nrow(params)))
  
  a0 <- params[1,]
  
  
  
  if(!is.null(X)){
    a <- params[2:(length(X)+1),]
    Xlinpred <- sweep(as.matrix(data[,X])%*%a, 2, a0, "+")
  } else{
    Xlinpred <- matrix(a0, nrow=nrow(data), ncol=Z, byrow = TRUE)
  }
  
  betamat <- params[(length(X)+2):(deg*length(M)+deg+length(X)+1), , drop=FALSE]
  
  if(length(D)>1){
    myarray <- array(NA_real_, dim = c(nrow(data), length(D), Z))
    
    for (ii in 1:(Z-1)) {
      myarray[,,ii] <- exposureRR(
        params = betamat[, ii],
        D = D,
        M = M,
        data = data,
        doseRRmod = doseRRmod,
        deg = deg
      )
    }
    myarray[,,Z] <- 1
    
    if(any(myarray<0)) stop("RR <0, please check supplied transformation")
    
    Xlinpred_array <- array(exp(pmin(Xlinpred, 7e1)), dim = c(nrow(data), 1, Z))
    Xlinpred_array <- Xlinpred_array[, rep(1, length(D)),, drop=FALSE]
    
    
    myarray <- myarray * Xlinpred_array # RR(d) * exp(X^t a)
    
    
    Asums <- colSums(aperm(myarray, c(3,1,2)))
    Asums_array <- array(Asums, dim=c(nrow(data), length(D), 1))
    Asums_array <- Asums_array[,, rep(1, Z), drop=FALSE]
    
    prob_array <- myarray/Asums_array # array with probabilities, indexed by individual x dose replicate x outcome category
    
    
    Ymat <- as.matrix(model.matrix(~data[,Y]-1))
    Y_array <- array(Ymat, dim=c(nrow(data), 1, Z))
    Y_array <- Y_array[,rep(1, length(D)), , drop=FALSE]
    
    myarray2 <- log(prob_array)*Y_array
    ls <- colSums(aperm(myarray2, c(1,3,2)), dims=2)
  } else{
    
    
    RRmat <- matrix(1, nrow=nrow(data), ncol=Z)
    
    for(ii in 1:(Z-1)){
      RRmat[,ii] <- exposureRR(
        params = betamat[,ii],
        D = D,
        M = M,
        data = data,
        doseRRmod = doseRRmod,
        deg = deg
      )
    }
    
    
    
    
    if(any(RRmat<0)) stop("RR <0, please check supplied transformation")
    
    RRmat <- RRmat * exp(pmin(Xlinpred, 7e1))
    probmat <- RRmat/rowSums(RRmat)[,drop=FALSE]
    
    
    
    Ymat <- diag(Z)[as.integer(data[,Y]), ]#as.matrix(model.matrix(~data[,Y]-1))
    ls <- sum(log(probmat)*Ymat)
    
    
  }
  
  
  if(ERC){ # Only one dose supplied for ERC
    
    A <- RRmat
    
    drdd <- apply(betamat, 2, function(betavec) dRRdD(params=betavec, D=D, M=M, data=data, doseRRmod=doseRRmod, deg=deg)$first)
    drdd2 <- apply(betamat, 2, function(betavec) dRRdD(params=betavec, D=D, M=M, data=data, doseRRmod=doseRRmod, deg=deg)$second)
    
    exp_X <- exp(pmin(Xlinpred, 7e1))
    rowSums_A <- rowSums(A)
    
    probvec <- rowSums(Ymat*probmat)#diag(tcrossprod(Ymat, probmat))
    
    drdd_expX <- drdd*exp_X
    Ymat_drdd_expX <- rowSums(Ymat*drdd_expX)#diag(tcrossprod(Ymat, drdd_expX))
    cross_drdd_exp_X <- rowSums(drdd*exp_X)#diag(tcrossprod(drdd, exp_X))
    
    dpdd <- (Ymat_drdd_expX - probvec * cross_drdd_exp_X ) / rowSums_A
    dpdd2 <- (rowSums(Ymat*exp_X*drdd2)-
                dpdd * cross_drdd_exp_X-
                probvec * rowSums(drdd2*exp_X))/rowSums_A+
      (cross_drdd_exp_X^2 * probvec - cross_drdd_exp_X * Ymat_drdd_expX) / rowSums_A^2
    
    mymat <- tcrossprod(dpdd/probvec)
    diag(mymat) <- dpdd2/probvec
    
    
    return(-1*(ls+log(max(1+.5*sum(mymat*Kmat), loglim))))
    
  } else{
    return(-1*ls)
  }
  
}



proflik <- memoise(function(parvalue, index, fun, inpar, optim.method=optim.method, ...){
  
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
    innerfit <- optim(par=inpar[-index],fn=innerfun, method=optim.method, ...) 
    if(optim.method=="Nelder-Mead"){
      count0 <- innerfit$counts
      innerfit <- optim(par=innerfit$par,fn=innerfun, method="BFGS", ...)
      innerfit$counts <- replace(count0, is.na(count0), 0) + replace(innerfit$counts, is.na(innerfit$counts), 0)
    }
    return(innerfit$value)
  }
  
})


profCIfun <- function(par,optval, index, fun, inpar, optim.method, ...){
  sapply(par, function(mypar){
    1-pchisq(2*(proflik(parvalue=mypar, index=index, fun=fun, inpar=inpar, optim.method=optim.method, ...)-optval),df=1)-.05
  })
}





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

ameras.fma <- function(family, dosevars, data, deg, transform=NULL,transform.jacobian=NULL, Y=NULL, M=NULL, X=NULL, offset=NULL, inpar=NULL, entry=NULL, exit=NULL, status=NULL, setnr=setnr, CI=NULL, unweighted=NULL, doseRRmod=NULL, MFMA=100000, optim.method="Nelder-Mead", ...){
  
  # Remove build warnings
  #HPDinterval <- as.mcmc <- NULL
  
  if(length(CI)==0) stop("No CI method specified")
  CI <- CI[CI %in% c("percentile","hpd")]
  if(length(CI)>1) stop("Provide one type of CI for FMA: one of percentile and hpd")
  if(length(CI)==0) stop("Incorrect CI method specified, should be one of percentile and hpd")
  
  
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
    
    names(coefs) <- names(sd) <- parnames
    
    if(CI=="percentile"){
      CIlower=apply(FMAsamples, 2, function(x) quantile(x, .025))
      CIupper=apply(FMAsamples, 2, function(x) quantile(x, .975))
      CIresult <- data.frame(lower=CIlower, upper=CIupper)
    } else if(CI=="hpd"){
      CIresult <- as.data.frame(HPDinterval(as.mcmc(FMAsamples)))
    }
    
    rownames(CIresult) <- parnames
    included.samples <- nrow(FMAsamples)
  } else {
    
    coefs <- sd <- CIlower <- CIupper <- NA * FMAfits[[1]]$coef
    CIresult <- data.frame(lower=CIlower, upper=CIupper)
    
    names(coefs) <- names(sd) <- parnames
    rownames(CIresult) <- parnames
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
    CI=CIresult,
    included.replicates=included.replicates,
    included.samples=included.samples,
    #includedfits=FMAfits,
    #allfits=allfits,
    runtime=runtime
  )
  
  return(out)
  
}

boundedSliceSampler <- nimbleFunction(
  contains = sampler_BASE,
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

ameras.bma <- function(family, dosevars, data, deg, Y=NULL, M=NULL, X=NULL, offset=NULL, entry=NULL, exit=NULL, status=NULL, setnr=setnr, CI=NULL, transform=NULL, inpar=NULL, doseRRmod=NULL, ERRprior="doubleexponential",prophaz_numints=10, nburnin=1000, niter=5000, included.replicates=1:length(dosevars), nchains=2, thin=10, optim.method="Nelder-Mead", ...){
  
  # Remove build warnings, local functions may need to be pulled out from this function
  # HPDinterval <- K <- Mlen <- Mmat <- N <- Xlen <- Xmat <- a <- as.mcmc <- 
  #   b <- bm <- bm1 <- bm2 <- buildMCMC <- nsets <- setmat <-
  #   col.ind <- compileNimble <- configureMCMC <- delta <- 
  #   dosemat <- equals <- h0 <- inprod <- nimbleModel <- runMCMC <- NULL
  
  
  
  ndoses <- length(dosevars)
  
  if(length(CI)==0) stop("No CI method specified")
  CI <- CI[CI %in% c("percentile","hpd")]
  if(length(CI)>1) stop("Provide one type of CI for BMA: one of percentile and hpd")
  if(length(CI)==0) stop("Incorrect CI method specified, should be one of percentile and hpd")
  
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
  
  
  if(CI=="percentile"){
    CIlower=apply(nimblesamples.stacked, 2, function(x) quantile(x, .025))
    CIupper=apply(nimblesamples.stacked, 2, function(x) quantile(x, .975))
    CIresult <- data.frame(lower=CIlower, upper=CIupper)
  } else if(CI=="hpd"){
    CIresult <- as.data.frame(HPDinterval(as.mcmc(nimblesamples.stacked)))
  }
  
  rownames(CIresult) <- parnames
  
  t1 <- proc.time()
  timedif <- t1-t0
  runtime <- paste(round(as.numeric(as.difftime(timedif["elapsed"], units="secs")),1), "seconds")
  
  prc_excluded <- round(100*(1-length(included.replicates)/ndoses), 1)
  if(length(included.replicates)/ndoses < .8) warning(paste0("WARNING: ",prc_excluded, "% of replicates excluded from model averaging. Try different bounds or starting values."))
  
  out <- list(coefficients=coef,
              sd=sd,
              CI=CIresult,
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


#-------------

coef.amerasfit <- function(object, ...) {
  object <- object[setdiff(names(object), "call")]
  
  bma <- ("BMA" %in% names(object))
  
  res <- as.data.frame(do.call("cbind",lapply(1:length(object), function(i){ 
    
    y <- object[[i]]
    
    if(bma){
      tmp <- c(y$coefficients,object$BMA$coefficients[-(1:length(y$coefficients))]*NA)
    } else{
      tmp <- y$coefficients
    }
    
    tmp
    
    
  })))
  names(res) <- names(object)
  res
}


summary.amerasfit <- function(object, ...) {
  
  object0 <- object[setdiff(names(object), "call")]
  
  bma <- ("BMA" %in% names(object0))
  
  summary_table <- do.call("rbind",lapply(1:length(object0), function(i){
    
    y <- object0[[i]]
    method <- names(object0)[i]
    
    coef <- y$coefficients
    
    se <- y$sd
    
    
    CI <- y$CI
    
    
    
    CI.lowerbound  <- CI.upperbound <- coef*NA
    CI.lowerbound[match(rownames(CI), names(coef))] <- CI$lower
    CI.upperbound[match(rownames(CI), names(coef))] <- CI$upper
    
    res <- data.frame(Method=method,
                      Term=names(coef),
                      Estimate = coef,
                      SE = se,
                      CI.lowerbound = CI.lowerbound,
                      CI.upperbound = CI.upperbound)
    if(bma){
      
      
      if(method=="BMA"){
        res <- cbind(res, y$Rhat)
      } else{
        res <- cbind(res, data.frame(Rhat=NA, n.eff=NA))
      }
    }
    rownames(res) <- NULL
    res
    
  })
  )
  
  runtime_table <- do.call("rbind",lapply(1:length(object0), function(i){
    
    y <- object0[[i]]
    method <- names(object0)[i]
    
    runtime <- as.numeric(strsplit(y$runtime, " seconds")[[1]])
    
    
    
    res <- data.frame(Method=method,
                      Runtime=runtime)
    rownames(res) <- NULL
    res
    
  })
  )
  
  total_runtime_seconds <- sum(sapply(object0, function(x) as.numeric(strsplit(x$runtime, " seconds")[[1]])))
  
  
  
  ans <- list(
    call = object$call,
    summary_table = summary_table,
    runtime_table = runtime_table,
    total_runtime_seconds = total_runtime_seconds
  )
  
  class(ans) <- "summary.amerasfit"
  return(ans)
}


print.summary.amerasfit <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("Call:\n")
  print(x$call)
  cat(paste0("\nTotal run time: ",x$total_runtime_seconds, " seconds\n\n"))
  
  cat("Runtime in seconds by method:\n\n")
  print(format(x$runtime_table, digits = digits, nsmall = 1), row.names = FALSE)
  
  cat("\nSummary of coefficients by method:\n\n")
  print(format(x$summary_table, digits = digits, nsmall = 2), row.names = FALSE)
  
  invisible(x)
}

traceplot <- function(object, ...) {
  UseMethod("traceplot")
}

traceplot.amerasfit <- function(object, iter=5000,  Rhat=TRUE, n.eff=TRUE, pdf=FALSE, ...){
  
  if("BMA" %in% names(object)){
    MCMCtrace(object$BMA$samples, iter=iter, Rhat=Rhat, n.eff=n.eff, pdf=pdf, ...)
    
  } else {
    stop("ERROR: traceplot() requires BMA output in the provided 'amerasfit' object.")
  }
  
}

