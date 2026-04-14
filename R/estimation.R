

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