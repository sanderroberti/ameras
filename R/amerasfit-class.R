new_amerasfit <- function(x = list()) {
  stopifnot(is.list(x))
  structure(
    x,
    class = "amerasfit"
  )
}



print.amerasfit <- function(x, digits = max(3, getOption("digits") - 3), ...) {

  object0 <- x[intersect(names(x), c("RC","ERC","MCML","FMA","BMA"))]
  
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
  
  
  cat("Call:\n")
  print(x$call)
  
  cat(paste0("\nNumber of individuals: ", x$num.individuals, "\n" ))
  cat(paste0("Number of dose replicates: ", x$num.replicates, "\n"))
  
  cat(paste0("\nTotal run time: ", total_runtime_seconds, " seconds\n\n"))
  
  cat("Runtime in seconds by method:\n\n")
  print(format(runtime_table, digits = digits, nsmall = 1), row.names = FALSE)
  
  
  
  
  invisible(x)
}



coef.amerasfit <- function(object, ...) {
  object <- object[intersect(names(object), c("RC","ERC","MCML","FMA","BMA"))]
  
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
  
  object0 <- object[intersect(names(object), c("RC","ERC","MCML","FMA","BMA"))]
  
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