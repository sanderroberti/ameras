
make_base_loglik_fn <- function(object, method_fit, data) {
  
  m    <- object$model
  ERC  <- method_fit$ERC
  if(is.null(ERC)) ERC <- FALSE
  Kmat <- if (ERC && m$family != "poisson") {
    cov(t(data[, m$dosevars, drop=FALSE]))
  } else {
    NULL
  }
  
  if (m$family == "gaussian") {
    function(params, D) {
      do.call(loglik.gaussian,
              c(list(params    = params,
                     D         = D,
                     X         = m$X,
                     Y         = m$Y,
                     M         = m$M,
                     data      = data,
                     deg       = m$deg,
                     ERC       = ERC,
                     Kmat      = Kmat,
                     loglim    = m$loglim,
                     transform = object$transform),
                object$other.args))
    }
  } else if (m$family == "poisson") {
    if (ERC) {
      function(params, D) {
        do.call(loglik.poisson.erc,
                c(list(params    = params,
                       D         = m$dosevars,
                       X         = m$X,
                       Y         = m$Y,
                       M         = m$M,
                       offset    = m$offset,
                       doseRRmod = m$doseRRmod,
                       data      = data,
                       deg       = m$deg,
                       loglim    = m$loglim,
                       transform = object$transform),
                  object$other.args))
      }
    } else {
      function(params, D) {
        do.call(loglik.poisson,
                c(list(params    = params,
                       D         = D,
                       X         = m$X,
                       Y         = m$Y,
                       M         = m$M,
                       offset    = m$offset,
                       doseRRmod = m$doseRRmod,
                       data      = data,
                       deg       = m$deg,
                       loglim    = m$loglim,
                       transform = object$transform),
                  object$other.args))
      }
    }
  } else if (m$family == "binomial") {
    function(params, D) {
      do.call(loglik.binomial,
              c(list(params    = params,
                     D         = D,
                     X         = m$X,
                     Y         = m$Y,
                     M         = m$M,
                     doseRRmod = m$doseRRmod,
                     data      = data,
                     deg       = m$deg,
                     ERC       = ERC,
                     Kmat      = Kmat,
                     loglim    = m$loglim,
                     transform = object$transform),
                object$other.args))
    }
  } else if (m$family == "prophaz") {
    function(params, D) {
      do.call(loglik.prophaz,
              c(list(params    = params,
                     D         = D,
                     X         = m$X,
                     status    = m$status,
                     entry     = m$entry,
                     exit      = m$exit,
                     M         = m$M,
                     doseRRmod = m$doseRRmod,
                     data      = data,
                     deg       = m$deg,
                     ERC       = ERC,
                     Kmat      = Kmat,
                     loglim    = m$loglim,
                     transform = object$transform),
                object$other.args))
    }
  } else if (m$family == "clogit") {
    designmat <- t(model.matrix(~as.factor(m$data[, m$setnr]) - 1))
    function(params, D) {
      do.call(loglik.clogit,
              c(list(params    = params,
                     D         = D,
                     X         = m$X,
                     status    = m$status,
                     designmat = designmat,
                     M         = m$M,
                     doseRRmod = m$doseRRmod,
                     data      = data,
                     deg       = m$deg,
                     ERC       = ERC,
                     Kmat      = Kmat,
                     loglim    = m$loglim,
                     transform = object$transform),
                object$other.args))
    }
  } else if (m$family == "multinomial") {
    function(params, D) {
      do.call(loglik.multinomial,
              c(list(params    = params,
                     D         = D,
                     X         = m$X,
                     Y         = m$Y,
                     M         = m$M,
                     doseRRmod = m$doseRRmod,
                     data      = data,
                     deg       = m$deg,
                     ERC       = ERC,
                     Kmat      = Kmat,
                     loglim    = m$loglim,
                     transform = object$transform),
                object$other.args))
    }
  } else {
    stop("Unknown family: ", m$family)
  }
}


make_loglik_fn <- function(object, method_name, method_fit, data) {
  
  m       <- object$model
  base_fn <- make_base_loglik_fn(object, method_fit, data)
  
  if (method_name == "MCML") {
    function(params) {
      logliks <- base_fn(params, D=m$dosevars)
      shifted <- pmax(pmin(-logliks - max(-logliks), 7e1), -7e1)
      -log(mean(exp(shifted))) - max(-logliks)
    }
  } else if (method_name %in% c("RC", "ERC")) {
    # Poisson ERC uses dosevars directly inside loglik.poisson.erc
    # so D argument is ignored in that case, but we pass it for consistency
    function(params) {
      base_fn(params, D="rcdose_ameras")
    }
  } else {
    stop("make_loglik_fn not applicable for method: ", method_name)
  }
}


resolve_data <- function(object, data=NULL) {
  
  if (!is.null(object$model$data)) {
    return(object$model$data)
  }
  
  if (is.null(data)) {
    stop(
      "Data not stored on object (keep.data=FALSE). ",
      "Please supply data argument to confint()."
    )
  }
  
  check_df(data)
  
  m            <- object$model
  needed       <- required_vars(m)
  missing_vars <- setdiff(needed, colnames(data))
  
  if (length(missing_vars)) {
    stop(
      "Supplied data is missing columns present during fitting: ",
      paste(missing_vars, collapse=", ")
    )
  }
  
  if (nrow(data) != object$num.rows) {
    warning(
      "Supplied data has ", nrow(data), " rows but original data had ",
      object$num.rows, " rows."
    )
  }
  
  if (!"rcdose_ameras" %in% colnames(data)) {
    data$rcdose_ameras <- rowMeans(data[, m$dosevars, drop=FALSE])
  }
  
  data
}







getCharVecStr <- function(x, sep=",") {
  
  ret <- paste0("'", x, "'")
  ret <- paste0(ret, collapse=sep)
  ret
  
}


# Function to check that an object is a string
isString <- function(obj) {
  
  if ((length(obj) == 1) && is.character(obj)) {
    ret <- TRUE
  } else {
    ret <- FALSE
  }
  
  ret
  
} # END: isString

check_pkgs <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Please install required packages: ", paste(missing, collapse = ", "))
  }
}

getVarNumbers <- function(vars, data) {
  
  if (!length(vars)) return(NULL)
  if (is.numeric(vars)) return(vars)
  cx  <- colnames(data)
  ret <- match(vars, cx)
  ret
}


