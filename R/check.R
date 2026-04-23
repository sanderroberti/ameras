
check_df <- function(x, nm="data") {
  
  if (!is.data.frame(x)) stop(paste0("ERROR: ", nm, " must be a data frame"))
  if (!nrow(x)) stop(paste0("ERROR: ", nm, " has no rows"))  
  if (!ncol(x)) stop(paste0("ERROR: ", nm, " has no columns"))
  NULL  
}

check_family <- function(x, nm="family") {
  valid <- c("gaussian", "binomial", "poisson", "prophaz", "multinomial", "clogit")
  check_char_vec(x, nm, valid=valid, def=NULL, len=1) 
}

check_doseRRmod <- function(x, nm="doseRRmod") {
  valid <- c("ERR","EXP","LINEXP")
  check_char_vec(x, nm, valid=valid, def=NULL, len=1) 
}

check_Y <- function(v, data, family) {
  nm <- "Y"
  check_vars(data, v, nm, minlen=1, maxlen=0)
  vec <- data[, v, drop=TRUE]
  if(family != "multinomial"){
    binary <- nonneg <- integer <- 0
    if (family %in% c("binomial", "prophaz", "clogit")) {
      binary <- 1
    } else if (family == "poisson") {
      nonneg  <- 1
      integer <- 1 
    }
    check_num_vec(vec, nm, binary=binary, nonneg=nonneg, integer=integer)
    NULL
  } else { # multinomial 
    check_factor_vec(vec, nm)
    NULL
  }
}

check_D <- function(vars, data, methods) {
  
  nm <- "dosevars"
  check_vars(data, vars, nm, minlen=1, maxlen=0)
  for (v in vars) {
    vec <- data[, v, drop=TRUE]
    nm2 <- paste0(nm, ":", v)
    check_num_vec(vec, nm2)
  }
  if(length(vars)==1 &  any(c("ERC", "MCML", "FMA", "BMA") %in% methods)) stop("Multiple exposure replicates required for ERC, MCML, FMA, and BMA. With a single exposure vector, use RC")
  NULL
}

check_M <- function(vars, data) {
  
  nm <- "M"
  check_vars(data, vars, nm, minlen=0, maxlen=0)
  for (v in vars) {
    vec <- data[, v, drop=TRUE]
    nm2 <- paste0(nm, ":", v)
    check_num_vec(vec, nm2, binary=1)
  }
  NULL
}

check_X <- function(vars, data) {
  
  nm <- "X"
  if (is.null(vars)) return(NULL)
  for (v in vars) {
    vec <- data[, v, drop=TRUE]
    nm2 <- paste0(nm, ":", v)
    check_num_vec(vec, nm2)
  }
  NULL
}

check_offset <- function(v, data) {
  if (!length(v)) return(NULL)
  nm <- "offset"
  check_vars(data, v, nm, minlen=0, maxlen=0)
  check_num_vec(data[, v, drop=TRUE], nm, nonneg=1)
  
  NULL
}

check_setnr <- function(v, data) {
  nm <- "setnr"
  check_vars(data, v, nm, minlen=1, maxlen=1)
  check_num_vec(data[, v, drop=TRUE], nm, nonneg=1, integer=1)
  
  nset_noncontributing <- sum(table(data[,v, drop=TRUE])==1)
  if(nset_noncontributing>0) warning(paste0("Data contains ", nset_noncontributing, " matched sets of size 1 that do not contribute to model estimation"))
  
  NULL
}

check_entry_exit <- function(entry, exit, data) {
  
  nm1 <- "entry"
  nm2 <- "exit"
  check_vars(data, entry, nm1, minlen=0, maxlen=0)
  check_vars(data, exit,  nm2, minlen=1, maxlen=1)
  vec2 <- data[, exit, drop=TRUE]
  check_num_vec(vec2, nm2)
  if (length(entry)) {
    vec1 <- data[, entry, drop=TRUE]
    check_num_vec(vec1, nm1)
    tmp <- entry > exit
    tmp[is.na(tmp)] <- FALSE
    if (any(tmp)) stop(paste0("ERROR: ", nm1, " > ", nm2, " for some values")) 
  }
  
  NULL
}

check_methods <- function(x) {
  
  nm    <- "methods"
  valid <- c("RC","ERC","MCML","FMA","BMA")
  ret   <- check_char_vec(x, nm, valid=valid, def="RC", len=0)
  ret   <- unique(ret)
  ret
}

check_deg <- function(x) {
  if (!length(x)) return(2)
  nm <- "deg"
  check_integer(x, nm, minlen=1, maxlen=1, min=1, max=2) 
  x
}

check_inpar <- function(x, family, M, X, deg, multinom_levels=0) {
  if (!length(x)) return(NULL)
  nm <- "inpar"
  if (family == "gaussian") {
    len <- 2+length(X)+length(M)*deg+deg
  } else if (family %in% c("binomial", "poisson")) {
    len <- 1+length(X)+length(M)*deg+deg 
  } else if (family %in% c("prophaz", "clogit")) {
    len <- length(X)+length(M)*deg+deg
  } else if (family=="multinomial"){
    len <- (multinom_levels-1) * (1+length(X)+length(M)*deg+deg) 
  } else {
    stop("ERROR")
  }
  check_num_vec(x, nm, binary=0, nonneg=0, integer=0, len=len)
  x
}


check_factor_vec <- function(x, nm, len=0) {
  
  if (!is.factor(x)) stop(paste0("ERROR: ", nm, " must be numeric"))
  if (len && (len != length(x))) {
    stop(paste0("ERROR: ", nm, " must be a numeric vector of length ", len))
  }
  
  if (length(levels(x))<3) stop(paste0("ERROR: ", nm, " must have at least 3 levels"))
  
  if (length(levels(x)) > length(unique(x))) stop(paste0("ERROR: ", nm, " contains unused levels"))
  
  NULL
}

check_num_vec <- function(x, nm, binary=0, nonneg=0, integer=0, len=0) {
  
  if (!is.numeric(x)) stop(paste0("ERROR: ", nm, " must be numeric"))
  if (len && (len != length(x))) {
    stop(paste0("ERROR: ", nm, " must be a numeric vector of length ", len))
  }
  tmp <- !is.finite(x)
  if (any(tmp)) stop(paste0("ERROR: ", nm, " must contain finite values"))
  if (binary) {
    tmp <- !(x %in% 0:1)
    if (any(tmp)) stop(paste0("ERROR: ", nm, " must contain binary (0 - 1) values"))
  }
  if (nonneg) {
    tmp <- x < 0
    if (any(tmp)) stop(paste0("ERROR: ", nm, " must contain non-negative values"))
  }
  if (integer) check_integer(x, nm, minlen=0, maxlen=0, min=NULL, max=NULL)  
  
  NULL
}

check_vars <- function(data, vars, nm, minlen=0, maxlen=0) {
  
  nv <- length(vars)
  if (minlen && (minlen == maxlen) && (nv != minlen)) stop(paste0("ERROR: ", nm, " must have length ", minlen))
  if (nv < minlen) stop(paste0("ERROR: ", nm, " must have length >= ", minlen))
  if (!nv) return(NULL)
  if (!is.vector(vars)) stop(paste0("ERROR: ", nm, " must be a vector of indices or variable names"))
  
  nc <- ncol(data)
  cx <- colnames(data)
  if (is.numeric(vars)) {
    check_integer(vars, nm, minlen=minlen, maxlen=nc, min=1, max=nc)
  } else if (is.character(vars)) {
    check_char_vec(vars, nm, valid=cx, def=NULL) 
  } else {
    stop(paste0("ERROR: ", nm, " must be a vector of indices or variable names"))
  }
  if (any(duplicated(vars))) stop(paste0("ERROR: ", nm, " contains duplicated values"))
  
  unique(vars)
  
}


check_char_vec <- function(x, nm, valid=NULL, def=NULL, len=0) {
  
  n <- length(x)
  if (len && (n != len)) stop(paste0("ERROR: ", nm,  " must have length ", len))
  if (!n) return(def)
  if (!is.character(x)) stop(paste0("ERROR: ", nm,  " must be character"))
  if (length(valid)) {
    tmp <- !(x %in% valid)
    if (any(tmp)) {
      err <- getCharVecStr(x[tmp])
      stop(paste0("ERROR ", nm, " contains invalid values: ", err))
    }  
  }
  x
}


required_vars <- function(m) {
  
  vars <- c(m$dosevars, m$X, m$M)
  
  if (m$family %in% c("gaussian", "binomial", "poisson", "multinomial")) {
    vars <- c(vars, m$Y)
  }
  
  if (m$family == "poisson" && !is.null(m$offset)) {
    vars <- c(vars, m$offset)
  }
  
  if (m$family == "prophaz") {
    vars <- c(vars, m$status, m$exit)
    if (!is.null(m$entry)) vars <- c(vars, m$entry)
  }
  
  if (m$family == "clogit") {
    vars <- c(vars, m$status, m$setnr)
  }
  
  vars[!is.null(vars)]
}

check_integer <- function(x, nm, minlen=1, maxlen=0, min=NULL, max=NULL) {
  
  n <- length(x)
  if (minlen && (minlen == maxlen) && (n != minlen)) stop(paste0("ERROR: ", nm, " must have length ", minlen))
  if (minlen && (n < minlen)) stop(paste0("ERROR: ", nm, " must have length >= ", minlen))
  if (!is.numeric(x)) stop(paste0("ERROR: ", nm, " must be integer"))
  if (any(!is.finite(x))) stop(paste0("ERROR: ", nm, " must be integer"))
  if (any(x != floor(x))) stop(paste0("ERROR: ", nm, " must be integer"))
  if (length(min) && any(x < min)) stop(paste0("ERROR: ", nm, " must be >= ", min))
  if (length(max) && any(x > max)) stop(paste0("ERROR: ", nm, " must be <= ", max)) 
  
  NULL
}



check_string <- function(obj, valid, parm) {
  
  # obj:   A character string (length 1)
  # valid: Character vector of valid values
  # parm:  The name of the argument being checked
  
  errFlag <- 0
  
  # Check for errors
  if (!isString(obj)) errFlag <- 1 
  if (!errFlag) {
    obj <- trimws(obj)
    if (!(obj %in% valid)) errFlag <- 1
  }
  if (errFlag) {
    msg <- getCharVecStr(valid)
    msg <- paste0("ERROR: ", parm, " contains the invalid values ", msg)
    stop(msg)
  }
  
  obj
  
} # END: check.string

getVarNumbers <- function(vars, data) {
  
  if (!length(vars)) return(NULL)
  if (is.numeric(vars)) return(vars)
  cx  <- colnames(data)
  ret <- match(vars, cx)
  ret
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