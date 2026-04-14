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


