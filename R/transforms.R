# The following 2 helper functions make sure no large environments get
# attached to the transform and transform.jacobian function

make_transform <- function(index.t = NULL, lowlimit = NULL) {
  e <- new.env(parent = asNamespace("ameras"))
  e$index.t <- index.t
  e$lowlimit <- lowlimit

  f <- function(params, boundcheck = FALSE, ...) {
    args <- list(params = params, boundcheck = boundcheck, ...)
    if (!is.null(index.t)) {
      args$index.t <- index.t
    }
    if (!is.null(lowlimit)) {
      args$lowlimit <- lowlimit
    }
    do.call(transform1, args)
  }
  environment(f) <- e
  f
}

make_transform.jacobian <- function(index.t = NULL, lowlimit = NULL) {
  e <- new.env(parent = asNamespace("ameras"))
  e$index.t <- index.t

  f <- function(params, ...) {
    args <- list(params = params, ...)
    if (!is.null(index.t)) {
      args$index.t <- index.t
    }
    do.call(transform1.jacobian, args)
  }
  environment(f) <- e
  f
}

transform1 <- function(
  params,
  index.t = 1:length(params),
  lowlimit = rep(0, length(index.t)),
  boundcheck = FALSE,
  boundtol = 1e-3,
  ...
) {
  # Transforms from reparametrized to original scale, i.e., xi to beta

  if (length(index.t) != length(lowlimit)) {
    stop("Length mismatch between index.t and lowlimit")
  }
  if (any(!(index.t %in% 1:length(params)))) {
    stop("Incorrect indices for transformation specified")
  }
  params[index.t] <- exp(params[index.t]) + lowlimit
  if (boundcheck) {
    if (any(params[index.t] - lowlimit < boundtol)) {
      warning(paste0(
        "WARNING: one or multiple parameter estimates within ",
        boundtol,
        " of lower bounds. Try different bounds or starting values."
      ))
    }
  }
  return(params)
}
transform1.inv <- function(
  params,
  index.t = 1:length(params),
  lowlimit = rep(0, length(index.t)),
  ...
) {
  # Transforms from original scale to reparametrized, i.e., beta to xi
  if (length(index.t) != length(lowlimit)) {
    stop("Length mismatch between index.t and lowlimit")
  }
  if (any(!(index.t %in% 1:length(params)))) {
    stop("Incorrect indices for transformation specified")
  }
  params[index.t] <- log(params[index.t] - lowlimit)
  return(params)
}
transform1.jacobian <- function(params, index.t = 1:length(params), ...) {
  # Transforms from reparametrized to original scale, i.e., xi to beta
  if (any(!(index.t %in% 1:length(params)))) {
    stop("Incorrect indices for transformation specified")
  }
  grad <- rep(1, length(params))
  grad[index.t] <- exp(params[index.t])
  if (length(params) > 1) {
    return(diag(grad))
  } else {
    return(matrix(grad))
  }
}
