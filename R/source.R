dldd_clogit <- function(designmat, RRs) {
  nr <- nrow(designmat)
  nc <- ncol(designmat)
  ret <- as.numeric(rep(-9999, nc))

  # designmat must be passed in as a vector with columns stacked
  tmp <- .C(
    "C_dldd_clogit",
    as.numeric(designmat),
    as.integer(nr),
    as.integer(nc),
    as.numeric(RRs),
    ret = ret,
    PACKAGE = "ameras"
  )
  ret <- tmp$ret
  ret
}

dldd_prophaz <- function(entry, exit, status, RRs) {
  n <- length(entry) # Number of individuals
  ret <- as.numeric(rep(-9999, n)) # Initialize return vector with placeholders

  # Ensure all inputs are numeric vectors (entry, exit, status, RRs)
  entry <- as.numeric(entry)
  exit <- as.numeric(exit)
  status <- as.integer(status)
  RRs <- as.numeric(RRs)

  # Call the C function using .C
  tmp <- .C(
    "C_dldd_prophaz",
    as.numeric(entry),
    as.numeric(exit),
    as.integer(status),
    as.numeric(RRs),
    as.integer(n),
    ret = ret,
    PACKAGE = "ameras"
  )

  # Return the result as a numeric vector
  ret <- tmp$ret
  return(ret)
}


compute_ERCmatrix_clogit <- function(designmat, RRs, drdd, drdd2, status) {
  nr <- nrow(designmat)
  nc <- ncol(designmat)
  ret <- as.numeric(rep(-9999, nc * nc))

  tmp <- .C(
    "C_compute_ERCmatrix_clogit",
    as.numeric(designmat),
    as.integer(nr),
    as.integer(nc),
    as.numeric(RRs),
    as.numeric(drdd),
    as.numeric(drdd2),
    as.integer(status),
    ret = ret,
    PACKAGE = "ameras"
  )
  ret <- matrix(tmp$ret, nrow = nc, ncol = nc, byrow = FALSE)
  ret
}


compute_ERCmatrix_prophaz <- function(entry, exit, status, RRs, drdd, drdd2) {
  n <- length(entry) # Number of individuals
  ret <- as.numeric(rep(-9999, n * n)) # Initialize return vector with placeholders

  # Ensure all inputs are numeric vectors (entry, exit, status, RRs)
  entry <- as.numeric(entry)
  exit <- as.numeric(exit)
  status <- as.integer(status)
  RRs <- as.numeric(RRs)
  drdd <- as.numeric(drdd)
  drdd2 <- as.numeric(drdd2)

  # Call the C function using .C
  tmp <- .C(
    "C_compute_ERCmatrix_prophaz",
    as.numeric(entry),
    as.numeric(exit),
    as.integer(status),
    as.numeric(RRs),
    as.numeric(drdd),
    as.numeric(drdd2),
    as.integer(n),
    ret = ret,
    PACKAGE = "ameras"
  )

  # Return the result as a numeric vector
  ret <- matrix(tmp$ret, nrow = n, ncol = n, byrow = FALSE)
  ret
  return(ret)
}

loglik_prophaz_rcpp <- function(
  exit_t,
  entry_t,
  RR_entry,
  RR_exit,
  status_ord,
  loglim = 1e-30
) {
  n <- length(exit_t)
  K <- ncol(RR_exit)

  # sanity checks (minimal)
  stopifnot(
    length(entry_t) == n,
    length(status_ord) == n,
    nrow(RR_entry) == n,
    nrow(RR_exit) == n,
    ncol(RR_entry) == K
  )

  ret <- double(K)

  tmp <- .C(
    "C_loglik_prophaz_rcpp",
    exit_t = as.double(exit_t),
    entry_t = as.double(entry_t),
    RR_entry = as.double(RR_entry),
    RR_exit = as.double(RR_exit),
    status_ord = as.integer(status_ord),
    n = as.integer(n),
    K = as.integer(K),
    loglim = as.double(loglim),
    ret = ret,
    PACKAGE = "ameras"
  )

  tmp$ret
}
