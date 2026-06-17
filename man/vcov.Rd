\name{vcov.amerasfit}
\alias{vcov.amerasfit}
\title{Extract Variance-Covariance Matrices from an amerasfit Object}

\description{
  Extracts the variance-covariance matrices of the parameter estimates
  from a fitted \code{amerasfit} object.
}

\usage{
\method{vcov}{amerasfit}(object, methods = c("RC", "ERC", "MCML", "FMA", "BMA"), \dots)
}

\arguments{
  \item{object}{
    A fitted model object of class \code{amerasfit}, as returned by
    \code{\link{ameras}}.
  }
  \item{methods}{
    Character vector specifying which estimation methods to extract
    variance-covariance matrices for. One or more of \code{"RC"},
    \code{"ERC"}, \code{"MCML"}, \code{"FMA"}, and \code{"BMA"}.
    Defaults to all. Methods that were not run are skipped; an error is
    raised if none of the requested methods are present.
  }
  \item{\dots}{
    Additional arguments, currently unused.
  }
}

\value{
  If a single method is requested or only one method was run, a named
  numeric matrix with rows and columns corresponding to model
  parameters. If multiple methods are requested and multiple methods
  were run, a named list of such matrices, one per method.

  Note that for FMA and BMA the variance-covariance matrix is based
  on the generated samples rather than a Hessian-based
  approximation, and may differ in interpretation from the RC, ERC,
  and MCML variance-covariance matrices. For BMA, the reliability of 
  the variance-covariance matrix depends on MCMC convergence as indicated by \code{Rhat}.
}

\seealso{
  \code{\link{ameras}} for model fitting,
  \code{\link[=coef.amerasfit]{coef}} for extracting coefficients,
  \code{\link{summary.amerasfit}} for a summary including standard
  errors,
  \code{\link[stats]{vcov}} for the generic function.
}

\examples{
\donttest{
data("data", package="ameras")
dosevars <- paste0("V", 1:10)

fit <- ameras(Y.binomial ~ dose(all_of(dosevars), model="ERR"),
              data=data, family="binomial", methods=c("RC", "ERC"))

## Extract vcov for all available methods
vcov(fit)

## Extract vcov for a single method, returns a matrix
vcov(fit, methods="RC")

## Extract vcov for multiple methods, returns a named list
vcov(fit, methods=c("RC", "ERC"))
}
}
