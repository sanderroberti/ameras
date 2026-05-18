\name{included_realizations}
\alias{included_realizations}
\alias{included_realizations.amerasfit}
\title{Extract Included Realizations from an amerasfit Object}

\description{
  Extracts the indices of dose realizations included in model averaging
  for FMA and/or BMA results of a fitted \code{amerasfit} object.
  Realizations are excluded from FMA if the optimization did not converge or
  the Hessian was not positive definite.
}

\usage{
included_realizations(x, \dots)

\method{included_realizations}{amerasfit}(x,
                              methods = c("FMA", "BMA"),
                              \dots)
}

\arguments{
  \item{x}{
    A fitted model object of class \code{amerasfit}, as returned by
    \code{\link{ameras}}, with FMA and/or BMA results present.
  }
  \item{methods}{
    Character vector specifying which method(s) to extract included
    realizations for. One or both of \code{"FMA"} and \code{"BMA"}.
    Defaults to both. Only methods that were run are returned.
  }
  \item{\dots}{
    Additional arguments, currently unused.
  }
}

\value{
  If a single method is requested or only one method was run, an
  integer vector of indices of the included dose realizations. If both
  FMA and BMA results are present and both are requested, a named list
  with elements \code{FMA} and \code{BMA}, each containing an integer
  vector of included realization indices.
}

\details{
  For FMA, realizations are excluded if the optimization did not
  converge or the Hessian was not invertible or not positive definite.
  If more than 20\% of realizations are excluded, a warning is issued
  during fitting.

  For BMA, all realizations are included by default. A subset of
  realizations can be specified via the \code{included.replicates.BMA}
  argument to \code{\link{ameras}}, for example by passing the indices
  returned by \code{included_realizations(fit, methods="FMA")} from a
  prior FMA fit.
}

\seealso{
  \code{\link{ameras}},
  \code{\link{Rhat.amerasfit}}
}

\examples{
\donttest{
data("data", package="ameras")
dosevars <- paste0("V", 1:10)

## FMA only
fit <- ameras(Y.binomial ~ dose(all_of(dosevars), model="ERR"),
              data=data, family="binomial", methods="FMA")
included_realizations(fit)
length(included_realizations(fit))

## Both FMA and BMA
fit2 <- ameras(Y.binomial ~ dose(all_of(dosevars), model="ERR"),
               data=data, family="binomial", methods=c("FMA", "BMA"))
included_realizations(fit2)
included_realizations(fit2, methods="FMA")
included_realizations(fit2, methods="BMA")
}
}