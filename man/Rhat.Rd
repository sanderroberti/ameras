\name{Rhat}
\alias{Rhat}
\alias{Rhat.amerasfit}
\title{Extract MCMC Convergence Diagnostics from an amerasfit Object}

\description{
  Extracts the Gelman-Rubin convergence diagnostics from the BMA
  results of a fitted \code{amerasfit} object.
}

\usage{
Rhat(x, \dots)

\method{Rhat}{amerasfit}(x, \dots)
}

\arguments{
  \item{x}{
    A fitted model object of class \code{amerasfit}, as returned by
    \code{\link{ameras}}, with BMA results present.
  }
  \item{\dots}{
    Additional arguments, currently unused.
  }
}

\value{
  A data frame with columns \code{Rhat} and \code{n.eff}, containing
  the Gelman-Rubin statistic and effective sample size for each
  parameter. Values of \code{Rhat} substantially above 1.05 indicate
  potential convergence problems, in which case longer chains via
  \code{niter.BMA} and \code{nburnin.BMA} are recommended. For
  single-chain BMA fits, convergence diagnostics cannot be estimated
  and these columns are returned as \code{NA}.
}

\seealso{
  \code{\link{ameras}},
  \code{\link{summary.amerasfit}},
  \code{\link{included_realizations.amerasfit}}
}

\examples{
\donttest{
data("data", package="ameras")
dosevars <- paste0("V", 1:10)
fit <- ameras(Y.binomial ~ dose(all_of(dosevars), model="ERR"),
              data=data, family="binomial", methods="BMA")
Rhat(fit)
}
}
