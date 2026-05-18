\name{CI}
\alias{CI}
\alias{CI.amerasfit}
\title{Extract Confidence Intervals from an amerasfit Object}

\description{
  Extracts confidence intervals from a fitted \code{amerasfit} object
  after they have been computed via
  \code{\link[=confint.amerasfit]{confint}}. This accessor is primarily
  intended for programmatic use, such as computing coverage in
  simulation studies, where printing is not desired.
}

\usage{
CI(x, \dots)

\method{CI}{amerasfit}(x, method = c("RC", "ERC", "MCML", "FMA", "BMA"),
   parm = NULL, \dots)
}

\arguments{
  \item{x}{
    A fitted model object of class \code{amerasfit}, as returned by
    \code{\link{ameras}}, with confidence intervals computed via
    \code{\link[=confint.amerasfit]{confint}}.
  }
  \item{method}{
    Character vector specifying which method(s) to extract confidence
    intervals for. One or more of \code{"RC"}, \code{"ERC"},
    \code{"MCML"}, \code{"FMA"}, and \code{"BMA"}. Defaults to all.
    Only methods that were run are returned.
  }
  \item{parm}{
    Either \code{NULL} to return all parameters, \code{"dose"} to
    return dose-related parameters only, \code{"all"} to explicitly
    return all parameters, or a character vector of specific parameter
    names. Defaults to \code{NULL}.
  }
  \item{\dots}{
    Additional arguments, currently unused.
  }
}

\value{
  If a single method is requested or only one method was run, a data
  frame with columns \code{lower} and \code{upper} and rows
  corresponding to parameters is returned. If multiple methods are requested and
  multiple methods were run, a named list of such data frames is returned, one per
  method.
}

\details{
  Confidence intervals must be computed via
  \code{\link[=confint.amerasfit]{confint}} before calling \code{CI}.
  If confidence intervals have not yet been computed, \code{NULL} is
  returned for all methods.
}

\seealso{
  \code{\link[=confint.amerasfit]{confint}} for computing confidence
  intervals,
  \code{\link{summary.amerasfit}} for a summary of the fitted model
  including confidence intervals,
  \code{\link{Rhat.amerasfit}},
  \code{\link{included_realizations.amerasfit}}
}

\examples{
data("data", package="ameras")
dosevars <- paste0("V", 1:10)

fit <- ameras(Y.binomial ~ dose(all_of(dosevars), model="ERR"),
              data=data, family="binomial", methods=c("RC", "FMA"))
fit <- confint(fit, method=c("wald.orig","percentile"))

## Extract all CIs for RC
CI(fit, method="RC")

## Extract dose parameters only for all methods
CI(fit, parm="dose")

## Extract specific parameter
CI(fit, method="RC", parm="dose")

}

