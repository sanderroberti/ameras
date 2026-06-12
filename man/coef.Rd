\name{coef.amerasfit}
\alias{coef.amerasfit}
\title{Estimated Coefficients for an amerasfit Object}

\description{
  Returns a data frame with all the parameters of a fitted
  \code{amerasfit} object. The resulting object has a column
  for every method supplied to \code{methods} when calling \code{ameras},
  with rows corresponding to parameters.
}

\usage{
\method{coef}{amerasfit}(object, \dots)
}

\arguments{
  \item{object}{
    A fitted model object of class \code{amerasfit}, as returned by
    \code{\link{ameras}}.
  }
  \item{\dots}{
    Additional arguments, currently unused.
  }
}

\value{
  Data frame with estimated model parameters. Column names correspond to the
  \code{methods} used in the \code{ameras} call, and row names correspond to parameter
  names.
}



\seealso{
  \code{\link{ameras}} for model fitting,
  \code{\link[=summary.amerasfit]{summary}} for a summary of the fitted model
  including confidence intervals if computed.
}

\examples{
\donttest{
data("data", package = "ameras")
dosevars <- paste0("V", 1:10)

## Fit the model
fit <- ameras(Y.binomial~dose(V1:V10, model="ERR"), data = data, family = "binomial", 
methods = c("RC","ERC"))
## Full matrix
coef(fit)

## Vector with RC parameters
coef(fit)$RC
}
}
