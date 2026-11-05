\name{coef.amerasfit}
\alias{coef.amerasfit}
\title{Estimated Coefficients for an amerasfit Object}

\description{
  Returns the estimated parameters of a fitted \code{amerasfit} object. For
  an object containing one fitted method, the result is a named numeric
  vector. For an object containing multiple methods, the result is a data
  frame with one column per method and rows corresponding to parameters.
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
  A named numeric vector when one method was fitted. When multiple methods
  were fitted, a data frame whose column names correspond to the methods and
  whose row names correspond to parameter names.
}



\seealso{
  \code{\link{ameras}} for model fitting,
  \code{\link[=summary.amerasfit]{summary}} for a summary of the fitted model
  including confidence intervals if computed.
}

\examples{
\donttest{% slow
dosevars <- paste0("V", 1:10)

## Fit the model
fit <- ameras(Y.binomial~dose(V1:V10, model="ERR"), data = ameras_simulated, family = "binomial",
methods = c("RC","ERC"))
## Data frame for a fit containing multiple methods
coef(fit)

## Vector with RC parameters
coef(fit)$RC
}
}
