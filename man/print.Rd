\name{print.amerasfit}
\alias{print.amerasfit}
\title{Simple summary for an amerasfit object}

\description{
  Prints a simple summary of a fitted \code{amerasfit} object.
}

\usage{
\method{print}{amerasfit}(x, digits = max(3, getOption("digits") - 3), \dots)
}

\arguments{
  \item{x}{
    A fitted model object of class \code{amerasfit}, as returned by
    \code{\link{ameras}}.
  }
  \item{digits}{
    Number of significant digits to be printed. Default is `max(3, getOption("digits") - 3)`
  }
  \item{\dots}{
    Additional arguments, currently unused.
  }
}

\value{
  Prints the `ameras` call, number of rows and dose replicates in the data, runtime,
  and model coefficients.
}



\seealso{
  \code{\link{ameras}} for model fitting,
  \code{\link[=summary.amerasfit]{summary}} for a more detailed summary of the fitted models
  including confidence intervals if computed.
}

\examples{
\donttest{
data("data", package = "ameras")
dosevars <- paste0("V", 1:10)

## Fit the model
fit <- ameras(data = data, family = "binomial", Y = "Y.binomial",
              dosevars = dosevars, methods = c("RC","ERC"), doseRRmod = "ERR")

## Default print
fit
print(fit)

## More digits
print(fit, digits=5)
}
}
