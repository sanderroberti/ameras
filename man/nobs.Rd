\name{nobs.amerasfit}
\alias{nobs.amerasfit}
\title{Number of Observations Used to Fit an ameras Model}

\description{
  Returns the number of observations used to fit an \code{amerasfit} object
  after missing-value handling and any family-specific row exclusions.
}

\usage{
\method{nobs}{amerasfit}(object, \dots)
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
  The number of observations used for fitting. This excludes rows omitted by
  \code{na.action}, uninformative matched-set rows for conditional logistic
  regression, and zero-follow-up rows for proportional hazards regression
  with entry times.
}

\seealso{
  \code{\link{ameras}} for model fitting,
  \code{\link[=summary.amerasfit]{summary}} for detailed row-exclusion counts,
  and \code{\link[stats]{nobs}} for the generic function.
}

\examples{
example_data <- ameras_simulated[1:100, ]
example_data$X1[1] <- NA

fit <- ameras(
  Y.gaussian ~ dose(V1:V10) + X1,
  data = example_data,
  family = "gaussian",
  methods = "RC",
  na.action = na.omit
)

nobs(fit)
}
