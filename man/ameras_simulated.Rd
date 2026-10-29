
\name{ameras_simulated}
\alias{ameras_simulated}
\docType{data}

\title{Simulated Data for ameras Examples}

\description{
A simulated dataset containing 10 exposure realizations, covariates, effect
modifiers, and outcomes for all six model families supported by
\pkg{ameras}. The outcomes were generated under linear excess relative risk
models. The dataset is intended for examples and demonstrations, not as a
benchmark for recovering known simulation parameters.
}

\format{
A data frame with 3,000 rows and 22 variables:
\describe{
  \item{Y.gaussian}{A numeric continuous outcome.}
  \item{Y.binomial}{A binary outcome coded as 0 or 1.}
  \item{Y.poisson}{A nonnegative integer count outcome.}
  \item{time}{Observed exit or follow-up time for proportional hazards models.}
  \item{status}{Event indicator for proportional hazards models, coded as 0 or 1.}
  \item{setnr}{Matched-set identifier for conditional logistic regression.}
  \item{Y.clogit}{Case indicator for conditional logistic regression, coded as 0 or 1.}
  \item{Y.multinomial}{A factor outcome with levels \code{"1"}, \code{"2"}, and \code{"3"}.}
  \item{X1, X2}{Binary covariates coded as 0 or 1.}
  \item{M1, M2}{Binary effect modifiers coded as 0 or 1.}
  \item{V1, V2, V3, V4, V5, V6, V7, V8, V9, V10}{Numeric simulated
  exposure realizations.}
}
}

\source{
Simulated by the package authors for illustrating the functionality of
\pkg{ameras}. The dataset contains no real participant data.
}

\examples{
# The dataset is available directly because ameras uses lazy data.
head(ameras_simulated)
}

\keyword{datasets}
