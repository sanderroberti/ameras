\name{convergence}
\alias{convergence}
\alias{convergence.amerasfit}
\title{Optimizer Convergence Diagnostics}

\description{
Compute or extract numerical gradient diagnostics for fitted \code{ameras}
models.
}

\usage{
convergence(object, ...)

\method{convergence}{amerasfit}(
  object,
  methods = c("RC", "ERC", "MCML"),
  data = NULL,
  recompute = FALSE,
  warn = FALSE,
  ...
)
}

\arguments{
\item{object}{an object of class \code{amerasfit}.}
\item{methods}{character vector of fitted methods for which diagnostics should
be returned. Gradient diagnostics are available for \code{"RC"}, \code{"ERC"},
and \code{"MCML"}.}
\item{data}{optional data frame used to reconstruct the likelihood when
diagnostics must be computed and the fit was created with
\code{keep.data=FALSE}.}
\item{recompute}{logical, whether to recompute the numerical gradient even when
stored diagnostics are available.}
\item{warn}{logical, whether to emit the same warning used during fitting when
\code{optim()} reported convergence but the numerical gradient remains large.}
\item{\dots}{currently unused.}
}

\value{
A data frame with one row per requested fitted method and columns:
\item{method}{method name.}
\item{convergence}{optimizer convergence code, where 0 indicates convergence.}
\item{gradient.rms}{root mean square of the numerical gradient at the fitted
parameter vector.}
\item{gradient.scaled.rms}{\code{gradient.rms} divided by
\code{max(1, abs(objective value))}.}
\item{gradient.max.abs}{maximum absolute numerical gradient component.}
\item{gradient.large}{logical indicator for the ameras diagnostic threshold
used to warn about potentially unreliable optimizer convergence.}
}

\details{
For objects fitted by recent versions of \pkg{ameras}, gradient diagnostics may
already be stored in the fitted method's \code{optim} component. For older
objects, or when \code{recompute=TRUE}, the likelihood is reconstructed from the
stored model metadata, as in profile-likelihood confidence interval
calculation, and the gradient is computed numerically at the stored optimizer
parameter vector.

These diagnostics do not change the fitted object. They are intended to help
identify fits where \code{optim()} returned convergence code 0 but the objective
is not locally flat enough to be fully reassuring. Poorly scaled continuous
covariates, such as uncentered calendar year variables, can contribute to this
kind of numerical behavior.

Gradient diagnostics are not returned for FMA or BMA. FMA combines
realization-specific optimizer fits and does not currently store the raw
realization-specific optima in the final object. BMA convergence is assessed
using MCMC diagnostics such as \code{\link{Rhat}}.
}

\examples{
\donttest{
  data(data, package="ameras")
  fit <- ameras(Y.gaussian ~ dose(V1:V10) + X1 + X2,
                data=data, family="gaussian", methods="RC")
  convergence(fit)
}
}

\seealso{
  \code{\link{ameras}}, \code{\link{Rhat}}
}
