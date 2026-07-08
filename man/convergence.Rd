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
\item{data}{optional data frame required when diagnostics must be recomputed
for a fit created with \code{keep.data=FALSE}.}
\item{recompute}{logical, whether to recompute the numerical gradient even when
stored diagnostics are available.}
\item{warn}{logical, whether to emit the same warning used during fitting when
\code{optim()} reported convergence but the optimizer diagnostics suggest that
the solution may not be fully stationary.}
\item{\dots}{currently unused.}
}

\value{
A data frame with one row per requested fitted method and columns:
\item{method}{method name.}
\item{optim.convergence}{optimizer convergence code, where 0 indicates
convergence.}
\item{gradient.rms}{root mean square of the numerical gradient at the fitted
parameter vector.}
\item{gradient.rms.scaled}{\code{gradient.rms} divided by
\code{max(1, abs(objective value))}.}
\item{newton.improvement}{approximate objective improvement still
available from one Newton step, computed as
\eqn{0.5 g^T H^{-1} g}.}
\item{newton.improvement.relative}{\code{newton.improvement}
divided by \code{max(1, abs(objective value))}. This makes the Newton-step
diagnostic less sensitive to the scale of the likelihood.}
\item{convergence.warning}{logical indicator for whether the ameras diagnostic
threshold would warn about potentially unreliable optimizer convergence.}
}

\details{
By default, \code{convergence()} returns stored diagnostics when available and
recomputes them only when needed or requested. If diagnostics must be recomputed
for a fit created with \code{keep.data=FALSE}, the original analysis data must
be supplied through \code{data}.

These diagnostics do not change the fitted object. They are intended to help
identify fits where \code{optim()} returned convergence code 0 but the objective
is not locally flat enough to be fully reassuring. The
\code{convergence.warning} column is based on \code{optim.convergence} and,
when available, \code{newton.improvement} and
\code{newton.improvement.relative}. Both Newton-improvement diagnostics must
exceed their warning thresholds: \code{newton.improvement > 1e-2} and
\code{newton.improvement.relative > 1e-6}. This is less sensitive to likelihood
scale than the raw gradient norm, because raw gradients can be large for
likelihoods with many observations even when the fitted coefficients are
effectively at the optimum. If the Hessian cannot be used for this calculation,
ameras falls back to \code{gradient.rms > 1e-3} and
\code{gradient.rms.scaled > 1e-4}.
Poorly scaled continuous covariates, such as uncentered calendar year variables,
can contribute to this kind of numerical behavior.

Gradient diagnostics are not returned for FMA or BMA. For BMA, use MCMC
diagnostics such as \code{\link{Rhat}}.
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
