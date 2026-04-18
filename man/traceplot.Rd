
\name{traceplot}
\alias{traceplot}
\alias{traceplot.amerasfit}
\title{Traceplots for MCMC samples}
\description{
Produce MCMC traceplots for \code{amerasfit} objects.
}
\usage{
traceplot(object, ...)

\method{traceplot}{amerasfit}(object, iter = 5000, Rhat = TRUE, n.eff = TRUE, pdf = FALSE, ...)
}
\arguments{
  \item{object}{a \code{amerasfit} object containing BMA output to be plotted}
  \item{iter}{number of iterations to include in the traceplot (defaults to last 5000)}
  \item{Rhat}{logical; whether to include R-hat diagnostics in the plot (default TRUE)}
  \item{n.eff}{logical; whether to include effective sample size in the plot (default TRUE)}
  \item{pdf}{logical; whether to save the output as a PDF (default FALSE)}
  \item{...}{additional arguments passed to \code{MCMCtrace}}
}
\details{
Wrapper for \code{MCMCvis::MCMCtrace} to produce MCMC diagnostic plots. See \code{?MCMCtrace} for more plotting options that can be provided through \code{...}.
}
\value{
Traceplots and posterior density plots.
}
\seealso{
\code{\link[MCMCvis]{MCMCtrace}}
}
\examples{
\donttest{
   data(data, package="ameras")
   fit <- ameras(Y.gaussian~dose(V1:V10), data, methods="BMA")
   traceplot(fit)
}
}