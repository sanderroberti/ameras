
\name{ecdfplot}
\alias{ecdfplot}
\title{Visualize multiple dose realizations}
\description{
Create a descriptive figure to visualize the distribution of dose and its uncertainty.
}
\usage{
ecdfplot(data, dosevars, xlab="Dose", 
  ylab="Cumulative distribution", log.xaxis=TRUE)
}
\arguments{
  \item{data}{data frame containing columns with dose vectors}
  \item{dosevars}{names or column indices of dose vectors.}
  \item{xlab}{label for the x-axis, default \code{"Dose"}}
  \item{ylab}{label for the y-axis, default \code{"Cumulative distribution"}}
  \item{log.xaxis}{logical, whether to use a log-scale for the dose axis (default \code{TRUE}). When \code{log.xaxis=TRUE}, any zeros in the doses are excluded while plotting.}
}
\details{
In the left panel, the empirical cumulative distribution function (ECDF) is
plotted for each dose realization. In other words, each curve shows one distribution
of dose across individuals. The spread within individual curves reflects the dose range
across individuals, while the spread between curves reflects between-realization variation
on the cohort level.

In the right panel, ECDFs are plotted for each individual,
showing distributions within individuals. A wide spread within individual curves
is indicative of large within-individual variation, while the spread between
curves reflects between-individual variation.

When using a log-scale for the x-axis, any zero dose values are excluded before plotting.
}

\examples{
\dontrun{
   data(data, package="ameras")
   ecdfplot(data,  dosevars=paste0("V", 1:10))
}
}