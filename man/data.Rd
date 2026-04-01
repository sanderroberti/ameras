
\name{data}
\alias{data}

\title{
  Example data
}
\description{
 Data includes outcomes of all six supported types in the appropriately named columns. For proportional hazards regression, the observed exit time is \code{time} and event status is \code{event}. For conditional logistic regression, the matched set variable is \code{setnr}. The data has 10 exposure replicates in columns \code{V1}-\code{V10}.
}
%\details{}
\examples{

 data(data, package="ameras")

 # Display a few rows of the data
 data[1:5, ]

}
\keyword{data}
%\seealso{  }

