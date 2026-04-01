
\title{Inverse of exponential parameter transformation}
\alias{transform1.inv}
\name{transform1.inv}

\description{
Inverse of \code{transform1} for the purpose of deriving initial values.
}

\usage{
transform1.inv(params, index.t=1:length(params), lowlimit=rep(0,length(index.t)), ... )
}
\arguments{
  \item{params}{full input parameter vector}
  \item{index.t}{indices of parameters to be transformed (default all)}
  \item{lowlimit}{lower limits to be applied (default zero), where the k-th component of \code{lowlimit} is applied to the k-th index in \code{index.t} }
  \item{...}{not used}
}
\value{
Transformed parameter vector.
}



\examples{
params <- c(.1, .5, 1) # Desired initial values on original scale
transform1.inv(params, lowlimit=c(0, -1, 1)) # Initial values to use on transformed scale
}
