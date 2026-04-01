
\name{transform1}
\alias{transform1}
\title{Exponential parameter transformation}

\description{
Applies exponential transformation \eqn{f(\theta_i)=\exp(\theta_i)+L_i} to one or multiple components of parameter vector \eqn{\bm \theta}, where \eqn{L_i} are lower limits that can be different for each component
}

\usage{
transform1(params, index.t=1:length(params), lowlimit=rep(0,length(index.t)), 
  boundcheck=FALSE, boundtol=1e-3, ... )
}
\arguments{
  \item{params}{full input parameter vector}
  \item{index.t}{indices of parameters to be transformed (default all)}
  \item{lowlimit}{lower limits to be applied (default zero), where the k-th component of \code{lowlimit} is applied to the k-th index in \code{index.t} }
  \item{boundcheck}{whether to produce a warning when any of the transformed parameters are within \code{boundtol} of \code{lowlimit}}
  \item{boundtol}{tolerance for producing a warning for reaching the boundary}
  \item{...}{not used}
}
\value{
Transformed parameter vector.
}



\examples{
params <- c(.1, .5, 1)
transform1(params, lowlimit=c(0, -1, 1))
}
