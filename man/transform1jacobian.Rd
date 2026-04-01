
\name{transform1.jacobian}
\alias{transform1.jacobian}
\title{Jacobian of the exponential parameter transformation}

\description{
Computes the Jacobian matrix of \link{transform1}. Note that lower limits do not need to be specified as the Jacobian is independent of those
}

\usage{
transform1.jacobian(params, index.t=1:length(params), ... )
}
\arguments{
  \item{params}{input parameter vector (before transformation) to evaluate the Jacobian at}
  \item{index.t}{indices of parameters to be transformed (default all)}
  \item{...}{not used}
}
\value{
Jacobian matrix.
}



\examples{
params <- c(.1, .5, 1)
transform1.jacobian(params)
}
