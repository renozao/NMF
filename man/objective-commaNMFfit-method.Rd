\docType{methods}
\name{objective,NMFfit-method}
\alias{objective,NMFfit-method}
\title{Returns the objective function associated with the algorithm that computed the
fitted NMF model \code{object}, or the objective value with respect to a given
target matrix \code{y} if it is supplied.}
\usage{
  \S4method{objective}{NMFfit}(object, y)
}
\arguments{
  \item{y}{optional target matrix used to compute the
  objective value.}

  \item{object}{an object computed using some algorithm, or
  that describes an algorithm itself.}
}
\description{
  Returns the objective function associated with the
  algorithm that computed the fitted NMF model
  \code{object}, or the objective value with respect to a
  given target matrix \code{y} if it is supplied.
}
\keyword{methods}

