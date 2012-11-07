\name{pmax.inplace}
\alias{neq.constraints.inplace}
\alias{pmax.inplace}
\title{Updating Objects In Place}
\usage{
  pmax.inplace(x, lim, skip = NULL)

  neq.constraints.inplace(x, constraints, ratio = NULL,
    value = NULL, copy = FALSE)
}
\arguments{
  \item{x}{an object to update in place.}

  \item{lim}{lower threshold value}

  \item{skip}{indexes to skip}

  \item{constraints}{constraint specification.}

  \item{ratio}{fixed ratio on which the constraint
  applies.}

  \item{value}{fixed value to enforce.}

  \item{copy}{a logical that indicates if \code{x} should
  be updated in place or not.}
}
\description{
  These functions modify objects (mainly matrix objects) in
  place, i.e. they act directly on the C pointer. Due to
  their side-effect, they are not meant to be called by the
  end-user.

  \code{neq.constraints.inplace} apply unequality
  constraints in place.
}
\details{
  \code{pmax.inplace} is a version of \code{\link{pmax}}
  that updates its first argument.
}
\keyword{internal}

