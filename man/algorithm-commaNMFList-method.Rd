\docType{methods}
\name{algorithm,NMFList-method}
\alias{algorithm,NMFList-method}
\title{Returns the method names used to compute the NMF fits in the list.
It returns \code{NULL} if the list is empty.}
\usage{
  \S4method{algorithm}{NMFList}(object, string = FALSE,
    unique = TRUE)
}
\arguments{
  \item{string}{a logical that indicate whether the names should be collapsed
into a comma-separated string.}

  \item{unique}{a logical that indicates whether the result should contain the
set of method names, removing duplicated names.
This argument is forced to \code{TRUE} when \code{string=TRUE}.}

  \item{object}{an object computed using some algorithm, or that describes an algorithm
itself.}
}
\description{
  Returns the method names used to compute the NMF fits in the list.
It returns \code{NULL} if the list is empty.
}
\keyword{methods}

