\docType{methods}
\name{runtime,NMFList-method}
\alias{runtime,NMFList-method}
\title{Returns the CPU time required to compute all NMF fits in the list.
It returns \code{NULL} if the list is empty.
If no timing data are available, the sequential time is returned.}
\usage{
  \S4method{runtime}{NMFList}(object, all = FALSE)
}
\arguments{
  \item{all}{logical that indicates if the CPU time of each fit should be
returned (\code{TRUE}) or only the total CPU time used to compute all
the fits in \code{object}.}

  \item{object}{an object computed using some algorithm, or that describes an algorithm
itself.}
}
\description{
  Returns the CPU time required to compute all NMF fits in the list.
It returns \code{NULL} if the list is empty.
If no timing data are available, the sequential time is returned.
}
\keyword{methods}

