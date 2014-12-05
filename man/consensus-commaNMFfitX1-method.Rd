\docType{methods}
\name{consensus,NMFfitX1-method}
\alias{consensus,NMFfitX1-method}
\title{Returns the consensus matrix computed while performing all NMF runs,
amongst which \code{object} was selected as the best fit.}
\usage{
  \S4method{consensus}{NMFfitX1}(object, no.attrib = FALSE)
}
\arguments{
  \item{object}{an object with a suitable \code{\link{predict}} method.}

  \item{no.attrib}{a logical that indicates if attributes containing information
about the NMF model should be attached to the result (\code{TRUE}) or not
(\code{FALSE}).}
}
\description{
  The result is the matrix stored in slot \sQuote{consensus}.
This method returns \code{NULL} if the consensus matrix is empty.
}
\keyword{methods}

