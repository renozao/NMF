\docType{methods}
\name{consensus,NMFfitXn-method}
\alias{consensus,NMFfitXn-method}
\alias{plot.NMF.consensus}
\title{Computes the consensus matrix of the set of fits stored in \code{object}, as
the mean connectivity matrix across runs.}
\usage{
  \S4method{consensus}{NMFfitXn}(object, ...,
    no.attrib = FALSE)
}
\arguments{
  \item{object}{an object with a suitable \code{\link{predict}} method.}

  \item{...}{extra arguments to allow extension.
They are passed to \code{\link{predict}}, except for the \code{vector} and
\code{factor} methods.}

  \item{no.attrib}{a logical that indicates if attributes containing information
about the NMF model should be attached to the result (\code{TRUE}) or not
(\code{FALSE}).}
}
\description{
  This method returns \code{NULL} on an empty object.
The result is a matrix with several attributes attached, that are used by
plotting functions such as \code{\link{consensusmap}} to annotate the plots.
}
\keyword{methods}

