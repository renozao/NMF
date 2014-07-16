\docType{methods}
\name{.basis<-,NMFarray,matrix-method}
\alias{.basis<-,NMFarray,matrix-method}
\title{Replaces a slice of the basis array.}
\usage{
  \S4method{.basis}{NMFarray,matrix}(object, ...,
    slice = 1)<-value
}
\arguments{
  \item{object}{an object from which to extract the factor
  matrices, typically an object of class
  \code{\linkS4class{NMF}}.}

  \item{...}{extra arguments to allow extension and passed
  to the low-level access functions \code{.coef} and
  \code{.basis}.

  Note that these throw an error if used in replacement
  functions \code{}.}

  \item{value}{replacement value}
}
\description{
  Replaces a slice of the basis array.
}
\keyword{methods}

