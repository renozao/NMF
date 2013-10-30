\docType{methods}
\name{bterms<-}
\alias{bterms<-}
\alias{bterms<--methods}
\alias{bterms<-,NMFstd-method}
\alias{cterms<-}
\alias{cterms<--methods}
\alias{cterms<-,NMFstd-method}
\title{Fixed Terms in NMF Models}
\usage{
  bterms(object)<-value

  cterms(object)<-value
}
\arguments{
  \item{object}{NMF object to be updated.}

  \item{value}{specification of the replacement value for
  fixed-terms.}
}
\description{
  These functions are for internal use and should not be
  called by the end-user.

  \code{cterms<-} sets fixed coefficient terms or indexes
  and should only be called on a newly created NMF object,
  i.e. in the constructor/factory generic
  \code{\link{nmfModel}}.
}
\details{
  They use \code{\link{model.matrix}(~ -1 + ., data=value)}
  to generate suitable term matrices.
}
\section{Methods}{
  \describe{

  \item{bterms<-}{\code{signature(object = "NMFstd")}:
  Default method tries to coerce \code{value} into a
  \code{data.frame} with \code{\link{as.data.frame}}. }

  \item{cterms<-}{\code{signature(object = "NMFstd")}:
  Default method tries to coerce \code{value} into a
  \code{data.frame} with \code{\link{as.data.frame}}. }

  }
}
\keyword{internal}
\keyword{methods}

