\docType{methods}
\name{bterms<-}
\alias{bterms<-}
\alias{bterms<-,NMFstd-method}
\alias{bterms<--methods}
\alias{cterms<-}
\alias{cterms<-,NMFstd-method}
\alias{cterms<--methods}
\title{Annotations in NMF Models}
\usage{
  bterms(object) <- value

  cterms(object) <- value
}
\arguments{
  \item{object}{NMF object to be updated.}

  \item{value}{specification of the replacement value for
  fixed-terms.}
}
\description{
  NMF models may contain annotations for columns/rows
  and/or rows/features, in a similar way gene expression
  data are annotated \code{\linkS4class{ExpressionSet}}
  objects in Bioconductor.

  \code{cterms<-} sets fixed coefficient terms or indexes
  and should only be called on a newly created NMF object,
  i.e. in the constructor/factory generic
  \code{\link{nmfModel}}.
}
\details{
  Fixed Terms in NMF Models

  These functions are for internal use and should not be
  called by the end-user.

  They use \code{\link{molde.matrix}(~ -1 + ., data=value)}
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

