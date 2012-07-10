\docType{methods}
\name{NMF-deprecated}
\alias{NMF-deprecated}
\alias{join}
\alias{metaHeatmap}
\alias{metaHeatmap,NMF-method}
\alias{metaHeatmap,NMFfitX-method}
\alias{metaHeatmap,matrix-method}
\alias{metaHeatmap-methods}
\title{Deprecated Functions in the Package NMF}
\usage{
  metaHeatmap(object, ...)

  join(...)
}
\arguments{
  \item{object}{an R object}

  \item{...}{extra arguments}
}
\description{
  Deprecated Functions in the Package NMF
}
\section{Methods}{
  \describe{

  \item{metaHeatmap}{\code{signature(object = "matrix")}:
  Defunct method substituted by \code{\link{aheatmap}}. }

  \item{metaHeatmap}{\code{signature(object = "NMF")}:
  Deprecated method that is substituted by
  \code{\link{coefmap}} and \code{\link{basismap}}. }

  \item{metaHeatmap}{\code{signature(object = "NMFfitX")}:
  Deprecated method subsituted by
  \code{\link{consensusmap}}. }

  }
}
\keyword{methods}

