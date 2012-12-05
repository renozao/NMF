\docType{methods}
\name{NMF-Defunct}
\alias{metaHeatmap}
\alias{metaHeatmap,matrix-method}
\alias{metaHeatmap-methods}
\alias{metaHeatmap,NMFfitX-method}
\alias{metaHeatmap,NMF-method}
\alias{NMF-Defunct}
\title{Defunct Functions and Classes in the NMF Package}
\usage{
  metaHeatmap(object, ...)
}
\arguments{
  \item{object}{an R object}

  \item{...}{other arguments}
}
\description{
  Defunct Functions and Classes in the NMF Package
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

