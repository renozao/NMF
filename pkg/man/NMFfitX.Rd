\docType{methods}
\name{NMFfitX}
\alias{NMFfitX}
\alias{NMFfitX,list-method}
\alias{NMFfitX-methods}
\alias{NMFfitX,NMFfit-method}
\alias{NMFfitX,NMFfitX-method}
\title{Factory Method for Multiple NMF Run Objects}
\usage{
  NMFfitX(object, ...)

  \S4method{NMFfitX}{list}(object, ..., .merge = FALSE)
}
\arguments{
  \item{object}{an object from which is created an
  \code{NMFfitX} object}

  \item{...}{extra arguments used to pass values for slots}

  \item{.merge}{a logical that indicates if the fits should
  be aggregated, only keeping the best fit, and return an
  \code{NMFfitX1} object. If \code{FALSE}, an
  \code{NMFfitXn} object containing the data of all the
  fits is returned.}
}
\description{
  Factory Method for Multiple NMF Run Objects
}
\section{Methods}{
  \describe{

  \item{NMFfitX}{\code{signature(object = "list")}: Create
  an \code{NMFfitX} object from a list of fits. }

  \item{NMFfitX}{\code{signature(object = "NMFfit")}:
  Creates an \code{NMFfitX1} object from a single fit. This
  is used in \code{\link{nmf}} when only the best fit is
  kept in memory or on disk. }

  \item{NMFfitX}{\code{signature(object = "NMFfitX")}:
  Provides a way to aggregate \code{NMFfitXn} objects into
  an \code{NMFfitX1} object. }

  }
}
\keyword{internal}
\keyword{methods}

