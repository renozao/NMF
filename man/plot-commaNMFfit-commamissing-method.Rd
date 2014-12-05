\docType{methods}
\name{plot,NMFfit,missing-method}
\alias{plot,NMFfit,missing-method}
\title{Plots the residual track computed at regular interval during the fit of
the NMF model \code{x}.}
\usage{
  \S4method{plot}{NMFfit,missing}(x, y, skip = -1, ...)
}
\arguments{
  \item{x}{an NMF object as returned by \code{\link{nmf}}.}

  \item{y}{other object (missing)}

  \item{skip}{an integer that indicates the number of points to skip/remove from the beginning
of the curve.
If \code{skip=1L} (default) only the initial residual -- that is computed before any iteration, is
skipped, if present in the track (it associated with iteration 0).}

  \item{...}{extra arguments passed to base graphics \code{\link[graphics]{plot}} function.}
}
\description{
  Plots the residual track computed at regular interval during the fit of
the NMF model \code{x}.
}
\keyword{methods}

