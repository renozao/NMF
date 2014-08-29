\docType{methods}
\name{residuals}
\alias{hasTrack}
\alias{residuals}
\alias{residuals<-}
\alias{residuals<--methods}
\alias{residuals-methods}
\alias{residuals<-,NMFfit-method}
\alias{residuals,NMFfit-method}
\alias{residuals,NMFfitX-method}
\alias{trackError}
\title{Residuals in NMF Models}
\usage{
  residuals(object, ...)

  \S4method{residuals}{NMFfit}(object, track = FALSE,
    niter = NULL, ...)

  residuals(object, ...)<-value

  \S4method{residuals}{NMFfit}(object, ..., niter = NULL,
    track = FALSE)<-value

  hasTrack(object, niter = NULL)

  trackError(object, value, niter, force = FALSE)
}
\arguments{
  \item{object}{an \code{NMFfit} object as fitted by function \code{\link{nmf}},
in single run mode.}

  \item{...}{extra parameters (not used)}

  \item{track}{a logical that indicates if the complete track of residuals
should be returned (if it has been computed during the fit), or only the last
value.}

  \item{niter}{specifies the iteration number for which one wants
to get/set/test a residual value. This argument is used only if not \code{NULL}}

  \item{value}{residual value}

  \item{force}{logical that indicates if the value should be added to the track
even if there already is a value for this iteration number or if the iteration
does not conform to the tracking interval \code{nmf.getOption('track.interval')}.}
}
\value{
  \code{residuals} returns a single numeric value if \code{track=FALSE}
or a numeric vector containing the residual values at some iterations.
The names correspond to the iterations at which the residuals were computed.
}
\description{
  The package NMF defines methods for the function
  \code{\link[stats]{residuals}} that returns the final
  residuals of an NMF fit or the track of the residuals
  along the fit process, computed according to the
  objective function associated with the algorithm that
  fitted the model.

  \code{residuals<-} sets the value of the last residuals,
  or, optionally, of the complete residual track.

  Tells if an \code{NMFfit} object contains a recorded
  residual track.

  \code{trackError} adds a residual value to the track of
  residuals.
}
\details{
  When called with \code{track=TRUE}, the whole residuals
  track is returned, if available. Note that method
  \code{\link{nmf}} does not compute the residuals track,
  unless explicitly required.

  It is a S4 methods defined for the associated generic
  functions from package \code{stats} (See
  \link[stats]{residuals}).
}
\note{
  Stricly speaking, the method \code{residuals,NMFfit} does
  not fulfill its contract as defined by the package
  \code{stats}, but rather acts as function
  \code{deviance}. The might be changed in a later release
  to make it behave as it should.
}
\section{Methods}{
\describe{

\item{residuals}{\code{signature(object = "NMFfit")}: Returns the residuals -- track -- between the target matrix and the NMF
fit \code{object}.
}

\item{residuals}{\code{signature(object = "NMFfitX")}: Returns the residuals achieved by the best fit object, i.e. the lowest
residual approximation error achieved across all NMF runs.
}

}

}
\seealso{
  Other stats: \code{\link{deviance}}, \code{\link{deviance,NMF-method}}, \code{\link{nmfDistance}}
}
\keyword{methods}

