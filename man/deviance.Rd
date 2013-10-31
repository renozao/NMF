\docType{methods}
\name{deviance}
\alias{deviance}
\alias{deviance-methods}
\alias{deviance,NMFfit-method}
\alias{deviance,NMFfitX-method}
\alias{deviance,NMF-method}
\alias{deviance,NMFStrategy-method}
\alias{nmfDistance}
\title{Distances and Objective Functions}
\usage{
  deviance(object, ...)

  \S4method{deviance}{NMF}(object, y,
    method = c("", "KL", "euclidean"), ...)

  nmfDistance(method = c("", "KL", "euclidean"))

  \S4method{deviance}{NMFfit}(object, y, method, ...)

  \S4method{deviance}{NMFStrategy}(object, x, y, ...)
}
\arguments{
  \item{y}{a matrix compatible with the NMF model
  \code{object}, i.e. \code{y} must have the same dimension
  as \code{fitted(object)}.}

  \item{method}{a character string or a function with
  signature \code{(x="NMF", y="matrix", ...)} that
  implements a distance measure between an NMF model
  \code{x} and a target matrix \code{y}, i.e. an objective
  function to use to compute the deviance. In
  \code{deviance}, it is passed to \code{nmfDistance} to
  get the function that effectively computes the deviance.}

  \item{...}{extra parameters passed to the objective
  function.}

  \item{x}{an NMF model that estimates \code{y}.}

  \item{object}{an object for which the deviance is
  desired.}
}
\value{
  \code{deviance} returns a nonnegative numerical value

  \code{nmfDistance} returns a function with least two
  arguments: an NMF model and a matrix.
}
\description{
  The NMF package defines methods for the generic
  \code{deviance} from the package \code{stats}, to compute
  approximation errors between NMF models and matrices,
  using a variety of objective functions.

  \code{nmfDistance} returns a function that computes the
  distance between an NMF model and a compatible matrix.
}
\section{Methods}{
  \describe{

  \item{deviance}{\code{signature(object = "NMF")}:
  Computes the distance between a matrix and the estimate
  of an \code{NMF} model. }

  \item{deviance}{\code{signature(object = "NMFfit")}:
  Returns the deviance of a fitted NMF model.

  This method returns the final residual value if the
  target matrix \code{y} is not supplied, or the
  approximation error between the fitted NMF model stored
  in \code{object} and \code{y}. In this case, the
  computation is performed using the objective function
  \code{method} if not missing, or the objective of the
  algorithm that fitted the model (stored in slot
  \code{'distance'}).

  If not computed by the NMF algorithm itself, the value is
  automatically computed at the end of the fitting process
  by the function \code{\link{nmf}}, using the objective
  function associated with the NMF algorithm, so that it
  should always be available. }

  \item{deviance}{\code{signature(object = "NMFfitX")}:
  Returns the deviance achieved by the best fit object,
  i.e. the lowest deviance achieved across all NMF runs. }

  \item{deviance}{\code{signature(object = "NMFStrategy")}:
  Computes the value of the objective function between the
  estimate \code{x} and the target \code{y}. }

  }
}
\seealso{
  Other stats: \code{\link{deviance,NMF-method}},
  \code{\link{hasTrack}}, \code{\link{residuals}},
  \code{\link{residuals<-}}, \code{\link{trackError}}
}
\keyword{methods}

