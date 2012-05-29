\docType{methods}
\name{NMFStrategy}
\alias{NMFStrategy}
\alias{NMFStrategy,character,character-method}
\alias{NMFStrategy,character,function-method}
\alias{NMFStrategy,character,missing-method}
\alias{NMFStrategy,character,NMFStrategy-method}
\alias{NMFStrategy-methods}
\alias{NMFStrategy,missing,character-method}
\alias{NMFStrategy,NMFStrategy,missing-method}
\alias{run,NMFStrategyFunction,matrix,NMFfit-method}
\alias{run,NMFStrategyIterative,matrix,NMFfit-method}
\alias{run,NMFStrategyIterativeX,matrix,NMFfit-method}
\alias{run,NMFStrategy,matrix,NMFfit-method}
\alias{run,NMFStrategy,matrix,NMF-method}
\title{Factory Method for NMFStrategy Objects}
\usage{
  NMFStrategy(name, method, ...)

  \S4method{run}{NMFStrategy,matrix,NMFfit}(object, y, x,
    ...)

  \S4method{run}{NMFStrategy,matrix,NMF}(object, y, x, ...)

  \S4method{run}{NMFStrategyFunction,matrix,NMFfit}(object,
    y, x, ...)

  \S4method{run}{NMFStrategyIterative,matrix,NMFfit}(object,
    y, x, .stop = NULL, maxIter = 2000, ...)

  \S4method{run}{NMFStrategyIterativeX,matrix,NMFfit}(object,
    y, x, maxIter, ...)
}
\arguments{
  \item{name}{name of an NMF algorithm.}

  \item{method}{definition of the algorithm}

  \item{...}{extra arguments passed to \code{\link{new}}.}

  \item{.stop}{specification of a stopping criterion, that
  is used instead of the one associated to the NMF
  algorithm. It may be specified as: \itemize{ \item the
  access key of a registered stopping criterion; \item a
  single integer that specifies the exact number of
  iterations to perform, which will be honoured unless a
  lower value is explicitly passed in argument
  \code{maxIter}. \item a single numeric value that
  specifies the stationnarity threshold for the objective
  function, used in with \code{\link{nmf.stop.stationary}};
  \item a function with signature
  \code{(object="NMFStrategy", i="integer", y="matrix",
  x="NMF", ...)}, where \code{object} is the
  \code{NMFStrategy} object that describes the algorithm
  being run, \code{i} is the current iteration, \code{y} is
  the target matrix and \code{x} is the current value of
  the NMF model. }}

  \item{maxIter}{maximum number of iterations to perform.}

  \item{object}{an object computed using some algorithm, or
  that describes an algorithm itself.}

  \item{y}{data object, e.g. a target matrix}

  \item{x}{a model object used as a starting point by the
  algorithm, e.g. a non-empty NMF model.}
}
\description{
  Creates NMFStrategy objects that wraps implementation of
  NMF algorithms into a unified interface.
}
\section{Methods}{
  \describe{

  \item{NMFStrategy}{\code{signature(name = "character",
  method = "function")}: Creates an
  \code{NMFStrategyFunction} object that wraps the function
  \code{method} into a unified interface.

  \code{method} must be a function with signature
  \code{(y="matrix", x="NMFfit", ...)}, and return an
  object of class \code{\linkS4class{NMFfit}}. }

  \item{NMFStrategy}{\code{signature(name = "character",
  method = "NMFStrategy")}: Creates an \code{NMFStrategy}
  object based on a template object (Constructor-Copy). }

  \item{NMFStrategy}{\code{signature(name = "NMFStrategy",
  method = "missing")}: Creates an \code{NMFStrategy} based
  on a template object (Constructor-Copy) using a randomly
  generated name. }

  \item{NMFStrategy}{\code{signature(name = "character",
  method = "character")}: Creates an \code{NMFStrategy}
  based on a registered NMF algorithm, that is used as a
  template. }

  \item{NMFStrategy}{\code{signature(name = "missing",
  method = "character")}: Creates an \code{NMFStrategy}
  based on a registered NMF algorithm (Constructor-Copy)
  using a randomly generated name. }

  \item{NMFStrategy}{\code{signature(name = "character",
  method = "missing")}: Creates an NMFStrategy, determining
  its type from the extra arguments passed in \code{...}:
  if there is an argument named \code{Update} then an
  \code{NMFStrategyIterative} is created, or if there is an
  argument named \code{algorithm} then an
  \code{NMFStrategyFunction} is created. Calls other than
  these generates an error. }

  \item{run}{\code{signature(object = "NMFStrategy", y =
  "matrix", x = "NMFfit")}: Pure virtual method defined for
  all NMF algorithms to ensure that a method \code{run} is
  defined by sub-classes of \code{NMFStrategy}.

  It throws an error if called directly. }

  \item{run}{\code{signature(object = "NMFStrategy", y =
  "matrix", x = "NMF")}: Method to run an NMF algorithm
  directly starting from a given NMF model. }

  \item{run}{\code{signature(object =
  "NMFStrategyFunction", y = "matrix", x = "NMFfit")}: Runs
  the NMF algorithms implemented by the single R function
  -- and stored in slot \code{'algorithm'} of
  \code{object}, on the data object \code{y}, using
  \code{x} as starting point. It is equivalent to calling
  \code{object@algorithm(y, x, ...)}.

  This method is usually not called directly, but only via
  the function \code{\link{nmf}}, which takes care of many
  other details such as seeding the computation, handling
  RNG settings, or setting up parallelisation. }

  \item{run}{\code{signature(object =
  "NMFStrategyIterative", y = "matrix", x = "NMFfit")}:
  Runs an NMF iterative algorithm on a target matrix
  \code{y}. }

  }
}
\keyword{methods}

