\docType{class}
\name{NMFStrategyFunction-class}
\alias{NMFStrategyFunction}
\alias{NMFStrategyFunction-class}
\title{Interface for Single Function NMF Strategies}
\description{
  This class implements the virtual interface
  \code{\link{NMFStrategy}} for NMF algorithms that are
  implemented by a single workhorse R function.
}
\section{Slots}{
  \describe{

  \item{algorithm}{a function that implements an NMF
  algorithm. It must have signature \code{(y='matrix',
  x='NMFfit')}, where \code{y} is the target matrix to
  approximate and \code{x} is the NMF model assumed to be
  seeded with an appropriate initial value -- as it is done
  internally by function \code{\link{nmf}}.

  Note that argument names currently do not matter, but it
  is recommended to name them as specified above.}

  }
}
\section{Methods}{
  \describe{

  \item{algorithm}{\code{signature(object =
  "NMFStrategyFunction")}: Returns the single R function
  that implements the NMF algorithm -- as stored in slot
  \code{algorithm}. }

  \item{algorithm<-}{\code{signature(object =
  "NMFStrategyFunction", value = "function")}: Sets the
  function that implements the NMF algorithm, stored in
  slot \code{algorithm}. }

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

  }
}

