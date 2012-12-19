\docType{class}
\name{NMFStrategyOctave-class}
\alias{NMFStrategyOctave-class}
\alias{show,NMFStrategyOctave-method}
\title{S4 Interface for Octave-Matlab NMF Algorithms}
\usage{
  \S4method{show}{NMFStrategyOctave}(object)
}
\arguments{
  \item{object}{Any R object}
}
\description{
  This class implements the virtual interface
  \code{\linkS4class{NMFStrategy}} for NMF algorithms that
  are implemented in Octave/Matlab, and provided as a set
  of .m files or as plain code.
}
\details{
  The \code{run} method for this class runs the algorithms
  via the \code{\link[RcppOctave]{RcppOctave}} package.
}
\section{Slots}{
  \describe{

  \item{algorithm}{character string that gives the name of
  the main Octave/Matlab function that implements the
  algorithm. The function must take at least two arguments:
  the target matrix and the initial NMF model, converted
  into an Octave list object, with elements corresponding
  to slots of the corresponding S4 class.}

  \item{mfiles}{character vector that contains a set of
  path to .m files. These files are (re-)sourced every time
  the strategy is called, and must be present at runtime in
  the current directory or in a directory from Octave
  path.}

  }
}
\section{Methods}{
  \describe{

  \item{algorithm}{\code{signature(object =
  "NMFStrategyOctave")}: Returns the name of the
  Octave/Matlab function that implements the NMF algorithm
  -- as stored in slot \code{algorithm}. }

  \item{algorithm<-}{\code{signature(object =
  "NMFStrategyOctave", value = "character")}: Sets the name
  of the Octave/Matlab function that implements the NMF
  algorithm. It is stored in slot \code{algorithm}. }

  \item{run}{\code{signature(object = "NMFStrategyOctave",
  y = "matrix", x = "NMFfit")}: Initialize method for
  \code{\linkS4class{NMFStrategyOctave}}. }

  }
}
\keyword{methods}

