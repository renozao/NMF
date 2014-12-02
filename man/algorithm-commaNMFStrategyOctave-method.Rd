\docType{methods}
\name{algorithm,NMFStrategyOctave-method}
\alias{algorithm,NMFStrategyOctave-method}
\title{Returns the name of the Octave/Matlab function that implements the NMF algorithm -- as stored in
slot \code{algorithm}.}
\usage{
  \S4method{algorithm}{NMFStrategyOctave}(object,
    load = FALSE)
}
\arguments{
  \item{load}{logical that indicates if the algorithm should be loaded as an
R function.}

  \item{object}{an object computed using some algorithm, or that describes an algorithm
itself.}
}
\description{
  Returns the name of the Octave/Matlab function that implements the NMF algorithm -- as stored in
slot \code{algorithm}.
}
\keyword{methods}

