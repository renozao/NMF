\docType{class}
\name{NMFStrategy-class}
\alias{is.mixed}
\alias{NMFStrategy-class}
\alias{objective<-,NMFStrategy,character-method}
\alias{objective<-,NMFStrategy,function-method}
\alias{objective,NMFStrategy-method}
\alias{show,NMFStrategy-method}
\title{Virtual Interface for NMF Algorithms}
\usage{
  \S4method{show}{NMFStrategy}(object)

  \S4method{objective}{NMFStrategy}(object)

  \S4method{objective}{NMFStrategy,character}(object) <-
    value

  \S4method{objective}{NMFStrategy,function}(object) <-
    value

  is.mixed(object)
}
\arguments{
  \item{object}{Any R object}

  \item{value}{replacement value}
}
\description{
  This class partially implements the generic interface
  defined for general algorithms defined in the \pkg{NMF}
  package (see \code{\link{algorithmic-NMF}}).

  \code{is.mixed} tells if an NMF algorithm works on
  mixed-sign data.
}
\section{Slots}{
  \describe{

  \item{objective}{the objective function associated with
  the algorithm (Frobenius, Kullback-Leibler, etc...).  It
  is either an access key of a registered objective
  function or a function definition.  In the latter case,
  the given function must have the following signature
  \code{(x="NMF", y="matrix")} and return a nonnegative
  real value.}

  \item{model}{a character string giving either the
  (sub)class name of the NMF-class instance used and
  returned by the strategy, or a function name.}

  \item{mixed}{a logical that indicates if the algorithm
  works on mixed-sign data.}

  }
}
\section{Methods}{
  \describe{

  \item{canFit}{\code{signature(x = "NMFStrategy", y =
  "character")}: Tells if an NMF algorithm can fit a given
  class of NMF models }

  \item{canFit}{\code{signature(x = "NMFStrategy", y =
  "NMF")}: Tells if an NMF algorithm can fit the same class
  of models as \code{y} }

  \item{deviance}{\code{signature(object = "NMFStrategy")}:
  Computes the value of the objective function between the
  estimate \code{x} and the target \code{y}. }

  \item{modelname}{\code{signature(object =
  "NMFStrategy")}: Returns the model(s) that an NMF
  algorithm can fit. }

  \item{NMFStrategy}{\code{signature(name = "NMFStrategy",
  method = "missing")}: Creates an \code{NMFStrategy} based
  on a template object (Constructor-Copy), in particular it
  uses the \strong{same} name. }

  \item{objective}{\code{signature(object =
  "NMFStrategy")}: Gets the objective function associated
  with an NMF algorithm.

  It is used in
  \code{\link[=deviance,NMFStrategy-method]{deviance}} to
  compute the objective value for an NMF model with respect
  to a given target matrix. }

  \item{objective}{\code{signature(object =
  "NMFStrategy")}: Gets the objective function associated
  with an NMF algorithm.

  It is used in
  \code{\link[=deviance,NMFStrategy-method]{deviance}} to
  compute the objective value for an NMF model with respect
  to a given target matrix. }

  \item{objective<-}{\code{signature(object =
  "NMFStrategy", value = "character")}: Sets the objective
  function associated with an NMF algorithm, with a
  character string that must be a registered objective
  function. }

  \item{objective<-}{\code{signature(object =
  "NMFStrategy", value = "character")}: Sets the objective
  function associated with an NMF algorithm, with a
  character string that must be a registered objective
  function. }

  \item{objective<-}{\code{signature(object =
  "NMFStrategy", value = "function")}: Sets the objective
  function associated with an NMF algorithm, with a
  function that computes the approximation error between an
  NMF model and a target matrix. }

  \item{objective<-}{\code{signature(object =
  "NMFStrategy", value = "function")}: Sets the objective
  function associated with an NMF algorithm, with a
  function that computes the approximation error between an
  NMF model and a target matrix. }

  \item{run}{\code{signature(object = "NMFStrategy", y =
  "matrix", x = "NMFfit")}: Pure virtual method defined for
  all NMF algorithms to ensure that a method \code{run} is
  defined by sub-classes of \code{NMFStrategy}.

  It throws an error if called directly. }

  \item{run}{\code{signature(object = "NMFStrategy", y =
  "matrix", x = "NMF")}: Method to run an NMF algorithm
  directly starting from a given NMF model. }

  }
}
\keyword{internal}
\keyword{methods}

