\docType{class}
\name{NMFStrategy-class}
\alias{NMFStrategy-class}
\alias{is.mixed}
\alias{name}
\alias{name,NMFSeed-method}
\alias{name,NMFStrategy-method}
\alias{name-methods}
\alias{name<-}
\alias{name<-,NMFSeed,character-method}
\alias{name<-,NMFStrategy,character-method}
\alias{name<--methods}
\alias{objective,NMFStrategy-method}
\alias{objective<-,NMFStrategy,character-method}
\alias{objective<-,NMFStrategy,function-method}
\alias{show,NMFStrategy-method}
\title{Virtual Interface for NMF Algorithms}
\usage{
  \S4method{show}{NMFStrategy}(object)

  name(object, ...)

  \S4method{name}{NMFStrategy}(object, all = FALSE)

  name(object, ...) <- value

  \S4method{objective}{NMFStrategy}(object, y, x, ...)

  \S4method{objective}{NMFStrategy,character}(object) <-
    value

  \S4method{objective}{NMFStrategy,function}(object) <-
    value

  is.mixed(object)
}
\arguments{
  \item{object}{an NMFStrategy object}

  \item{...}{extra arguments to allow extension}

  \item{value}{replacement value}

  \item{all}{a logical that indicates if all the names
  associated with a strategy should be return
  (\code{TRUE}), or only the first (primary) one
  (\code{FALSE}).}

  \item{y}{data object, e.g. a target matrix}

  \item{x}{a model object used as a starting point by the
  algorithm, e.g. a non-empty NMF model.}
}
\description{
  This class partially implements the generic interface
  defined for general algorithms defined in the package NMF
  (see \code{\link{algorithmic-NMF}}).

  \code{name} and \code{name<-} gets and sets the name
  associated with an object. In the case of
  \code{NMFStrategy} object it is the the name of the
  algorithm.

  \code{is.mixed} tells if an NMF algorithm works on
  mixed-sign data.
}
\section{Slots}{
  \describe{

  \item{name}{character string giving the name of the
  algorithm}

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

  \item{NMFStrategy}{\code{signature(name = "NMFStrategy",
  method = "missing")}: Creates an \code{NMFStrategy} based
  on a template object (Constructor-Copy) using a randomly
  generated name. }

  \item{canFit}{\code{signature(x = "NMFStrategy", y =
  "character")}: Tells if an NMF algorithm can fit a given
  class of NMF models }

  \item{canFit}{\code{signature(x = "NMFStrategy", y =
  "NMF")}: Tells if an NMF algorithm can fit the same class
  of models as \code{y} }

  \item{modelname}{\code{signature(object =
  "NMFStrategy")}: Returns the model(s) that an NMF
  algorithm can fit. }

  \item{name}{\code{signature(object = "NMFStrategy")}:
  Returns the name of an algorithm }

  \item{name}{\code{signature(object = "NMFStrategy")}:
  Returns the name of an algorithm }

  \item{name}{\code{signature(object = "NMFSeed")}: Returns
  the name of the seeding method described by
  \code{object}. }

  \item{name<-}{\code{signature(object = "NMFStrategy",
  value = "character")}: Sets the name(s) of an NMF
  algorithm }

  \item{name<-}{\code{signature(object = "NMFStrategy",
  value = "character")}: Sets the name(s) of an NMF
  algorithm }

  \item{name<-}{\code{signature(object = "NMFSeed", value =
  "character")}: Sets the name of the seeding method
  described by \code{object}.

  This method is not meant to be called directly }

  \item{objective}{\code{signature(object =
  "NMFStrategy")}: Gets the objective function associated
  with an NMF algorithm, or computes the objective value
  for an NMF model with respect to a given target matrix. }

  \item{objective}{\code{signature(object =
  "NMFStrategy")}: Gets the objective function associated
  with an NMF algorithm, or computes the objective value
  for an NMF model with respect to a given target matrix. }

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

