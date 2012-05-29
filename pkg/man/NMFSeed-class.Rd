\docType{class}
\name{NMFSeed-class}
\alias{NMFSeed}
\alias{NMFSeed-class}
\title{Base class that defines the interface for NMF seeding methods.}
\description{
  This class implements a simple wrapper strategy object
  that defines a unified interface to seeding methods, that
  are used to initialise NMF models before fitting them
  with any NMF algorithm.
}
\section{Slots}{
  \describe{

  \item{name}{character string giving the name of the
  seeding strategy}

  \item{method}{workhorse function that implements the
  seeding strategy. It must have signature
  \code{(object="NMF", x="matrix", ...)} and initialise the
  NMF model \code{object} with suitable values for fitting
  the target matrix \code{x}.}

  }
}
\section{Methods}{
  \describe{

  \item{algorithm}{\code{signature(object = "NMFSeed")}:
  Returns the workhorse function of the seeding method
  described by \code{object}. }

  \item{algorithm<-}{\code{signature(object = "NMFSeed",
  value = "function")}: Sets the workhorse function of the
  seeding method described by \code{object}. }

  \item{name}{\code{signature(object = "NMFSeed")}: Returns
  the name of the seeding method described by
  \code{object}. }

  \item{name<-}{\code{signature(object = "NMFSeed", value =
  "character")}: Sets the name of the seeding method
  described by \code{object}.

  This method is not meant to be called directly }

  \item{show}{\code{signature( = "NMFSeed")}: Show method
  for objects of class \code{NMFSeed} }

  }
}

