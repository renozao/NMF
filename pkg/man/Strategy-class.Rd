\docType{class}
\name{Strategy-class}
\alias{Strategy-class}
\alias{name}
\alias{name,Strategy-method}
\alias{name-methods}
\alias{name<-}
\alias{name<-,Strategy,character-method}
\alias{name<--methods}
\title{Generic Strategy Class}
\usage{
  name(object, ...)

  \S4method{name}{Strategy}(object, all = FALSE)

  name(object, ...) <- value
}
\arguments{
  \item{object}{an R object with a defined \code{name}
  method}

  \item{...}{extra arguments to allow extension}

  \item{value}{replacement value}

  \item{all}{a logical that indicates if all the names
  associated with a strategy should be return
  (\code{TRUE}), or only the first (primary) one
  (\code{FALSE}).}
}
\description{
  This class defines a common interface for generic
  algorithm strategies (eg.,
  \code{\linkS4class{NMFStrategy}}).

  \code{name} and \code{name<-} gets and sets the name
  associated with an object. In the case of \code{Strategy}
  objects it is the the name of the algorithm.
}
\section{Slots}{
  \describe{

  \item{name}{character string giving the name of the
  algorithm}

  \item{package}{name of the package that defined the
  strategy.}

  }
}
\section{Methods}{
  \describe{

  \item{name}{\code{signature(object = "Strategy")}:
  Returns the name of an algorithm }

  \item{name}{\code{signature(object = "Strategy")}:
  Returns the name of an algorithm }

  \item{name<-}{\code{signature(object = "Strategy", value
  = "character")}: Sets the name(s) of an NMF algorithm }

  \item{name<-}{\code{signature(object = "Strategy", value
  = "character")}: Sets the name(s) of an NMF algorithm }

  }
}
\keyword{internal}
\keyword{methods}

