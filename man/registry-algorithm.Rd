\name{methods-NMF}
\alias{existsNMFMethod}
\alias{getNMFMethod}
\alias{methods-NMF}
\alias{removeNMFMethod}
\alias{selectNMFMethod}
\title{Registry for NMF Algorithms}
\usage{
  selectNMFMethod(name, model, load = FALSE, exact = FALSE,
    all = FALSE, quiet = FALSE)

  getNMFMethod(...)

  existsNMFMethod(name, exact = TRUE)

  removeNMFMethod(name, ...)
}
\arguments{
  \item{name}{name of a registered NMF algorithm}

  \item{model}{class name of an NMF model, i.e. a class that inherits from class
\code{\linkS4class{NMF}}.}

  \item{load}{a logical that indicates if the selected algorithms should be loaded
into \code{NMFStrategy} objects}

  \item{all}{a logical that indicates if all algorithms that can fit \code{model}
should be returned or only the default or first found.}

  \item{quiet}{a logical that indicates if the operation should be performed quietly,
without throwing errors or warnings.}

  \item{...}{extra arguments passed to \code{\link[pkgmaker]{pkgreg_fetch}}
or \code{\link[pkgmaker]{pkgreg_remove}}.}

  \item{exact}{a logical that indicates if the access key should be matched
exactly (\code{TRUE}) or partially (\code{FALSE}).}
}
\value{
  \code{selectNMFMethod} returns a character vector or \code{NMFStrategy} objects,
or NULL if no suitable algorithm was found.
}
\description{
  Registry for NMF Algorithms

  \code{selectNMFMethod} tries to select an appropriate NMF
  algorithm that is able to fit a given the NMF model.

  \code{getNMFMethod} retrieves NMF algorithm objects from
  the registry.

  \code{existsNMFMethod} tells if an NMF algorithm is
  registered under the

  \code{removeNMFMethod} removes an NMF algorithm from the
  registry.
}

