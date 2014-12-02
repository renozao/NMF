\docType{methods}
\name{NMFSeed}
\alias{NMFSeed}
\alias{NMFSeed,character-method}
\alias{NMFSeed-methods}
\alias{NMFSeed,NMFSeed-method}
\alias{removeNMFSeed}
\alias{setNMFSeed}
\title{\code{NMFSeed} is a constructor method that instantiate
\code{\linkS4class{NMFSeed}} objects.}
\usage{
  NMFSeed(key, method, ...)

  setNMFSeed(..., overwrite = isLoadingNamespace(),
    verbose = TRUE)

  removeNMFSeed(name, ...)
}
\arguments{
  \item{key}{access key as a single character string}

  \item{method}{specification of the seeding method, as a function that takes
at least the following arguments:
\describe{
\item{object}{uninitialised/empty NMF model, i.e. that it has 0 rows and
columns, but has already the rank requested in the call to \code{\link{nmf}}
or \code{\link{seed}}.}
\item{x}{target matrix}
\item{...}{extra arguments}
}}

  \item{...}{arguments passed to \code{NMFSeed} and used to initialise slots
in the \code{\linkS4class{NMFSeed}} object, or to \code{\link[pkgmaker]{pkgreg_remove}}.}

  \item{name}{name of the seeding method.}

  \item{overwrite}{logical that indicates if any existing NMF method with the
same name should be overwritten (\code{TRUE}) or not (\code{FALSE}),
in which case an error is thrown.}

  \item{verbose}{a logical that indicates if information about the registration
should be printed (\code{TRUE}) or not (\code{FALSE}).}
}
\description{
  \code{NMFSeed} is a constructor method that instantiate
\code{\linkS4class{NMFSeed}} objects.

NMF seeding methods are registered via the function \code{setNMFSeed}, which
stores them as \code{\linkS4class{NMFSeed}} objects in a dedicated registry.

\code{removeNMFSeed} removes an NMF seeding method from the registry.
}
\section{Methods}{
\describe{

\item{NMFSeed}{\code{signature(key = "character")}: Default method simply calls \code{\link{new}} with the same arguments.
}

\item{NMFSeed}{\code{signature(key = "NMFSeed")}: Creates an \code{NMFSeed} based on a template object (Constructor-Copy),
in particular it uses the \strong{same} name.
}

}

}
\keyword{methods}

