\docType{methods}
\name{registerDoBackend}
\alias{ForeachBackend}
\alias{ForeachBackend,ANY-method}
\alias{ForeachBackend,character-method}
\alias{ForeachBackend,cluster-method}
\alias{ForeachBackend,doMPI_backend-method}
\alias{ForeachBackend,doParallel_backend-method}
\alias{ForeachBackend-methods}
\alias{ForeachBackend,missing-method}
\alias{ForeachBackend,mpicluster-method}
\alias{ForeachBackend,NULL-method}
\alias{ForeachBackend,numeric-method}
\alias{getDoBackend}
\alias{getDoParHosts}
\alias{getDoParHosts,ANY-method}
\alias{getDoParHosts-methods}
\alias{getDoParNHosts}
\alias{register}
\alias{registerDoBackend}
\alias{setDoBackend}
\title{Utilities and Extensions for Foreach Loops}
\usage{
  registerDoBackend(object, ...)

  getDoBackend()

  setDoBackend(data)

  register(x, ...)

  ForeachBackend(object, ...)

  \S4method{ForeachBackend}{doParallel_backend}(object, cl)

  \S4method{ForeachBackend}{doMPI_backend}(object, cl)

  getDoParHosts(object, ...)

  getDoParNHosts(object)
}
\arguments{
  \item{object}{specification of a foreach backend, e.g.
  \sQuote{SEQ}, \sQuote{PAR} (for doParallel),
  \sQuote{MPI}, etc\ldots}

  \item{...}{extra arguments passed to the backend own
  registration function.}

  \item{data}{internal data of a foreach \%dopar\%
  backend.}

  \item{x}{specification of a foreach backend}

  \item{cl}{cluster specification: a cluster object or a
  numeric that indicates the number of nodes to use.}
}
\description{
  \code{registerDoBackend} is a unified register function
  for foreach backends.

  \code{getDoBackend} returns the internal data of the
  currently registered foreach \%dopar\% backend.

  \code{setDoBackend} is identical to
  \code{\link[foreach]{setDoPar}}, but returns the internal
  of the previously registered backend.

  \code{register} is a generic function that register
  objects. It is used to as a unified interface to register
  foreach backends.

  \code{ForeachBackend} is a factory method for foreach
  backend objects.

  \code{getDoParHosts} is a generic function that returns
  the hostname of the worker nodes used by a backend.

  \code{getDoParNHosts} returns the number of hosts used by
  a backend.
}
\section{Methods}{
  \describe{

  \item{ForeachBackend}{\code{signature(object = "ANY")}:
  Default method defined to throw an informative error
  message, when no other method was found. }

  \item{ForeachBackend}{\code{signature(object =
  "character")}: Creates a foreach backend object based on
  its name. }

  \item{ForeachBackend}{\code{signature(object =
  "missing")}: Creates a foreach backend object for the
  currently registered backend. }

  \item{ForeachBackend}{\code{signature(object = "NULL")}:
  Dummy method that returns \code{NULL}, defined for
  correct dispatch. }

  \item{ForeachBackend}{\code{signature(object =
  "cluster")}: Creates a doParallel foreach backend that
  uses the cluster described in \code{object}. }

  \item{ForeachBackend}{\code{signature(object =
  "numeric")}: Creates a doParallel foreach backend with
  \code{object} processes. }

  \item{ForeachBackend}{\code{signature(object =
  "doParallel_backend")}: doParallel-specific backend
  factory }

  \item{ForeachBackend}{\code{signature(object =
  "mpicluster")}: Creates a doMPI foreach backend that uses
  the MPI cluster described in \code{object}. }

  \item{ForeachBackend}{\code{signature(object =
  "doMPI_backend")}: doMPI-specific backend factory }

  \item{getDoParHosts}{\code{signature(object = "ANY")}:
  Default method that tries to heuristaically infer the
  number of hosts and in last resort temporarly register
  the backend and performs a foreach loop, to retrieve the
  nodename from each worker. }

  }
}
\keyword{internal}
\keyword{methods}

