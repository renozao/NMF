\name{parallel-NMF}
\alias{gVariable}
\alias{hostfile}
\alias{parallel-NMF}
\alias{ts_eval}
\alias{ts_tempfile}
\title{Utilities for Parallel Computations}
\usage{
  ts_eval(mutex = synchronicity::boost.mutex(),
    verbose = FALSE)

  ts_tempfile(pattern = "file", ..., host = TRUE,
    pid = TRUE)

  hostfile(pattern = "file", tmpdir = tempdir(),
    fileext = "", host = TRUE, pid = TRUE)

  gVariable(init, shared = FALSE)
}
\arguments{
  \item{mutex}{a mutex or a mutex descriptor.
If missing, a new mutex is created via the function \code{\link[synchronicity]{boost.mutex}}.}

  \item{verbose}{a logical that indicates if messages should be printed when
locking and unlocking the mutex.}

  \item{...}{extra arguments passed to \code{\link[base]{tempfile}}.}

  \item{host}{logical that indicates if the host machine name should
be appear in the filename.}

  \item{pid}{logical that indicates if the current process id
be appear in the filename.}

  \item{init}{initial value}

  \item{shared}{a logical that indicates if the variable should be stored in shared
memory or in a local environment.}

  \item{pattern}{a non-empty character vector giving the initial part
    of the name.}

  \item{tmpdir}{a non-empty character vector giving the directory name}

  \item{fileext}{a non-empty character vector giving the file extension}
}
\description{
  Utilities for Parallel Computations

  \code{ts_eval} generates a thread safe version of
  \code{\link{eval}}. It uses boost mutexes provided by the
  \code{\link[synchronicity]{synchronicity}} package. The
  generated function has arguments \code{expr} and
  \code{envir}, which are passed to \code{\link{eval}}.

  \code{ts_tempfile} generates a \emph{unique} temporary
  filename that includes the name of the host machine
  and/or the caller's process id, so that it is thread
  safe.

  \code{hostfile} generates a temporary filename composed
  with the name of the host machine and/or the current
  process id.

  \code{gVariable} generates a function that access a
  global static variable, possibly in shared memory (only
  for numeric matrix-coercible data in this case). It is
  used primarily in parallel computations, to preserve data
  accross computations that are performed by the same
  process.
}

