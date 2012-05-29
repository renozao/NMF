\name{setupBackend}
\alias{setupBackend}
\alias{setupLibPaths}
\alias{setupRNG}
\alias{setupSharedMemory}
\alias{setupTempDirectory}
\title{Computational Setup Functions}
\usage{
  setupBackend(spec, optional, backend, verbose = FALSE)

  setupSharedMemory(verbose)

  setupTempDirectory(verbose)

  setupLibPaths(pkg = "NMF")

  setupRNG(seed, n, verbose = FALSE)
}
\arguments{
  \item{spec}{target parallel specification: either
  \code{TRUE} or \code{FALSE}, or a single numeric value
  that specifies the number of cores to setup.}

  \item{optional}{a logical that indicates if the
  specification must be fully satisfied, throwing an error
  if it is not, or if one can switch back to sequential,
  only outputting a verbose message.}

  \item{backend}{value from argument \code{.pbackend} of
  \code{nmf}.}

  \item{verbose}{logical or integer level of verbosity for
  message outputs.}

  \item{pkg}{package name whose path should be exported the
  workers.}

  \item{seed}{initial RNG seed specification}

  \item{n}{number of RNG seeds to generate}
}
\description{
  Functions used internally to setup the computational
  environment.

  \code{setupBackend} sets up the foreach backend.

  \code{setupSharedMemory} checks if one can use the
  packages \emph{bigmemory} and \emph{sychronicity} to
  speed-up parallel computations when not keeping all the
  fits. When both these packages are available, only one
  result per host is written on disk, with its achieved
  deviance stored in shared memory, that is accessible to
  all cores on a same host. It returns \code{TRUE} if both
  packages are available and NMF option \code{'shared'} is
  toggled on.

  \code{setupTempDirectory} creates a temporary directory
  to store the best fits computed on each host. It ensures
  each worker process has access to it.

  \code{setupLibPaths} add the path to the NMF package to
  each workers' libPaths.

  \code{setupRNG} sets the RNG for use by the function nmf.
  It returns the old RNG as an rstream object or the result
  of set.seed if the RNG is not changed due to one of the
  following reason: - the settings are not compatible with
  rstream
}
\keyword{internals}

