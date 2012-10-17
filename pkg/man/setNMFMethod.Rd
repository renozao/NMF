\name{setNMFMethod}
\alias{nmfRegisterAlgorithm}
\alias{setNMFMethod}
\title{Registering NMF Algorithms}
\usage{
  setNMFMethod(name, method, ..., overwrite = FALSE,
    verbose = nmf.getOption("verbose"))

  nmfRegisterAlgorithm(name, method, ...,
    overwrite = FALSE, verbose = nmf.getOption("verbose"))
}
\arguments{
  \item{overwrite}{logical that indicates if any existing
  NMF method with the same name should be overwritten
  (\code{TRUE}) or not (\code{FALSE}), in which case an
  error is thrown.}

  \item{verbose}{a logical that indicates if information
  about the registration should be printed (\code{TRUE}) or
  not (\code{FALSE}).}

  \item{...}{arguments passed to the factory function
  \code{\link{NMFStrategy}}, which instantiate the
  \code{\linkS4class{NMFStrategy}} object that is stored in
  registry.}

  \item{name}{name/key of an NMF algorithm.}

  \item{method}{definition of the algorithm}
}
\description{
  Adds a new algorithm to the registry of algorithms that
  perform Nonnegative Matrix Factorization.

  \code{nmfRegisterAlgorithm} is an alias to
  \code{setNMFMethod} for backward compatibility.
}

