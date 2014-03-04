\name{setNMFMethod}
\alias{nmfRegisterAlgorithm}
\alias{setNMFMethod}
\title{Registering NMF Algorithms}
\usage{
  setNMFMethod(name, method, ...,
    overwrite = isLoadingNamespace(), verbose = TRUE)

  nmfRegisterAlgorithm(name, method, ...,
    overwrite = isLoadingNamespace(), verbose = TRUE)
}
\arguments{
  \item{...}{arguments passed to the factory function
  \code{\link{NMFStrategy}}, which instantiate the
  \code{\linkS4class{NMFStrategy}} object that is stored in
  registry.}

  \item{overwrite}{logical that indicates if any existing
  NMF method with the same name should be overwritten
  (\code{TRUE}) or not (\code{FALSE}), in which case an
  error is thrown.}

  \item{verbose}{a logical that indicates if information
  about the registration should be printed (\code{TRUE}) or
  not (\code{FALSE}).}

  \item{name}{name/key of an NMF algorithm.}

  \item{method}{definition of the algorithm}
}
\description{
  Adds a new algorithm to the registry of algorithms that
  perform Nonnegative Matrix Factorization.

  \code{nmfRegisterAlgorithm} is an alias to
  \code{setNMFMethod} for backward compatibility.
}
\examples{
\dontshow{# roxygen generated flag
options(R_CHECK_RUNNING_EXAMPLES_=TRUE)
}

# define/regsiter a new -- dummy -- NMF algorithm with the minimum arguments
# y: target matrix
# x: initial NMF model (i.e. the seed)
# NB: this algorithm simply return the seed unchanged
setNMFMethod('mynmf', function(y, x, ...){ x })

# check algorithm on toy data
res <- nmfCheck('mynmf')
# the NMF seed is not changed
stopifnot( nmf.equal(res, nmfCheck('mynmf', seed=res)) )
}

