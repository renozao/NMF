\name{nmfCheck}
\alias{nmfCheck}
\title{Checking NMF Algorithm}
\usage{
  nmfCheck(method = NULL, rank = max(ncol(x)/5, 3),
    x = NULL, seed = 1234, ...)
}
\arguments{
  \item{method}{name of the NMF algorithm to be tested.}

  \item{rank}{rank of the factorization}

  \item{x}{target data. If \code{NULL}, a random 20 x 10
  matrix is generated}

  \item{seed}{specifies a seed or seeding method for the
  computation.}

  \item{...}{other arguments passed to the call to
  \code{\link{nmf}}.}
}
\value{
  the result of the NMF fit invisibly.
}
\description{
  \code{nmfCheck} enables to quickly check that a given NMF
  algorithm runs properly, by applying it to some small
  random data.
}
\examples{
\dontshow{# roxygen generated flag
options(R_CHECK_RUNNING_EXAMPLES_=TRUE)
}

# test default algorithm
nmfCheck()

# test 'lee' algorithm
nmfCheck('lee')
}

