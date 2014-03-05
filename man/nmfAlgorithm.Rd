\name{nmfAlgorithm}
\alias{nmfAlgorithm}
\title{Listing and Retrieving NMF Algorithms}
\usage{
  nmfAlgorithm(name = NULL, version = NULL, all = FALSE,
    ...)
}
\arguments{
  \item{name}{Access key. If not missing, it must be a
  single character string that is partially matched against
  the available algorithms in the registry. In this case,
  if \code{all=FALSE} (default), then the algorithm is
  returned as an \code{NMFStrategy} object that can be
  directly passed to \code{\link{nmf}}. An error is thrown
  if no matching algorithm is found.

  If missing or \code{NULL}, then access keys of algorithms
  -- that match the criteria \code{version}, are returned.
  This argument is assumed to be regular expression if
  \code{all=TRUE} or \code{version} is not \code{NULL}.}

  \item{version}{version of the algorithm(s) to retrieve.
  Currently only value \code{'R'} is supported, which
  searched for plain R implementations.}

  \item{all}{a logical that indicates if all algorithm keys
  should be returned, including the ones from alternative
  algorithm versions (e.g. plain R implementations of
  algorithms, for which a version based on optimised C
  updates is used by default).}

  \item{...}{extra arguments passed to
  \code{\link{getNMFMethod}} when \code{name} is not
  \code{NULL} and \code{all=FALSE}. It is not used
  otherwise.}
}
\value{
  an \code{\linkS4class{NMFStrategy}} object if \code{name}
  is not \code{NULL} and \code{all=FALSE}, or a named
  character vector that contains the access keys of the
  matching algorithms. The names correspond to the access
  key of the primary algorithm: e.g. algorithm \sQuote{lee}
  has two registered versions, one plain R
  (\sQuote{.R#lee}) and the other uses optimised C updates
  (\sQuote{lee}), which will all get named \sQuote{lee}.
}
\description{
  \code{nmfAlgorithm} lists access keys or retrieves NMF
  algorithms that are stored in registry. It allows to list
}
\examples{
\dontshow{# roxygen generated flag
options(R_CHECK_RUNNING_EXAMPLES_=TRUE)
}

# list all main algorithms
nmfAlgorithm()
# list all versions of algorithms
nmfAlgorithm(all=TRUE)
# list all plain R versions
nmfAlgorithm(version='R')
}
\seealso{
  Other regalgo: \code{\link{canFit}}
}

