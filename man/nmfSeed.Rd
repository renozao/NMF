\name{nmfSeed}
\alias{existsNMFSeed}
\alias{getNMFSeed}
\alias{nmfSeed}
\title{Seeding Strategies for NMF Algorithms}
\usage{
  nmfSeed(name = NULL, ...)

  getNMFSeed(name = NULL, ...)

  existsNMFSeed(name, exact = TRUE)
}
\arguments{
  \item{name}{access key of a seeding method stored in
  registry. If missing, \code{nmfSeed} returns the list of
  all available seeding methods.}

  \item{...}{extra arguments used for internal calls}

  \item{exact}{a logical that indicates if the access key
  should be matched exactly or partially.}
}
\description{
  \code{nmfSeed} lists and retrieves NMF seeding methods.

  \code{getNMFSeed} is an alias for \code{nmfSeed}.

  \code{existsNMFSeed} tells if a given seeding method
  exists in the registry.
}
\details{
  Currently the internal registry contains the following
  seeding methods, which may be specified to the function
  \code{\link{nmf}} via its argument \code{seed} using
  their access keys:

  \describe{ \item{random}{ The entries of each factors are
  drawn from a uniform distribution over \eqn{[0, max(x)]},
  where $x$ is the target matrix.} \item{nndsvd}{
  Nonnegative Double Singular Value Decomposition.

  The basic algorithm contains no randomization and is
  based on two SVD processes, one approximating the data
  matrix, the other approximating positive sections of the
  resulting partial SVD factors utilising an algebraic
  property of unit rank matrices.

  It is well suited to initialise NMF algorithms with
  sparse factors. Simple practical variants of the
  algorithm allows to generate dense factors.

  \strong{Reference:} \cite{Boutsidis et al. (2008)}}
  \item{ica}{ Uses the result of an Independent Component
  Analysis (ICA) (from the \code{fastICA} package). Only
  the positive part of the result are used to initialise
  the factors.} \item{none}{ Fixed seed.

  This method allows the user to manually provide initial
  values for both matrix factors.} }
}
\examples{
\dontshow{# roxygen generated flag
options(R_CHECK_RUNNING_EXAMPLES_=TRUE)
}

# list all registered seeding methods
nmfSeed()
# retrieve one of the methods
nmfSeed('ica')
}
\references{
  Boutsidis C and Gallopoulos E (2008). "SVD based
  initialization: A head start for nonnegative matrix
  factorization." _Pattern Recognition_, *41*(4), pp.
  1350-1362. ISSN 00313203, <URL:
  http://dx.doi.org/10.1016/j.patcog.2007.09.010>, <URL:
  http://linkinghub.elsevier.com/retrieve/pii/S0031320307004359>.
}

