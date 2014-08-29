\docType{methods}
\name{getRNG1}
\alias{.getRNG}
\alias{getRNG1}
\alias{getRNG1-methods}
\alias{getRNG1,NMFfitX1-method}
\alias{getRNG1,NMFfitX-method}
\alias{getRNG1,NMFfitXn-method}
\alias{.getRNG-methods}
\alias{.getRNG,NMFfitXn-method}
\title{Extracting RNG Data from NMF Objects}
\usage{
  getRNG1(object, ...)

  .getRNG(object, ...)
}
\arguments{
  \item{object}{an R object from which RNG settings can be
  extracted, e.g. an integer vector containing a suitable
  value for \code{.Random.seed} or embedded RNG data, e.g.,
  in S3/S4 slot \code{rng} or \code{rng$noise}.}

  \item{...}{extra arguments to allow extension and passed
  to a suitable S4 method \code{.getRNG} or
  \code{.setRNG}.}
}
\description{
  The \code{\link{nmf}} function returns objects that
  contain embedded RNG data, that can be used to exactly
  reproduce any computation. These data can be extracted
  using dedicated methods for the S4 generics
  \code{\link[rngtools]{getRNG}} and
  \code{\link[rngtools]{getRNG1}}.
}
\section{Methods}{
\describe{

\item{.getRNG}{\code{signature(object = "NMFfitXn")}: Returns the RNG settings used for the best fit.

This method throws an error if the object is empty.
}


\item{getRNG1}{\code{signature(object = "NMFfitX")}: Returns the RNG settings used for the first NMF run of multiple NMF runs.
}

\item{getRNG1}{\code{signature(object = "NMFfitX1")}: Returns the RNG settings used to compute the first of all NMF runs, amongst
which \code{object} was selected as the best fit.
}

\item{getRNG1}{\code{signature(object = "NMFfitXn")}: Returns the RNG settings used for the first run.

This method throws an error if the object is empty.
}

}

}
\examples{
\dontshow{# roxygen generated flag
options(R_CHECK_RUNNING_EXAMPLES_=TRUE)
}

# For multiple NMF runs, the RNG settings used for the first run is also stored
V <- rmatrix(20,10)
res <- nmf(V, 3, nrun=4)
# RNG used for the best fit
getRNG(res)
# RNG used for the first of all fits
getRNG1(res)
# they may differ if the best fit is not the first one
rng.equal(res, getRNG1(res))
}
\keyword{methods}

