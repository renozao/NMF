\docType{methods}
\name{RNG-NMF}
\alias{getRNG}
\alias{.getRNG}
\alias{getRNG1}
\alias{getRNG1,ANY-method}
\alias{getRNG1-methods}
\alias{getRNG1,NMFfitX1-method}
\alias{getRNG1,NMFfitX-method}
\alias{getRNG1,NMFfitXn-method}
\alias{.getRNG,ANY-method}
\alias{.getRNG,list-method}
\alias{.getRNG-methods}
\alias{.getRNG,missing-method}
\alias{.getRNG,NMFfitXn-method}
\alias{.getRNG,numeric-method}
\alias{nextRNG}
\alias{rng1.equal}
\alias{RNGdigest}
\alias{rng.equal}
\alias{RNGinfo}
\alias{RNG-NMF}
\alias{RNGrecovery}
\alias{RNGtype}
\alias{setRNG}
\alias{.setRNG}
\alias{.setRNG,character-method}
\alias{.setRNG-methods}
\alias{.setRNG,numeric-method}
\title{RNG Settings}
\usage{
  RNGrecovery()

  RNGtype(object)

  RNGinfo(object = getRNG(), indent = "")

  getRNG(object, ...)

  .getRNG(object, ...)

  getRNG1(object, ...)

  nextRNG(object, ...)

  setRNG(object, ..., verbose = FALSE)

  .setRNG(object, ...)

  RNGdigest(x)

  rng.equal(x, y)

  rng1.equal(x, y)
}
\arguments{
  \item{indent}{character string to use as indentation
  prefix in the output from \code{RNGinfo}.}

  \item{object}{an R object from which RNG settings can be
  extracted, e.g. an integer vector containing a suitable
  value for \code{.Random.seed} or an object returned by
  the function \code{\link{nmf}}.}

  \item{...}{extra arguments passed to a suitable method
  \code{.getRNG}.}

  \item{object}{an R object from which RNG settings can be
  extracted, e.g. an integer vector containing a suitable
  value for \code{.Random.seed} or an object returned by
  the function \code{\link{nmf}}.}

  \item{...}{extra arguments to allow extension}

  \item{verbose}{a logical that indicates if the new RNG
  settings should be displayed.}

  \item{x}{objects from which RNG settings are extracted}

  \item{y}{object from which RNG settings are extracted}
}
\value{
  the RNG settings as a single integer vector as
  \code{\link{.Random.seed}} or \code{NULL} if no RNG data
  was found.

  the RNG settings as a single integer vector as
  \code{\link{.Random.seed}} or \code{NULL} if no RNG data
  was found.

  \code{getRNG}, \code{getRNG1}, \code{nextRNG} and
  \code{setRNG} return an integer vector of length greater
  than 3 (see \code{\link{.Random.seed}}.

  \code{setRNG} invisibly returns the old RNG settings as
  they were before changing them.

  \code{rng.equal} and \code{rng.equal1} return a
  \code{TRUE} or \code{FALSE}.
}
\description{
  The functions documented here provide a unified interface
  to work with RNG settings.

  \code{RNGrecovery} allows to recover from a broken state
  of \code{.Random.seed}.

  \code{RNGtype} returns the same type of values as
  \code{RNGkind()}, except that it can extract the RNG
  settings from an object. If \code{object} is missing it
  returns the kinds of the current RNG settings, i.e. it is
  identical to \code{RNGkind()}.

  \code{RNGinfo} shows displays human readable information
  about RNG settings. If \code{object} is missing it
  displays information about the current RNG.

  The functions documented here provide a unified interface
  to work with RNG settings.

  \code{getRNG} is an S4 generic that returns the Random
  Number Generator (RNG) settings used for computing an
  object. For example, in the case of results from multiple
  NMF runs, it returns the RNG settings used to compute the
  best fit.

  \code{getRNG1} is defined to provide separate access to
  the RNG settings as they were at the very beginning of a
  whole computation, which might differ from the RNG
  settings returned by \code{getRNG}, that allows to
  reproduce the result only.

  \code{nextRNG} returns the RNG settings as they would be
  after seeding with \code{seed}.

  \code{setRNG} tries to extract RNG settings from
  \code{object} and use a suitable \code{.setRNG} method to
  set these settings. All arguments are passed to the next
  call to \code{.setRNG}.

  \code{.setRNG} is an S4 generic that sets the current RNG
  settings, from a variety of format. Its methods define
  the workhorse functions that are called by \code{setRNG}.

  \code{RNGdigest} computes a hash from the RNG settings
  associated with an object.

  \code{rng.equal} and \code{rng1.equal} return \code{TRUE}
  the RNG settings associated with two objects are
  identical, and \code{FALSE} otherwise. The comparison is
  made between the hashes returned by \code{RNGdigest}.

  The function \code{rng1.equal} tests whether two objects
  have identical \strong{initial} RNG settings.
}
\details{
  \code{getRNG} is an S4 generic that returns the Random
  Number Generator (RNG) settings used for computing an
  object. For example, in the case of results from multiple
  NMF runs, it returns the RNG settings used to compute the
  best fit.

  Think of a sequence of separate computations, from which
  only one result is used for the result (e.g. the one that
  maximise a likelihood): \code{getRNG1} would return the
  RNG settings to reproduce the complete sequence of
  computations, while \code{getRNG} would return the RNG
  settings necessary to reproduce only the computation
  whose result has maximum likelihood.
}
\section{Methods}{
  \describe{

  \item{.getRNG}{\code{signature(object = "missing")}:
  Returns the current RNG settings. }

  \item{.getRNG}{\code{signature(object = "ANY")}: Default
  method that tries to extract RNG information from
  \code{object}, by looking sequentially to a slot named
  \code{'rng'}, a slot named \code{'rng.seed'} or an
  attribute names \code{'rng'}.

  It returns \code{NULL} if no RNG data was found. }

  \item{.getRNG}{\code{signature(object = "list")}: Method
  for S3 objects, that aims at reproducing the behaviour of
  the function \code{getRNG} of the package \code{getRNG}.

  It sequentially looks for RNG data in elements
  \code{'rng'}, \code{noise$rng} if element \code{'noise'}
  exists and is a \code{list}, or in attribute
  \code{'rng'}. }

  \item{.getRNG}{\code{signature(object = "numeric")}:
  Methods for numeric object, which returns the object
  itself, if it has more than one element, coerced into an
  integer vector if necessary, as it is assumed to already
  represent a value for \code{\link{.Random.seed}}.

  Or if \code{object} has a single element, the value of
  \code{.Random.seed} as it would be after calling
  \code{set.seed(object, ...)} In this case, all arguments
  in \code{...} are passed to \code{\link{set.seed}}. }

  \item{.getRNG}{\code{signature(object = "NMFfitXn")}:
  Returns the RNG settings used for the best fit.

  This method throws an error if the object is empty. }

  \item{getRNG1}{\code{signature(object = "ANY")}: Default
  method that is identical to \code{getRNG(object, ...)}. }

  \item{getRNG1}{\code{signature(object = "NMFfitX")}:
  Returns the RNG settings used for the first NMF run of
  multiple NMF runs. }

  \item{getRNG1}{\code{signature(object = "NMFfitX1")}:
  Returns the RNG settings used to compute the first of all
  NMF runs, amongst which \code{object} was selected as the
  best fit. }

  \item{getRNG1}{\code{signature(object = "NMFfitXn")}:
  Returns the RNG settings used for the first run.

  This method throws an error if the object is empty. }

  \item{.setRNG}{\code{signature(object = "character")}:
  Sets the RNG to kind \code{object}, assuming is a valid
  RNG kind: it is equivalent to \code{RNGkind(object, ...}.
  All arguments in \code{...} are passed to
  \code{\link{RNGkind}}. }

  \item{.setRNG}{\code{signature(object = "numeric")}: Sets
  the RNG settings using \code{object} directly the new
  value for \code{.Random.seed} or to initialise it with
  \code{\link{set.seed}}. }

  }
}
\examples{
# get current RNG settings
head(getRNG())
obj <- list(x=10, rng=123)
setRNG(obj)
rng <- getRNG()
runif(10)
set.seed(123)
rng.equal(rng)
# set RNG kind
old <- setRNG('Marsaglia')
# restore
setRNG(old)
# directly set .Random.seed
rng <- getRNG()
r <- runif(10)
setRNG(rng)
rng.equal(rng)

# initialise from a single number (<=> set.seed)
setRNG(123)
rng <- getRNG()
runif(10)
set.seed(123)
rng.equal(rng)
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

