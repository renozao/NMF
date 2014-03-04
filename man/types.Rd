\name{is.nmf}
\alias{hasBasis}
\alias{hasCoef}
\alias{is.empty.nmf}
\alias{is.nmf}
\alias{isNMFfit}
\alias{is.partial.nmf}
\title{Testing NMF Objects}
\usage{
  is.nmf(x)

  is.empty.nmf(x, ...)

  hasBasis(x)

  hasCoef(x)

  is.partial.nmf(x)

  isNMFfit(object, recursive = TRUE)
}
\arguments{
  \item{x}{an R object. See section \emph{Details}, for how
  each function uses this argument.}

  \item{...}{extra parameters to allow extension or passed
  to subsequent calls}

  \item{object}{any R object.}

  \item{recursive}{if \code{TRUE} and \code{object} is a
  plain list then \code{isNMFfit} tests each element of the
  list. Note that the recursive test only applies in the
  case of lists that are not themselves NMFfit objects,
  like \code{NMFfitXn} objects for which the result of
  \code{isNMFfit} will always be \code{TRUE}, although they
  are list objects (a single logical value).}
}
\value{
  \code{isNMFfit} returns a \code{logical} vector (or a
  list if \code{object} is a list of list) of the same
  length as \code{object}.
}
\description{
  The functions documented here tests different
  characteristics of NMF objects.

  \code{is.nmf} tests if an object is an NMF model or a
  class that extends the class NMF.

  \code{hasBasis} tests whether an objects contains a basis
  matrix -- returned by a suitable method \code{basis} --
  with at least one row.

  \code{hasBasis} tests whether an objects contains a
  coefficient matrix -- returned by a suitable method
  \code{coef} -- with at least one column.

  \code{is.partial.nmf} tests whether an NMF model object
  contains either an empty basis or coefficient matrix. It
  is a shorcut for \code{!hasCoef(x) || !hasBasis(x)}.
}
\details{
  \code{is.nmf} tests if \code{object} is the name of a
  class (if a \code{character} string), or inherits from a
  class, that extends \code{\linkS4class{NMF}}.

  \code{is.empty.nmf} returns \code{TRUE} if the basis and
  coefficient matrices of \code{x} have respectively zero
  rows and zero columns. It returns \code{FALSE} otherwise.

  In particular, this means that an empty model can still
  have a non-zero number of basis components, i.e. a
  factorization rank that is not null. This happens, for
  example, in the case of NMF models created calling the
  factory method \code{\link{nmfModel}} with a value only
  for the factorization rank.

  \emph{isNMFfit} checks if \code{object} inherits from
  class \code{\linkS4class{NMFfit}} or
  \code{\linkS4class{NMFfitX}}, which are the two types of
  objects returned by the function \code{\link{nmf}}. If
  \code{object} is a plain \code{list} and
  \code{recursive=TRUE}, then the test is performed on each
  element of the list, and the return value is a logical
  vector (or a list if \code{object} is a list of list) of
  the same length as \code{object}.
}
\note{
  The function \code{is.nmf} does some extra work with the
  namespace as this function needs to return correct
  results even when called in \code{.onLoad}. See
  discussion on r-devel:
  \url{https://stat.ethz.ch/pipermail/r-devel/2011-June/061357.html}
}
\examples{
\dontshow{# roxygen generated flag
options(R_CHECK_RUNNING_EXAMPLES_=TRUE)
}

#----------
# is.nmf
#----------
# test if an object is an NMF model, i.e. that it implements the NMF interface
is.nmf(1:4)
is.nmf('NMFstd')
is.nmf('NMFblah')
is.nmf( nmfModel(3) )
is.nmf( nmf(rmatrix(20,10), 3) )

#----------
# is.empty.nmf
#----------
# empty model
is.empty.nmf( nmfModel(3) )
# non empty models
is.empty.nmf( nmfModel(3, 10, 0) )
is.empty.nmf( rnmf(3, 10, 5) )

#----------
# isNMFfit
#----------
## Testing results of fits
# generate a random
V <- rmatrix(20, 10)

# single run -- using very low value for maxIter to speed up the example
res <- nmf(V, 3, maxIter=3L)
isNMFfit(res)

# multiple runs - keeping single fit
resm <- nmf(V, 3, nrun=3, maxIter=3L)
isNMFfit(resm)

# multiple runs - keeping all fits
resM <- nmf(V, 3, nrun=3, .opt='k', maxIter=3L)
isNMFfit(resM)

# with a list of results
isNMFfit(list(res, resm, resM, 'not a result'))
isNMFfit(list(res, list(resm, resM), 'not a result')) # list of list
isNMFfit(list(res, resm, resM, 'not a result'), recursive=FALSE)
}
\seealso{
  \code{\linkS4class{NMFfit}},
  \code{\linkS4class{NMFfitX}},
  \code{\linkS4class{NMFfitXn}}
}

