\docType{methods}
\name{rnmf}
\alias{rnmf}
\alias{rnmf,ANY,data.frame-method}
\alias{rnmf,ANY,matrix-method}
\alias{rnmf,NMF,missing-method}
\alias{rnmf,NMF,numeric-method}
\alias{rnmf,NMFOffset,numeric-method}
\alias{rnmf,formula,ANY-method}
\alias{rnmf,missing,missing-method}
\alias{rnmf,numeric,missing-method}
\alias{rnmf,numeric,numeric-method}
\alias{rnmf-methods}
\title{Generating Random NMF Models}
\usage{
  rnmf(x, target, ...)

  \S4method{rnmf}{NMF,numeric}(x, target, ncol = NULL,
    keep.names = TRUE, dist = runif)

  \S4method{rnmf}{ANY,matrix}(x, target, ...,
    dist = list(max = max(max(target, na.rm = TRUE), 1)),
    use.dimnames = TRUE)

  \S4method{rnmf}{numeric,missing}(x, target, ..., W, H,
    dist = runif)

  \S4method{rnmf}{missing,missing}(x, target, ..., W, H)

  \S4method{rnmf}{numeric,numeric}(x, target, ncol = NULL,
    ..., dist = runif)

  \S4method{rnmf}{formula,ANY}(x, target, ...,
    dist = runif)
}
\arguments{
  \item{x}{an object that determines the rank, dimension
  and/or class of the generated NMF model, e.g. a numeric
  value or an object that inherits from class
  \code{\linkS4class{NMF}}. See the description of the
  specific methods for more details on the supported
  types.}

  \item{target}{optional specification of target
  dimensions. See section \emph{Methods} for how this
  parameter is used by the different methods.}

  \item{...}{extra arguments to allow extensions and passed
  to the next method eventually down to
  \code{\link{nmfModel}}, where they are used to initialise
  slots that are specific to the instantiating NMF model.}

  \item{ncol}{single numeric value that specifies the
  number of columns of the coefficient matrix. Only used
  when \code{target} is a single numeric value.}

  \item{keep.names}{a logical that indicates if the
  dimension names of the original NMF object \code{x}
  should be conserved (\code{TRUE}) or discarded
  (\code{FALSE}).}

  \item{dist}{specification of the random distribution to
  use to draw the entries of the basis and coefficient
  matrices. It may be specified as: \itemize{

  \item a \code{function} which must be a distribution
  function such as e.g. \code{\link{runif}} that is used to
  draw the entries of both the basis and coefficient
  matrices. It is passed in the \code{dist} argument of
  \code{\link{rmatrix}}.

  \item a \code{list} of arguments that are passed
  internally to \code{\link{rmatrix}}, via
  \code{do.call('rmatrix', dist)}.

  \item a \code{character} string that is partially matched
  to \sQuote{basis} or \sQuote{coef}, that specifies which
  matrix in should be drawn randomly, the other remaining
  as in \code{x} -- unchanged.

  \item a \code{list} with elements \sQuote{basis} and/or
  \sQuote{coef}, which specify the \code{dist} argument
  separately for the basis and coefficient matrix
  respectively.

  These elements may be either a distribution function, or
  a list of arguments that are passed internally to
  \code{\link{rmatrix}}, via \code{do.call('rmatrix',
  dist$basis)} or \code{do.call('rmatrix', dist$coef)}. }}

  \item{use.dimnames}{a logical that indicates whether the
  dimnames of the target matrix should be set on the
  returned NMF model.}

  \item{W}{value for the basis matrix. \code{data.frame}
  objects are converted into matrices with
  \code{\link{as.matrix}}.}

  \item{H}{value for the mixture coefficient matrix
  \code{data.frame} objects are converted into matrices
  with \code{\link{as.matrix}}.}
}
\value{
  An NMF model, i.e. an object that inherits from class
  \code{\linkS4class{NMF}}.
}
\description{
  Generates NMF models with random values drawn from a
  uniform distribution. It returns an NMF model with basis
  and mixture coefficient matrices filled with random
  values. The main purpose of the function \code{rnmf} is
  to provide a common interface to generate random seeds
  used by the \code{\link{nmf}} function.
}
\details{
  If necessary, extensions of the standard NMF model or
  custom models must define a method
  "rnmf,<NMF.MODEL.CLASS>,numeric" for initialising their
  specific slots other than the basis and mixture
  coefficient matrices. In order to benefit from the
  complete built-in interface, the overloading methods
  should call the generic version using function
  \code{\link{callNextMethod}}, prior to set the values of
  the specific slots. See for example the method
  \code{\link[=rnmf,NMFOffset,numeric-method]{rnmf}}
  defined for \code{\linkS4class{NMFOffset}} models:
  \code{showMethods(rnmf, class='NMFOffset',
  include=TRUE))}.

  For convenience, shortcut methods for working on
  \code{data.frame} objects directly are implemented.
  However, note that conversion of a \code{data.frame} into
  a \code{matrix} object may take some non-negligible time,
  for large datasets. If using this method or other
  NMF-related methods several times, consider converting
  your data \code{data.frame} object into a matrix once for
  good, when first loaded.
}
\section{Methods}{
  \describe{

  \item{rnmf}{\code{signature(x = "NMFOffset", target =
  "numeric")}: Generates a random NMF model with offset,
  from class \code{NMFOffset}.

  The offset values are drawn from a uniform distribution
  between 0 and the maximum entry of the basis and
  coefficient matrices, which are drawn by the next
  suitable \code{\link{rnmf}} method, which is the
  workhorse method \code{rnmf,NMF,numeric}. }

  \item{rnmf}{\code{signature(x = "NMF", target =
  "numeric")}: Generates a random NMF model of the same
  class and rank as another NMF model.

  This is the workhorse method that is eventually called by
  all other methods. It generates an NMF model of the same
  class and rank as \code{x}, compatible with the
  dimensions specified in \code{target}, that can be a
  single or 2-length numeric vector, to specify a square or
  rectangular target matrix respectively.

  The second dimension can also be passed via argument
  \code{ncol}, so that calling \code{rnmf(x, 20, 10, ...)}
  is equivalent to \code{rnmf(x, c(20, 10), ...)}, but
  easier to write.

  The entries are uniformly drawn between \code{0} and
  \code{max} (optionally specified in \code{...}) that
  defaults to 1.

  By default the dimnames of \code{x} are set on the
  returned NMF model. This behaviour is disabled with
  argument \code{keep.names=FALSE}. See
  \code{\link{nmfModel}}. }

  \item{rnmf}{\code{signature(x = "ANY", target =
  "matrix")}: Generates a random NMF model compatible and
  consistent with a target matrix.

  The entries are uniformly drawn between \code{0} and
  \code{max(target)}. It is more or less a shortcut for:
  \samp{ rnmf(x, dim(target), max=max(target), ...)}

  It returns an NMF model of the same class as \code{x}. }

  \item{rnmf}{\code{signature(x = "ANY", target =
  "data.frame")}: Shortcut for \code{rnmf(x,
  as.matrix(target))}. }

  \item{rnmf}{\code{signature(x = "NMF", target =
  "missing")}: Generates a random NMF model of the same
  dimension as another NMF model.

  It is a shortcut for \code{rnmf(x, nrow(x), ncol(x),
  ...)}, which returns a random NMF model of the same class
  and dimensions as \code{x}. }

  \item{rnmf}{\code{signature(x = "numeric", target =
  "missing")}: Generates a random NMF model of a given
  rank, with known basis and/or coefficient matrices.

  This methods allow to easily generate partially random
  NMF model, where one or both factors are known. Although
  the later case might seems strange, it makes sense for
  NMF models that have fit extra data, other than the basis
  and coefficient matrices, that are drawn by an
  \code{rnmf} method defined for their own class, which
  should internally call \code{rnmf,NMF,numeric} and let it
  draw the basis and coefficient matrices. (e.g. see
  \code{\linkS4class{NMFOffset}} and
  \code{\link{rnmf,NMFOffset,numeric-method}}).

  Depending on whether arguments \code{W} and/or \code{H}
  are missing, this method interprets \code{x} differently:
  \itemize{

  \item \code{W} provided, \code{H} missing: \code{x} is
  taken as the number of columns that must be drawn to
  build a random coefficient matrix (i.e. the number of
  columns in the target matrix).

  \item \code{W} is missing, \code{H} is provided: \code{x}
  is taken as the number of rows that must be drawn to
  build a random basis matrix (i.e. the number of rows in
  the target matrix).

  \item both \code{W} and \code{H} are provided: \code{x}
  is taken as the target rank of the model to generate.
  \item Having both \code{W} and \code{H} missing produces
  an error, as the dimension of the model cannot be
  determined in this case. }

  The matrices \code{W} and \code{H} are reduced if
  necessary and possible to be consistent with this value
  of the rank, by the internal call to
  \code{\link{nmfModel}}.

  All arguments in \code{...} are passed to the function
  \code{\link{nmfModel}} which is used to build an initial
  NMF model, that is in turn passed to
  \code{rnmf,NMF,numeric} with \code{dist=list(coef=dist)}
  or \code{dist=list(basis=dist)} when suitable. The type
  of NMF model to generate can therefore be specified in
  argument \code{model} (see \code{\link{nmfModel}} for
  other possible arguments).

  The returned NMF model, has a basis matrix equal to
  \code{W} (if not missing) and a coefficient matrix equal
  to \code{H} (if not missing), or drawn according to the
  specification provided in argument \code{dist} (see
  method \code{rnmf,NMF,numeric} for details on the
  supported values for \code{dist}). }

  \item{rnmf}{\code{signature(x = "missing", target =
  "missing")}: Generates a random NMF model with known
  basis and coefficient matrices.

  This method is a shortcut for calling
  \code{rnmf,numeric,missing} with a suitable value for
  \code{x} (the rank), when both factors are known:
  code{rnmf(min(ncol(W), nrow(H)), ..., W=W, H=H)}.

  Arguments \code{W} and \code{H} are required. Note that
  calling this method only makes sense for NMF models that
  contains data to fit other than the basis and coefficient
  matrices, e.g. \code{\linkS4class{NMFOffset}}. }

  \item{rnmf}{\code{signature(x = "numeric", target =
  "numeric")}: Generates a random standard NMF model of
  given dimensions.

  This is a shortcut for \code{rnmf(nmfModel(x, target,
  ncol, ...)), dist=dist)}. It generates a standard NMF
  model compatible with the dimensions passed in
  \code{target}, that can be a single or 2-length numeric
  vector, to specify a square or rectangular target matrix
  respectively. See \code{\link{nmfModel}}. }

  \item{rnmf}{\code{signature(x = "formula", target =
  "ANY")}: Generate a random formula-based NMF model, using
  the method \code{\link{nmfModel,formula,ANY-method}}. }

  }
}
\examples{
# random NMF model with offset
x <- rnmf(2, 3, model='NMFOffset')
x
offset(x)
# from a matrix
x <- rnmf(2, rmatrix(5,3, max=10), model='NMFOffset')
offset(x)
## random NMF of same class and rank as another model

x <- nmfModel(3, 10, 5)
x
rnmf(x, 20) # square
rnmf(x, 20, 13)
rnmf(x, c(20, 13))

# using another distribution
rnmf(x, 20, dist=rnorm)

# other than standard model
y <- rnmf(3, 50, 10, model='NMFns')
y
\dontshow{ stopifnot( identical(dim(y), c(50L,10L,3L)) ) }
\dontshow{ stopifnot( is(y, 'NMFns') ) }
# random NMF compatible with a target matrix
x <- nmfModel(3, 10, 5)
y <- rmatrix(20, 13)
rnmf(x, y) # rank of x
rnmf(2, y) # rank 2
## random NMF from another model

a <- nmfModel(3, 100, 20)
b <- rnmf(a)
\dontshow{ stopifnot( !nmf.equal(a,b) ) }
# random NMF model with known basis matrix
x <- rnmf(5, W=matrix(1:18, 6)) # 6 x 5 model with rank=3
basis(x) # fixed
coef(x) # random

# random NMF model with known coefficient matrix
x <- rnmf(5, H=matrix(1:18, 3)) # 5 x 6 model with rank=3
basis(x) # random
coef(x) # fixed

# random model other than standard NMF
x <- rnmf(5, H=matrix(1:18, 3), model='NMFOffset')
basis(x) # random
coef(x) # fixed
offset(x) # random
# random model other than standard NMF
x <- rnmf(W=matrix(1:18, 6), H=matrix(21:38, 3), model='NMFOffset')
basis(x) # fixed
coef(x) # fixed
offset(x) # random
## random standard NMF of given dimensions

# generate a random NMF model with rank 3 that fits a 100x20 matrix
rnmf(3, 100, 20)
\dontshow{ stopifnot( identical(dim(rnmf(3, 100, 20)), c(100L,20L,3L)) ) }
# generate a random NMF model with rank 3 that fits a 100x100 matrix
rnmf(3, 100)
\dontshow{ stopifnot( identical(dim(rnmf(3, 100)), c(100L,100L,3L)) ) }
}
\seealso{
  \code{\link{rmatrix}}

  Other NMF-interface: \code{\link{$,NMF-method}},
  \code{\link{$<-,NMF-method}},
  \code{\link{.DollarNames,NMF-method}},
  \code{\link{.basis}}, \code{\link{.basis<-}},
  \code{\link{.coef}}, \code{\link{.coef<-}},
  \code{\link{NMF-class}}, \code{\link{basis}},
  \code{\link{basis<-}}, \code{\link{coef}},
  \code{\link{coef<-}}, \code{\link{coefficients}},
  \code{\link{loadings,NMF-method}}, \code{\link{misc}},
  \code{\link{nmfModel}}, \code{\link{nmfModels}},
  \code{\link{scoef}}
}
\keyword{methods}

