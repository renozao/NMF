\docType{methods}
\name{nmfModel}
\alias{nmfModel}
\alias{nmfModel,matrix,ANY-method}
\alias{nmfModel,matrix,matrix-method}
\alias{nmfModel-methods}
\alias{nmfModel,missing,ANY-method}
\alias{nmfModel,missing,missing-method}
\alias{nmfModel,NULL,ANY-method}
\alias{nmfModel,numeric,matrix-method}
\alias{nmfModel,numeric,missing-method}
\alias{nmfModel,numeric,numeric-method}
\alias{nmfModels}
\title{Factory Methods NMF Models}
\usage{
  nmfModel(rank, target = 0L, ...)

  \S4method{nmfModel}{numeric,numeric}(rank, target,
    ncol = NULL, model = "NMFstd", W, H, ...,
    force.dim = TRUE, order.basis = TRUE)

  \S4method{nmfModel}{numeric,matrix}(rank, target, ...,
    use.names = TRUE)

  nmfModels(builtin.only = FALSE)
}
\arguments{
  \item{rank}{specification of the target factorization
  rank (i.e. the number of components).}

  \item{target}{an object that specifies the dimension of
  the estimated target matrix.}

  \item{...}{extra arguments to allow extension, that are
  passed down to the workhorse method
  \code{nmfModel,numeric.numeric}, where they are used to
  initialise slots specific to the instantiating NMF model
  class.}

  \item{ncol}{a numeric value that specifies the number of
  columns of the target matrix, fitted the NMF model. It is
  used only if not missing and when argument \code{target}
  is a single numeric value.}

  \item{model}{the class of the object to be created. It
  must be a valid class name that inherits from class
  \code{NMF}. Default is the standard NMF model
  \code{\linkS4class{NMFstd}}.}

  \item{W}{value for the basis matrix}

  \item{H}{value for the mixture coefficient matrix}

  \item{force.dim}{logical that indicates whether the
  method should try lowering the rank or shrinking
  dimensions of the input matrices to make them compatible}

  \item{order.basis}{logical that indicates whether the
  basis components should reorder the rows of the mixture
  coefficient matrix to match the order of the basis
  components, based on their respective names. It is only
  used if the basis and coefficient matrices have common
  unique column and row names respectively.}

  \item{use.names}{a logical that indicates whether the
  dimension names of the target matrix should be set on the
  returned NMF model.}

  \item{builtin.only}{logical that indicates whether only
  built-in NMF models, i.e. defined within the NMF package,
  should be listed.}
}
\value{
  an object that inherits from class
  \code{\linkS4class{NMF}}.

  a list
}
\description{
  \code{nmfModel} is a S4 generic function which provides a
  convenient way to build NMF models. It implements a
  unified interface for creating \code{NMF} objects from
  any NMF models, which is designed to resolve potential
  dimensions inconsistencies.

  \code{nmfModels} lists all available NMF models currently
  defined that can be used to create NMF objects, i.e. --
  more or less -- all S4 classes that inherit from class
  \code{\linkS4class{NMF}}.
}
\details{
  All \code{nmfModel} methods return an object that
  inherits from class \code{NMF}, that is suitable for
  seeding NMF algorithms via arguments \code{rank} or
  \code{seed} of the \code{\link{nmf}} method, in which
  case the factorisation rank is implicitly set by the
  number of basis components in the seeding model (see
  \code{\link{nmf}}).
}
\section{Methods}{
  \describe{

  \item{nmfModel}{\code{signature(rank = "numeric", target
  = "numeric")}: Main factory method for NMF models

  This method is the workhorse method that is eventually
  called by all other methods. See section \emph{Main
  factory method} for more details. }

  \item{nmfModel}{\code{signature(rank = "numeric", target
  = "missing")}: Creates an empty NMF model of a given
  rank.

  This call is equivalent to \code{nmfModel(rank, 0L,
  ...)}, which creates \emph{empty} \code{NMF} object with
  a basis and mixture coefficient matrix of dimension 0 x
  \code{rank} and \code{rank} x 0 respectively. }

  \item{nmfModel}{\code{signature(rank = "missing", target
  = "ANY")}: Creates an empty NMF model of null rank and a
  given dimension.

  This call is equivalent to \code{nmfModel(0, target,
  ...)}. }

  \item{nmfModel}{\code{signature(rank = "NULL", target =
  "ANY")}: Creates an empty NMF model of null rank and
  given dimension.

  This call is equivalent to \code{nmfModel(0, target,
  ...)}, and is meant for internal usage only. }

  \item{nmfModel}{\code{signature(rank = "missing", target
  = "missing")}: Creates an empty NMF model or from
  existing factors

  This method is equivalent to \code{nmfModel(0, 0, ...,
  force.dim=FALSE)}. This means that the dimensions of the
  NMF model will be taken from the optional basis and
  mixture coefficient arguments \code{W} and \code{H}. An
  error is thrown if their dimensions are not compatible.

  Hence, this method may be used to generate an NMF model
  from existing factor matrices, by providing the named
  arguments \code{W} and/or \code{H}:

  \code{nmfModel(W=w)} or \code{nmfModel(H=h)} or
  \code{nmfModel(W=w, H=h)}

  Note that this may be achieved using the more convenient
  interface is provided by the method
  \code{nmfModel,matrix,matrix} (see its dedicated
  description).

  See the description of the appropriate method below. }

  \item{nmfModel}{\code{signature(rank = "numeric", target
  = "matrix")}: Creates an NMF model compatible with a
  target matrix.

  This call is equivalent to \code{nmfModel(rank,
  dim(target), ...)}. That is that the returned NMF object
  fits a target matrix of the same dimension as
  \code{target}.

  Only the dimensions of \code{target} are used to
  construct the \code{NMF} object. The matrix slots are
  filled with \code{NA} values if these are not specified
  in arguments \code{W} and/or \code{H}. However, dimension
  names are set on the return NMF model if present in
  \code{target} and argument \code{use.names=TRUE}. }

  \item{nmfModel}{\code{signature(rank = "matrix", target =
  "matrix")}: Creates an NMF model based on two existing
  factors.

  This method is equivalent to \code{nmfModel(0, 0, W=rank,
  H=target..., force.dim=FALSE)}. This allows for a natural
  shortcut for wrapping existing \strong{compatible}
  matrices into NMF models: \samp{nmfModel(w, h)}

  Note that an error is thrown if their dimensions are not
  compatible. }

  \item{nmfModel}{\code{signature(rank = "matrix", target =
  "ANY")}: Creates an NMF model with arguments \code{rank}
  and \code{target} swapped.

  This call is equivalent to \code{nmfModel(rank=target,
  target=rank, ...)}. This allows to call the
  \code{nmfModel} function with arguments \code{rank} and
  \code{target} swapped. It exists for convenience:
  \itemize{ \item allows typing \code{nmfModel(V)} instead
  of \code{nmfModel(target=V)} to create a model compatible
  with a given matrix \code{V} (i.e. of dimension
  \code{nrow(V), 0, ncol(V)}) \item one can pass the
  arguments in any order (the one that comes to the user's
  mind first) and it still works as expected. } }

  }
}
\section{Main factory method}{
  The main factory engine of NMF models is implemented by
  the method with signature \code{numeric, numeric}. Other
  factory methods provide convenient ways of creating NMF
  models from e.g. a given target matrix or known
  basis/coef matrices (see section \emph{Other Factory
  Methods}).

  This method creates an object of class \code{model},
  using the extra arguments in \code{...} to initialise
  slots that are specific to the given model.

  All NMF models implement get/set methods to access the
  matrix factors (see \code{\link{basis}}), which are
  called to initialise them from arguments \code{W} and
  \code{H}. These argument names derive from the definition
  of all built-in models that inherit derive from class
  \code{\linkS4class{NMFstd}}, which has two slots, \var{W}
  and \var{H}, to hold the two factors -- following the
  notations used in \cite{Lee et al. (1999)}.

  If argument \code{target} is missing, the method creates
  a standard NMF model of dimension 0x\code{rank}x0. That
  is that the basis and mixture coefficient matrices,
  \var{W} and \var{H}, have dimension 0x\code{rank} and
  \code{rank}x0 respectively.

  If target dimensions are also provided in argument
  \code{target} as a 2-length vector, then the method
  creates an \code{NMF} object compatible to fit a target
  matrix of dimension \code{target[1]}x\code{target[2]}.
  That is that the basis and mixture coefficient matrices,
  \var{W} and \var{H}, have dimension
  \code{target[1]}x\code{rank} and
  \code{rank}x\code{target[2]} respectively. The target
  dimensions can also be specified using both arguments
  \code{target} and \code{ncol} to define the number of
  rows and the number of columns of the target matrix
  respectively. If no other argument is provided, these
  matrices are filled with NAs.

  If arguments \code{W} and/or \code{H} are provided, the
  method creates a NMF model where the basis and mixture
  coefficient matrices, \var{W} and \var{H}, are
  initialised using the values of \code{W} and/or \code{H}.

  The dimensions given by \code{target}, \code{W} and
  \code{H}, must be compatible. However if
  \code{force.dim=TRUE}, the method will reduce the
  dimensions to the achieve dimension compatibility
  whenever possible.

  When \code{W} and \code{H} are both provided, the
  \code{NMF} object created is suitable to seed a NMF
  algorithm in a call to the \code{\link{nmf}} method. Note
  that in this case the factorisation rank is implicitly
  set by the number of basis components in the seed.
}
\examples{
# data
n <- 20; r <- 3; p <- 10
V <- rmatrix(n, p) # some target matrix

# create a r-ranked NMF model with a given target dimensions n x p as a 2-length vector
nmfModel(r, c(n,p)) # directly
nmfModel(r, dim(V)) # or from an existing matrix <=> nmfModel(r, V)
# or alternatively passing each dimension separately
nmfModel(r, n, p)

# trying to create a NMF object based on incompatible matrices generates an error
w <- rmatrix(n, r)
h <- rmatrix(r+1, p)
try( new('NMFstd', W=w, H=h) )
try( nmfModel(w, h) )
try( nmfModel(r+1, W=w, H=h) )
# The factory method can be force the model to match some target dimensions
# but warnings are thrown
nmfModel(r, W=w, H=h)
nmfModel(r, n-1, W=w, H=h)
## Empty model of given rank
nmfModel(3)
nmfModel(target=10) #square
nmfModel(target=c(10, 5))
# Build an empty NMF model
nmfModel()

# create a NMF object based on one random matrix: the missing matrix is deduced
# Note this only works when using factory method NMF
n <- 50; r <- 3;
w <- rmatrix(n, r)
nmfModel(W=w)

# create a NMF object based on random (compatible) matrices
p <- 20
h <- rmatrix(r, p)
nmfModel(H=h)

# specifies two compatible matrices
nmfModel(W=w, H=h)
# error if not compatible
try( nmfModel(W=w, H=h[-1,]) )
# create a r-ranked NMF model compatible with a given target matrix
obj <- nmfModel(r, V)
all(is.na(basis(obj)))
## From two existing factors

# allows a convenient call without argument names
w <- rmatrix(n, 3); h <- rmatrix(3, p)
nmfModel(w, h)

# Specify the type of NMF model (e.g. 'NMFns' for non-smooth NMF)
mod <- nmfModel(w, h, model='NMFns')
mod

# One can use such an NMF model as a seed when fitting a target matrix with nmf()
V <- rmatrix(mod)
res <- nmf(V, mod)
nmf.equal(res, nmf(V, mod))

# NB: when called only with such a seed, the rank and the NMF algorithm
# are selected based on the input NMF model.
# e.g. here rank was 3 and the algorithm "nsNMF" is used, because it is the default
# algorithm to fit "NMFns" models (See ?nmf).
## swapped arguments `rank` and `target`
V <- rmatrix(20, 10)
nmfModel(V) # equivalent to nmfModel(target=V)
nmfModel(V, 3) # equivalent to nmfModel(3, V)
# show all the NMF models available (i.e. the classes that inherit from class NMF)
nmfModels()
# show all the built-in NMF models available
nmfModels(builtin.only=TRUE)
}
\references{
  Lee DD and Seung HS (1999). "Learning the parts of
  objects by non-negative matrix factorization." _Nature_,
  *401*(6755), pp. 788-91. ISSN 0028-0836, <URL:
  http://dx.doi.org/10.1038/44565>, <URL:
  http://www.ncbi.nlm.nih.gov/pubmed/10548103>.
}
\seealso{
  \code{\link{is.empty.nmf}}

  Other NMF-interface: \code{\link{basis}},
  \code{\link{basis<-}}, \code{\link{coef}},
  \code{\link{coef<-}}, \code{\link{coefficients}},
  \code{\link{loadings,NMF-method}}, \code{\link{misc}},
  \code{\link{NMF-class}}, \code{\link{$<-,NMF-method}},
  \code{\link{$,NMF-method}}, \code{\link{rnmf}},
  \code{\link{scoef}}
}
\keyword{methods}

