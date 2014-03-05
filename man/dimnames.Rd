\docType{methods}
\name{basisnames}
\alias{basisnames}
\alias{basisnames<-}
\alias{basisnames<-,ANY-method}
\alias{basisnames,ANY-method}
\alias{basisnames<--methods}
\alias{basisnames-methods}
\alias{dimnames-NMF}
\alias{dimnames<-,NMF-method}
\alias{dimnames,NMF-method}
\title{Dimension names for NMF objects}
\usage{
  basisnames(x, ...)

  basisnames(x, ...)<-value

  \S4method{dimnames}{NMF}(x)

  \S4method{dimnames}{NMF}(x)<-value
}
\arguments{
  \item{x}{an object with suitable \code{basis} and
  \code{coef} methods, such as an object that inherit from
  \code{\linkS4class{NMF}}.}

  \item{...}{extra argument to allow extension.}

  \item{value}{a character vector, or \code{NULL} or, in
  the case of \code{dimnames<-}, a list 2 or 3-length list
  of character vectors.}
}
\description{
  The methods \code{dimnames}, \code{rownames},
  \code{colnames} and \code{basisnames} and their
  respective replacement form allow to get and set the
  dimension names of the matrix factors in a NMF model.

  \code{dimnames} returns all the dimension names in a
  single list. Its replacement form \code{dimnames<-}
  allows to set all dimension names at once.

  \code{rownames}, \code{colnames} and \code{basisnames}
  provide separate access to each of these dimension names
  respectively. Their respective replacement form allow to
  set each dimension names separately.
}
\details{
  The function \code{basisnames} is a new S4 generic
  defined in the package NMF, that returns the names of the
  basis components of an object. Its default method should
  work for any object, that has a suitable \code{basis}
  method defined for its class.

  The method \code{dimnames} is implemented for the base
  generic \code{\link{dimnames}}, which make the base
  function \code{\link{rownames}} and
  \code{\link{colnames}} work directly.

  Overall, these methods behave as their equivalent on
  \code{matrix} objects. The function \code{basisnames<-}
  ensures that the dimension names are handled in a
  consistent way on both factors, enforcing the names on
  both matrix factors simultaneously.

  The function \code{basisnames<-} is a new S4 generic
  defined in the package NMF, that sets the names of the
  basis components of an object. Its default method should
  work for any object, that has suitable \code{basis<-} and
  \code{coef<-} methods method defined for its class.
}
\section{Methods}{
  \describe{

  \item{basisnames}{\code{signature(x = "ANY")}: Default
  method which returns the column names of the basis matrix
  extracted from \code{x}, using the \code{basis} method.

  For NMF objects these also correspond to the row names of
  the coefficient matrix. }

  \item{basisnames<-}{\code{signature(x = "ANY")}: Default
  method which sets, respectively, the row and the column
  names of the basis matrix and coefficient matrix of
  \code{x} to \code{value}. }

  \item{dimnames}{\code{signature(x = "NMF")}: Returns the
  dimension names of the NMF model \code{x}.

  It returns either NULL if no dimnames are set on the
  object, or a 3-length list containing the row names of
  the basis matrix, the column names of the mixture
  coefficient matrix, and the column names of the basis
  matrix (i.e. the names of the basis components). }

  \item{dimnames<-}{\code{signature(x = "NMF")}: sets the
  dimension names of the NMF model \code{x}.

  \code{value} can be \code{NULL} which resets all
  dimension names, or a 1, 2 or 3-length list providing
  names at least for the rows of the basis matrix.

  The optional second element of \code{value} (NULL if
  absent) is used to set the column names of the
  coefficient matrix. The optional third element of
  \code{value} (NULL if absent) is used to set both the
  column names of the basis matrix and the row names of the
  coefficient matrix. }

  }
}
\examples{
\dontshow{# roxygen generated flag
options(R_CHECK_RUNNING_EXAMPLES_=TRUE)
}

# create a random NMF object
a <- rnmf(2, 5, 3)

# set dimensions
dims <- list( features=paste('f', 1:nrow(a), sep='')
				, samples=paste('s', 1:ncol(a), sep='')
				, basis=paste('b', 1:nbasis(a), sep='') )
dimnames(a) <- dims
dimnames(a)
basis(a)
coef(a)

# access the dimensions separately
rownames(a)
colnames(a)
basisnames(a)

# set only the first dimension (rows of basis): the other two dimnames are set to NULL
dimnames(a) <- dims[1]
dimnames(a)
basis(a)
coef(a)

# set only the two first dimensions (rows and columns of basis and coef respectively):
# the basisnames are set to NULL
dimnames(a) <- dims[1:2]
dimnames(a)
basis(a)

# reset the dimensions
dimnames(a) <- NULL
dimnames(a)
basis(a)
coef(a)

# set each dimensions separately
rownames(a) <- paste('X', 1:nrow(a), sep='') # only affect rows of basis
basis(a)

colnames(a) <- paste('Y', 1:ncol(a), sep='') # only affect columns of coef
coef(a)

basisnames(a) <- paste('Z', 1:nbasis(a), sep='') # affect both basis and coef matrices
basis(a)
coef(a)
}
\keyword{methods}

