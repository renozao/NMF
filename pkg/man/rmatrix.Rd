\name{rmatrix: Generating random matrices}
\alias{rmatrix}

\alias{rmatrix-methods}
\alias{rmatrix,numeric-method}
\alias{rmatrix,matrix-method}
\alias{rmatrix,NMF-method}

\title{ Creates a Random Matrix Using Any Given Distribution Function}
\description{
This function provides a short-cut to create a random matrix whose entries are 
drawn any given random distribution, as soon as this one is implemented as a random variate 
generation function similar to: \code{runif}, \code{rnorm}, etc ...  

It essentially wraps the following common call:

\code{matrix(dist(nrow*ncol, ...), nrow, ncol)}

into the following shorter call -- that should also be less prone to errors:

\code{rmatrix(nrow, ncol, dist, ...)}

}

\usage{

\S4method{rmatrix}{numeric}(x, y, dist=runif, byrow = FALSE, dimnames = NULL, ...)

\S4method{rmatrix}{matrix}(x, ...)

\S4method{rmatrix}{NMF}(x, ...)

}

\arguments{

  \item{x}{ a \code{numeric} value giving the number of rows of the result matrix,
  a matrix, or an \code{\linkS4class{NMF}} model.}
  
  \item{y}{ a \code{numeric} value giving the number of columns of the result matrix.
  If \var{y} is missing, then the created matrix is square (i.e. \code{y=x}).}
  
  \item{dist}{ a random variate generation function (e.g. \code{\link{runif}}, \code{\link{rnorm}}, etc...) 
  from which to draw the matrix entries.
  It must be a function whose first parameter is the number of values \code{n} to be drawn, 
  and return a \code{numeric} vector of length \code{n}.}
  
  \item{byrow}{logical. If \code{FALSE} (the default) the matrix is filled by
          columns, otherwise the matrix is filled by rows. See \code{\link{matrix}}.}
  
  \item{dimnames}{A \code{dimnames} attribute for the matrix: \code{NULL} or a \code{list} of
          length 2 giving the row and column names respectively.  An
          empty list is treated as \code{NULL}, and a list of length one as
          row names.  The list can be named, and the list names will be
          used as names for the dimensions. See \code{\link{matrix}}.}
  
  \item{\dots}{
  	\itemize{
	\item If \code{x} is a \code{numeric}: extra parameters passed to function \code{dist}.  	
	\item If \code{x} is a \code{matrix}: extra parameters passed to function the internal call 
	\code{rmatrix(nrow(x), ncol(x), ...)}.
  	}
  }
}

\value{ 

	a matrix whose entries are drawn from distribution \code{dist} and of 
	dimension \code{x} X \code{y} if \var{x} is a \code{numeric}, or of dimension 
	\code{nrow(x)} X \code{ncol(y)} if \var{x} is a \code{matrix} or an 
	\code{NMF} model.

}

\author{ Renaud Gaujoux  }

\seealso{ \code{\link{runif}}, \code{\link{rnorm}} or any other \code{\link[=distributions]{rxxx}} similar random variate generation function.}

\examples{

## Generate a random matrix of a given size
rmatrix(5, 3)
\dontshow{ stopifnot( identical(dim(rmatrix(5, 3)), c(5L,3L)) ) }

## Generate a random matrix of the same dimension of a template matrix
a <- matrix(1, 3, 4)
rmatrix(a)
\dontshow{ stopifnot( identical(dim(rmatrix(a)), c(3L,4L)) ) }

## Generate a random matrix of the dimension of the target matrix of an NMF model 
a <- nmfModel(2, 10, 5)
rmatrix(a)
\dontshow{ stopifnot( identical(dim(rmatrix(a)), c(10L,5L)) ) }

## Specificy the distribution to use

# the default is uniform
a <- rmatrix(1000, 50)
\dontrun{ hist(a) }

# use normal ditribution
a <- rmatrix(1000, 50, rnorm)
\dontrun{ hist(a) }

# extra arguments can be passed to the random variate generation function 
a <- rmatrix(1000, 50, rnorm, mean=2, sd=0.5)
\dontrun{ hist(a) }

}
\keyword{ distribution }
