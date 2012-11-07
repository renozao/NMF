\name{dimnames-methods}
\docType{methods}
\alias{dimnames-NMF}

\alias{rownames-NMF}
\alias{colnames-NMF}

\alias{basisnames}
\alias{basisnames,NMF-method}
\alias{basisnames<-}
\alias{basisnames<-,NMF-method}
\alias{basisnames<-.NMF}

\alias{dimnames,NMF-method}
\alias{dimnames<-,NMF-method}
\alias{dimnames<-.NMF}

\title{ Dimension names for NMF objects }

\description{
The methods \code{dimnames}, \code{rownames}, \code{colnames} and \code{basisnames} 
and their respective replacement form allow to get and set the dimension names of 
the matrix factors in a NMF model.

They behave as their equivalent on \code{matrix} objects, and ensure that the dimension names 
are handled in a consistent way on both factors -- especially \code{basisnames<-} which affects 
both matrix factors simultaneously.

The methods \code{dimnames} and \code{basisnames} are implemented as S4 methods, while the 
methods \code{\link{rownames}} and \code{\link{colnames}} are the default ones that make use of the result from 
\code{dimnames}.
}

\usage{

\S4method{basisnames}{NMF}(x)
\method{basisnames}{NMF}(x, \dots) <- value

\S4method{dimnames}{NMF}(x)
\method{dimnames}{NMF}(x) <- value

## This allows to call: 
# rownames(x)
# rownames(x) <- value
# and:
# colnames(x)
# colnames(x) <- value


}
\arguments{

	\item{x}{ an object of class \code{\linkS4class{NMF}}.}
	
	\item{value}{ a character vector, or \code{NULL} or, in the case of \code{dimnames<-}, 
	a list 2 or 3-length list of character vectors.
	See section \emph{Details} for more details.}
	
	\item{...}{ extra argument to pass to internal methods. Not used.}
}

\details{
\describe{

    \item{basisnames, basisnames<-}{: 
    returns (resp. simultaneously sets) the names of the columns of the matrix of basis vectors 
    and the rows of the mixture coefficient matrix.    
    }
    
    \item{colnames, colnames<-}{: 
    returns/sets the names of the columns of the mixture coefficient matrix.
    Note that the standard arguments \code{do.NULL}, \code{prefix} as described in \code{\link{colnames}} 
    should not be used (it will always return a character vector of length 1, 
    which is likely to be incorrect).    
    }
    
    \item{rownames, rownames<-}{: 
    returns/sets the names of the columns of the basis vector matrix.
    Note that the standard arguments \code{do.NULL}, \code{prefix} as described in \code{\link{rownames}} 
    must not be used (it will always return a character vector of length 1, 
    which is likely to be incorrect).    
    }
    
    \item{dimnames}{\code{signature(x = "NMF")}: returns the dimension names of the 
    NMF model \code{x}.
    It returns either NULL if no dimnames are set on the object, or a 3-length 
    list containing the row names of the basis matrix, the column names of the 
    mixture coefficient matrix, and the column names of the basis matrix 
    (i.e. the basis vector names). 
    }
    
    \item{dimnames<-}{\code{signature(x = "NMF", value)}: sets the dimension names of the 
    NMF model \code{x}. 
    \code{value} can be \code{NULL} which resets all dimension names, 
    or a 1, 2 or 3-length list providing names at least for the rows of the 
    basis vector matrix. 
    The optional second element of \code{value} (NULL if absent) is used to 
    set the column names of the mixture coefficient matrix.
    The optional third element of \code{value} (NULL if absent) is used to set 
    both the column names of the basis vector matrix and the row names of the 
    mixture coefficient matrix.
    }
    
    
}}
\examples{

# create a random NMF object
a <- rnmf(2, 5, 3)

# set dimensions
dims <- list( features=paste('f', 1:nrow(a), sep=''), samples=paste('s', 1:ncol(a), sep=''), basis=paste('b', 1:nbasis(a), sep='') )
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
