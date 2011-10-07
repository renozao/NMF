\name{dim-methods: Dimensions for NMF objects}
\docType{methods}
\alias{nbasis}
\alias{nbasis,matrix-method}
\alias{nbasis,NMF-method}
\alias{nbasis,NMFfitXn-method}
\alias{dim,NMF-method}
\alias{dim,NMFfitXn-method}

\title{ Dimension for NMF objects }

\description{
The methods \code{dim}, \code{nrow}, \code{ncol} and \code{nbasis} return the 
different dimensions of an NMF model.

The package NMF defines a S4 method \code{dim} for objects of class 
\code{\linkS4class{NMF}} (resp. \code{\linkS4class{NMFfitXn}}).
This allows the base methods \code{\link{nrow}} and \code{\link{ncol}} to be 
used, to get the number of rows and columns of the target matrix estimated by an 
NMF model (resp. fit).

The method \code{nbasis} is a new S4 generic that returns the third element of 
the vector returned \code{dim}.
}

\usage{

\S4method{dim}{NMF}(x)
## This allows to call: 
# nrow(x)
# ncol(x)

\S4method{nbasis}{NMF}(x)

}
\arguments{

	\item{x}{ an object of class \code{\linkS4class{NMF}} 
	or \code{\linkS4class{NMFfitXn}}.}
	
}

\details{

\describe{

    \item{dim}{
    returns a 3-length vector containing the dimensions of the target matrix 
    fitted by \code{x} together with the factorization rank of \code{x}. 
    e.g. It will return \code{c(2000, 30, 3)} for an \code{NMF} model that fits 
    a 2000 x 30 target matrix using 3 basis vectors.
    
    For \code{\linkS4class{NMFfitXn}} objects it returns the common dimension of 
    the NMF problem fitted by all the runs -- based on the dimension of the first
    fit in the list.
    
    }
    
    \item{nrow}{: 
    returns the number of rows of the basis matrix.
    It corresponds to the number of rows of the fitted target matrix.
    }
    
    \item{ncol}{: 
    returns the number of columns of the mixture coefficient matrix.
    It corresponds to the number of columns of the fitted target matrix.
    }
    
    \item{nbasis}{: 
    returns the number of columns of the basis matrix.
    It corresponds to the rank of the factorization.
    
    When \code{x} is a matrix, then it returns the attribute 'nbasis', which is 
    for example attached to the consensus matrix returned by the method 
    \code{\link{consensus}}. This is used to keep track of data about the parent 
    fit and annotate plots.
    }
    
}}

\value{
	For \code{dim}, a 3-length integer vector.
	
	For \code{nrow}, \code{ncol} and \code{nbasis} a single integer.
}
\seealso{
	\code{\link[=dimnames,NMF-method]{dimnames}}
}
\examples{

# dimensions of an empty NMF model
dim( nmfModel() )

# dimensions of a random NMF model
x <- nmfModel(3, 100, 20)
dim(x)
nrow(x)
ncol(x)
nbasis(x)

# dimensions of a single NMF fit
V <- rmatrix(100, 20)
x <- nmf(V, 3)
dim(x)
nrow(x)
ncol(x)
nbasis(x)

# dimensions of a multiple NMF fit
x <- nmf(V, 3, nrun=3)
dim(x)
nrow(x)
ncol(x)
nbasis(x)

# dimensions of a multiple NMF fit (when keeping of the fits)
x <- nmf(V, 3, nrun=3, .opt='k')
dim(x)
nrow(x)
ncol(x)
nbasis(x)


}
\keyword{methods}
