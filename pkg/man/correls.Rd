\name{Correlations between basis/profiles}
\docType{methods}

\alias{basiscor}
\alias{basiscor-methods}
\alias{basiscor,NMF,NMF-method}
\alias{basiscor,NMF,matrix-method}
\alias{basiscor,matrix,NMF-method}
\alias{basiscor,NMF,missing-method}

\alias{profcor}
\alias{profcor-methods}
\alias{profcor,NMF,NMF-method}
\alias{profcor,NMF,matrix-method}
\alias{profcor,matrix,NMF-method}
\alias{profcor,NMF,missing-method}

\title{Compute the Correlation Between NMF Basis and Mixture Coefficient Matrices}

\description{
Functions \code{basiscor} (resp. \code{profcor}) computes the correlation matrix 
between the basis vectors (resp. the basis profiles) of two NMF models, or of 
an NMF model and a given compatible matrix. 

The computation uses the base function \code{\link{cor}}.
}

\usage{

\S4method{basiscor}{NMF,NMF}(x, y, ...)
\S4method{basiscor}{NMF,matrix}(x, y, ...)

\S4method{profcor}{NMF,NMF}(x, y, ...)
\S4method{profcor}{NMF,matrix}(x, y, ...)

}

\arguments{

	\item{x}{ A \code{matrix} or an object that inherits from class \code{\linkS4class{NMF}}.}
	
	\item{y}{ A \code{matrix} or an object that inherits from class \code{\linkS4class{NMF}}.
	If \var{y} is missing then the correlations are computed using \var{y=x}.
	}
	
	\item{...}{ extra arguments passed to function \code{\link{cor}}. }
}

\details{
\describe{

\item{basiscor}{ computes the correlation matrix between the basis vectors 
(i.e. the \emph{columns} of the basis vector matrix) or the \code{columns} (for 
matrix arguments) of \var{x} and \var{y}. 

The arguments' dimensions must be compatible, i.e same number of rows and 
basis vectors (or columns for matrix arguments).

If \var{y} is missing, then the correlations are computed between \var{x} and \code{y=x}. 
}


\item{profcor}{: computes the correlation matrix between the basis profiles 
(i.e. the \emph{rows} of the mixture coefficient matrix) or the rows (for 
matrix arguments) of \var{x} and \var{y}. 

The arguments' dimensions must be compatible, i.e same number of columns and 
basis vectors (or rows for matrix arguments).

If \var{y} is missing, then the correlations are computed between \var{x} and \code{y=x}. 
}

}
}

\seealso{ \code{\link{basis}}, \code{\link[=coef]{coef,NMF-method}}}

\author{Renaud Gaujoux}

\examples{

# generate two random NMF model
a <- rnmf(100, 3, 20)
b <- rnmf(100, 3, 20)

# Compute auto-correlations
basiscor(a)
# Compute correlations with b
basiscor(a, b)

# try to recover the underlying NMF model 'a'
res <- nmf(fitted(a), 3)

# Compute correlations with the true model
basiscor(a, res)
profcor(a, res)

# Compute correlations with a random compatible matrix
W <- rmatrix(basis(a))
basiscor(a, W)
\dontshow{ basiscor(W, a) }

H <- rmatrix(coef(a))
profcor(a, H)
\dontshow{ profcor(H, a) }

}

\keyword{methods}
\keyword{plot}

