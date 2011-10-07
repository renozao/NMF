\name{basis/coef-methods: Accessing NMF matrix factors}
\docType{methods}
\alias{basis-methods}
\alias{basis}
\alias{basis,NMF-method}
\alias{basis<-}
\alias{basis<--methods}
\alias{basis<-,NMF,matrix-method}

\alias{coef}
\alias{coefficients}
\alias{coef,NMF-method}
\alias{coef<--methods}
\alias{coef<-}
\alias{coef<-,NMF,matrix-method}
\alias{coefficients,NMF-method}
\alias{coefficients<-,NMF,matrix-method}

\alias{scoef}
\alias{scoef-methods}
\alias{scoef,NMF-method}

\title{ Get/Set the matrix factors in a NMF model }
\description{
 \code{basis} and \code{basis<-} are S4 generic functions which respectively 
 extract and set the matrix of basis vectors (i.e. the first matrix factor) 
 of a NMF model. 
 For example, in the case of the standard NMF model \eqn{V \equiv WH}, 
 method \code{basis} will return matrix \eqn{W}.
 
 \code{coef} and \code{coef<-} are S4 methods defined for the associated generic 
 functions from package \code{stats} (See \link[stats]{coef}). 
 They respectively extract and set the matrix of mixture coefficients 
 (i.e. the second matrix factor) of a NMF model. 
 For example, in the case of the standard NMF model \eqn{V \equiv WH}, 
 method \code{coef} will return matrix \eqn{H}.
 
 Methods \code{coefficients} and \code{coefficients<-} are simple aliases for methods 
 \code{coef} and \code{coef<-} respectively.
 
 \code{scoef} returns the mixture coefficient matrix of a NMF model with the 
 columns scaled so that they sum up to one.
 
}
\section{Methods}{
\describe{

	\item{basis}{\code{signature(object = "NMF")}:
	Extracts the matrix of basis vectors from NMF model \code{object}.
	
	 Note that it is implemented as a pure virtual method, that must be 
	 overloaded by sub-classes that implement concrete NMF models. It throws an 
	 error if directly called.
	 See \code{\linkS4class{NMF}} for more details.
	 }
	
	\item{basis<-}{\code{signature(object = "NMF", value = "matrix"}:
	 Sets the matrix of basis vectors from NMF model \code{object}.
	
	 Note that it is implemented as a pure virtual method, that must be 
	 overloaded by sub-classes that implement concrete NMF models. It throws an 
	 error if directly called.
	 See \code{\linkS4class{NMF}} for more details.
	}
	

	\item{coef}{\code{signature(object = "NMF")}:
	Extracts the matrix of mixture coefficients from NMF model \code{object}.
	
	 Note that it is implemented as a pure virtual method, that must be 
	 overloaded by sub-classes that implement concrete NMF models. It throws an 
	 error if directly called.
	 See \code{\linkS4class{NMF}} for more details.
	 }
	
	\item{coef<-}{\code{signature(object = "NMF", value = "matrix"}:
	 Sets the matrix of mixture coefficients from NMF model \code{object}.
	
	 Note that it is implemented as a pure virtual method, that must be 
	 overloaded by sub-classes that implement concrete NMF models. It throws an 
	 error if directly called.
	 See \code{\linkS4class{NMF}} for more details.
	}
	
	\item{scoef}{ returns the mixture coefficient matrix of a NMF model with the 
	columns scaled so that they sum up to one.}
}
}
\keyword{methods}
\seealso{ NMF, NMFstd, NMFfit}
