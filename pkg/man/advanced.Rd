\name{advanced}
\alias{NMF-advanced}
\alias{isNMFfit}
\alias{is.nmf}

\title{ Advanced usage of package NMF }
\description{
	The functions documented here provide advanced functionalities useful when developing 
	within the framework implemented in the NMF package. 
	
}
\usage{

isNMFfit(object, recursive=TRUE)

is.nmf(object)

}

\arguments{

	\item{object}{any R object.}
	
	\item{recursive}{ if \code{TRUE} and \code{object} is a list then the check is performed on each 
	element of the list. Note that the recursivity only applies in the case of lists that are not 
	themselves NMFfit objects, unlike \code{NMFfitXn} objects for which the result of \code{isNMFfit} 
	will always be \code{TRUE} (a single logical value).  
	 
	}
		
}
\details{

\describe{

	\item{isNMFfit}{ tells if an object results from an NMF fit. That is it checks if \code{object} 
	inherits from class \code{\linkS4class{NMFfit}} or form class \code{\linkS4class{NMFfitX}}, which 
	are returned by the function \code{\link{nmf}}.
	If \code{object} is a \code{list} and \code{recursive=TRUE}, then the check is performed on each 
	element of the list, and the return value is a vector (or a list if \code{object} is a list of list) 
	of the same length as \code{object}. 
	}
	
	\item{is.nmf}{ tests if an object is an NMF model or a class that extends 
	the class NMF.
	}	
}
}
\value{
	For \code{isNMFfit}, a \code{logical} vector (or a list if \code{object} is a list of list) 
	of the same length as \code{object}.  
	    
}

\author{ Renaud Gaujoux  }

\seealso{
	\code{\linkS4class{NMFfit}}, \code{\linkS4class{NMFfitX}}
}

\examples{

	# generate a random 50 x 10 matrix
	V <- rmatrix(50, 10)
	
	# single run
	res <- nmf(V, 3)
	isNMFfit(res)
	
	# multiple runs - keeping single fit
	resm <- nmf(V, 3, nrun=3)
	isNMFfit(resm)
	
	# multiple runs - keeping all fits
	resM <- nmf(V, 3, nrun=3, .opt='k') 
	isNMFfit(resM)
	
	# with a list of results
	isNMFfit(list(res, resm, resM, 'not a result'))
	isNMFfit(list(res, list(resm, resM), 'not a result')) # list of list
	isNMFfit(list(res, resm, resM, 'not a result'), recursive=FALSE)
	
}
