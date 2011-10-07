\name{nmf.equal-methods}

\docType{methods}
\alias{nmf.equal}
\alias{nmf.equal-methods}
\alias{nmf.equal,NMFfit,NMF-method}
\alias{nmf.equal,NMFfit,NMFfit-method}
\alias{nmf.equal,NMFfitX,NMF-method}
\alias{nmf.equal,NMFfitX1,NMFfitX1-method}
\alias{nmf.equal,list,list-method}
\alias{nmf.equal,NMF,NMF-method}
\alias{nmf.equal,NMF,NMFfit-method}
\alias{nmf.equal,NMF,NMFfitX-method}

\title{Methods to Compare Two NMF Models}
\description{
 The method \code{nmf.equal} is used to compare two NMF models. The comparison is 
 performed entry-wise on both the basis and the mixture coefficient matrices.
  
 The two compared objects can be of class \code{\linkS4class{NMF}}, \code{\linkS4class{NMFfit}}, 
 \code{\linkS4class{NMFfitX}} or \code{list}, allowing the direct comparison of 
 models, fits or list of fits.
 In the case of \code{\linkS4class{NMFfitXn}} or \code{list} objects, either only the 
 best fits (default) or each fit separately are compared.
 
 }

\details{
The comparison is performed using either the function \code{\link{identical}} or 
the function \code{\link{all.equal}} on the \code{\linkS4class{NMF}} objects 
stored in \code{x} and \code{y}.

Currently, comparing lists of different length will return \code{FALSE}.
}

\usage{

\S4method{nmf.equal}{NMF,NMF}(x, y, identical=TRUE, ...)
\S4method{nmf.equal}{list,list}(x, y, all=FALSE, vector=FALSE, ...)

}

\arguments{

	\item{x}{ an object that inherits from class \code{\linkS4class{NMF}}, 
	\code{\linkS4class{NMFfit}}, \code{\linkS4class{NMFfitX}} or a list of these.	
	}
	
	\item{y}{ an object that inherits from class \code{\linkS4class{NMF}}, 
	\code{\linkS4class{NMFfit}}, \code{\linkS4class{NMFfitX}} or a list of these.	
	}	
	
	\item{identical}{ a single \code{logical} (default to \code{TRUE}) that specifies if the comparison 
	should be made using the function \code{\link{identical}} (i.e. exact equality) or the function 
	\code{\link{all.equal}} that allows for some numerical tolerance (controlled by argument \code{tolerance}).
	}
	
	\item{all}{If \code{all=FALSE} (default) the comparison is done between the 
	best fits of each list, i.e. the NMF model with the lowest deviance.
	Otherwise, the comparison is made pairwise.}
	
	\item{vector}{a logical used when \code{all=TRUE} that specifies if all the 
	pairwise comparisons must be returned as a vector. If \code{FALSE} (default), 
	the function will return \code{TRUE} only if all the fits in the lists 
	are identical/equals.
	}
	
	\item{...}{ extra arguments that are passed to \code{\link{all.equal}} 
	(e.g. \code{tolerance}) when \code{identical=FALSE}, or to the internal 
	call to \code{nmf.equal} when doing pairwise comparison.}
}

\value{

	When \code{identical=TRUE} this function returns \code{TRUE} if the two 
	\code{NMF} models are identical, and \code{FALSE} otherwise.
	
	When \code{identical=FALSE} this function returns the result of 
	\code{all.equal} applied to the two \code{NMF} models. 
	That is it returns either \code{TRUE} if the objects are considered equals 
	(given the numerical tolerance), or a \code{character} vector giving a 
	summary of the observed differences. 
	The numerical tolerance and comparison parameters can be passed through \code{...}.
	
	When doing pairwise comparison (\code{all=TRUE}) and \code{vector=FALSE}, 
	the function returns \code{TRUE} only if all the fits in the lists are 
	identical/equals.
	When \code{vector=TRUE} the returned value is a logical vector of the same 
	length as the input lists.
}

\seealso{ \code{\link{identical}}, \code{\link{all.equal}} } 
\keyword{methods}

\examples{

# create a random NMF model
a <- rnmf(100, 3, 20)

# compare a with itself
nmf.equal(a, a)
nmf.equal(a, a, identical=FALSE)
\dontshow{
stopifnot( nmf.equal(a, a) )
stopifnot( nmf.equal(a, a, identical=FALSE) )
}

# compare with a slightly different model
b <- a
basis(b) <- basis(b) + .Machine$double.eps^0.6
nmf.equal(a, b)
\dontshow{ stopifnot( !nmf.equal(a, b) ) }
nmf.equal(a, b, identical=FALSE) 
\dontshow{ stopifnot( nmf.equal(a, b, identical=FALSE) ) }

# compare with a model that is more different
b <- a
basis(b) <- basis(b) + .Machine$double.eps^0.4
nmf.equal(a, b)
\dontshow{ stopifnot( !nmf.equal(a, b) ) }
nmf.equal(a, b, identical=FALSE)
\dontshow{ stopifnot( !isTRUE(nmf.equal(a, b, identical=FALSE)) ) }

# the same works on results from function 'nmf'
V <- rmatrix(100, 20)
a <- nmf(V, 3)

# compare the result NMFfit object with the fitted NMF object
nmf.equal(a, fit(a))

# comparing with another run is unlikely to give the same result
b <- nmf(V, 3)
nmf.equal(a, b)
\dontshow{ stopifnot( !nmf.equal(a, b) ) }
# even when allowing for non identical objects
nmf.equal(a, b, FALSE)
\dontshow{ stopifnot( !isTRUE(nmf.equal(a, b, FALSE)) ) }
# but using the same random seed returns identical results
b <- nmf(V, 3, seed=getRNG(a))
nmf.equal(a, b)
\dontshow{ stopifnot( nmf.equal(a, b) ) }

# compare for multiple runs
a <- nmf(V, 3, nrun=3)

# nmf.equal compares the best fit
nmf.equal(a, fit(a))
\dontshow{ stopifnot( nmf.equal(a, fit(a)) ) }
# compare with another (random) NMF model
nmf.equal(a, rnmf(a))
# running another set of fits is unlikely to give an identical fit
b <- nmf(V, 3, nrun=3)
nmf.equal(a, b)

# but again: using the same random seed returns an identical fit
b <- nmf(V, 3, nrun=3, seed=getRNG1(a))
nmf.equal(a, b)
\dontshow{ stopifnot( nmf.equal(a, b) ) }


# compare all the runs from multiple runs
a <- nmf(V, 3, nrun=3, .opt='k')

# the comparison is done between the best fit
nmf.equal(a, fit(a))
\dontshow{ stopifnot( nmf.equal(a, fit(a)) ) }
# running another set of fits is unlikely to give identical results (for any of the runs)
b <- nmf(V, 3, nrun=3, .opt='k')
nmf.equal(a, b)
nmf.equal(a, b, all=TRUE)
nmf.equal(a, b, all=TRUE, vector=TRUE)

# but again: using the same random seed returns identical results (for each run)
b <- nmf(V, 3, nrun=3, seed=getRNG1(a), .opt='k')
nmf.equal(a, b)
nmf.equal(a, b, all=TRUE)
nmf.equal(a, b, all=TRUE, vector=TRUE)
\dontshow{ 
stopifnot( nmf.equal(a, b) )
stopifnot( nmf.equal(a, b, all=TRUE) )
stopifnot( all(nmf.equal(a, b, all=TRUE, vector=TRUE)) )
}

}