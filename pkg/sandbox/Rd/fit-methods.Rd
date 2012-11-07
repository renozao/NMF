\name{fit-methods: Extracting NMF models}
\docType{methods}

\alias{fit}
\alias{fit-methods}
\alias{fit,NMFfit-method}
\alias{fit,NMFfitX-method}

\alias{minfit}
\alias{minfit-methods}
\alias{minfit,NMFfit-method}
\alias{minfit,NMFfitX-method}

\title{ Extract the Best Model and Result Object From NMF Runs}
\description{
 Methods \code{fit} and \code{minfit} extract the best NMF model and fit objects respectively,  
 from the result of function \code{\link{nmf}}. That is the run that achieves the best approximation 
 error across all runs performed in the call to \code{nmf}.
 
 }
 
 \details{
 
 Method \code{fit} returns the NMF model itself as an object of class \code{\linkS4class{NMF}}, while 
 method \code{minfit} returns the result of the best run of function \code{nmf}. It is an object 
 of class \code{\linkS4class{NMFfit}} that, besides the NMF model, contains data about the run 
 (runtime, RNG seed, etc...). 
 
 They define a single interface to access the NMF model or fit result for objects of class 
 \code{NMFfit}, \code{NMFfitX1} and \code{NMFfitXn}, obtained from single runs, non-conservative multiple 
 runs (i.e. that only return the best fit), and conservative multiple runs (i.e. that return the list of 
 all the fits).
 
 IMPORTANT: note that the behaviour of \code{fit,NMFfitX1} and \code{fit,NMFfitXn} changed in version 0.5.3, 
 as these functions returned the best fit as an \code{NMFfit} object, which is what function \code{minfit} now does.
}

\usage{

\S4method{fit}{NMFfit}(object)

\S4method{minfit}{NMFfitX}(object)

}
\arguments{

\item{object}{ an object that inherits from class \code{\linkS4class{NMFfit}} or \code{\linkS4class{NMFfitX}}} 

\item{...}{ extra arguments to allow future extension.}

}

\value{

\itemize{

\item \code{fit} An object of class \code{\linkS4class{NMF}} that holds the NMF model that achieves the best 
approximation error.

\item \code{minfit}: An object of class \code{\linkS4class{NMFfit}} that contains data about the best run in addition 
to the actual NMF model.
 
}
	
}

\seealso{ \code{\link{residuals}}, \code{\link{nmf}} }

\keyword{methods}

\examples{

# generate a random target matrix
V <- rmatrix(100, 20)

# fit a single NMF model
res <- nmf(V, 3)

# extract best NMF model
fit(res)
# extract best run 
minfit(res)
# in the case of single runs it is the result itself
identical(minfit(res), res)

# perform non-conservative multiple runs
res <- nmf(V, 3, nrun=3)
# extract best NMF model
fit(res)
# extract best run 
minfit(res)

# perform conservative multiple runs
res <- nmf(V, 3, nrun=3, .opt='k')
# extract best NMF model
fit(res)
# extract best run 
minfit(res)


}
