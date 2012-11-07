\name{Utilities for Random Number Generators (RNG)}

\docType{methods}

\alias{RNG-NMF}

\alias{getRNG}
\alias{getRNG-methods}
\alias{getRNG,missing-method}
\alias{getRNG,list-method}
\alias{getRNG,ANY-method}
\alias{getRNG,NMFfitXn-method}
%\alias{getRNG,rstream-method}
\alias{getRNG,integer-method}

\alias{getRNG1}
\alias{getRNG1-methods}
\alias{getRNG1,NMFfit-method}
\alias{getRNG1,NMFfitX-method}
\alias{getRNG1,NMFfitX1-method}
\alias{getRNG1,NMFfitXn-method}

\alias{RNGdigest}
\alias{rng.equal}
\alias{rng1.equal}

\title{Get the Random Number Generator State}
\description{
 The function \code{getRNG} is a S4 generic function that returns either the 
 current RNG settings if called with no argument, or the RNG settings used to 
 generate an object.
 In the case of results from multiple NMF runs, it returns the RNG settings 
 used to compute the best fit.
 
 The function \code{getRNG1} is a S4 generic function that returns the RNG 
 settings used for the first of the multiple NMF runs performed to generate an 
 \code{\linkS4class{NMFfitX}} object.
 This is useful to reproduce or compare results from different multiple runs. 
 
 The function \code{RNGdigest} returns a hash that uniquely identifies given 
 RNG settings. 
 %The hash can be -- and is -- used to compare two \code{\link{rstream}}
 %objects and test if they represent the same random streams, in the same state.
 The hash is computed based on the result of \code{getRNG}.
  
 The function \code{rng.equal} tests whether two -- embedded -- RNG objects 
 represent the same random stream, in the same state. It does so by comparing 
 the hashes returned by \code{RNGdigest}.
 
 The function \code{rng1.equal} tests whether two results from multiple NMF runs 
 used identical initial RNG settings.
 }
 
\details{
The current methods for \code{getRNG} extracts the RNG settings from slot \code{rng} 
for S4 objects, or from element \code{rng} or \code{noise$rng} for a \code{list} 
(wich includes S3 objects), or from attribute \code{rng} otherwise.

}
\usage{

\S4method{getRNG}{missing}(object)
\S4method{getRNG}{ANY}(object)
\S4method{getRNG1}{NMFfitX}(object)

RNGdigest(x)
rng.equal(x, y)
rng1.equal(x, y)

}

\arguments{

\item{object}{ either \code{missing}, a \code{list}, or any object, usually with 
relevant stored RNG information, such as \code{\linkS4class{NMFfit}} objects.
If missing, the returned value is based on the current RNG settings.
 
For \code{getRNG1}, an \code{NMFfitX} object.}

%\item{packed}{a logical that specifies if the returned \code{rstream} object 
%should be packed or not.}

\item{x}{ any object handled by \code{getRNG}.
For \code{rng1.equal}, \code{x} must be an \code{\linkS4class{NMFfitX}} object, such as 
returned by multiple NMF runs.}

\item{y}{ any object handled by \code{getRNG}. If \code{y} is missing, then 
the current RNG settings are used.
For \code{rng1.equal}, \code{y} must be an \code{\linkS4class{NMFfitX}} object, 
such as returned by multiple NMF runs.}

}

\value{

\code{getRNG} returns an integer vector (see \code{\link{.Random.seed}}) 
%For \code{getRNG}, the class of the returned value depends on the type of the 
%extracted object (see details).
%For \code{NMFfit} objects, the result is an \code{\linkS4class{rstream}} object 
%that contains the RNG settings as before the computation starts, i.e. before the 
%seeding step.

\code{RNGdigest} returns the RNG hash as a single character string.

\code{rng.equal} and \code{rng1.equal} returns \code{TRUE} if the two 
-- embedded -- RNG objects are identical, \code{FALSE} otherwise.
See details for the differences between these two functions.
 
}

\keyword{methods}

\seealso{ \code{\link{.Random.seed}}, \code{\link{set.seed}}. }

\examples{

# current random seed (by default it is a 626-length numeric vector)
getRNG()

# fit an NMF model on a random target matrix
V <- rmatrix(100,20)
s <- getRNG()
res <- nmf(V, 3)

# the random seed changed since by default nmf use a randomly generating starting point 
rng.equal(s)
# but the starting RNG settings are stored in the object and can be accessed by getRNG
getRNG(res)
# the setting used are the one in use before running NMF  
rng.equal(s, res)
\dontshow{ stopifnot(rng.equal(s, res)) }

# show the digest version
RNGdigest(res)

# For multiple runs, the RNG settings used for the first run is also stored
res <- nmf(V, 3, nrun=4)
# RNG used for the best fit
getRNG(res)
# RNG used for the first fit
getRNG1(res)
# they may differ if the best fit is not the first one
rng.equal(res, getRNG1(res))


}
