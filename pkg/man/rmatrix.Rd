\name{rmatrix}
\alias{rmatrix}

\title{ Generates a Random Matrix Using Any Given Distribution Function}
\description{
This function provides a short-cut to generate a random matrix whose entries are 
drawn any given random distribution, as soon as this one is implemented as an \R function 
similar to: \code{runif}, \code{rnorm}, etc ...  

It basically wraps the following call:

\code{matrix(dist(nrow*ncol, ...), nrow, ncol)}

in a quicker and intuitive call:

\code{rmatrix(nrow, ncol, dist, ...)}

}
\usage{

rmatrix(nrow, ncol, dist = runif, byrow = FALSE, dimnames = NULL, ...)

}

\arguments{

  \item{nrow}{ the desired number of rows of the generated matrix.}
  
  \item{ncol}{ the desired number of columns of the generated matrix.}
  
  \item{dist}{ the distribution function from which to draw the matrix entries.
  It must be a -- distribution -- function whose first parameter is the number of draws 
  \code{n} to be done, and return a vector of length \code{n}.}
  
  \item{byrow}{logical. If \code{FALSE} (the default) the matrix is filled by
          columns, otherwise the matrix is filled by rows. See \code{\link{matrix}}.}
  
  \item{dimnames}{A \code{dimnames} attribute for the matrix: \code{NULL} or a \code{list} of
          length 2 giving the row and column names respectively.  An
          empty list is treated as \code{NULL}, and a list of length one as
          row names.  The list can be named, and the list names will be
          used as names for the dimensions. See \code{\link{matrix}}.}
  
  \item{\dots}{ any extra parameters to pass to function \code{dist}.}
}

\value{ a matrix of dimension \code{nrow} x \code{ncol} whose entries are drawn from 
distribution \code{dist}}

\author{ Renaud Gaujoux \email{renaud@cbio.uct.ac.za} }

\seealso{ \code{\link{runif}}, \code{\link{rnorm}} or any other \code{rABCDE} similar distribution functions.}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(nrow, ncol, dist=runif, ...){
	# check that 'dist' is a function.
	if( !is.function(dist) )
		stop("NMF::rmatrix - invalid value for argument 'dist': must be a function [class(dist)='", class(dist), "'].")
	
	# build the random matrix using the distribution function
	matrix(dist(nrow*ncol, ...), nrow, ncol)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
