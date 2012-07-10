\docType{methods}
\name{rmatrix}
\alias{rmatrix}
\alias{rmatrix,ANY-method}
\alias{rmatrix,NMF-method}
\alias{rmatrix,numeric-method}
\alias{rmatrix-methods}
\title{Generating Random Matrices}
\usage{
  rmatrix(x, ...)

  \S4method{rmatrix}{numeric}(x, y = NULL, dist = runif,
    byrow = FALSE, dimnames = NULL, ...)
}
\arguments{
  \item{x}{object from which to generate a random matrix}

  \item{y}{optional specification of number of columns}

  \item{dist}{a random distribution function (see details
  of method \code{rmatrix,numeric})}

  \item{byrow}{a logical passed in the internal call to the
  function \code{\link{matrix}}}

  \item{dimnames}{\code{NULL} or a \code{list} passed in
  the internal call to the function \code{\link{matrix}}}

  \item{...}{extra arguments passed to the distribution
  function \code{dist}.}
}
\description{
  The S4 generic \code{rmatrix} generates a random matrix
  from a given object. Methods are provided to generate
  matrices with entries drawn from any given random
  distribution function, e.g. \code{\link{runif}} or
  \code{\link{rnorm}}.
}
\section{Methods}{
  \describe{

  \item{rmatrix}{\code{signature(x = "numeric")}: Generates
  a random matrix of given dimensions, whose entries are
  drawn using the distribution function \code{dist}.

  This is the workhorse method that is eventually called by
  all other methods. It returns a matrix with: \itemize{
  \item \code{x} rows and \code{y} columns if \code{y} is
  not missing and not \code{NULL}; \item dimension
  \code{x[1]} x \code{x[2]} if \code{x} has at least two
  elements; \item dimension \code{x} (i.e. a square matrix)
  otherwise. }

  The default is to draw its entries from the standard
  uniform distribution using the base function
  \code{\link{runif}}, but any other function that
  generates random numeric vectors of a given length may be
  specified in argument \code{dist}. All arguments in
  \code{...} are passed to the function specified in
  \code{dist}.

  The only requirement is that the function in \code{dist}
  is of the following form:

  \samp{ function(n, ...){ # return vector of length n ...
  }}

  This is the case of all base random draw function such as
  \code{\link{rnorm}}, \code{\link{rgamma}}, etc\ldots }

  \item{rmatrix}{\code{signature(x = "ANY")}: Default
  method which calls \code{rmatrix,vector} on the
  dimensions of \code{x} that is assumed to be returned by
  a suitable \code{dim} method: it is equivalent to
  \code{rmatrix(dim(x), y=NULL, ...)}. }

  \item{rmatrix}{\code{signature(x = "NMF")}: Returns the
  target matrix estimate of the NMF model \code{x},
  perturbated by adding a random matrix generated using the
  default method of \code{rmatrix}: it is a equivalent to
  \code{fitted(x) + rmatrix(fitted(x), ...)}.

  This method can be used to generate random target
  matrices that depart from a known NMF model to a
  controlled extend. This is useful to test the robustness
  of NMF algorithms to the presence of certain types of
  noise in the data. }

  }
}
\examples{
## Generate a random matrix of a given size
rmatrix(5, 3)
\dontshow{ stopifnot( identical(dim(rmatrix(5, 3)), c(5L,3L)) ) }

## Generate a random matrix of the same dimension of a template matrix
a <- matrix(1, 3, 4)
rmatrix(a)
\dontshow{ stopifnot( identical(dim(rmatrix(a)), c(3L,4L)) ) }

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
# random matrix of the same dimension as another matrix
x <- matrix(3,4)
dim(rmatrix(x))
# generate noisy fitted target from an NMF model (the true model)
gr <- as.numeric(mapply(rep, 1:3, 3))
h <- outer(1:3, gr, '==') + 0
x <- rnmf(10, H=h)
y <- rmatrix(x)
\dontrun{
# show heatmap of the noisy target matrix: block patterns should be clear
aheatmap(y)
}
\dontshow{ stopifnot( identical(dim(y), dim(x)[1:2]) ) }

# test NMF algorithm on noisy data
# add some noise to the true model (drawn from uniform [0,1])
res <- nmf(rmatrix(x), 3)
summary(res)

# add more noise to the true model (drawn from uniform [0,10])
res <- nmf(rmatrix(x, max=10), 3)
summary(res)
}
\keyword{methods}

