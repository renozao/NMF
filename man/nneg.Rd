\docType{methods}
\name{nneg}
\alias{nneg}
\alias{nneg,matrix-method}
\alias{nneg-methods}
\alias{nneg,NMF-method}
\alias{posneg}
\alias{rposneg}
\alias{rposneg,matrix-method}
\alias{rposneg-methods}
\alias{rposneg,NMF-method}
\title{Transforming from Mixed-sign to Nonnegative Data}
\usage{
  nneg(object, ...)

  \S4method{nneg}{matrix}(object,
    method = c("pmax", "posneg", "absolute", "min"),
    threshold = 0, shift = TRUE)

  posneg(...)

  rposneg(object, ...)

  \S4method{rposneg}{matrix}(object, unstack = TRUE)
}
\arguments{
  \item{object}{The data object to transform}

  \item{...}{extra arguments to allow extension or passed down to \code{nneg,matrix}
or \code{rposneg,matrix} in subsequent calls.}

  \item{method}{Name of the transformation method to use, that is partially
matched against the following possible methods:
\describe{
\item{pmax}{Each entry is constrained to be above threshold \code{threshold}.}

\item{posneg}{The matrix is split into its "positive" and "negative" parts,
with the entries of each part constrained to be above threshold \code{threshold}.
The result consists in these two parts stacked in rows (i.e. \code{\link{rbind}}-ed)
into a single matrix, which has double the number of rows of the input
matrix \code{object}.}

\item{absolute}{The absolute value of each entry is constrained to be above
threshold \code{threshold}.}

\item{min}{Global shift by adding the minimum entry to each entry, only if
it is negative, and then apply threshold.
}

}}

  \item{threshold}{Nonnegative lower threshold value (single numeric).
See argument \code{shit} for details on how the threshold is used and affects
the result.}

  \item{shift}{a logical indicating whether the entries below the threshold
value \code{threshold} should be forced (shifted) to 0 (default) or to
the threshold value itself.
In other words, if \code{shift=TRUE} (default) all entries in
the result matrix are either 0 or strictly greater than \code{threshold}.
They are all greater or equal than \code{threshold} otherwise.}

  \item{unstack}{Logical indicating whether the positive and negative parts
should be unstacked and combined into a matrix as \code{pos - neg}, which contains
half the number of rows of \code{object} (default), or left
stacked as \code{[pos; -neg]}.}
}
\value{
  an object of the same class as argument \code{object}.

an object of the same type of \code{object}
}
\description{
  \code{nneg}  is a generic function to transform a data objects that
contains negative values into a similar object that only contains
values that are nonnegative or greater than a given threshold.

\code{posneg} is a shortcut for \code{nneg(..., method='posneg')}, to split
mixed-sign data into its positive and negative part.
See description for method \code{"posneg"}, in \code{\link{nneg}}.

\code{rposneg} performs the "reverse" transformation of the \code{\link{posneg}} function.
}
\section{Methods}{
\describe{

\item{nneg}{\code{signature(object = "matrix")}: Transforms a mixed-sign matrix into a nonnegative matrix, optionally apply a
lower threshold.
This is the workhorse method, that is eventually called by all other
methods defined in the \code{\link{NMF}} package.
}

\item{nneg}{\code{signature(object = "NMF")}: Apply \code{nneg} to the basis matrix of an \code{\link{NMF}}
object (i.e. \code{basis(object)}).
All extra arguments in \code{...} are passed to the method \code{nneg,matrix}.
}


\item{rposneg}{\code{signature(object = "NMF")}: Apply \code{rposneg} to the basis matrix of an \code{\link{NMF}} object.
}

}

}
\examples{
\dontshow{# roxygen generated flag
options(R_CHECK_RUNNING_EXAMPLES_=TRUE)
}

#----------
# nneg,matrix-method
#----------
# random mixed sign data (normal distribution)
set.seed(1)
x <- rmatrix(5,5, rnorm, mean=0, sd=5)
x

# pmax (default)
nneg(x)
# using a threshold
nneg(x, threshold=2)
# without shifting the entries lower than threshold
nneg(x, threshold=2, shift=FALSE)

# posneg: split positive and negative part
nneg(x, method='posneg')
nneg(x, method='pos', threshold=2)

# absolute
nneg(x, method='absolute')
nneg(x, method='abs', threshold=2)

# min
nneg(x, method='min')
nneg(x, method='min', threshold=2)

#----------
# nneg,NMF-method
#----------
# random
M <- nmfModel(x, rmatrix(ncol(x), 3))
nnM <- nneg(M)
basis(nnM)
# mixture coefficients are not affected
identical( coef(M), coef(nnM) )

#----------
# posneg
#----------
# shortcut for the "posneg" transformation
posneg(x)
posneg(x, 2)

#----------
# rposneg,matrix-method
#----------
# random mixed sign data (normal distribution)
set.seed(1)
x <- rmatrix(5,5, rnorm, mean=0, sd=5)
x

# posneg-transform: split positive and negative part
y <- posneg(x)
dim(y)
# posneg-reverse
z <- rposneg(y)
identical(x, z)
rposneg(y, unstack=FALSE)

# But posneg-transformation with a non zero threshold is not reversible
y1 <- posneg(x, 1)
identical(rposneg(y1), x)

#----------
# rposneg,NMF-method
#----------
# random mixed signed NMF model
M <- nmfModel(rmatrix(10, 3, rnorm), rmatrix(3, 4))
# split positive and negative part
nnM <- posneg(M)
M2 <- rposneg(nnM)
identical(M, M2)
}
\seealso{
  \code{\link{pmax}}

Other transforms: \code{\link{t.NMF}}
}
\keyword{methods}

