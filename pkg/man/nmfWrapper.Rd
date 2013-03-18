\name{nmfWrapper}
\alias{nmfWrapper}
\title{Wrapping NMF Algorithms}
\usage{
  nmfWrapper(method, ..., .FIXED = FALSE)
}
\arguments{
  \item{method}{Name of the NMF algorithm to be wrapped. It
  should be the name of a registered algorithm as returned
  by \code{\link{nmfAlgorithm}}, or an NMF algorithm object
  (i.e. an instance of \code{\linkS4class{NMFStrategy}}).}

  \item{...}{extra named arguments that define default
  values for any arguments of \code{\link{nmf}} or the
  algorithm itself.}

  \item{.FIXED}{a logical that indicates if the default
  arguments defined in \code{...} must be considered as
  fixed, i.e. that they are forced to have the defined
  values and cannot be used in a call to the wrapper
  function, in which case, a warning about discarding them
  is thrown if they are used. Non fixed arguments may have
  their value changed at call time, in which case it is
  honoured and passed to the \code{nmf} call.

  \code{.FIXED} may also be a character vector that
  specifies which argument amongst \code{...} should be
  considered as fixed.}
}
\value{
  a function with argument \code{...} and a set of default
  arguments defined in \code{...} in the call to
  \code{nmfWrapper}.
}
\description{
  This function creates a wrapper function for calling the
  function \code{\link{nmf}} with a given NMF algorithm.
}
\examples{
# wrap Lee & Seung algorithm into a function
lee <- nmfWrapper('lee', seed=12345)
args(lee)

# test on random data
x <- rmatrix(100,20)
res <- nmf(x, 3, 'lee', seed=12345)
res2 <- lee(x, 3)
nmf.equal(res, res2)
res3 <- lee(x, 3, seed=123)
nmf.equal(res, res3)

\dontshow{
stopifnot(nmf.equal(res, res2))
stopifnot( !nmf.equal(res, res3))
}

# argument 'method' has no effect
res4 <- lee(x, 3, method='brunet')
nmf.equal(res, res4)

\dontshow{
stopifnot(nmf.equal(res, res4))
}
}
\seealso{
  \code{\link{nmfAlgorithm}}, \code{\link{nmf}}
}
\keyword{internal}

