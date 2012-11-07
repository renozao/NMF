\name{nmfWrapper}
\alias{nmfWrapper}
\title{Wrapping NMF Algorithms}
\usage{
  nmfWrapper(method, ...)
}
\arguments{
  \item{method}{Name of the NMF algorithm to be wrapped. It
  should be the name of a registered algorithm as returned
  by \code{\link{nmfAlgorithm}}, or an NMF algorithm object
  (i.e. an instance of \code{\linkS4class{NMFStrategy}}).}

  \item{...}{extra named arguments that define default
  values for any arguments of \code{\link{nmf}} or the
  algorithm itself.}
}
\value{
  a function with attribute \code{'algorithm'} set to the
  value of \code{key}.
}
\description{
  This function creates a wrapper function for calling the
  function \code{\link{nmf}} with a given NMF algorithm.
}
\examples{
# wrap Lee & Seung algorithm into a function
lee <- nmfWrapper('lee')
lee

# test on random data
x <- rmatrix(100,20)
res <- nmf(x, 3, 'lee', seed=12345)
res2 <- lee(x, 3, seed=12345)
nmf.equal(res, res2)

\dontshow{ stopifnot(nmf.equal(res, res2)) }
}
\seealso{
  \code{\link{nmfAlgorithm}}, \code{\link{nmf}}
}
\keyword{internal}

