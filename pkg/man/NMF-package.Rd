\docType{package}
\name{NMF-package}
\alias{NMF}
\alias{NMF-package}
\title{Algorithms and framework for Nonnegative Matrix Factorization (NMF).}
\description{
  This package provides a framework to perform Non-negative
  Matrix Factorization (NMF). It implements a set of
  already published algorithms and seeding methods, and
  provides a framework to test, develop and plug new/custom
  algorithms. Most of the built-in algorithms have been
  optimized in C++, and the main interface function
  provides an easy way of performing parallel computations
  on multicore machines.
}
\details{
  \code{\link{nmf}} Run a given NMF algorithm
}
\examples{
# generate a synthetic dataset with known classes
n <- 50; counts <- c(5, 5, 8);
V <- syntheticNMF(n, counts)

# perform a 3-rank NMF using the default algorithm
res <- nmf(V, 3)

basismap(res)
coefmap(res)
}
\author{
  Renaud Gaujoux \email{renaud@cbio.uct.ac.za}
}
\references{
  \url{http://www.r-project.org/}
}
\seealso{
  \code{\link{nmf}}
}
\keyword{package}

