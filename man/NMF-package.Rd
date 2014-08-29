\docType{package}
\name{NMF-package}
\alias{NMF}
\alias{NMF-package}
\title{Algorithms and framework for Nonnegative Matrix Factorization (NMF)}
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
  \tabular{ll}{ Package: \tab NMF\cr Type: \tab Package\cr
  Version: \tab 0.21.3\cr Date: \tab 2013-09-17\cr License:
  \tab GPL (>=2)\cr LazyLoad: \tab yes\cr }
}
\examples{
\dontshow{# roxygen generated flag
options(R_CHECK_RUNNING_EXAMPLES_=TRUE)
}

# generate a synthetic dataset with known classes
n <- 50; counts <- c(5, 5, 8);
V <- syntheticNMF(n, counts)

# perform a 3-rank NMF using the default algorithm
res <- nmf(V, 3)

basismap(res)
coefmap(res)
}
\author{
  Renaud Gaujoux

Cathal Seoighe

Maintainer: Renaud Gaujoux <renaud@tx.technion.ac.il>
}
\references{
  \url{http://cran.r-project.org/}

  \url{http://renozao.github.io/NMF}
}
\seealso{
  \code{\link{nmf}}
}
\keyword{package}

