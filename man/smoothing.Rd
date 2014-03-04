\name{smoothing}
\alias{smoothing}
\title{Smoothing Matrix in Nonsmooth NMF Models}
\usage{
  smoothing(x, theta = x@theta, ...)
}
\arguments{
  \item{x}{a object of class \code{NMFns}.}

  \item{theta}{the smoothing parameter (numeric) between 0
  and 1.}

  \item{...}{extra arguments to allow extension (not used)}
}
\value{
  if \code{x} estimates a \eqn{r}-rank NMF, then the result
  is a \eqn{r \times r} square matrix.
}
\description{
  The function \code{smoothing} builds a smoothing matrix
  for using in Nonsmooth NMF models.
}
\details{
  For a \eqn{r}-rank NMF, the smoothing matrix of parameter
  \eqn{\theta} is built as follows: \deqn{S = (1-\theta)I +
  \frac{\theta}{r} 11^T ,} where \eqn{I} is the identity
  matrix and \eqn{1} is a vector of ones (cf.
  \code{\link{NMFns-class}} for more details).
}
\examples{
\dontshow{# roxygen generated flag
options(R_CHECK_RUNNING_EXAMPLES_=TRUE)
}

x <- nmfModel(3, model='NMFns')
smoothing(x)
smoothing(x, 0.1)
}

