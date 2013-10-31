\docType{methods}
\name{sparseness}
\alias{sparseness}
\alias{sparseness,matrix-method}
\alias{sparseness-methods}
\alias{sparseness,NMF-method}
\alias{sparseness,numeric-method}
\title{Sparseness}
\usage{
  sparseness(x, ...)
}
\arguments{
  \item{x}{an object whose sparseness is computed.}

  \item{...}{extra arguments to allow extension}
}
\value{
  usually a single numeric value -- in [0,1], or a numeric
  vector. See each method for more details.
}
\description{
  Generic function that computes the \emph{sparseness} of
  an object, as defined by \cite{Hoyer (2004)}. The
  sparseness quantifies how much energy of a vector is
  packed into only few components.
}
\details{
  In \cite{Hoyer (2004)}, the sparseness is defined for a
  real vector \eqn{x} as: \deqn{Sparseness(x) =
  \frac{\sqrt{n} - \frac{\sum |x_i|}{\sqrt{\sum
  x_i^2}}}{\sqrt{n}-1}}{ (srqt(n) - ||x||_1 / ||x||_2) /
  (sqrt(n) - 1)}

  , where \eqn{n} is the length of \eqn{x}.

  The sparseness is a real number in \eqn{[0,1]}. It is
  equal to 1 if and only if \code{x} contains a single
  nonzero component, and is equal to 0 if and only if all
  components of \code{x} are equal. It interpolates
  smoothly between these two extreme values. The closer to
  1 is the sparseness the sparser is the vector.

  The basic definition is for a \code{numeric} vector, and
  is extended for matrices as the mean sparseness of its
  column vectors.
}
\section{Methods}{
  \describe{

  \item{sparseness}{\code{signature(x = "numeric")}: Base
  method that computes the sparseness of a numeric vector.

  It returns a single numeric value, computed following the
  definition given in section \emph{Description}. }

  \item{sparseness}{\code{signature(x = "matrix")}:
  Computes the sparseness of a matrix as the mean
  sparseness of its column vectors. It returns a single
  numeric value. }

  \item{sparseness}{\code{signature(x = "NMF")}: Compute
  the sparseness of an object of class \code{NMF}, as the
  sparseness of the basis and coefficient matrices computed
  separately.

  It returns the two values in a numeric vector with names
  \sQuote{basis} and \sQuote{coef}. }

  }
}
\references{
  Hoyer P (2004). "Non-negative matrix factorization with
  sparseness constraints." _The Journal of Machine Learning
  Research_, *5*, pp. 1457-1469. <URL:
  http://portal.acm.org/citation.cfm?id=1044709>.
}
\seealso{
  Other assess: \code{\link{entropy}}, \code{\link{purity}}
}
\keyword{methods}

