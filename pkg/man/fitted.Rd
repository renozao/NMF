\docType{methods}
\name{fitted}
\alias{fitted}
\alias{fitted-methods}
\alias{fitted,NMFfit-method}
\alias{fitted,NMF-method}
\alias{fitted,NMFns-method}
\alias{fitted,NMFOffset-method}
\alias{fitted,NMFstd-method}
\title{Fitted Matrix in NMF Models}
\usage{
  fitted(object, ...)

  \S4method{fitted}{NMFstd}(object, W, H, ...)

  \S4method{fitted}{NMFOffset}(object, W, H,
    offset = object@offset)

  \S4method{fitted}{NMFns}(object, W, H, S, ...)
}
\arguments{
  \item{object}{an object that inherit from class
  \code{NMF}}

  \item{...}{extra arguments to allow extension}

  \item{W}{a matrix to use in the computation as the basis
  matrix in place of \code{basis(object)}. It must be
  compatible with the coefficient matrix used in the
  computation (i.e. number of columns in \code{W} = number
  of rows in \code{H}).}

  \item{H}{a matrix to use in the computation as the
  coefficient matrix in place of \code{coef(object)}. It
  must be compatible with the basis matrix used in the
  computation (i.e. number of rows in \code{H} = number of
  columns in \code{W}).}

  \item{offset}{offset vector}

  \item{S}{smoothing matrix to use instead of
  \code{smoothing(object)} It must be a square matrix
  compatible with the basis and coefficient matrices used
  in the computation.}
}
\value{
  the target matrix estimate as fitted by the model
  \code{object}
}
\description{
  Computes the estimated target matrix based on a given
  \emph{NMF} model. The estimation depends on the
  underlying NMF model. For example in the standard model
  \eqn{V \equiv W H}{V ~ W H}, the target matrix is
  estimated by the matrix product \eqn{W H}. In other
  models, the estimate may depend on extra
  parameters/matrix (cf. Non-smooth NMF in
  \code{\link{NMFns-class}}).
}
\details{
  This function is a S4 generic function imported from
  \link[stats]{fitted} in the package \emph{stats}. It is
  implemented as a pure virtual method for objects of class
  \code{NMF}, meaning that concrete NMF models must provide
  a definition for their corresponding class (i.e.
  sub-classes of class \code{NMF}). See
  \code{\linkS4class{NMF}} for more details.
}
\section{Methods}{
  \describe{

  \item{fitted}{\code{signature(object = "NMF")}: Pure
  virtual method for objects of class
  \code{\linkS4class{NMF}}, that should be overloaded by
  sub-classes, and throws an error if called. }

  \item{fitted}{\code{signature(object = "NMFstd")}:
  Compute the target matrix estimate in \emph{standard NMF
  models}.

  The estimate matrix is computed as the product of the two
  matrix slots \code{W} and \code{H}: \deqn{\hat{V} = W
  H}{V ~ W H} }

  \item{fitted}{\code{signature(object = "NMFOffset")}:
  Computes the target matrix estimate for an NMFOffset
  object.

  The estimate is computed as: \deqn{ W H + offset } }

  \item{fitted}{\code{signature(object = "NMFns")}: Compute
  estimate for an NMFns object, according to the Nonsmooth
  NMF model (cf. \code{\link{NMFns-class}}).

  Extra arguments in \code{...} are passed to method
  \code{smoothing}, and are typically used to pass a value
  for \code{theta}, which is used to compute the smoothing
  matrix instead of the one stored in \code{object}. }

  \item{fitted}{\code{signature(object = "NMFfit")}:
  Computes and return the estimated target matrix from an
  NMF model fitted with function \code{\link{nmf}}.

  It is a shortcut for \code{fitted(fit(object), ...)},
  dispatching the call to the \code{fitted} method of the
  actual NMF model. }

  }
}
\examples{
# random standard NMF model
x <- rnmf(3, 10, 5)
all.equal(fitted(x), basis(x) \%*\% coef(x))
}
\keyword{methods}

