\docType{methods}
\name{basis}
\alias{.basis}
\alias{.basis,NMF-method}
\alias{.basis,NMFfit-method}
\alias{.basis,NMFstd-method}
\alias{.basis-methods}
\alias{.basis<-}
\alias{.basis<-,NMF,matrix-method}
\alias{.basis<-,NMFfit,matrix-method}
\alias{.basis<-,NMFstd,matrix-method}
\alias{.basis<--methods}
\alias{.coef}
\alias{.coef,NMF-method}
\alias{.coef,NMFfit-method}
\alias{.coef,NMFstd-method}
\alias{.coef-methods}
\alias{.coef<-}
\alias{.coef<-,NMF,matrix-method}
\alias{.coef<-,NMFfit,matrix-method}
\alias{.coef<-,NMFstd,matrix-method}
\alias{.coef<--methods}
\alias{basis}
\alias{basis,NMF-method}
\alias{basis,NMFfitXn-method}
\alias{basis-methods}
\alias{basis<-}
\alias{basis<-,ANY-method}
\alias{basis<--methods}
\alias{coef}
\alias{coef,NMF-method}
\alias{coef,NMFfitXn-method}
\alias{coef-methods}
\alias{coef<-}
\alias{coef<-,ANY-method}
\alias{coef<--methods}
\alias{coefficients}
\alias{coefficients,NMF-method}
\alias{coefficients-methods}
\alias{loadings,NMF-method}
\alias{scoef}
\alias{scoef,NMF-method}
\alias{scoef-methods}
\title{Accessing NMF Factors}
\usage{
  basis(object, ...)

  \S4method{basis}{NMF}(object, all = TRUE, ...)

  .basis(object, ...)

  basis(object) <- value

  .basis(object) <- value

  \S4method{loadings}{NMF}(x)

  coef(object, ...)

  \S4method{coef}{NMF}(object, all = TRUE, ...)

  .coef(object, ...)

  coef(object) <- value

  .coef(object) <- value

  coefficients(object, ...)

  \S4method{coefficients}{NMF}(object, all = TRUE, ...)

  scoef(object, ...)

  \S4method{scoef}{NMF}(object, scale = 1)
}
\arguments{
  \item{object}{an object from which to extract the factor
  matrices, typically an object of class
  \code{\linkS4class{NMF}}.}

  \item{...}{extra arguments to allow extension and passed
  to the low-level access functions \code{.coef} and
  \code{.basis}.}

  \item{all}{a logical that indicates whether the complete
  matrix factor should be returned (\code{TRUE}) or only
  the non-fixed part. This is relevant only for
  formula-based NMF models that include fixed basis or
  coefficient terms.}

  \item{value}{replacement value}

  \item{scale}{scaling factor, which indicates to the value
  the columns of the coefficient matrix should sum up to.}

  \item{x}{an object of class \code{"\link{factanal}"} or
  \code{"\link{princomp}"} or the \code{loadings} component
  of such an object.}
}
\description{
  \code{basis} and \code{basis<-} are S4 generic functions
  which respectively extract and set the matrix of basis
  components of an NMF model (i.e. the first matrix
  factor).

  The methods \code{.basis}, \code{.coef} and their
  replacement versions are implemented as pure virtual
  methods for the interface class \code{NMF}, meaning that
  concrete NMF models must provide a definition for their
  corresponding class (i.e. sub-classes of class
  \code{NMF}). See \code{\linkS4class{NMF}} for more
  details.

  \code{coef} and \code{coef<-} respectively extract and
  set the coefficient matrix of an NMF model (i.e. the
  second matrix factor). For example, in the case of the
  standard NMF model \eqn{V \equiv WH}{V ~ W H}, the method
  \code{coef} will return the matrix \eqn{H}.

  \code{.coef} and \code{.coef<-} are low-level S4 generics
  that simply return/set coefficient data in an object,
  leaving some common processing to be performed in
  \code{coef} and \code{coef<-}.

  Methods \code{coefficients} and \code{coefficients<-} are
  simple aliases for methods \code{coef} and \code{coef<-}
  respectively.

  \code{scoef} is similar to \code{coef}, but returns the
  mixture coefficient matrix of an NMF model, with the
  columns scaled so that they sum up to a given value (1 by
  default).
}
\details{
  For example, in the case of the standard NMF model \eqn{V
  \equiv W H}{V ~ W H}, the method \code{basis} will return
  the matrix \eqn{W}.

  \code{basis} and \code{basis<-} are defined for the top
  virtual class \code{\linkS4class{NMF}} only, and rely
  internally on the low-level S4 generics \code{.basis} and
  \code{.basis<-} respectively that effectively extract/set
  the coefficient data. These data are post/pre-processed,
  e.g., to extract/set only their non-fixed terms or check
  dimension compatibility.

  \code{coef} and \code{coef<-} are S4 methods defined for
  the corresponding generic functions from package
  \code{stats} (See \link[stats]{coef}). Similarly to
  \code{basis} and \code{basis<-}, they are defined for the
  top virtual class \code{\linkS4class{NMF}} only, and rely
  internally on the S4 generics \code{.coef} and
  \code{.coef<-} respectively that effectively extract/set
  the coefficient data. These data are post/pre-processed,
  e.g., to extract/set only their non-fixed terms or check
  dimension compatibility.
}
\section{Methods}{
  \describe{

  \item{.basis}{\code{signature(object = "NMF")}: Pure
  virtual method for objects of class
  \code{\linkS4class{NMF}}, that should be overloaded by
  sub-classes, and throws an error if called. }

  \item{.basis}{\code{signature(object = "NMFstd")}: Get
  the basis matrix in standard NMF models

  This function returns slot \code{W} of \code{object}. }

  \item{.basis}{\code{signature(object = "NMFfit")}:
  Returns the basis matrix from an NMF model fitted with
  function \code{\link{nmf}}.

  It is a shortcut for \code{.basis(fit(object), ...)},
  dispatching the call to the \code{.basis} method of the
  actual NMF model. }

  \item{.basis<-}{\code{signature(object = "NMF", value =
  "matrix")}: Pure virtual method for objects of class
  \code{\linkS4class{NMF}}, that should be overloaded by
  sub-classes, and throws an error if called. }

  \item{.basis<-}{\code{signature(object = "NMFstd", value
  = "matrix")}: Set the basis matrix in standard NMF models

  This function sets slot \code{W} of \code{object}. }

  \item{.basis<-}{\code{signature(object = "NMFfit", value
  = "matrix")}: Sets the the basis matrix of an NMF model
  fitted with function \code{\link{nmf}}.

  It is a shortcut for \code{.basis(fit(object)) <- value},
  dispatching the call to the \code{.basis<-} method of the
  actual NMF model. It is not meant to be used by the user,
  except when developing NMF algorithms, to update the
  basis matrix of the seed object before returning it. }

  \item{.coef}{\code{signature(object = "NMF")}: Pure
  virtual method for objects of class
  \code{\linkS4class{NMF}}, that should be overloaded by
  sub-classes, and throws an error if called. }

  \item{.coef}{\code{signature(object = "NMFstd")}: Get the
  mixture coefficient matrix in standard NMF models

  This function returns slot \code{H} of \code{object}. }

  \item{.coef}{\code{signature(object = "NMFfit")}: Returns
  the the coefficient matrix from an NMF model fitted with
  function \code{\link{nmf}}.

  It is a shortcut for \code{.coef(fit(object), ...)},
  dispatching the call to the \code{.coef} method of the
  actual NMF model. }

  \item{.coef<-}{\code{signature(object = "NMF", value =
  "matrix")}: Pure virtual method for objects of class
  \code{\linkS4class{NMF}}, that should be overloaded by
  sub-classes, and throws an error if called. }

  \item{.coef<-}{\code{signature(object = "NMFstd", value =
  "matrix")}: Set the mixture coefficient matrix in
  standard NMF models

  This function sets slot \code{H} of \code{object}. }

  \item{.coef<-}{\code{signature(object = "NMFfit", value =
  "matrix")}: Sets the the coefficient matrix of an NMF
  model fitted with function \code{\link{nmf}}.

  It is a shortcut for \code{.coef(fit(object)) <- value},
  dispatching the call to the \code{.coef<-} method of the
  actual NMF model. It is not meant to be used by the user,
  except when developing NMF algorithms, to update the
  coefficient matrix in the seed object before returning
  it. }

  \item{basis}{\code{signature(object = "NMFfitXn")}:
  Returns the basis matrix of the best fit amongst all the
  fits stored in \code{object}. It is a shortcut for
  \code{basis(fit(object))}. }

  \item{basis<-}{\code{signature(object = "ANY")}: Default
  methods that calls \code{.basis<-} and check the validity
  of the updated object. }

  \item{coef}{\code{signature(object = "NMFfitXn")}:
  Returns the coefficient matrix of the best fit amongst
  all the fits stored in \code{object}. It is a shortcut
  for \code{coef(fit(object))}. }

  \item{coef<-}{\code{signature(object = "ANY")}: Default
  methods that calls \code{.coef<-} and check the validity
  of the updated object. }

  \item{coefficients}{\code{signature(object = "NMF")}:
  Alias to \code{coef,NMF}, therefore also pure virtual. }

  \item{loadings}{\code{signature(x = "NMF")}: Method
  loadings for NMF Models

  The method \code{loadings} is identical to \code{basis},
  but do not accept any extra argument.

  The method \code{loadings} is provided to standardise the
  NMF interface against the one defined in the
  \code{\link{stats}} package, and emphasises the
  similarities between NMF and PCA or factorial analysis
  (see \code{\link{loadings}}). }

  }
}
\examples{
# Scaled coefficient matrix
x <- rnmf(3, 10, 5)
scoef(x)
scoef(x, 100)
# random standard NMF model
x <- rnmf(3, 10, 5)
basis(x)
coef(x)

# set matrix factors
basis(x) <- matrix(1, nrow(x), nbasis(x))
coef(x) <- matrix(1, nbasis(x), ncol(x))
# set random factors
basis(x) <- rmatrix(basis(x))
coef(x) <- rmatrix(coef(x))

# incompatible matrices generate an error:
try( coef(x) <- matrix(1, nbasis(x)-1, nrow(x)) )
# but the low-level method allow it
.coef(x) <- matrix(1, nbasis(x)-1, nrow(x))
try( validObject(x) )
}
\seealso{
  Other NMF-interface: \code{\link{$,NMF-method}},
  \code{\link{$<-,NMF-method}},
  \code{\link{.DollarNames,NMF-method}},
  \code{\link{NMF-class}}, \code{\link{misc}},
  \code{\link{nmfModel}}, \code{\link{nmfModels}},
  \code{\link{rnmf}}
}
\keyword{methods}

