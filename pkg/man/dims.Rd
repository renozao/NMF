\docType{methods}
\name{nbasis}
\alias{dim,NMF-method}
\alias{dim,NMFfitXn-method}
\alias{dim-NMF}
\alias{nbasis}
\alias{nbasis,ANY-method}
\alias{nbasis,NMFfitXn-method}
\alias{nbasis,matrix-method}
\alias{nbasis-methods}
\title{Dimension of NMF Objects}
\usage{
  nbasis(x, ...)

  \S4method{dim}{NMF}(x)

  \S4method{dim}{NMFfitXn}(x)
}
\arguments{
  \item{x}{an object with suitable \code{basis} and
  \code{coef} methods, such as an object that inherit from
  \code{\linkS4class{NMF}}.}

  \item{...}{extra arguments to allow extension.}
}
\value{
  a single integer value or, for \code{dim}, a length-3
  integer vector, e.g. \code{c(2000, 30, 3)} for an
  \code{NMF} model that fits a 2000 x 30 matrix using 3
  basis components.
}
\description{
  The methods \code{dim}, \code{nrow}, \code{ncol} and
  \code{nbasis} return the different dimensions associated
  with an NMF model.

  \code{dim} returns all dimensions in a length-3 integer
  vector: the number of row and columns of the estimated
  target matrix, as well as the factorization rank (i.e.
  the number of basis components).

  \code{nrow}, \code{ncol} and \code{nbasis} provide
  separate access to each of these dimensions respectively.
}
\details{
  The NMF package does not implement specific functions
  \code{nrow} and \code{ncol}, but rather the S4 method
  \code{dim} for objects of class \code{\linkS4class{NMF}}.
  This allows the base methods \code{\link{nrow}} and
  \code{\link{ncol}} to directly work with such objects, to
  get the number of rows and columns of the target matrix
  estimated by an NMF model.

  The function \code{nbasis} is a new S4 generic defined in
  the package NMF, that returns the number of basis
  components of an object. Its default method should work
  for any object, that has a suitable \code{basis} method
  defined for its class.
}
\section{Methods}{
  \describe{

  \item{dim}{\code{signature(x = "NMF")}: method for NMF
  objects for the base generic \code{\link{dim}}. It
  returns all dimensions in a length-3 integer vector: the
  number of row and columns of the estimated target matrix,
  as well as the factorization rank (i.e. the number of
  basis components). }

  \item{dim}{\code{signature(x = "NMFfitXn")}: Returns the
  dimension common to all fits.

  Since all fits have the same dimensions, it returns the
  dimension of the first fit. This method returns
  \code{NULL} if the object is empty. }

  \item{nbasis}{\code{signature(x = "ANY")}: Default method
  which returns the number of columns of the basis matrix
  extracted from \code{x} using a suitable method
  \code{basis}.

  For NMF models, this also corresponds to the number of
  rows in the coefficient matrix. }

  \item{nbasis}{\code{signature(x = "matrix")}: For a
  matrix, the attribute 'nbasis' is returned, e.g. as
  attached to the consensus matrix returned by the method
  \code{\link{consensus}}.

  This is used internally to keep track of data about the
  parent fit and annotate plots. }

  \item{nbasis}{\code{signature(x = "NMFfitXn")}: Returns
  the number of basis components common to all fits.

  Since all fits have been computed using the same rank, it
  returns the factorization rank of the first fit. This
  method returns \code{NULL} if the object is empty. }

  }
}
\keyword{methods}

