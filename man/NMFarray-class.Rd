\docType{class}
\name{NMFarray-class}
\alias{NMFarray-class}
\title{TRANSFAC Model - Multidimensional Nonnegative Matrix Factorization}
\description{
  This class provides data structure for NMF array models
  (e.g., TRANSFAC models), where the basis and/or the
  coefficient data is stored in a multidimensional array.
}
\section{Methods}{
  \describe{

  \item{.basis}{\code{signature(object = "NMFarray")}: This
  function returns slot \code{W} of \code{object}. }

  \item{.basis<-}{\code{signature(object = "NMFarray",
  value = "matrix")}: Replaces a slice of the basis array.
  }

  \item{.basis<-}{\code{signature(object = "NMFarray",
  value = "array")}: Set the basis array in NMF array
  models

  This function sets slot \code{W} of \code{object}. }

  \item{.coef}{\code{signature(object = "NMFarray")}: This
  function returns slot \code{H} of \code{object}. }

  \item{.coef<-}{\code{signature(object = "NMFarray", value
  = "array")}: Set the mixture coefficient matrix in
  standard NMF models

  This function sets slot \code{H} of \code{object}. }

  \item{.coef<-}{\code{signature(object = "NMFarray", value
  = "matrix")}: Replaces a slice of the coefficent array. }

  \item{fitted}{\code{signature(object = "NMFarray")}:
  Computes the target matrix estimate for an NMF array
  model. }

  \item{show}{\code{signature(object = "NMFarray")}: Show
  method for objects of class \code{NMFarray} }

  }
}
\examples{
\dontshow{# roxygen generated flag
options(R_CHECK_RUNNING_EXAMPLES_=TRUE)
}

# create a completely empty NMF object
new('NMFarray')

# create a NMF object based on random (compatible) arrays
n <- 4; r <- 3; p <- 5; q <- 2
w <- array(seq(n*r*q), dim = c(n, r, q))
h <- rmatrix(r, p)
new('NMFarray', W=w, H=h)
}
\seealso{
  Other NMF-model:
  \code{\link{initialize,NMFOffset-method}},
  \code{\link{NMFns-class}}, \code{\link{NMFOffset-class}},
  \code{\link{NMFstd-class}}
}

