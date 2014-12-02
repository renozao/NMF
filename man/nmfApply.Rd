\name{nmfApply}
\alias{nmfApply}
\title{Apply Function for NMF Objects}
\usage{
  nmfApply(X, MARGIN, FUN, ..., simplify = TRUE,
    USE.NAMES = TRUE)
}
\arguments{
  \item{X}{an object that has suitable \code{\link{basis}} and \code{coef} methods,
e.g. an NMF model.}

  \item{MARGIN}{a single numeric (integer) value that specifies over which margin(s)
the function \code{FUN} is applied.
See section \emph{Details} for a list of possible values.}

  \item{FUN}{a function to apply over the specified margins.}

  \item{...}{extra arguments passed to \code{FUN}}

  \item{simplify}{a logical only used when \code{MARGIN=3}, that indicates if \code{sapply}
should try to simplify result if possible.
Since this argument follows \sQuote{...} its name cannot be abbreviated.}

  \item{USE.NAMES}{a logical only used when \code{MARGIN=3}, that indicates if \code{sapply}
should use the names of the basis components to name the results if present.
Since this argument follows \sQuote{...} its name cannot be abbreviated.}
}
\value{
  a vector or a list.
See \code{\link[base]{apply}} and \code{\link[base]{sapply}} for more details on
the output format.
}
\description{
  The function \code{nmfApply} provides exteneded \code{apply}-like
functionality for objects of class \code{NMF}.
It enables to easily apply a function over different margins of
NMF models.
}
\details{
  The function \code{FUN} is applied via a call to \code{\link{apply}}
or \code{\link{sapply}} according to the value of argument \code{MARGIN}
as follows:

\describe{
\item{MARGIN=1}{ apply \code{FUN} to each \emph{row} of the basis matrix:
\code{apply(basis(X), 1L, FUN, ...)}.}

\item{MARGIN=2}{ apply \code{FUN} to each \emph{column} of the coefficient matrix:
\code{apply(coef(X), 2L, FUN, ...)}.}

\item{MARGIN=3}{ apply \code{FUN} to each \emph{pair} of associated basis component
and basis profile:
more or less \code{sapply(seq(nbasis(X)), function(i, ...) FUN(basis(X)[,i], coef(X)[i, ], ...), ...)}.

In this case \code{FUN} must be have at least two arguments, to which are passed
each basis components and basis profiles respectively -- as numeric vectors.}

\item{MARGIN=4}{ apply \code{FUN} to each \emph{column} of the basis matrix, i.e. to each
basis component:
\code{apply(basis(X), 2L, FUN, ...)}.}

\item{MARGIN=5}{ apply \code{FUN} to each \emph{row} of the coefficient matrix:
\code{apply(coef(X), 1L, FUN, ...)}.}

}
}

