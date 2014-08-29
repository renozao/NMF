\name{t.NMF}
\alias{t.NMF}
\title{Transformation NMF Model Objects}
\usage{
  \method{t}{NMF} (x)
}
\arguments{
  \item{x}{NMF model object.}
}
\description{
  \code{t} transpose an NMF model, by transposing and
  swapping its basis and coefficient matrices:
  \eqn{t([W,H]) = [t(H), t(W)]}.
}
\details{
  The function \code{t} is a generic defined in the
  \pkg{base} package. The method \code{t.NMF} defines the
  trasnformation for the general NMF interface. This method
  may need to be overloaded for NMF models, whose structure
  requires specific handling.
}
\examples{
\dontshow{# roxygen generated flag
options(R_CHECK_RUNNING_EXAMPLES_=TRUE)
}

x <- rnmf(3, 100, 20)
x
# transpose
y <- t(x)
y

# factors are swapped-transposed
stopifnot( identical(basis(y), t(coef(x))) )
stopifnot( identical(coef(y), t(basis(x))) )
}
\seealso{
  Other transforms: \code{\link{nneg}}, \code{\link{posneg}}, \code{\link{rposneg}}
}

