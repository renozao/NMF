\docType{methods}
\name{basiscor}
\alias{basiscor}
\alias{basiscor,matrix,NMF-method}
\alias{basiscor-methods}
\alias{basiscor,NMF,matrix-method}
\alias{basiscor,NMF,missing-method}
\alias{basiscor,NMF,NMF-method}
\alias{profcor}
\alias{profcor,matrix,NMF-method}
\alias{profcor-methods}
\alias{profcor,NMF,matrix-method}
\alias{profcor,NMF,missing-method}
\alias{profcor,NMF,NMF-method}
\title{Correlations in NMF Models}
\usage{
  basiscor(x, y, ...)

  profcor(x, y, ...)
}
\arguments{
  \item{x}{a matrix or an object with suitable methods \code{\link{basis}}
or \code{\link{coef}}.}

  \item{y}{a matrix or an object with suitable methods \code{\link{basis}}
or \code{\link{coef}}, and dimensions compatible with \code{x}.
If missing the correlations are computed between \code{x} and \code{y=x}.}

  \item{...}{extra arguments passed to \code{\link{cor}}.}
}
\description{
  \code{basiscor} computes the correlation matrix between basis vectors, i.e.
the \emph{columns} of its basis matrix -- which is the model's first matrix factor.

\code{profcor} computes the correlation matrix between basis profiles,
i.e. the \emph{rows} of the coefficient matrix -- which is the model's second
matrix factor.
}
\details{
  Each generic has methods defined for computing correlations between NMF models
and/or compatible matrices.
The computation is performed by the base function \code{\link{cor}}.
}
\section{Methods}{
\describe{

\item{basiscor}{\code{signature(x = "NMF", y = "matrix")}: Computes the correlations between the basis vectors of \code{x} and
the columns of \code{y}.
}

\item{basiscor}{\code{signature(x = "matrix", y = "NMF")}: Computes the correlations between the columns of \code{x}
and the the basis vectors of \code{y}.
}

\item{basiscor}{\code{signature(x = "NMF", y = "NMF")}: Computes the correlations between the basis vectors of \code{x} and \code{y}.
}

\item{basiscor}{\code{signature(x = "NMF", y = "missing")}: Computes the correlations between the basis vectors of \code{x}.
}


\item{profcor}{\code{signature(x = "NMF", y = "matrix")}: Computes the correlations between the basis profiles of \code{x} and
the rows of \code{y}.
}

\item{profcor}{\code{signature(x = "matrix", y = "NMF")}: Computes the correlations between the rows of \code{x} and the basis
profiles of \code{y}.
}

\item{profcor}{\code{signature(x = "NMF", y = "NMF")}: Computes the correlations between the basis profiles of \code{x} and \code{y}.
}

\item{profcor}{\code{signature(x = "NMF", y = "missing")}: Computes the correlations between the basis profiles of \code{x}.
}

}

}
\examples{
\dontshow{# roxygen generated flag
options(R_CHECK_RUNNING_EXAMPLES_=TRUE)
}

# generate two random NMF models
a <- rnmf(3, 100, 20)
b <- rnmf(3, 100, 20)

# Compute auto-correlations
basiscor(a)
profcor(a)
# Compute correlations with b
basiscor(a, b)
profcor(a, b)

# try to recover the underlying NMF model 'a' from noisy data
res <- nmf(fitted(a) + rmatrix(a), 3)

# Compute correlations with the true model
basiscor(a, res)
profcor(a, res)

# Compute correlations with a random compatible matrix
W <- rmatrix(basis(a))
basiscor(a, W)
identical(basiscor(a, W), basiscor(W, a))

H <- rmatrix(coef(a))
profcor(a, H)
identical(profcor(a, H), profcor(H, a))
}
\keyword{methods}

