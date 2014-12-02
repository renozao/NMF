\docType{class}
\name{NMFstd-class}
\alias{NMFstd-class}
\title{NMF Model - Standard model}
\description{
  This class implements the standard model of Nonnegative Matrix
Factorization.
It provides a general structure and generic functions to manage
factorizations that follow the standard NMF model, as defined by
\cite{Lee et al. (2001)}.
}
\details{
  Let \eqn{V} be a \eqn{n \times m} non-negative matrix and \eqn{r} a positive
integer.  In its standard form (see references below), a NMF of \eqn{V} is
commonly defined as a pair of matrices \eqn{(W, H)} such that:

\deqn{V \equiv W H,}

where:
\itemize{
\item \eqn{W} and \eqn{H} are \eqn{n \times r} and \eqn{r
\times m} matrices respectively with non-negative entries;
\item \eqn{\equiv} is to be understood with respect to some loss function.
Common choices of loss functions are based on Frobenius norm or Kullback-Leibler
divergence.
}

Integer \eqn{r} is called the \emph{factorization rank}.
Depending on the context of application of NMF, the columns of \eqn{W}
and \eqn{H} are given different names:
\describe{
\item{columns of \code{W}}{basis vector, metagenes, factors, source, image basis}
\item{columns of \code{H}}{mixture coefficients, metagene sample expression profiles, weights}
\item{rows of \code{H}}{basis profiles, metagene expression profiles}
}

NMF approaches have been successfully applied to several fields.
The package NMF was implemented trying to use names as generic as possible
for objects and methods.

The following terminology is used:
\describe{
\item{samples}{the columns of the target matrix \eqn{V}}
\item{features}{the rows of the target matrix \eqn{V}}
\item{basis matrix}{the first matrix factor \eqn{W}}
\item{basis vectors}{the columns of first matrix factor \eqn{W}}
\item{mixture matrix}{the second matrix factor \eqn{H}} \item{mixtures
coefficients}{the columns of second matrix factor \eqn{H}}
}

However, because the package NMF was primarily implemented to work with gene
expression microarray data, it also provides a layer to easily and
intuitively work with objects from the Bioconductor base framework.
See \link{bioc-NMF} for more details.
}
\section{Slots}{
\describe{

\item{W}{A \code{matrix} that contains the basis matrix, i.e. the \emph{first}
matrix factor of the factorisation}

\item{H}{A \code{matrix} that contains the coefficient matrix, i.e. the
\emph{second} matrix factor of the factorisation}

\item{bterms}{a \code{data.frame} that contains the primary data that
define fixed basis terms. See \code{\link{bterms}}.}

\item{ibterms}{integer vector that contains the indexes of the basis components
that are fixed, i.e. for which only the coefficient are estimated.

IMPORTANT: This slot is set on construction of an NMF model via
\code{\link[=nmfModel,formula,ANY-method]{nmfModel}} and is not recommended to
not be subsequently changed by the end-user.}

\item{cterms}{a \code{data.frame} that contains the primary data that
define fixed coefficient terms. See \code{\link{cterms}}.}

\item{icterms}{integer vector that contains the indexes of the basis components
that have fixed coefficients, i.e. for which only the basis vectors are estimated.

IMPORTANT: This slot is set on construction of an NMF model via
\code{\link[=nmfModel,formula,ANY-method]{nmfModel}} and is not recommended to
not be subsequently changed by the end-user.}

}

}
\section{Methods}{
\describe{

\item{.basis}{\code{signature(object = "NMFstd")}: Get the basis matrix in standard NMF models

This function returns slot \code{W} of \code{object}.
}


\item{.basis<-}{\code{signature(object = "NMFstd", value = "array")}: Set the basis matrix in standard NMF models

This function sets slot \code{W} of \code{object}.
}

\item{.basis<-}{\code{signature(object = "NMFstd", value = "matrix")}: Replaces a slice of the basis array.
}


\item{bterms<-}{\code{signature(object = "NMFstd")}: Default method tries to coerce \code{value} into a \code{data.frame} with
\code{\link{as.data.frame}}.
}


\item{.coef}{\code{signature(object = "NMFstd")}: Get the mixture coefficient matrix in standard NMF models

This function returns slot \code{H} of \code{object}.
}


\item{.coef<-}{\code{signature(object = "NMFstd", value = "array")}: Set the mixture coefficient matrix in standard NMF models

This function sets slot \code{H} of \code{object}.
}

\item{.coef<-}{\code{signature(object = "NMFstd", value = "matrix")}: Replaces a slice of the coefficent array.
}


\item{cterms<-}{\code{signature(object = "NMFstd")}: Default method tries to coerce \code{value} into a \code{data.frame} with
\code{\link{as.data.frame}}.
}


\item{fitted}{\code{signature(object = "NMFstd")}: Compute the target matrix estimate in \emph{standard NMF models}.

The estimate matrix is computed as the product of the two matrix slots
\code{W} and \code{H}:
\deqn{\hat{V} = W H}{V ~ W H}
}


\item{ibterms}{\code{signature(object = "NMFstd")}: Method for standard NMF models, which returns the integer vector that is
stored in slot \code{ibterms} when a formula-based NMF model is instantiated.
}


\item{icterms}{\code{signature(object = "NMFstd")}: Method for standard NMF models, which returns the integer vector that is
stored in slot \code{icterms} when a formula-based NMF model is instantiated.
}

}

}
\examples{
\dontshow{# roxygen generated flag
options(R_CHECK_RUNNING_EXAMPLES_=TRUE)
}

# create a completely empty NMFstd object
new('NMFstd')

# create a NMF object based on one random matrix: the missing matrix is deduced
# Note this only works when using factory method NMF
n <- 50; r <- 3;
w <- rmatrix(n, r)
nmfModel(W=w)

# create a NMF object based on random (compatible) matrices
p <- 20
h <- rmatrix(r, p)
nmfModel(W=w, H=h)

# create a NMF object based on incompatible matrices: generate an error
h <- rmatrix(r+1, p)
try( new('NMFstd', W=w, H=h) )
try( nmfModel(w, h) )

# Giving target dimensions to the factory method allow for coping with dimension
# incompatibilty (a warning is thrown in such case)
nmfModel(r, W=w, H=h)

# create a NMF array object based on random (compatible) arrays
# extra dimension (levels)
q <- 2
w <- array(seq(n*r*q), dim = c(n, r, q))
h <- rmatrix(r, p)
nmfModel(W = w, H = h)
}
\references{
  Lee DD and Seung H (2001). "Algorithms for non-negative matrix factorization." _Advances in neural information processing
systems_. <URL: http://scholar.google.com/scholar?q=intitle:Algorithms+for+non-negative+matrix+factorization\#0>.
}
\seealso{
  Other NMF-model: \code{\link{initialize,NMFOffset-method}}, \code{\link{NMFns-class}}, \code{\link{NMFOffset-class}}
}

