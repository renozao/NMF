\docType{methods}
\name{[,NMF-method}
\alias{[,NMF-method}
\title{Sub-setting NMF Objects}
\usage{
  \S4method{[}{NMF}(x, i, j, ..., drop = FALSE)
}
\arguments{
  \item{i}{index used to subset on the \strong{rows} of the basis matrix (i.e.
the features).
It can be a \code{numeric}, \code{logical}, or \code{character} vector
(whose elements must match the row names of \code{x}).
In the case of a \code{logical} vector the entries are recycled if necessary.}

  \item{j}{index used to subset on the \strong{columns} of the mixture
coefficient matrix (i.e. the samples).
It can be a \code{numeric}, \code{logical}, or \code{character} vector
(whose elements must match the column names of \code{x}).
In the case of a \code{logical} vector the entries are recycled if necessary.}

  \item{...}{used to specify a third index to subset on the basis components,
i.e. on both the columns and rows of the basis matrix and mixture
coefficient respectively.
It can be a \code{numeric}, \code{logical}, or \code{character} vector
(whose elements must match the basis names of \code{x}).
In the case of a \code{logical} vector the entries are recycled if necessary.

Note that only the first extra subset index is used.
A warning is thrown if more than one extra argument is passed in \code{...}.}

  \item{drop}{single \code{logical} value used to drop the \code{NMF-class}
wrapping and only return subsets of one of the factor matrices (see \emph{Details})}

  \item{x}{
    object from which to extract element(s) or in which to replace element(s).
  }
}
\description{
  This method provides a convenient way of sub-setting objects of class \code{NMF},
using a matrix-like syntax.

It allows to consistently subset one or both matrix factors in the NMF model, as well
as retrieving part of the basis components or part of the mixture coefficients with
a reduced amount of code.
}
\details{
  The returned value depends on the number of subset index passed and the
value of argument \code{drop}:

\itemize{ \item No index as in \code{x[]} or \code{x[,]}: the value is the
object \code{x} unchanged.

\item One single index as in \code{x[i]}: the value is the complete NMF
model composed of the selected basis components, subset by \code{i},
except if argument \code{drop=TRUE}, or if it is missing and \code{i} is of length 1.
Then only the basis matrix is returned with dropped dimensions:
\code{x[i, drop=TRUE]} <=> \code{drop(basis(x)[, i])}.

This means for example that \code{x[1L]} is the first basis vector,
and \code{x[1:3, drop = TRUE]} is the matrix composed of the 3 first basis vectors -- in columns.

Note that in version <= 0.18.3, the call \code{x[i, drop = TRUE.or.FALSE]} was equivalent to
\code{basis(x)[, i, drop=TRUE.or.FALSE]}.

\item More than one index with \code{drop=FALSE} (default) as in
\code{x[i,j]}, \code{x[i,]}, \code{x[,j]}, \code{x[i,j,k]}, \code{x[i,,k]},
etc...: the value is a \code{NMF} object whose basis and/or mixture
coefficient matrices have been subset accordingly. The third index \code{k}
affects simultaneously the columns of the basis matrix AND the rows of the
mixture coefficient matrix. In this case argument \code{drop} is not used.

\item More than one index with \code{drop=TRUE} and \code{i} xor \code{j}
missing: the value returned is the matrix that is the more affected by the
subset index. That is that \code{x[i, , drop=TRUE]} and \code{x[i, , k,
drop=TRUE]} return the basis matrix subset by \code{[i,]} and \code{[i,k]}
respectively, while \code{x[, j, drop=TRUE]} and \code{x[, j, k, drop=TRUE]}
return the mixture coefficient matrix subset by \code{[,j]} and \code{[k,j]}
respectively.

}
}
\examples{
\dontshow{# roxygen generated flag
options(R_CHECK_RUNNING_EXAMPLES_=TRUE)
}

# create a dummy NMF object that highlight the different way of subsetting
a <- nmfModel(W=outer(seq(1,5),10^(0:2)), H=outer(10^(0:2),seq(-1,-10)))
basisnames(a) <- paste('b', 1:nbasis(a), sep='')
rownames(a) <- paste('f', 1:nrow(a), sep='')
colnames(a) <- paste('s', 1:ncol(a), sep='')

# or alternatively:
# dimnames(a) <- list( features=paste('f', 1:nrow(a), sep='')
#					, samples=paste('s', 1:ncol(a), sep='')
#					, basis=paste('b', 1:nbasis(a)) )

# look at the resulting NMF object
a
basis(a)
coef(a)

# extract basis components
a[1]
a[1, drop=FALSE] # not dropping matrix dimension
a[2:3]

# subset on the features
a[1,]
a[2:4,]
# dropping the NMF-class wrapping => return subset basis matrix
a[2:4,, drop=TRUE]

# subset on the samples
a[,1]
a[,2:4]
# dropping the NMF-class wrapping => return subset coef matrix
a[,2:4, drop=TRUE]

# subset on the basis => subsets simultaneously basis and coef matrix
a[,,1]
a[,,2:3]
a[4:5,,2:3]
a[4:5,,2:3, drop=TRUE] # return subset basis matrix
a[,4:5,2:3, drop=TRUE] # return subset coef matrix

# 'drop' has no effect here
a[,,2:3, drop=TRUE]
}
\keyword{methods}

