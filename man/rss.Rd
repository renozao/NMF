\docType{methods}
\name{rss}
\alias{evar}
\alias{evar,ANY-method}
\alias{evar-methods}
\alias{rss}
\alias{rss,ANY-method}
\alias{rss,matrix-method}
\alias{rss-methods}
\title{Residual Sum of Squares and Explained Variance}
\usage{
  rss(object, ...)

  \S4method{rss}{matrix}(object, target)

  evar(object, ...)

  \S4method{evar}{ANY}(object, target, ...)
}
\arguments{
  \item{object}{an R object with a suitable \code{\link{fitted}}, \code{rss} or
\code{evar} method.}

  \item{...}{extra arguments to allow extension, e.g. passed to \code{rss}
in \code{evar} calls.}

  \item{target}{target matrix}
}
\value{
  a single numeric value
}
\description{
  \code{rss} and \code{evar} are S4 generic functions that respectively computes
the Residual Sum of Squares (RSS) and explained variance achieved by a model.

The explained variance for a target \eqn{V} is computed as:
\deqn{evar = 1 - \frac{RSS}{\sum_{i,j} v_{ij}^2} }{evar = 1 - RSS/sum v_{ij}^2},
}
\details{
  where RSS is the residual sum of squares.

The explained variance is usefull to compare the performance of different
models and their ability to accurately reproduce the original target matrix.
Note, however, that a possible caveat is that some models explicitly aim at
minimizing the RSS (i.e. maximizing the explained variance), while others do not.
}
\section{Methods}{
\describe{

\item{evar}{\code{signature(object = "ANY")}: Default method for \code{evar}.

It requires a suitable \code{rss} method to be defined
for \code{object}, as it internally calls \code{rss(object, target, ...)}.
}


\item{rss}{\code{signature(object = "matrix")}: Computes the RSS between a target matrix and its estimate \code{object},
which must be a matrix of the same dimensions as \code{target}.

The RSS between a target matrix \eqn{V} and its estimate \eqn{v} is computed as:
\deqn{RSS = \sum_{i,j} (v_{ij} - V_{ij})^2}

Internally, the computation is performed using an optimised C++ implementation,
that is light in memory usage.
}

\item{rss}{\code{signature(object = "ANY")}: Residual sum of square between a given target matrix and a model that has a
suitable \code{\link{fitted}} method.
It is equivalent to \code{rss(fitted(object), ...)}

In the context of NMF, \cite{Hutchins et al. (2008)} used the variation of the RSS
in combination with the algorithm from \cite{Lee et al. (1999)} to estimate the
correct number of basis vectors.
The optimal rank is chosen where the graph of the RSS first shows an inflexion
point, i.e. using a screeplot-type criterium.
See section \emph{Rank estimation} in \code{\link{nmf}}.

Note that this way of estimation may not be suitable for all models.
Indeed, if the NMF optimisation problem is not based on the Frobenius norm,
the RSS is not directly linked to the quality of approximation of the NMF model.
However, it is often the case that it still decreases with the rank.
}

}

}
\examples{
\dontshow{# roxygen generated flag
options(R_CHECK_RUNNING_EXAMPLES_=TRUE)
}

#----------
# rss,matrix-method
#----------
# RSS bewteeen random matrices
x <- rmatrix(20,10, max=50)
y <- rmatrix(20,10, max=50)
rss(x, y)
rss(x, x + rmatrix(x, max=0.1))

#----------
# rss,ANY-method
#----------
# RSS between an NMF model and a target matrix
x <- rmatrix(20, 10)
y <- rnmf(3, x) # random compatible model
rss(y, x)

# fit a model with nmf(): one should do better
y2 <- nmf(x, 3) # default minimizes the KL-divergence
rss(y2, x)
y2 <- nmf(x, 3, 'lee') # 'lee' minimizes the RSS
rss(y2, x)
}
\references{
  Hutchins LN, Murphy SM, Singh P and Graber JH (2008). "Position-dependent motif characterization using non-negative matrix
factorization." _Bioinformatics (Oxford, England)_, *24*(23), pp. 2684-90. ISSN 1367-4811, <URL:
http://dx.doi.org/10.1093/bioinformatics/btn526>, <URL: http://www.ncbi.nlm.nih.gov/pubmed/18852176>.

Lee DD and Seung HS (1999). "Learning the parts of objects by non-negative matrix factorization." _Nature_, *401*(6755), pp.
788-91. ISSN 0028-0836, <URL: http://dx.doi.org/10.1038/44565>, <URL: http://www.ncbi.nlm.nih.gov/pubmed/10548103>.
}
\keyword{methods}

