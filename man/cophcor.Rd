\docType{methods}
\name{cophcor}
\alias{cophcor}
\alias{cophcor,matrix-method}
\alias{cophcor-methods}
\alias{cophcor,NMFfitX-method}
\title{Cophenetic Correlation Coefficient}
\usage{
  cophcor(object, ...)

  \S4method{cophcor}{matrix}(object, linkage = "average")
}
\arguments{
  \item{object}{an object from which is extracted a consensus matrix.}

  \item{...}{extra arguments to allow extension and passed to subsequent calls.}

  \item{linkage}{linkage method used in the hierarchical clustering.
It is passed to \code{\link{hclust}}.}
}
\description{
  The function \code{cophcor} computes the cophenetic correlation coefficient
from consensus matrix \code{object}, e.g. as obtained from multiple NMF runs.
}
\details{
  The cophenetic correlation coeffificient is based on the consensus matrix
(i.e. the average of connectivity matrices) and was proposed by
\cite{Brunet et al. (2004)} to measure the stability of the clusters obtained from NMF.

It is defined as the Pearson correlation between the samples' distances
induced by the consensus matrix (seen as a similarity matrix) and their
cophenetic distances from a hierachical clustering based on these very
distances (by default an average linkage is used).
See \cite{Brunet et al. (2004)}.
}
\section{Methods}{
\describe{

\item{cophcor}{\code{signature(object = "matrix")}: Workhorse method for matrices.
}

\item{cophcor}{\code{signature(object = "NMFfitX")}: Computes the cophenetic correlation coefficient on the consensus matrix
of \code{object}.
All arguments in \code{...} are passed to the method \code{cophcor,matrix}.
}

}

}
\references{
  Brunet J, Tamayo P, Golub TR and Mesirov JP (2004). "Metagenes and molecular pattern discovery using matrix factorization."
_Proceedings of the National Academy of Sciences of the United States of America_, *101*(12), pp. 4164-9. ISSN 0027-8424,
<URL: http://dx.doi.org/10.1073/pnas.0308531101>, <URL: http://www.ncbi.nlm.nih.gov/pubmed/15016911>.
}
\seealso{
  \code{\link{cophenetic}}
}
\keyword{methods}

