\docType{methods}
\name{connectivity}
\alias{connectivity}
\alias{connectivity,ANY-method}
\alias{connectivity,factor-method}
\alias{connectivity-methods}
\alias{connectivity,NMF-method}
\alias{connectivity,numeric-method}
\alias{consensus}
\alias{consensus-methods}
\alias{consensus,NMFfitX-method}
\alias{consensus,NMF-method}
\title{Clustering Connectivity and Consensus Matrices}
\usage{
  connectivity(object, ...)

  \S4method{connectivity}{NMF}(object, no.attrib = FALSE)

  consensus(object, ...)
}
\arguments{
  \item{object}{an object with a suitable \code{\link{predict}} method.}

  \item{...}{extra arguments to allow extension.
They are passed to \code{\link{predict}}, except for the \code{vector} and
\code{factor} methods.}

  \item{no.attrib}{a logical that indicates if attributes containing information
about the NMF model should be attached to the result (\code{TRUE}) or not
(\code{FALSE}).}
}
\value{
  a square matrix of dimension the number of samples in the model, full
of 0s or 1s.
}
\description{
  \code{connectivity} is an S4 generic that computes the
  connectivity matrix based on the clustering of samples
  obtained from a model's \code{\link{predict}} method.

  The consensus matrix has been proposed by \cite{Brunet et
  al. (2004)} to help visualising and measuring the
  stability of the clusters obtained by NMF approaches. For
  objects of class \code{NMF} (e.g. results of a single NMF
  run, or NMF models), the consensus matrix reduces to the
  connectivity matrix.
}
\details{
  The connectivity matrix of a given partition of a set of
  samples (e.g. given as a cluster membership index) is the
  matrix \eqn{C} containing only 0 or 1 entries such that:
  \deqn{C_{ij} = \left\{\begin{array}{l} 1\mbox{ if sample
  }i\mbox{ belongs to the same cluster as sample }j\\
  0\mbox{ otherwise} \end{array}\right..}{ C_{ij} = 1 if
  sample i belongs to the same cluster as sample j, 0
  otherwise}
}
\section{Methods}{
\describe{

\item{connectivity}{\code{signature(object = "ANY")}: Default method which computes the connectivity matrix
using the result of \code{predict(x, ...)} as cluster membership index.
}

\item{connectivity}{\code{signature(object = "factor")}: Computes the connectivity matrix using \code{x} as cluster membership index.
}

\item{connectivity}{\code{signature(object = "numeric")}: Equivalent to \code{connectivity(as.factor(x))}.
}

\item{connectivity}{\code{signature(object = "NMF")}: Computes the connectivity matrix for an NMF model, for which cluster
membership is given by the most contributing basis component in each sample.
See \code{\link{predict,NMF-method}}.
}


\item{consensus}{\code{signature(object = "NMFfitX")}: Pure virtual method defined to ensure \code{consensus} is defined for sub-classes of \code{NMFfitX}.
It throws an error if called.
}

\item{consensus}{\code{signature(object = "NMF")}: This method is provided for completeness and is identical to
\code{\link{connectivity}}, and returns the connectivity matrix,
which, in the case of a single NMF model, is also the consensus matrix.
}

\item{consensus}{\code{signature(object = "NMFfitX1")}: The result is the matrix stored in slot \sQuote{consensus}.
This method returns \code{NULL} if the consensus matrix is empty.

See \code{\link{consensus,NMFfitX1-method}} for more details.
}

\item{consensus}{\code{signature(object = "NMFfitXn")}: This method returns \code{NULL} on an empty object.
The result is a matrix with several attributes attached, that are used by
plotting functions such as \code{\link{consensusmap}} to annotate the plots.

See \code{\link{consensus,NMFfitXn-method}} for more details.
}

}

}
\examples{
\dontshow{# roxygen generated flag
options(R_CHECK_RUNNING_EXAMPLES_=TRUE)
}

#----------
# connectivity,ANY-method
#----------
# clustering of random data
h <- hclust(dist(rmatrix(10,20)))
connectivity(stats::cutree(h, 2))

#----------
# connectivity,factor-method
#----------
connectivity(gl(2, 4))
}
\references{
  Brunet J, Tamayo P, Golub TR and Mesirov JP (2004).
  "Metagenes and molecular pattern discovery using matrix
  factorization." _Proceedings of the National Academy of
  Sciences of the United States of America_, *101*(12), pp.
  4164-9. ISSN 0027-8424, <URL:
  http://dx.doi.org/10.1073/pnas.0308531101>, <URL:
  http://www.ncbi.nlm.nih.gov/pubmed/15016911>.
}
\seealso{
  \code{\link{predict}}
}
\keyword{methods}

