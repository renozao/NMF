\docType{methods}
\name{predict}
\alias{predict}
\alias{predict-methods}
\alias{predict,NMFfitX-method}
\alias{predict,NMF-method}
\title{Clustering and Prediction}
\usage{
  predict(object, ...)

  \S4method{predict}{NMF}(object,
    what = c("columns", "rows", "samples", "features"),
    prob = FALSE, dmatrix = FALSE)

  \S4method{predict}{NMFfitX}(object,
    what = c("columns", "rows", "samples", "features", "consensus", "chc"),
    dmatrix = FALSE, ...)
}
\arguments{
  \item{object}{an NMF model}

  \item{what}{a character string that indicates the type of
  cluster membership should be returned: \sQuote{columns}
  or \sQuote{rows} for clustering the colmuns or the rows
  of the target matrix respectively. The values
  \sQuote{samples} and \sQuote{features} are aliases for
  \sQuote{colmuns} and \sQuote{rows} respectively.}

  \item{prob}{logical that indicates if the relative
  contributions of/to the dominant basis component should
  be computed and returned. See \emph{Details}.}

  \item{dmatrix}{logical that indicates if a dissimiliarity
  matrix should be attached to the result. This is notably
  used internally when computing NMF clustering
  silhouettes.}

  \item{...}{additional arguments affecting the predictions
  produced.}
}
\description{
  The methods \code{predict} for NMF models return the
  cluster membership of each sample or each feature.
  Currently the classification/prediction of new data is
  not implemented.
}
\details{
  The cluster membership is computed as the index of the
  dominant basis component for each sample
  (\code{what='samples' or 'columns'}) or each feature
  (\code{what='features' or 'rows'}), based on their
  corresponding entries in the coefficient matrix or basis
  matrix respectively.

  For example, if \code{what='samples'}, then the dominant
  basis component is computed for each column of the
  coefficient matrix as the row index of the maximum within
  the column.

  If argument \code{prob=FALSE} (default), the result is a
  \code{factor}. Otherwise a list with two elements is
  returned: element \code{predict} contains the cluster
  membership index (as a \code{factor}) and element
  \code{prob} contains the relative contribution of the
  dominant component to each sample (resp. the relative
  contribution of each feature to the dominant basis
  component):

  \itemize{ \item Samples: \deqn{p_j = x_{k_0} / \sum_k
  x_k}{p(j) = x(k0) / sum_k x(k)}, for each sample
  \eqn{1\leq j \leq p}, where \eqn{x_k}{x(k)} is the
  contribution of the \eqn{k}-th basis component to
  \eqn{j}-th sample (i.e. \code{H[k ,j]}), and
  \eqn{x_{k_0}}{x(k0)} is the maximum of these
  contributions.

  \item Features: \deqn{p_i = y_{k_0} / \sum_k y_k}{p(i) =
  y(k0) / sum_k y(k)}, for each feature \eqn{1\leq i \leq
  p}, where \eqn{y_k}{y(k)} is the contribution of the
  \eqn{k}-th basis component to \eqn{i}-th feature (i.e.
  \code{W[i, k]}), and \eqn{y_{k_0}}{y(k0)} is the maximum
  of these contributions.

  }
}
\section{Methods}{
  \describe{

  \item{predict}{\code{signature(object = "NMF")}: Default
  method for NMF models }

  \item{predict}{\code{signature(object = "NMFfitX")}:
  Returns the cluster membership index from an NMF model
  fitted with multiple runs.

  Besides the type of clustering available for any NMF
  models (\code{'columns', 'rows', 'samples', 'features'}),
  this method can return the cluster membership index based
  on the consensus matrix, computed from the multiple NMF
  runs.

  Argument \code{what} accepts the following extra types:
  \describe{ \item{\code{'chc'}}{ returns the cluster
  membership based on the hierarchical clustering of the
  consensus matrix, as performed by
  \code{\link{consensushc}}.} \item{\code{'consensus'}}{
  same as \code{'chc'} but the levels of the membership
  index are re-labeled to match the order of the clusters
  as they would be displayed on the associated dendrogram,
  as re-ordered on the default annotation track in
  consensus heatmap produced by
  \code{\link{consensusmap}}.} } }

  }
}
\examples{
# random target matrix
v <- rmatrix(20, 10)
# fit an NMF model
x <- nmf(v, 5)

# predicted column and row clusters
predict(x)
predict(x, 'rows')

# with relative contributions of each basis component
predict(x, prob=TRUE)
predict(x, 'rows', prob=TRUE)
}
\references{
  Brunet J, Tamayo P, Golub TR and Mesirov JP (2004).
  "Metagenes and molecular pattern discovery using matrix
  factorization." _Proceedings of the National Academy of
  Sciences of the United States of America_, *101*(12), pp.
  4164-9. ISSN 0027-8424, <URL:
  http://dx.doi.org/10.1073/pnas.0308531101>, <URL:
  http://www.ncbi.nlm.nih.gov/pubmed/15016911>.

  Pascual-Montano A, Carazo JM, Kochi K, Lehmann D and
  Pascual-marqui RD (2006). "Nonsmooth nonnegative matrix
  factorization (nsNMF)." _IEEE Trans. Pattern Anal. Mach.
  Intell_, *28*, pp. 403-415.
}
\keyword{methods}

