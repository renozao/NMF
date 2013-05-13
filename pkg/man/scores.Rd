\docType{methods}
\name{featureScore}
\alias{extractFeatures}
\alias{extractFeatures,matrix-method}
\alias{extractFeatures-methods}
\alias{extractFeatures,NMF-method}
\alias{featureScore}
\alias{featureScore,matrix-method}
\alias{featureScore-methods}
\alias{featureScore,NMF-method}
\title{Feature Selection in NMF Models}
\usage{
  featureScore(object, ...)

  \S4method{featureScore}{matrix}(object,
    method = c("kim", "max"))

  extractFeatures(object, ...)

  \S4method{extractFeatures}{matrix}(object,
    method = c("kim", "max"),
    format = c("list", "combine", "subset"), nodups = TRUE)
}
\arguments{
  \item{object}{an object from which scores/features are
  computed/extracted}

  \item{...}{extra arguments to allow extension}

  \item{method}{scoring or selection method. It specifies
  the name of one of the method described in sections
  \emph{Feature scores} and \emph{Feature selection}.

  Additionally for \code{extractFeatures}, it may be an
  integer vector that indicates the number of top most
  contributing features to extract from each column of
  \code{object}, when ordered in decreasing order, or a
  numeric value between 0 and 1 that indicates the minimum
  relative basis contribution above which a feature is
  selected (i.e. basis contribution threshold). In the case
  of a single numeric value (integer or percentage), it is
  used for all columns.

  Note that \code{extractFeatures(x, 1)} means relative
  contribution threshold of 100\%, to select the top
  contributing features one must explicitly specify an
  integer value as in \code{extractFeatures(x, 1L)}.
  However, if all elements in methods are > 1, they are
  automatically treated as if they were integers:
  \code{extractFeatures(x, 2)} means the top-2 most
  contributing features in each component.}

  \item{format}{output format. The following values are
  accepted: \describe{ \item{\sQuote{list}}{(default)
  returns a list with one element per column in
  \code{object}, each containing the indexes of the
  selected features, as an integer vector. If \code{object}
  has row names, these are used to name each index vector.
  Components for which no feature were selected are
  assigned a \code{NA} value.}

  \item{\sQuote{combine}}{ returns all indexes in a single
  vector. Duplicated indexes are made unique if
  \code{nodups=TRUE} (default).}

  \item{\sQuote{subset}}{ returns an object of the same
  class as \code{object}, but subset with the selected
  indexes, so that it contains data only from
  basis-specific features.} }}

  \item{nodups}{logical that indicates if duplicated
  indexes, i.e. features selected on multiple basis
  components (which should in theory not happen), should be
  only appear once in the result. Only used when
  \code{format='combine'}.}
}
\value{
  \code{featureScore} returns a numeric vector of the
  length the number of rows in \code{object} (i.e. one
  score per feature).

  \code{extractFeatures} returns the selected features as a
  list of indexes, a single integer vector or an object of
  the same class as \code{object} that only contains the
  selected features.
}
\description{
  The function \code{featureScore} implements different
  methods to computes basis-specificity scores for each
  feature in the data.

  The function \code{extractFeatures} implements different
  methods to select the most basis-specific features of
  each basis component.
}
\details{
  One of the properties of Nonnegative Matrix Factorization
  is that is tend to produce sparse representation of the
  observed data, leading to a natural application to
  bi-clustering, that characterises groups of samples by a
  small number of features.

  In NMF models, samples are grouped according to the basis
  components that contributes the most to each sample, i.e.
  the basis components that have the greatest coefficient
  in each column of the coefficient matrix (see
  \code{\link{predict,NMF-method}}). Each group of samples
  is then characterised by a set of features selected based
  on basis-specifity scores that are computed on the basis
  matrix.
}
\section{Methods}{
  \describe{

  \item{extractFeatures}{\code{signature(object =
  "matrix")}: Select features on a given matrix, that
  contains the basis component in columns. }

  \item{extractFeatures}{\code{signature(object = "NMF")}:
  Select basis-specific features from an NMF model, by
  applying the method \code{extractFeatures,matrix} to its
  basis matrix. }

  \item{featureScore}{\code{signature(object = "matrix")}:
  Computes feature scores on a given matrix, that contains
  the basis component in columns. }

  \item{featureScore}{\code{signature(object = "NMF")}:
  Computes feature scores on the basis matrix of an NMF
  model. }

  }
}
\section{Feature scores}{
  The function \code{featureScore} can compute
  basis-specificity scores using the following methods:

  \describe{

  \item{\sQuote{kim}}{ Method defined by \cite{Kim et al. (2007)}.

  The score for feature \eqn{i} is defined as: \deqn{S_i =
  1 + \frac{1}{\log_2 k} \sum_{q=1}^k p(i,q) \log_2
  p(i,q)}{ S_i = 1 + 1/log2(k) sum_q [ p(i,q) log2( p(i,q)
  ) ] },

  where \eqn{p(i,q)} is the probability that the \eqn{i}-th
  feature contributes to basis \eqn{q}: \deqn{p(i,q) =
  \frac{W(i,q)}{\sum_{r=1}^k W(i,r)} }{ p(i,q) = W(i,q) /
  (sum_r W(i,r)) }

  The feature scores are real values within the range
  [0,1]. The higher the feature score the more
  basis-specific the corresponding feature. }

  \item{\sQuote{max}}{Method defined by
  \cite{Carmona-Saez et al. (2006)}.

  The feature scores are defined as the row maximums. }

  }
}

\section{Feature selection}{
  The function \code{extractFeatures} can select features
  using the following methods: \describe{
  \item{\sQuote{kim}}{ uses \cite{Kim et al. (2007)} scoring schema
  and feature selection method.

  The features are first scored using the function
  \code{featureScore} with method \sQuote{kim}. Then only
  the features that fulfil both following criteria are
  retained:

  \itemize{ \item score greater than \eqn{\hat{\mu} + 3
  \hat{\sigma}}, where \eqn{\hat{\mu}} and
  \eqn{\hat{\sigma}} are the median and the median absolute
  deviation (MAD) of the scores respectively;

  \item the maximum contribution to a basis component is
  greater than the median of all contributions (i.e. of all
  elements of W). }

  }

  \item{\sQuote{max}}{ uses the selection method used in
  the \code{bioNMF} software package and described in
  \cite{Carmona-Saez et al. (2006)}.

  For each basis component, the features are first sorted
  by decreasing contribution. Then, one selects only the
  first consecutive features whose highest contribution in
  the basis matrix is effectively on the considered basis.
  }

  }
}
\examples{
# random NMF model
x <- rnmf(3, 50,20)

# probably no feature is selected
extractFeatures(x)
# extract top 5 for each basis
extractFeatures(x, 5L)
# extract features that have a relative basis contribution above a threshold
extractFeatures(x, 0.5)
# ambiguity?
extractFeatures(x, 1) # means relative contribution above 100\%
extractFeatures(x, 1L) # means top contributing feature in each component
}
\references{
  Carmona-Saez P, Pascual-Marqui RD, Tirado F, Carazo JM
  and Pascual-Montano A (2006). "Biclustering of gene
  expression data by Non-smooth Non-negative Matrix
  Factorization." _BMC bioinformatics_, *7*, pp. 78. ISSN
  1471-2105, <URL:
  http://dx.doi.org/10.1186/1471-2105-7-78>, <URL:
  http://www.ncbi.nlm.nih.gov/pubmed/16503973>.

  Kim H and Park H (2007). "Sparse non-negative matrix
  factorizations via alternating non-negativity-constrained
  least squares for microarray data analysis."
  _Bioinformatics (Oxford, England)_, *23*(12), pp.
  1495-502. ISSN 1460-2059, <URL:
  http://dx.doi.org/10.1093/bioinformatics/btm134>, <URL:
  http://www.ncbi.nlm.nih.gov/pubmed/17483501>.
}
\keyword{methods}

