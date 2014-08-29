\name{nmfAlgorithm.SNMF_R}
\alias{nmfAlgorithm.SNMF_L}
\alias{nmfAlgorithm.SNMF_R}
\alias{SNMF/L-nmf}
\alias{SNMF/R-nmf}
\title{NMF Algorithm - Sparse NMF via Alternating NNLS}
\usage{
  nmfAlgorithm.SNMF_R(..., maxIter = 20000L, eta = -1,
    beta = 0.01, bi_conv = c(0, 10), eps_conv = 1e-04)

  nmfAlgorithm.SNMF_L(..., maxIter = 20000L, eta = -1,
    beta = 0.01, bi_conv = c(0, 10), eps_conv = 1e-04)
}
\arguments{
  \item{maxIter}{maximum number of iterations.}

  \item{eta}{parameter to suppress/bound the L2-norm of \code{W} and in
\code{H} in \sQuote{SNMF/R} and \sQuote{SNMF/L} respectively.

If \code{eta < 0}, then it is set to the maximum value in the target matrix is used.}

  \item{beta}{regularisation parameter for sparsity control, which
balances the trade-off between the accuracy of the approximation and the
sparseness of \code{H} and \code{W} in \sQuote{SNMF/R} and \sQuote{SNMF/L} respectively.

Larger beta generates higher sparseness on \code{H} (resp. \code{W}).
Too large beta is not recommended.}

  \item{bi_conv}{parameter of the biclustering convergence test.
It must be a size 2 numeric vector \code{bi_conv=c(wminchange, iconv)},
with:
\describe{
\item{\code{wminchange}:}{the minimal allowance of change in row-clusters.}
\item{\code{iconv}:}{ decide convergence if row-clusters
(within the allowance of \code{wminchange})
and column-clusters have not changed for \code{iconv} convergence checks.}
}

Convergence checks are performed every 5 iterations.}

  \item{eps_conv}{threshold for the KKT convergence test.}

  \item{...}{extra argument not used.}
}
\description{
  NMF algorithms proposed by \cite{Kim et al. (2007)} that
  enforces sparsity constraint on the basis matrix
  (algorithm \sQuote{SNMF/L}) or the mixture coefficient
  matrix (algorithm \sQuote{SNMF/R}).
}
\details{
  The algorithm \sQuote{SNMF/R} solves the following NMF
  optimization problem on a given target matrix \eqn{A} of
  dimension \eqn{n \times p}{n x p}: \deqn{
  \begin{array}{ll} & \min_{W,H} \frac{1}{2} \left(|| A -
  WH ||_F^2 + \eta ||W||_F^2 + \beta (\sum_{j=1}^p
  ||H_{.j}||_1^2)\right)\\ s.t. & W\geq 0, H\geq 0
  \end{array} }{ min_{W,H} 1/2 (|| A - WH ||_F^2 + eta
  ||W||_F^2 + beta (sum_j ||H[,j]||_1^2))

  s.t. W>=0, H>=0 }

  The algorithm \sQuote{SNMF/L} solves a similar problem on
  the transposed target matrix \eqn{A}, where \eqn{H} and
  \eqn{W} swap roles, i.e. with sparsity constraints
  applied to \code{W}.
}
\references{
  Kim H and Park H (2007). "Sparse non-negative matrix
  factorizations via alternating non-negativity-constrained
  least squares for microarray data analysis."
  _Bioinformatics (Oxford, England)_, *23*(12), pp.
  1495-502. ISSN 1460-2059, <URL:
  http://dx.doi.org/10.1093/bioinformatics/btm134>, <URL:
  http://www.ncbi.nlm.nih.gov/pubmed/17483501>.
}

