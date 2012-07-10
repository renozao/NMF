\docType{methods}
\name{dispersion}
\alias{dispersion}
\alias{dispersion,NMFfitX-method}
\alias{dispersion,matrix-method}
\alias{dispersion-methods}
\title{Dispersion of a Matrix}
\usage{
  dispersion(object, ...)
}
\arguments{
  \item{object}{an object from which the dispersion is
  computed}

  \item{...}{extra arguments to allow extension}
}
\description{
  Computes the dispersion coefficient of a -- consensus --
  matrix \code{object}, generally obtained from multiple
  NMF runs.
}
\details{
  The dispersion coefficient is based on the consensus
  matrix (i.e. the average of connectivity matrices) and
  was proposed by \cite{Kim et al. (2007)} to measure the
  reproducibility of the clusters obtained from NMF.

  It is defined as: \deqn{\rho = \sum_{i,j=1}^n 4 (C_{ij} -
  \frac{1}{2})^2 , } where \eqn{n} is the total number of
  samples.

  By construction, \eqn{0 \leq \rho \leq 1} and \eqn{\rho =
  1} only for a perfect consensus matrix, where all entries
  0 or 1. A perfect consensus matrix is obtained only when
  all the connectivity matrices are the same, meaning that
  the algorithm gave the same clusters at each run. See
  \cite{Kim et al. (2007)}.
}
\section{Methods}{
  \describe{

  \item{dispersion}{\code{signature(object = "matrix")}:
  Workhorse method that computes the dispersion on a given
  matrix. }

  \item{dispersion}{\code{signature(object = "NMFfitX")}:
  Computes the dispersion on the consensus matrix obtained
  from multiple NMF runs. }

  }
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
\keyword{methods}

