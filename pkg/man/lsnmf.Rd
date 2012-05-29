\name{nmf.update.lsnmf}
\alias{lsnmf-nmf}
\alias{nmf.update.lsnmf}
\alias{wrss}
\title{Multiplicative Updates for LS-NMF}
\usage{
  nmf.update.lsnmf(i, X, object, weight, eps = 10^-9, ...)

  wrss(object, X, weight)
}
\arguments{
  \item{i}{current iteration}

  \item{X}{target matrix}

  \item{object}{current NMF model}

  \item{weight}{value for \eqn{\Sigma}{S}, i.e. the weights
  that are applied to each entry in \code{X} by \code{X *
  weight} (= entry wise product). Weights are usually
  specified as a matrix of the same dimension as \code{X}
  (e.g. uncertainty estimates for each measurement), but
  may also be passed as a vector, in which case the
  standard rules for entry wise product between matrices
  and vectors apply (e.g. recylcing elements).}

  \item{eps}{small number passed to the standard
  euclidean-based NMF updates (see
  \code{\link{nmf.update.euclidean}}).}

  \item{...}{extra arguments (not used)}
}
\value{
  updated object \code{object}
}
\description{
  Implementation of the updates for the LS-NMF algorithm
  from \cite{Wang et al. (2006)}.

  \code{wrss} implements the objective function used by the
  LS-NMF algorithm.
}
\references{
  Wang G, Kossenkov AV and Ochs MF (2006). "LS-NMF: a
  modified non-negative matrix factorization algorithm
  utilizing uncertainty estimates." _BMC bioinformatics_,
  *7*, pp. 175. ISSN 1471-2105, <URL:
  http://dx.doi.org/10.1186/1471-2105-7-175>, <URL:
  http://www.ncbi.nlm.nih.gov/pubmed/16569230>.
}

