\name{nmf_update.lsnmf}
\alias{lsNMF-nmf}
\alias{nmfAlgorithm.lsNMF}
\alias{nmf_update.lsnmf}
\alias{wrss}
\title{Multiplicative Updates for LS-NMF}
\usage{
  nmf_update.lsnmf(i, X, object, weight, eps = 10^-9, ...)

  wrss(object, X, weight)

  nmfAlgorithm.lsNMF(..., .stop = NULL,
    maxIter = nmf.getOption("maxIter") \%||\% 2000, weight,
    eps = 10^-9, stationary.th = .Machine$double.eps,
    check.interval = 5 * check.niter, check.niter = 10L)
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
  \code{\link{nmf_update.euclidean}}).}

  \item{...}{extra arguments (not used)}

  \item{.stop}{specification of a stopping criterion, that
  is used instead of the one associated to the NMF
  algorithm. It may be specified as: \itemize{ \item the
  access key of a registered stopping criterion; \item a
  single integer that specifies the exact number of
  iterations to perform, which will be honoured unless a
  lower value is explicitly passed in argument
  \code{maxIter}. \item a single numeric value that
  specifies the stationnarity threshold for the objective
  function, used in with \code{\link{nmf.stop.stationary}};
  \item a function with signature
  \code{(object="NMFStrategy", i="integer", y="matrix",
  x="NMF", ...)}, where \code{object} is the
  \code{NMFStrategy} object that describes the algorithm
  being run, \code{i} is the current iteration, \code{y} is
  the target matrix and \code{x} is the current value of
  the NMF model. }}

  \item{maxIter}{maximum number of iterations to perform.}

  \item{stationary.th}{maximum absolute value of the
  gradient, for the objective function to be considered
  stationary.}

  \item{check.interval}{interval (in number of iterations)
  on which the stopping criterion is computed.}

  \item{check.niter}{number of successive iteration used to
  compute the stationnary criterion.}
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

