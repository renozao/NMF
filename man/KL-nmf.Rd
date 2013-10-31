\name{nmf_update.brunet_R}
\alias{brunet_M-nmf}
\alias{brunet-nmf}
\alias{brunet_R-nmf}
\alias{KL-nmf}
\alias{nmfAlgorithm.brunet}
\alias{nmfAlgorithm.brunet_M}
\alias{nmfAlgorithm.brunet_R}
\alias{nmfAlgorithm.KL}
\alias{nmf_update.brunet}
\alias{nmf_update.brunet_R}
\title{NMF Algorithm/Updates for Kullback-Leibler Divergence}
\source{
  Original MATLAB files and references can be found at:

  \url{http://www.broadinstitute.org/mpr/publications/projects/NMF/nmf.m}

  \url{http://www.broadinstitute.org/publications/broad872}

  Original license terms:

  This software and its documentation are copyright 2004 by
  the Broad Institute/Massachusetts Institute of
  Technology. All rights are reserved. This software is
  supplied without any warranty or guaranteed support
  whatsoever. Neither the Broad Institute nor MIT can not
  be responsible for its use, misuse, or functionality.
}
\usage{
  nmf_update.brunet_R(i, v, x, eps = .Machine$double.eps,
    ...)

  nmf_update.brunet(i, v, x, copy = FALSE,
    eps = .Machine$double.eps, ...)

  nmfAlgorithm.brunet_R(..., .stop = NULL, maxIter = 2000,
    eps = .Machine$double.eps, stopconv = 40,
    check.interval = 10)

  nmfAlgorithm.brunet(..., .stop = NULL, maxIter = 2000,
    copy = FALSE, eps = .Machine$double.eps, stopconv = 40,
    check.interval = 10)

  nmfAlgorithm.KL(..., .stop = NULL, maxIter = 2000,
    copy = FALSE, eps = .Machine$double.eps,
    stationary.th = .Machine$double.eps,
    check.interval = 5 * check.niter, check.niter = 10L)

  nmfAlgorithm.brunet_M(..., object, y, x)
}
\arguments{
  \item{i}{current iteration number.}

  \item{v}{target matrix.}

  \item{x}{current NMF model, as an
  \code{\linkS4class{NMF}} object.}

  \item{eps}{small numeric value used to ensure numeric
  stability, by shifting up entries from zero to this fixed
  value.}

  \item{...}{extra arguments. These are generally not used
  and present only to allow other arguments from the main
  call to be passed to the initialisation and stopping
  criterion functions (slots \code{onInit} and \code{Stop}
  respectively).}

  \item{copy}{logical that indicates if the update should
  be made on the original matrix directly (\code{FALSE}) or
  on a copy (\code{TRUE} - default). With \code{copy=FALSE}
  the memory footprint is very small, and some speed-up may
  be achieved in the case of big matrices. However, greater
  care should be taken due the side effect. We recommend
  that only experienced users use \code{copy=TRUE}.}

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

  \item{object}{an object computed using some algorithm, or
  that describes an algorithm itself.}

  \item{y}{data object, e.g. a target matrix}

  \item{stopconv}{number of iterations intervals over which
  the connectivity matrix must not change for stationarity
  to be achieved.}

  \item{check.interval}{interval (in number of iterations)
  on which the stopping criterion is computed.}

  \item{stationary.th}{maximum absolute value of the
  gradient, for the objective function to be considered
  stationary.}

  \item{check.niter}{number of successive iteration used to
  compute the stationnary criterion.}
}
\description{
  The built-in NMF algorithms described here minimise the
  Kullback-Leibler divergence (KL) between an NMF model and
  a target matrix. They use the updates for the basis and
  coefficient matrices (\eqn{W} and \eqn{H}) defined by
  \cite{Brunet et al. (2004)}, which are essentially those
  from \cite{Lee et al. (2001)}, with an stabilisation step
  that shift up all entries from zero every 10 iterations,
  to a very small positive value.

  \code{nmf_update.brunet} implements in C++ an optimised
  version of the single update step.

  Algorithms \sQuote{brunet} and \sQuote{.R#brunet} provide
  the complete NMF algorithm from \cite{Brunet et al.
  (2004)}, using the C++-optimised and pure R updates
  \code{\link{nmf_update.brunet}} and
  \code{\link{nmf_update.brunet_R}} respectively.

  Algorithm \sQuote{KL} provides an NMF algorithm based on
  the C++-optimised version of the updates from
  \cite{Brunet et al. (2004)}, which uses the stationarity
  of the objective value as a stopping criterion
  \code{\link{nmf.stop.stationary}}, instead of the
  stationarity of the connectivity matrix
  \code{\link{nmf.stop.connectivity}} as used by
  \sQuote{brunet}.

  \samp{ library(RcppOctave)
  file.show(system.mfile('brunet.m', package='NMF')) }
}
\details{
  \code{nmf_update.brunet_R} implements in pure R a single
  update step, i.e. it updates both matrices.
}
\author{
  Original implementation in MATLAB: Jean-Philippe Brunet
  \email{brunet@broad.mit.edu}

  Port to R and optimisation in C++: Renaud Gaujoux
}
\references{
  Brunet J, Tamayo P, Golub TR and Mesirov JP (2004).
  "Metagenes and molecular pattern discovery using matrix
  factorization." _Proceedings of the National Academy of
  Sciences of the United States of America_, *101*(12), pp.
  4164-9. ISSN 0027-8424, <URL:
  http://dx.doi.org/10.1073/pnas.0308531101>, <URL:
  http://www.ncbi.nlm.nih.gov/pubmed/15016911>.

  Lee DD and Seung H (2001). "Algorithms for non-negative
  matrix factorization." _Advances in neural information
  processing systems_. <URL:
  http://scholar.google.com/scholar?q=intitle:Algorithms+for+non-negative+matrix+factorization\#0>.
}

