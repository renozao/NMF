\name{nmf_update.lee_R}
\alias{Frobenius-nmf}
\alias{lee-nmf}
\alias{nmf_update.lee}
\alias{nmf_update.lee_R}
\alias{Rlee-nmf}
\title{NMF Algorithm/Updates for Frobenius Norm}
\usage{
  nmf_update.lee_R(i, v, x, rescale = TRUE, eps = 10^-9,
    ...)

  nmf_update.lee(i, v, x, rescale = TRUE, copy = FALSE,
    eps = 10^-9, weight = NULL, ...)
}
\arguments{
  \item{rescale}{logical that indicates if the basis matrix
  \eqn{W} should be rescaled so that its columns sum up to
  one.}

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
  criterion functions.}

  \item{copy}{logical that indicates if the update should
  be made on the original matrix directly (\code{FALSE}) or
  on a copy (\code{TRUE} - default). With \code{copy=FALSE}
  the memory footprint is very small, and some speed-up may
  be achieved in the case of big matrices. However, greater
  care should be taken due the side effect. We recommend
  that only experienced users use \code{copy=TRUE}.}

  \item{weight}{numeric vector of sample weights, e.g.,
  used to normalise samples coming from multiple datasets.
  It must be of the same length as the number of
  samples/columns in \code{v} -- and \code{h}.}
}
\description{
  The built-in NMF algorithms \sQuote{Frobenius} and
  \sQuote{lee} minimise the Frobenius norm (Euclidean
  distance) between an NMF model and a target matrix. They
  use the updates for the basis and coefficient matrices
  (\eqn{W} and \eqn{H}) defined by \cite{Lee et al.
  (2001)}.

  \code{nmf_update.lee} implements in C++ an optimised
  version of the single update step.

  Algorithms \sQuote{lee} and \sQuote{.R#lee} provide the
  complete NMF algorithm from \cite{Lee et al. (2001)},
  using the C++-optimised and pure R updates
  \code{\link{nmf_update.lee}} and
  \code{\link{nmf_update.lee_R}} respectively.

  Algorithm \sQuote{Frobenius} provides an NMF algorithm
  based on the C++-optimised version of the updates from
  \cite{Lee et al. (2001)}, but uses the stopping criterion
  based on the stationarity of the objective value.
}
\details{
  The only difference between algorithms \sQuote{Frobenius}
  and \sQuote{lee} is that \sQuote{Frobenius} uses the
  stopping criterion based on the stationarity of the
  objective value \code{\link{nmf.stop.stationary}},
  whereas \sQuote{lee} uses the stationarity of the
  connectivity matrix \code{\link{nmf.stop.connectivity}}.

  \code{nmf_update.lee_R} implements in pure R a single
  update step, i.e. it updates both matrices.
}
\author{
  Original update definition: D D Lee and HS Seung

  Port to R and optimisation in C++: Renaud Gaujoux
}
\references{
  Lee DD and Seung H (2001). "Algorithms for non-negative
  matrix factorization." _Advances in neural information
  processing systems_. <URL:
  http://scholar.google.com/scholar?q=intitle:Algorithms+for+non-negative+matrix+factorization\#0>.
}

