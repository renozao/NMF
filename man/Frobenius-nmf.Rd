\name{nmf_update.lee_R}
\alias{Frobenius-nmf}
\alias{lee-nmf}
\alias{lee_R-nmf}
\alias{nmfAlgorithm.Frobenius}
\alias{nmfAlgorithm.lee}
\alias{nmfAlgorithm.lee_R}
\alias{nmf_update.lee}
\alias{nmf_update.lee_R}
\title{NMF Algorithm/Updates for Frobenius Norm}
\usage{
  nmf_update.lee_R(i, v, x, rescale = TRUE, eps = 10^-9,
    ...)

  nmf_update.lee(i, v, x, rescale = TRUE, copy = FALSE,
    eps = 10^-9, weight = NULL, ...)

  nmfAlgorithm.lee_R(..., .stop = NULL,
    maxIter = nmf.getOption("maxIter") \%||\% 2000,
    rescale = TRUE, eps = 10^-9, stopconv = 40,
    check.interval = 10)

  nmfAlgorithm.lee(..., .stop = NULL,
    maxIter = nmf.getOption("maxIter") \%||\% 2000,
    rescale = TRUE, copy = FALSE, eps = 10^-9,
    weight = NULL, stopconv = 40, check.interval = 10)

  nmfAlgorithm.Frobenius(..., .stop = NULL,
    maxIter = nmf.getOption("maxIter") \%||\% 2000,
    rescale = TRUE, copy = FALSE, eps = 10^-9,
    weight = NULL, stationary.th = .Machine$double.eps,
    check.interval = 5 * check.niter, check.niter = 10L)
}
\arguments{
  \item{rescale}{logical that indicates if the basis matrix \eqn{W} should be
rescaled so that its columns sum up to one.}

  \item{i}{current iteration number.}

  \item{v}{target matrix.}

  \item{x}{current NMF model, as an \code{\linkS4class{NMF}} object.}

  \item{eps}{small numeric value used to ensure numeric stability, by shifting up
entries from zero to this fixed value.}

  \item{...}{extra arguments. These are generally not used and present
only to allow other arguments from the main call to be passed to the
initialisation and stopping criterion functions (slots \code{onInit} and
\code{Stop} respectively).}

  \item{copy}{logical that indicates if the update should be made on the original
matrix directly (\code{FALSE}) or on a copy (\code{TRUE} - default).
With \code{copy=FALSE} the memory footprint is very small, and some speed-up may be
achieved in the case of big matrices.
However, greater care should be taken due the side effect.
We recommend that only experienced users use \code{copy=TRUE}.}

  \item{.stop}{specification of a stopping criterion, that is used instead of the
one associated to the NMF algorithm.
It may be specified as:
\itemize{
\item the access key of a registered stopping criterion;
\item a single integer that specifies the exact number of iterations to perform, which will
be honoured unless a lower value is explicitly passed in argument \code{maxIter}.
\item a single numeric value that specifies the stationnarity threshold for the
objective function, used in with \code{\link{nmf.stop.stationary}};
\item a function with signature \code{(object="NMFStrategy", i="integer", y="matrix", x="NMF", ...)},
where \code{object} is the \code{NMFStrategy} object that describes the algorithm being run,
\code{i} is the current iteration, \code{y} is the target matrix and \code{x} is the current value of
the NMF model.
}}

  \item{maxIter}{maximum number of iterations to perform.}

  \item{stopconv}{number of iterations intervals over which the connectivity
matrix must not change for stationarity to be achieved.}

  \item{check.interval}{interval (in number of iterations) on which the stopping
criterion is computed.}

  \item{stationary.th}{maximum absolute value of the gradient, for the objective
function to be considered stationary.}

  \item{check.niter}{number of successive iteration used to compute the stationnary
criterion.}

  \item{weight}{numeric vector of sample weights, e.g., used to normalise samples
coming from multiple datasets.
It must be of the same length as the number of samples/columns in \code{v}
-- and \code{h}.}
}
\description{
  The built-in NMF algorithms described here minimise
the Frobenius norm (Euclidean distance) between an NMF model and a target matrix.
They use the updates for the basis and coefficient matrices (\eqn{W} and \eqn{H})
defined by \cite{Lee et al. (2001)}.

\code{nmf_update.lee} implements in C++ an optimised version of the single update step.

Algorithms \sQuote{lee} and \sQuote{.R#lee} provide the complete NMF algorithm from \cite{Lee et al. (2001)},
using the C++-optimised and pure R updates \code{\link{nmf_update.lee}} and \code{\link{nmf_update.lee_R}}
respectively.

Algorithm \sQuote{Frobenius} provides an NMF algorithm based on the C++-optimised version of
the updates from \cite{Lee et al. (2001)}, which uses the stationarity of the objective value
as a stopping criterion \code{\link{nmf.stop.stationary}}, instead of the
stationarity of the connectivity matrix \code{\link{nmf.stop.connectivity}} as used by
\sQuote{lee}.
}
\details{
  \code{nmf_update.lee_R} implements in pure R a single update step, i.e. it updates
both matrices.
}
\author{
  Original update definition: D D Lee and HS Seung

Port to R and optimisation in C++: Renaud Gaujoux
}
\references{
  Lee DD and Seung H (2001). "Algorithms for non-negative matrix factorization." _Advances in neural information processing
systems_. <URL: http://scholar.google.com/scholar?q=intitle:Algorithms+for+non-negative+matrix+factorization\#0>.
}

