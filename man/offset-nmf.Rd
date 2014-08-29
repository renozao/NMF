\name{nmf_update.euclidean_offset.h}
\alias{nmfAlgorithm.offset}
\alias{nmfAlgorithm.offset_R}
\alias{nmf_update.euclidean_offset.h}
\alias{nmf_update.euclidean_offset.w}
\alias{nmf_update.offset}
\alias{nmf_update.offset_R}
\alias{offset_R-nmf}
\title{NMF Multiplicative Update for NMF with Offset Models}
\usage{
  nmf_update.euclidean_offset.h(v, w, h, offset,
    eps = 10^-9, copy = TRUE)

  nmf_update.euclidean_offset.w(v, w, h, offset,
    eps = 10^-9, copy = TRUE)

  nmf_update.offset_R(i, v, x, eps = 10^-9, ...)

  nmf_update.offset(i, v, x, copy = FALSE, eps = 10^-9,
    ...)

  nmfAlgorithm.offset_R(..., .stop = NULL,
    maxIter = nmf.getOption("maxIter") \%||\% 2000,
    eps = 10^-9, stopconv = 40, check.interval = 10)

  nmfAlgorithm.offset(..., .stop = NULL,
    maxIter = nmf.getOption("maxIter") \%||\% 2000,
    copy = FALSE, eps = 10^-9, stopconv = 40,
    check.interval = 10)
}
\arguments{
  \item{offset}{current value of the offset/intercept vector.
It must be of length equal to the number of rows in the target matrix.}

  \item{v}{target matrix.}

  \item{eps}{small numeric value used to ensure numeric stability, by shifting up
entries from zero to this fixed value.}

  \item{copy}{logical that indicates if the update should be made on the original
matrix directly (\code{FALSE}) or on a copy (\code{TRUE} - default).
With \code{copy=FALSE} the memory footprint is very small, and some speed-up may be
achieved in the case of big matrices.
However, greater care should be taken due the side effect.
We recommend that only experienced users use \code{copy=TRUE}.}

  \item{i}{current iteration number.}

  \item{x}{current NMF model, as an \code{\linkS4class{NMF}} object.}

  \item{...}{extra arguments. These are generally not used and present
only to allow other arguments from the main call to be passed to the
initialisation and stopping criterion functions (slots \code{onInit} and
\code{Stop} respectively).}

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

  \item{w}{current basis matrix}

  \item{h}{current coefficient matrix}
}
\value{
  an \code{\linkS4class{NMFOffset}} model object.
}
\description{
  These update rules proposed by \cite{Badea (2008)} are
  modified version of the updates from \cite{Lee et al.
  (2001)}, that include an offset/intercept vector, which
  models a common baseline for each feature accross all
  samples: \deqn{V \approx W H + I}

  \code{nmf_update.euclidean_offset.h} and
  \code{nmf_update.euclidean_offset.w} compute the updated
  NMFOffset model, using the optimized \emph{C++}
  implementations.

  \code{nmf_update.offset_R} implements a complete single
  update step, using plain R updates.

  \code{nmf_update.offset} implements a complete single
  update step, using C++-optimised updates.

  Algorithms \sQuote{offset} and \sQuote{.R#offset} provide
  the complete NMF-with-offset algorithm from \cite{Badea
  (2008)}, using the C++-optimised and pure R updates
  \code{\link{nmf_update.offset}} and
  \code{\link{nmf_update.offset_R}} respectively.
}
\details{
  The associated model is defined as an
  \code{\linkS4class{NMFOffset}} object. The details of the
  multiplicative updates can be found in \cite{Badea
  (2008)}. Note that the updates are the ones defined for a
  single datasets, not the simultaneous NMF model, which is
  fit by algorithm \sQuote{siNMF} from formula-based NMF
  models.
}
\author{
  Original update definition: Liviu Badea

Port to R and optimisation in C++: Renaud Gaujoux
}
\references{
  Badea L (2008). "Extracting gene expression profiles
  common to colon and pancreatic adenocarcinoma using
  simultaneous nonnegative matrix factorization." _Pacific
  Symposium on Biocomputing. Pacific Symposium on
  Biocomputing_, *290*, pp. 267-78. ISSN 1793-5091, <URL:
  http://www.ncbi.nlm.nih.gov/pubmed/18229692>.

  Lee DD and Seung H (2001). "Algorithms for non-negative
  matrix factorization." _Advances in neural information
  processing systems_. <URL:
  http://scholar.google.com/scholar?q=intitle:Algorithms+for+non-negative+matrix+factorization\#0>.
}

