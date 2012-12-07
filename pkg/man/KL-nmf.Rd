\name{nmf_update.brunet_R}
\alias{brunet-nmf}
\alias{KL-nmf}
\alias{Mbrunet-nmf}
\alias{nmf_update.brunet}
\alias{nmf_update.brunet_R}
\alias{Rbrunet-nmf}
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
  criterion functions.}

  \item{copy}{logical that indicates if the update should
  be made on the original matrix directly (\code{FALSE}) or
  on a copy (\code{TRUE} - default). With \code{copy=FALSE}
  the memory footprint is very small, and some speed-up may
  be achieved in the case of big matrices. However, greater
  care should be taken due the side effect. We recommend
  that only experienced users use \code{copy=TRUE}.}
}
\description{
  The built-in NMF algorithms \sQuote{KL} and
  \sQuote{brunet} minimise the Kullback-Leibler divergence
  (KL) between an NMF model and a target matrix. They use
  the updates for the basis and coefficient matrices
  (\eqn{W} and \eqn{H}) defined by \cite{Brunet et al.
  (2004)}, which are essentially those from \cite{Lee et
  al. (2001)}, with an stabilisation step that shift up all
  entries from zero every 10 iterations, to a very small
  positive value.

  \code{nmf_update.brunet} implements in C++ an optimised
  version of the single update step.

  Algorithms \sQuote{brunet} and \sQuote{.R#brunet} provide
  the complete NMF algorithm from \cite{Brunet et al.
  (2004)}, using the C++-optimised and pure R updates
  \code{\link{nmf_update.brunet}} and
  \code{\link{nmf_update.brunet_R}} respectively.

  Algorithm \sQuote{KL} provides an NMF algorithm based on
  the C++-optimised version of the updates from
  \cite{Brunet et al. (2004)}, but using the stationarity
  of the objective value as a stopping criterion.

  \code{file.show(system.file('matlab', 'brunet.m',
  package='NMF'))}
}
\details{
  The only difference between algorithms \sQuote{KL} and
  \sQuote{brunet} is that \sQuote{KL} uses the stopping
  criterion based on the stationarity of the objective
  value \code{\link{nmf.stop.stationary}}, whereas
  \sQuote{brunet} uses the stationarity of the connectivity
  matrix \code{\link{nmf.stop.connectivity}}.

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

