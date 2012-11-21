\name{nmf_update.ns}
\alias{nmf_update.ns}
\alias{nmf_update.ns_R}
\alias{nsNMF-nmf}
\alias{RnsNMF-nmf}
\title{NMF Multiplicative Update for Nonsmooth Nonnegative Matrix Factorization (nsNMF).}
\usage{
  nmf_update.ns(i, v, x, copy = FALSE, ...)

  nmf_update.ns_R(i, v, x, ...)
}
\arguments{
  \item{i}{current iteration number.}

  \item{v}{target matrix.}

  \item{x}{current NMF model, as an
  \code{\linkS4class{NMF}} object.}

  \item{copy}{logical that indicates if the update should
  be made on the original matrix directly (\code{FALSE}) or
  on a copy (\code{TRUE} - default). With \code{copy=FALSE}
  the memory footprint is very small, and some speed-up may
  be achieved in the case of big matrices. However, greater
  care should be taken due the side effect. We recommend
  that only experienced users use \code{copy=TRUE}.}

  \item{...}{extra arguments. These are generally not used
  and present only to allow other arguments from the main
  call to be passed to the initialisation and stopping
  criterion functions.}
}
\value{
  an \code{\linkS4class{NMFns}} model object.
}
\description{
  These update rules, defined for the
  \code{\linkS4class{NMFns}} model \eqn{V \approx W S H}
  from \cite{Pascual-Montano et al. (2006)}, that
  introduces an intermediate smoothing matrix to enhance
  sparsity of the factors.

  \code{nmf_update.ns_R} implements the same updates in
  \emph{plain R}.

  Algorithms \sQuote{nsNMF} and \sQuote{.R#nsNMF} provide
  the complete NMF algorithm from \cite{Pascual-Montano et
  al. (2006)}, using the C++-optimised and plain R updates
  \code{\link{nmf_update.brunet}} and
  \code{\link{nmf_update.brunet_R}} respectively. The
  stopping criterion is based on the stationarity of the
  connectivity matrix.
}
\details{
  \code{nmf_update.ns} computes the updated nsNMF model. It
  uses the optimized \emph{C++} implementations
  \code{\link{nmf_update.KL.w}} and
  \code{\link{nmf_update.KL.h}} to update \eqn{W} and
  \eqn{H} respectively.

  The multiplicative updates are based on the updates
  proposed by \cite{Brunet et al. (2004)}, except that the
  NMF estimate \eqn{W H} is replaced by \eqn{W S H} and
  \eqn{W} (resp. \eqn{H}) is replaced by \eqn{W S} (resp.
  \eqn{S H}) in the update of \eqn{H} (resp. \eqn{W}).

  See \code{\link{nmf_update.KL}} for more details on the
  update formula.
}
\references{
  Pascual-Montano A, Carazo JM, Kochi K, Lehmann D and
  Pascual-marqui RD (2006). "Nonsmooth nonnegative matrix
  factorization (nsNMF)." _IEEE Trans. Pattern Anal. Mach.
  Intell_, *28*, pp. 403-415.

  Brunet J, Tamayo P, Golub TR and Mesirov JP (2004).
  "Metagenes and molecular pattern discovery using matrix
  factorization." _Proceedings of the National Academy of
  Sciences of the United States of America_, *101*(12), pp.
  4164-9. ISSN 0027-8424, <URL:
  http://dx.doi.org/10.1073/pnas.0308531101>, <URL:
  http://www.ncbi.nlm.nih.gov/pubmed/15016911>.
}

