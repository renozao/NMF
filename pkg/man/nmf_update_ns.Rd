\name{nmf_update.ns}
\alias{nmf_update.ns}
\alias{nmf_update.ns_R}
\title{NMF Multiplicative Update for Nonsmooth Nonnegative Matrix Factorization (nsNMF).}
\usage{
  nmf_update.ns(i, v, data, copy = FALSE, ...)

  nmf_update.ns_R(i, v, data, ...)
}
\arguments{
  \item{i}{current iteration}

  \item{v}{target matrix}

  \item{data}{current NMF model}

  \item{copy}{logical that indicates if the update should
  be made in place (\code{FALSE}) or on a copy of the
  current NMF model (\code{TRUE} - default).}

  \item{...}{extra arguments to cope with arguments that
  are not aimed at this function}
}
\value{
  an \code{\linkS4class{NMFns}} model object.
}
\description{
  These update rules, defined for the \code{nsNMF} model
  \eqn{V \approx W S H} from \cite{Pascual-Montano et al.
  (2006)}, that introduces an intermediate smoothing matrix
  to enhance sparsity of the factors.

  \code{nmf_update.ns_R} implements the same updates in
  \emph{plain R}.
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

