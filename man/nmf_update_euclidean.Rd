\name{nmf_update.euclidean.h}
\alias{nmf_update.euclidean}
\alias{nmf_update.euclidean.h}
\alias{nmf_update.euclidean.h_R}
\alias{nmf_update.euclidean.w}
\alias{nmf_update.euclidean.w_R}
\title{NMF Multiplicative Updates for Euclidean Distance}
\usage{
  nmf_update.euclidean.h(v, w, h, eps = 10^-9,
    nbterms = 0L, ncterms = 0L, copy = TRUE)

  nmf_update.euclidean.h_R(v, w, h, wh = NULL, eps = 10^-9)

  nmf_update.euclidean.w(v, w, h, eps = 10^-9,
    nbterms = 0L, ncterms = 0L, weight = NULL, copy = TRUE)

  nmf_update.euclidean.w_R(v, w, h, wh = NULL, eps = 10^-9)
}
\arguments{
  \item{eps}{small numeric value used to ensure numeric stability, by shifting up
entries from zero to this fixed value.}

  \item{wh}{already computed NMF estimate used to compute the denominator term.}

  \item{weight}{numeric vector of sample weights, e.g., used to normalise samples
coming from multiple datasets.
It must be of the same length as the number of samples/columns in \code{v}
-- and \code{h}.}

  \item{v}{target matrix}

  \item{w}{current basis matrix}

  \item{h}{current coefficient matrix}

  \item{nbterms}{number of fixed basis terms}

  \item{ncterms}{number of fixed coefficient terms}

  \item{copy}{logical that indicates if the update should be made on the original
matrix directly (\code{FALSE}) or on a copy (\code{TRUE} - default).
With \code{copy=FALSE} the memory footprint is very small, and some speed-up may be
achieved in the case of big matrices.
However, greater care should be taken due the side effect.
We recommend that only experienced users use \code{copy=TRUE}.}
}
\value{
  a matrix of the same dimension as the input matrix to update
(i.e. \code{w} or \code{h}).
If \code{copy=FALSE}, the returned matrix uses the same memory as the input object.
}
\description{
  Multiplicative updates from \cite{Lee et al. (2001)} for standard Nonnegative Matrix Factorization
models \eqn{V \approx W H}, where the distance between the target matrix and its NMF
estimate is measured by the -- euclidean -- Frobenius norm.

\code{nmf_update.euclidean.w} and \code{nmf_update.euclidean.h} compute the updated basis and coefficient
matrices respectively.
They use a \emph{C++} implementation which is optimised for speed and memory usage.

\code{nmf_update.euclidean.w_R} and \code{nmf_update.euclidean.h_R} implement the same updates
in \emph{plain R}.
}
\details{
  The coefficient matrix (\code{H}) is updated as follows:
\deqn{
H_{kj} \leftarrow \frac{\max(H_{kj} W^T V)_{kj}, \varepsilon) }{(W^T W H)_{kj} + \varepsilon}
}{
H_kj <- max(H_kj (W^T V)_kj, eps) / ( (W^T W H)_kj + eps )
}

These updates are used by the built-in NMF algorithms \code{\link[=Frobenius-nmf]{Frobenius}} and
\code{\link[=lee-nmf]{lee}}.

The basis matrix (\code{W}) is updated as follows:
\deqn{
W_ik \leftarrow \frac{\max(W_ik (V H^T)_ik, \varepsilon) }{ (W H H^T)_ik + \varepsilon}
}{
W_ik <- max(W_ik (V H^T)_ik, eps) / ( (W H H^T)_ik + eps )
}
}
\author{
  Update definitions by \cite{Lee2001}.

C++ optimised implementation by Renaud Gaujoux.
}
\references{
  Lee DD and Seung H (2001). "Algorithms for non-negative matrix factorization." _Advances in neural information processing
systems_. <URL: http://scholar.google.com/scholar?q=intitle:Algorithms+for+non-negative+matrix+factorization\#0>.
}

