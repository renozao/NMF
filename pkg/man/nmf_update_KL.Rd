\name{nmf_update.KL.h}
\alias{nmf_update.KL}
\alias{nmf_update.KL.h}
\alias{nmf_update.KL.h_R}
\alias{nmf_update.KL.w}
\alias{nmf_update.KL.w_R}
\title{NMF Multiplicative Updates for Kullback-Leibler Divergence}
\usage{
  nmf_update.KL.h(v, w, h, nbterms = 0L, ncterms = 0L,
    copy = TRUE)

  nmf_update.KL.h_R(v, w, h, wh = NULL)

  nmf_update.KL.w(v, w, h, nbterms = 0L, ncterms = 0L,
    copy = TRUE)

  nmf_update.KL.w_R(v, w, h, wh = NULL)
}
\arguments{
  \item{v}{target matrix}

  \item{w}{current basis matrix}

  \item{h}{current coefficient matrix}

  \item{nbterms}{number of fixed basis terms}

  \item{ncterms}{number of fixed coefficient terms}

  \item{copy}{logical that indicates if the update should
  be made on the original matrix directly (\code{FALSE}) or
  on a copy (\code{TRUE} - default). With \code{copy=FALSE}
  the memory footprint is very small, and some speed-up may
  be achieved in the case of big matrices. However, greater
  care should be taken due the side effect. We recommend
  that only experienced users use \code{copy=TRUE}.}

  \item{wh}{already computed NMF estimate used to compute
  the denominator term.}
}
\value{
  a matrix of the same dimension as the input matrix to
  update (i.e. \code{w} or \code{h}). If \code{copy=FALSE},
  the returned matrix uses the same memory as the input
  object.
}
\description{
  Multiplicative updates from \cite{Lee et al. (2001)} for
  standard Nonnegative Matrix Factorization models \eqn{V
  \approx W H}, where the distance between the target
  matrix and its NMF estimate is measured by the
  Kullback-Leibler divergence.

  \code{nmf_update.KL.w_R} and \code{nmf_update.KL.h_R}
  implement the same updates in \emph{plain R}.
}
\details{
  \code{nmf_update.KL.w} and \code{nmf_update.KL.h} compute
  the updated basis and coefficient matrices respectively.
  They use a \emph{C++} implementation which is optimised
  for speed and memory usage.

  The coefficient matrix (\code{H}) is updated as follows:
  \deqn{ H_{kj} \leftarrow H_{kj} \frac{\left( sum_i
  \frac{W_{ik} V_{ij}}{(WH)_{ij}} \right)}{ sum_i W_{ik} }.
  }{ H_kj <- H_kj ( sum_i [ W_ik V_ij / (WH)_ij ] ) / (
  sum_i W_ik ) }

  These updates are used in built-in NMF algorithms
  \code{\link[=KL-nmf]{KL}} and
  \code{\link[=brunet-nmf]{brunet}}.

  The basis matrix (\code{W}) is updated as follows: \deqn{
  W_{ik} \leftarrow W_{ik} \frac{ sum_j [\frac{H_{kj}
  A_{ij}}{(WH)_{ij}} ] }{sum_j H_{kj} } }{ W_ik <- W_ik (
  sum_u [H_kl A_il / (WH)_il ] ) / ( sum_l H_kl ) }
}
\author{
  Update definitions by \cite{Lee2001}.

  C++ optimised implementation by Renaud Gaujoux.
}
\references{
  Lee DD and Seung H (2001). "Algorithms for non-negative
  matrix factorization." _Advances in neural information
  processing systems_. <URL:
  http://scholar.google.com/scholar?q=intitle:Algorithms+for+non-negative+matrix+factorization\#0>.
}

