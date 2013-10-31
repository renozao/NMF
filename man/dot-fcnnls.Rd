\name{.fcnnls}
\alias{.fcnnls}
\title{Internal Routine for Fast Combinatorial Nonnegative Least-Squares}
\usage{
  .fcnnls(x, y, verbose = FALSE, pseudo = FALSE, eps = 0)
}
\arguments{
  \item{x}{the coefficient matrix}

  \item{y}{the target matrix to be approximated by \eqn{X
  K}.}

  \item{verbose}{logical that indicates if log messages
  should be shown.}

  \item{pseudo}{By default (\code{pseudo=FALSE}) the
  algorithm uses Gaussian elimination to solve the
  successive internal linear problems, using the
  \code{\link{solve}} function.  If \code{pseudo=TRUE} the
  algorithm uses Moore-Penrose generalized
  \code{\link[corpcor]{pseudoinverse}} from the
  \code{corpcor} package instead of \link{solve}.}

  \item{eps}{threshold for considering entries as
  nonnegative. This is an experimental parameter, and it is
  recommended to leave it at 0.}
}
\value{
  A list with the following elements:

  \item{coef}{the fitted coefficient matrix.}
  \item{Pset}{the set of passive constraints, as a logical
  matrix of the same size as \code{K} that indicates which
  element is positive.}
}
\description{
  This is the workhorse function for the higher-level
  function \code{\link{fcnnls}}, which implements the fast
  nonnegative least-square algorithm for multiple
  right-hand-sides from \cite{Van Benthem et al. (2004)} to
  solve the following problem:

  \deqn{ \begin{array}{l} \min \|Y - X K\|_F\\ \mbox{s.t. }
  K>=0 \end{array} }{min ||Y - X K||_F, s.t. K>=0}

  where \eqn{Y} and \eqn{X} are two real matrices of
  dimension \eqn{n \times p}{n x p} and \eqn{n \times r}{n
  x r} respectively, and \eqn{\|.\|_F}{|.|_F} is the
  Frobenius norm.

  The algorithm is very fast compared to other approaches,
  as it is optimised for handling multiple right-hand
  sides.
}
\references{
  Van Benthem M and Keenan MR (2004). "Fast algorithm for
  the solution of large-scale non-negativity-constrained
  least squares problems." _Journal of Chemometrics_,
  *18*(10), pp. 441-450. ISSN 0886-9383, <URL:
  http://dx.doi.org/10.1002/cem.889>, <URL:
  http://doi.wiley.com/10.1002/cem.889>.
}

