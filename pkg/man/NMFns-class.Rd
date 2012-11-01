\docType{class}
\name{NMFns-class}
\alias{NMFns-class}
\title{NMF Model - Nonsmooth Nonnegative Matrix Factorization}
\description{
  This class implements the \emph{Nonsmooth Nonnegative
  Matrix Factorization} (nsNMF) model, required by the
  Nonsmooth NMF algorithm.
}
\details{
  The Nonsmooth NMF algorithm is defined by
  \cite{Pascual-Montano et al. (2006)} as a modification of
  the standard divergence based NMF algorithm (see section
  Details and references below).  It aims at obtaining
  sparser factor matrices, by the introduction of a
  smoothing matrix.

  The Nonsmooth NMF algorithm is a modification of the
  standard divergence based NMF algorithm (see
  \code{\linkS4class{NMF}}). Given a non-negative \eqn{n
  \times p}{n x p} matrix \eqn{V} and a factorization rank
  \eqn{r}, it fits the following model:

  \deqn{V \equiv W S(\theta) H,}{V &#126; W S(theta) H,}
  where: \itemize{

  \item \eqn{W} and \eqn{H} are such as in the standard
  model, i.e. non-negative matrices of dimension \eqn{n
  \times r}{n x r} and \eqn{r \times p}{r x p}
  respectively;

  \item \eqn{S} is a \eqn{r \times r} square matrix whose
  entries depends on an extra parameter \eqn{0\leq \theta
  \leq 1} in the following way: \deqn{S = (1-\theta)I +
  \frac{\theta}{r} 11^T ,} where \eqn{I} is the identity
  matrix and \eqn{1} is a vector of ones.

  }

  The interpretation of S as a smoothing matrix can be
  explained as follows: Let \eqn{X} be a positive, nonzero,
  vector. Consider the transformed vector \eqn{Y = S X}. If
  \eqn{\theta = 0}, then \eqn{Y = X} and no smoothing on
  \eqn{X} has occurred.  However, as \eqn{\theta \to
  1}{theta tends to 1}, the vector \eqn{Y} tends to the
  constant vector with all elements almost equal to the
  average of the elements of \eqn{X}. This is the smoothest
  possible vector in the sense of non-sparseness because
  all entries are equal to the same nonzero value, instead
  of having some values close to zero and others clearly
  nonzero.
}
\section{Methods}{
  \describe{

  \item{fitted}{\code{signature(object = "NMFns")}: Compute
  estimate for an NMFns object, according to the Nonsmooth
  NMF model (cf. \code{\link{NMFns-class}}).

  Extra arguments in \code{...} are passed to method
  \code{smoothing}, and are typically used to pass a value
  for \code{theta}, which is used to compute the smoothing
  matrix instead of the one stored in \code{object}. }

  \item{show}{\code{signature(object = "NMFns")}: Show
  method for objects of class \code{NMFns} }

  }
}
\section{Creating objects from the Class}{
  Object of class \code{NMFns} can be created using the
  standard way with operator \code{\link{new}}

  However, as for all NMF model classes -- that extend
  class \code{\linkS4class{NMF}}, objects of class
  \code{NMFns} should be created using factory method
  \code{\link{nmfModel}} :

  \code{new('NMFns')}

  \code{nmfModel(model='NMFns')}

  \code{nmfModel(model='NMFns', W=w, theta=0.3}

  See \code{\link{nmfModel}} for more details on how to use
  the factory method.
}

\section{Algorithm}{
  The Nonsmooth NMF algorithm uses a modified version of
  the multiplicative update equations in Lee & Seung's
  method for Kullback-Leibler divergence minimization. The
  update equations are modified to take into account the --
  constant -- smoothing matrix. The modification reduces to
  using matrix \eqn{W S} instead of matrix \eqn{W} in the
  update of matrix \eqn{H}, and similarly using matrix
  \eqn{S H} instead of matrix \eqn{H} in the update of
  matrix \eqn{W}.

  After the matrix \eqn{W} has been updated, each of its
  columns is scaled so that it sums up to 1.
}
\examples{
# create a completely empty NMFns object
new('NMFns')

# create a NMF object based on random (compatible) matrices
n <- 50; r <- 3; p <- 20
w <- rmatrix(n, r)
h <- rmatrix(r, p)
nmfModel(model='NMFns', W=w, H=h)

# apply Nonsmooth NMF algorithm to a random target matrix
V <- rmatrix(n, p)
\dontrun{nmf(V, r, 'ns')}

# random nonsmooth NMF model
rnmf(3, 10, 5, model='NMFns', theta=0.3)
}
\references{
  Pascual-Montano A, Carazo JM, Kochi K, Lehmann D and
  Pascual-marqui RD (2006). "Nonsmooth nonnegative matrix
  factorization (nsNMF)." _IEEE Trans. Pattern Anal. Mach.
  Intell_, *28*, pp. 403-415.
}
\seealso{
  Other NMF-model:
  \code{\link{initialize,NMFOffset-method}},
  \code{\link{NMFOffset-class}}, \code{\link{NMFstd-class}}
}

