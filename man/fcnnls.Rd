\docType{methods}
\name{fcnnls}
\alias{fcnnls}
\alias{fcnnls,ANY,numeric-method}
\alias{fcnnls,matrix,matrix-method}
\alias{fcnnls-methods}
\alias{fcnnls,numeric,matrix-method}
\title{Fast Combinatorial Nonnegative Least-Square}
\usage{
  fcnnls(x, y, ...)

  \S4method{fcnnls}{matrix,matrix}(x, y, verbose = FALSE,
    pseudo = TRUE, ...)
}
\arguments{
  \item{...}{extra arguments passed to the internal function \code{.fcnnls}.
Currently not used.}

  \item{verbose}{toggle verbosity (default is \code{FALSE}).}

  \item{x}{the coefficient matrix}

  \item{y}{the target matrix to be approximated by \eqn{X K}.}

  \item{pseudo}{By default (\code{pseudo=FALSE}) the algorithm uses Gaussian
elimination to solve the successive internal linear problems, using the
\code{\link{solve}} function.  If \code{pseudo=TRUE} the algorithm uses
Moore-Penrose generalized \code{\link[corpcor]{pseudoinverse}} from the
\code{corpcor} package instead of \link{solve}.}
}
\value{
  A list containing the following components:

\item{x}{ the estimated optimal matrix \eqn{K}.} \item{fitted}{ the fitted
matrix \eqn{X K}.} \item{residuals}{ the residual matrix \eqn{Y - X K}.}
\item{deviance}{ the residual sum of squares between the fitted matrix
\eqn{X K} and the target matrix \eqn{Y}. That is the sum of the square
residuals.} \item{passive}{ a \eqn{r x p} logical matrix containing the
passive set, that is the set of entries in \eqn{K} that are not null (i.e.
strictly positive).} \item{pseudo}{ a logical that is \code{TRUE} if the
computation was performed using the pseudoinverse. See argument
\code{pseudo}.}
}
\description{
  This function solves the following nonnegative least square linear problem
using normal equations and the fast combinatorial strategy from \cite{Van Benthem et al. (2004)}:

\deqn{
 \begin{array}{l}
 \min \|Y - X K\|_F\\
 \mbox{s.t. } K>=0
 \end{array}
}{min ||Y - X K||_F, s.t. K>=0}

where \eqn{Y} and \eqn{X} are two real matrices of dimension \eqn{n \times p}{n x p}
and \eqn{n \times r}{n x r} respectively,
and \eqn{\|.\|_F}{|.|_F} is the Frobenius norm.

The algorithm is very fast compared to other approaches, as it is optimised
for handling multiple right-hand sides.
}
\details{
  Within the \code{NMF} package, this algorithm is used internally by the
SNMF/R(L) algorithm from \cite{Kim et al. (2007)} to solve general Nonnegative
Matrix Factorization (NMF) problems, using alternating nonnegative
constrained least-squares.
That is by iteratively and alternatively estimate each matrix factor.

The algorithm is an active/passive set method, which rearrange the
right-hand side to reduce the number of pseudo-inverse calculations.
It uses the unconstrained solution \eqn{K_u} obtained from the
unconstrained least squares problem,
i.e. \eqn{\min \|Y - X K\|_F^2}{min ||Y - X K||_F^2} , so as to determine
the initial passive sets.

The function \code{fcnnls} is provided separately so that it can be
used to solve other types of nonnegative least squares problem.
For faster computation, when multiple nonnegative least square fits
are needed, it is recommended to directly use the function \code{\link{.fcnnls}}.

The code of this function is a port from the original MATLAB code
provided by \cite{Kim et al. (2007)}.
}
\section{Methods}{
\describe{

\item{fcnnls}{\code{signature(x = "matrix", y = "matrix")}: This method wraps a call to the internal function \code{.fcnnls}, and
formats the results in a similar way as other lest-squares methods such
as \code{\link{lm}}.
}

\item{fcnnls}{\code{signature(x = "numeric", y = "matrix")}: Shortcut for \code{fcnnls(as.matrix(x), y, ...)}.
}

\item{fcnnls}{\code{signature(x = "ANY", y = "numeric")}: Shortcut for \code{fcnnls(x, as.matrix(y), ...)}.
}

}

}
\examples{
\dontshow{# roxygen generated flag
options(R_CHECK_RUNNING_EXAMPLES_=TRUE)
}

## Define a random nonnegative matrix matrix
n <- 200; p <- 20; r <- 3
V <- rmatrix(n, p)

## Compute the optimal matrix K for a given X matrix
X <- rmatrix(n, r)
res <- fcnnls(X, V)

## Compute the same thing using the Moore-Penrose generalized pseudoinverse
res <- fcnnls(X, V, pseudo=TRUE)

## It also works in the case of single vectors
y <- runif(n)
res <- fcnnls(X, y)
# or
res <- fcnnls(X[,1], y)
}
\author{
  Original MATLAB code : Van Benthem and Keenan

Adaption of MATLAB code for SNMF/R(L): H. Kim

Adaptation to the NMF package framework: Renaud Gaujoux
}
\references{
  Original MATLAB code from Van Benthem and Keenan, slightly modified by H.
Kim:\cr \url{http://www.cc.gatech.edu/~hpark/software/fcnnls.m}

Van Benthem M and Keenan MR (2004). "Fast algorithm for the solution of large-scale non-negativity-constrained least squares
problems." _Journal of Chemometrics_, *18*(10), pp. 441-450. ISSN 0886-9383, <URL: http://dx.doi.org/10.1002/cem.889>, <URL:
http://doi.wiley.com/10.1002/cem.889>.

Kim H and Park H (2007). "Sparse non-negative matrix factorizations via alternating non-negativity-constrained least squares
for microarray data analysis." _Bioinformatics (Oxford, England)_, *23*(12), pp. 1495-502. ISSN 1460-2059, <URL:
http://dx.doi.org/10.1093/bioinformatics/btm134>, <URL: http://www.ncbi.nlm.nih.gov/pubmed/17483501>.
}
\seealso{
  \code{\link{nmf}}
}
\keyword{methods}
\keyword{multivariate}
\keyword{optimize}
\keyword{regression}

