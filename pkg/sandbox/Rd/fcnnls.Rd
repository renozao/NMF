\name{fcnnls}
\docType{methods}
\alias{fcnnls}
\alias{fcnnls,matrix,matrix-method}
\alias{fcnnls,numeric,matrix-method}
\alias{fcnnls,ANY,numeric-method}

\title{Fast Combinatorial Non-Negative Least-Square}
\description{
This function solves the following non-negative least square linear problem using normal equations and the 
fast combinatorial strategy from Benthem and Keenan (2004):

\deqn{%
  \begin{array}{l}%
  \min \|Y - X K\|_F\\%
  \mbox{s.t. } K>=0%
  \end{array}%
 }{min ||Y - X K||_F, s.t. K>=0}
where \eqn{\|.\|_F} is the Frobenius norm.

The resulting algorithm is very fast to converge compared to other approaches. 

Within the \code{NMF} package, this algorithm is used internally by the SNMF/R(L) algorithm from Kim and Park (2007) 
to solve general Nonnegative Matrix Factorization (NMF) problems, using alternating non-negative constrained least-squares.
That is by iteratively and alternatively estimate each matrix factor (see section \emph{References}).

It is provided separately so that it can be used to solve other types of non-negative least squares problem.
For faster computation, please the internal -- non-exported -- function \code{NMF:::.fcnnls}
The code is a port from the original MATLAB code used in Kim and Park (2007) (see references). 
   
}
\usage{
\S4method{fcnnls}{matrix,matrix}(x, y, verbose=FALSE, pseudo=TRUE, ...)
}
\arguments{
  \item{x}{ the coefficient matrix}
  \item{y}{ the target matrix to be approximated by \eqn{X K}.}
  \item{verbose}{ toggle verbosity (default is \code{FALSE}).}
  \item{pseudo}{ By default (\code{pseudo=FALSE}) the algorithm uses Gaussian elimination to solve the 
  successive internal linear problems, using the \code{\link{solve}} function.
  If \code{pseudo=TRUE} the algorithm uses Moore-Penrose generalized \code{\link[corpcor]{pseudoinverse}} 
  from the \code{corpcor} package instead of \link{solve}.}
  \item{...}{ extra arguments passed to the internal function \code{.fcnnls}. Currently not used.}
}
\details{
  Given two real matrices \eqn{Y} and \eqn{X}, of dimension \eqn{n \times p}{n x p} and \eqn{n \times r}{n x r} 
  respectively, this algorithm solves for the optimal nonnegative matrix \eqn{K} (\eqn{r \times p}{r x p}) such that:
  \deqn{%
  \begin{array}{l}%
  \min \|Y - X K\|_F\\%
  \mbox{s.t. } K>=0%
  \end{array}%
  }{min ||Y - X K||_F, s.t. K>=0}
  where \eqn{\|.\|_F} is the Frobenius norm.

  It is based on the active/passive set method. It uses the unconstrained solution \eqn{K_u} obtained from the 
  unconstrained least squares problem, i.e. \eqn{\min \|Y - X K\|_F^2}{min ||Y - X K||_F^2} , so as to determine the initial passive sets.
  
  
}
\value{
  The returned value is a list containing the following components:
  \item{x}{ the estimated optimal matrix \eqn{K}.}
  \item{fitted}{ the fitted matrix \eqn{X K}.}
  \item{residuals}{ the residual matrix \eqn{Y - X K}.}
  \item{deviance}{ the residual sum of squares between the fitted matrix \eqn{X K} and the target matrix \eqn{Y}. That is the sum of the square residuals.}
  \item{passive}{ a \eqn{r x p} logical matrix containing the passive set, that is the set of entries in \eqn{K} that are not null (i.e. strictly positive).}
  \item{pseudo}{ a logical that is \code{TRUE} if the computation was performed using the pseudoinverse. See argument \code{pseudo}.}
}

\references{
	M. H. van Benthem and M. R. Keenan (2004).
	\emph{Fast algorithm for the solution of large-scale non-negativity-constrained least squares problems}.
	J. Chemometrics 2004, \bold{18}:441-450.	

	Kim, H. and Park, H. (2007).
	\emph{Sparse non-negative matrix factorizations via alternating non-negativity-constrained least squares for microarray data analysis}.
	Bioinformatics 2007; \bold{23(12)}:1495-502.
	
Original MATLAB code from Van Benthem and Keenan, slightly modified by H. Kim:\cr
\url{http://www.cc.gatech.edu/~hpark/software/fcnnls.m} 

	
}
\author{ Renaud Gaujoux  }

\seealso{ \code{\link{nmf}} }
\examples{
## Define a random non-negative matrix matrix
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
\keyword{ optimize }
\keyword{ multivariate }
\keyword{ regression }
