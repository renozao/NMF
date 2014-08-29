\name{NMFStop}
\alias{NMFStop}
\alias{nmf.stop.connectivity}
\alias{nmf.stop.iteration}
\alias{nmf.stop.stationary}
\alias{nmf.stop.threshold}
\alias{stop-NMF}
\title{Stopping Criteria for NMF Iterative Strategies}
\usage{
  NMFStop(s, check = TRUE)

  nmf.stop.iteration(n)

  nmf.stop.threshold(threshold)

  nmf.stop.stationary(object, i, y, x,
    stationary.th = .Machine$double.eps,
    check.interval = 5 * check.niter, check.niter = 10L,
    ...)

  nmf.stop.connectivity(object, i, y, x, stopconv = 40,
    check.interval = 10, ...)
}
\arguments{
  \item{s}{specification of the stopping criterion.
See section \emph{Details} for the supported formats and how they are processed.}

  \item{check}{logical that indicates if the validity of the stopping criterion
function should be checked before returning it.}

  \item{n}{maximum number of iteration to perform.}

  \item{threshold}{default stationarity threshold}

  \item{object}{an NMF strategy object}

  \item{i}{the current iteration}

  \item{y}{the target matrix}

  \item{x}{the current NMF model}

  \item{stationary.th}{maximum absolute value of the gradient, for the objective
function to be considered stationary.}

  \item{check.interval}{interval (in number of iterations) on which the stopping
criterion is computed.}

  \item{check.niter}{number of successive iteration used to compute the stationnary
criterion.}

  \item{...}{extra arguments passed to the function \code{\link{objective}},
which computes the objective value between \code{x} and \code{y}.}

  \item{stopconv}{number of iterations intervals over which the connectivity
matrix must not change for stationarity to be achieved.}
}
\value{
  a function that can be passed to argument \code{.stop} of function
\code{\link{nmf}}, which is typically used when the algorith is implemented as
an iterative strategy.

a function that can be used as a stopping criterion for NMF algorithms
defined as \code{\linkS4class{NMFStrategyIterative}} objects.
That is a function with arguments \code{(strategy, i, target, data, ...)}
that returns \code{TRUE} if the stopping criterion is satisfied -- which in
turn stops the iterative process, and \code{FALSE} otherwise.
}
\description{
  The function documented here implement
  stopping/convergence criteria commonly used in NMF
  algorithms.

  \code{NMFStop} acts as a factory method that creates
  stopping criterion functions from different types of
  values, which are subsequently used by
  \code{\linkS4class{NMFStrategyIterative}} objects to
  determine when to stop their iterative process.

  \code{nmf.stop.iteration} generates a function that
  implements the stopping criterion that limits the number
  of iterations to a maximum of \code{n}), i.e. that
  returns \code{TRUE} if \code{i>=n}, \code{FALSE}
  otherwise.

  \code{nmf.stop.threshold} generates a function that
  implements the stopping criterion that stops when a given
  stationarity threshold is achieved by successive
  iterations. The returned function is identical to
  \code{nmf.stop.stationary}, but with the default
  threshold set to \code{threshold}.

  More precisely, the objective function is computed over
  \eqn{n} successive iterations (specified in argument
  \code{check.niter}), every \code{check.interval}
  iterations. The criterion stops when the absolute
  difference between the maximum and the minimum objective
  values over these iterations is lower than a given
  threshold \eqn{\alpha} (specified in
  \code{stationary.th}):

  \code{nmf.stop.connectivity} implements the stopping
  criterion that is based on the stationarity of the
  connectivity matrix.
}
\details{
  \code{NMFStop} can take the following values: \describe{
  \item{function}{ is returned unchanged, except when it
  has no arguments, in which case it assumed to be a
  generator, which is immediately called and should return
  a function that implements the actual stopping
  criterion;} \item{integer}{ the value is used to create a
  stopping criterion that stops at that exact number of
  iterations via \code{nmf.stop.iteration};}
  \item{numeric}{ the value is used to create a stopping
  criterion that stops when at that stationary threshold
  via \code{nmf.stop.threshold};} \item{character}{ must be
  a single string which must be an access key for
  registered criteria (currently available:
  \dQuote{connectivity} and \dQuote{stationary}), or the
  name of a function in the global environment or the
  namespace of the loading package.} }

  \deqn{ \left| \frac{\max_{i- N_s + 1 \leq k \leq i} D_k -
  \min_{i - N_s +1 \leq k \leq i} D_k}{n} \right| \leq
  \alpha, }{ | [max( D(i- N_s + 1), ..., D(i) ) - min( D(i-
  N_s + 1), ..., D(i) )] / n | <= alpha }
}

