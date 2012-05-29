\docType{class}
\name{NMFStrategyIterative-class}
\alias{NMFStrategyIterative}
\alias{NMFStrategyIterative-class}
\title{Interface for Algorithms: Implementation for Iterative NMF Algorithms}
\description{
  This class provides a specific implementation for the
  function \code{run} -- completing the the interface class
  \code{\linkS4class{NMFStrategy}}, for NMF algorithms that
  conform to the following iterative schema:

  \itemize{ \item 1. Initialisation \item 2. Update the
  model at each iteration \item 3. Stop if some criterion
  is satisfied \item 4. Wrap up }

  This schema could possibly apply to all NMF algorithms,
  since these are essentially optimisation algorithms,
  almost all of which use iterative methods to approximate
  a solution of the optimisation problem. The main
  advantage is that it allows to implement updates and
  stopping criterion separately, and combine them in
  different ways. In particular, many NMF algorithms are
  based on multiplicative updates, following the approach
  from \cite{Lee et al. (2001)}, which are specially
  suitable to be cast into this simple schema.
}
\section{Slots}{
  \describe{

  \item{Update}{function that implement the update step,
  which computes new values for the model, based on its
  previous value. It is called at each iteration, until the
  stopping criterion is met or the maximum number of
  iteration is achieved.}

  \item{Stop}{function that implements the stopping
  criterion. It is called \emph{before} each Update step.}

  \item{WrapNMF}{function that wraps up the result into an
  NMF object. It is called just before returning the}

  }
}
\section{Methods}{
  \describe{

  \item{run}{\code{signature(object =
  "NMFStrategyIterative", y = "matrix", x = "NMFfit")}:
  Runs an NMF iterative algorithm on a target matrix
  \code{y}. }

  \item{show}{\code{signature( = "NMFStrategyIterative")}:
  Show method for objects of class
  \code{NMFStrategyIterative} }

  }
}
\references{
  Lee DD and Seung H (2001). "Algorithms for non-negative
  matrix factorization." _Advances in neural information
  processing systems_. <URL:
  http://scholar.google.com/scholar?q=intitle:Algorithms+for+non-negative+matrix+factorization\#0>.
}

