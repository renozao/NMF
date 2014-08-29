\name{advanced-NMF}
\alias{advanced-NMF}
\alias{which.best}
\title{Advanced Usage of the Package NMF}
\usage{
  which.best(object, FUN = deviance, ...)
}
\arguments{
  \item{object}{an NMF model fitted by multiple runs.}

  \item{FUN}{the function that computes the quantitative measure.}

  \item{...}{extra arguments passed to \code{FUN}.}
}
\description{
  The functions documented here provide advanced
  functionalities useful when developing within the
  framework implemented in the NMF package.

  \code{which.best} returns the index of the best fit in a
  list of NMF fit, according to some quantitative measure.
  The index of the fit with the lowest measure is returned.
}

