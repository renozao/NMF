\name{nmfFormals}
\alias{nmfArgs}
\alias{nmfFormals}
\title{Showing Arguments of NMF Algorithms}
\usage{
  nmfFormals(x, ...)

  nmfArgs(x)
}
\arguments{
  \item{x}{algorithm specification}

  \item{...}{extra argument to allow extension}
}
\description{
  This function returns the extra arguments that can be
  passed to a given NMF algorithm in call to
  \code{\link{nmf}}.

  \code{nmfArgs} is a shortcut for
  \code{args(nmfWrapper(x))}, to display the arguments of a
  given NMF algorithm.
}
\examples{
# show arguments of an NMF algorithm
nmfArgs('brunet')
nmfArgs('snmf/r')
}

