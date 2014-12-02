\name{utils-NMF}
\alias{str_args}
\alias{utils-NMF}
\title{Utility Function in the NMF Package}
\usage{
  str_args(x, exdent = 10L)
}
\arguments{
  \item{x}{a function}

  \item{exdent}{indentation for extra lines if the output takes more than one line.}
}
\description{
  Utility Function in the NMF Package

\code{str_args} formats the arguments of a function using \code{\link{args}},
but returns the output as a string.
}
\examples{
\dontshow{# roxygen generated flag
options(R_CHECK_RUNNING_EXAMPLES_=TRUE)
}

args(library)
str_args(library)
}

