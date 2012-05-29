\name{NMFStop}
\alias{NMFStop}
\alias{stop-nmf}
\alias{stop-NMF}
\title{Get Stopping Criterion}
\usage{
  NMFStop(val)
}
\arguments{
  \item{val}{access key that can be a character string, a
  single integer or numeric, or a function.}
}
\value{
  a function
}
\description{
  Convert a key into a stopping criterion (i.e. a
  \code{function}).
}
\details{
  Functions are returned unchanged, integer values are used
  to create a stopping criterion of that number of
  iterations, numeric values are used to create a stopping
  criterion of that stationary threshold.
}

