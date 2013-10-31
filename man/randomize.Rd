\name{randomize}
\alias{randomize}
\title{Randomizing Data}
\usage{
  randomize(x, ...)
}
\arguments{
  \item{x}{data to be permutated. It must be an object
  suitable to be passed to the function
  \code{\link{apply}}.}

  \item{...}{extra arguments passed to the function
  \code{\link{sample}}.}
}
\value{
  a matrix
}
\description{
  \code{randomize} permutates independently the entries in
  each column of a matrix-like object, to produce random
  data that can be used in permutation tests or bootstrap
  analysis.
}
\details{
  In the context of NMF, it may be used to generate random
  data, whose factorization serves as a reference for
  selecting a factorization rank, that does not overfit the
  data.
}
\examples{
x <- matrix(1:32, 4, 8)
randomize(x)
randomize(x)
}

