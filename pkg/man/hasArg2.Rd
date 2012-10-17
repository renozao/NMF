\name{hasArg2}
\alias{hasArg2}
\title{Checking for Missing Arguments}
\usage{
  hasArg2(name)
}
\arguments{
  \item{name}{a character string.}
}
\description{
  This function is identical to \code{\link{hasArg}},
  except that it accepts the argument name as a character
  string. This avoids to have a check NOTE about invisible
  binding variable.
}
\examples{
f <- function(...){ hasArg2('abc') }
f(a=1)
f(abc=1)
f(b=1)
}

