\name{staticVar}
\alias{staticVar}
\title{Get/Set a Static Variable in NMF Algorithms}
\usage{
  staticVar(name, value, init = FALSE)
}
\arguments{
  \item{name}{Name of the static variable (as a single
  character string)}

  \item{value}{New value of the static variable}

  \item{init}{a logical used when a \code{value} is
  provided, that specifies if the variable should be set to
  the new value only if it does not exist yet
  (\code{init=TRUE}).}
}
\value{
  The value of the static variable
}
\description{
  This function is used in iterative NMF algorithms to
  manage variables stored in a local workspace, that are
  accessible to all functions that define the iterative
  schema described in
  \code{\linkS4class{NMFStrategyIterative}}.

  It is specially useful for computing stopping criteria,
  which often require model data from different iterations.
}

