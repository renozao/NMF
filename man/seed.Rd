\docType{methods}
\name{seed}
\alias{seed}
\alias{seed,ANY,ANY,character-method}
\alias{seed,ANY,ANY,function-method}
\alias{seed,ANY,ANY,missing-method}
\alias{seed,ANY,ANY,NULL-method}
\alias{seed,ANY,ANY,numeric-method}
\alias{seed,ANY,list,NMFSeed-method}
\alias{seed,ANY,numeric,NMFSeed-method}
\alias{seed,matrix,NMF,NMFSeed-method}
\alias{seed-methods}
\title{Interface for NMF Seeding Methods}
\usage{
  seed(x, model, method, ...)

  \S4method{seed}{matrix,NMF,NMFSeed}(x, model, method,
    rng, ...)

  \S4method{seed}{ANY,ANY,function}(x, model, method, name,
    ...)
}
\arguments{
  \item{x}{target matrix one wants to approximate with NMF}

  \item{model}{specification of the NMF model, e.g., the factorization rank.}

  \item{method}{specification of a seeding method.
See each method for details on the supported formats.}

  \item{...}{extra to allow extensions and passed down to the actual seeding method.}

  \item{rng}{rng setting to use.
If not missing the RNG settings are set and restored on exit using
\code{\link{setRNG}}.

All arguments in \code{...} are passed to teh seeding strategy.}

  \item{name}{optional name of the seeding method for custom seeding strategies.}
}
\value{
  an \code{\linkS4class{NMFfit}} object.
}
\description{
  The function \code{seed} provides a single interface for calling all seeding
methods used to initialise NMF computations.
These methods at least set the basis and coefficient matrices of the initial
\code{object} to valid nonnegative matrices.
They will be used as a starting point by any NMF algorithm that accept
initialisation.

IMPORTANT: this interface is still considered experimental and is subject
to changes in future release.
}
\section{Methods}{
\describe{

\item{seed}{\code{signature(x = "matrix", model = "NMF", method = "NMFSeed")}: This is the workhorse method that seeds an NMF model object using a given
seeding strategy defined by an \code{NMFSeed} object, to fit a given
target matrix.
}

\item{seed}{\code{signature(x = "ANY", model = "ANY", method = "function")}: Seeds an NMF model using a custom seeding strategy, defined by a function.

\code{method} must have signature \code{(x='NMFfit', y='matrix', ...)}, where
\code{x} is the unseeded NMF model and \code{y} is the target matrix to fit.
It must return an \code{\linkS4class{NMF}} object, that contains the seeded
NMF model.
}

\item{seed}{\code{signature(x = "ANY", model = "ANY", method = "missing")}: Seeds the model with the default seeding method given by
\code{nmf.getOption('default.seed')}
}

\item{seed}{\code{signature(x = "ANY", model = "ANY", method = "NULL")}: Use NMF method \code{'none'}.
}

\item{seed}{\code{signature(x = "ANY", model = "ANY", method = "numeric")}: Use \code{method} to set the RNG with \code{\link{setRNG}} and use method
\dQuote{random} to seed the NMF model.

Note that in this case the RNG settings are not restored.
This is due to some internal technical reasons, and might change in future
releases.
}

\item{seed}{\code{signature(x = "ANY", model = "ANY", method = "character")}: Use the registered seeding method whose access key is \code{method}.
}

\item{seed}{\code{signature(x = "ANY", model = "list", method = "NMFSeed")}: Seed a model using the elements in \code{model} to instantiate it with
\code{\link{nmfModel}}.
}

\item{seed}{\code{signature(x = "ANY", model = "numeric", method = "NMFSeed")}: Seeds a standard NMF model (i.e. of class \code{\linkS4class{NMFstd}}) of rank
\code{model}.
}

}

}
\keyword{methods}

