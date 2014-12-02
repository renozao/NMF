\docType{methods}
\name{canFit}
\alias{canFit}
\alias{canFit,character,ANY-method}
\alias{canFit-methods}
\alias{canFit,NMFStrategy,character-method}
\alias{canFit,NMFStrategy,NMF-method}
\title{Testing Compatibility of Algorithm and Models}
\usage{
  canFit(x, y, ...)

  \S4method{canFit}{NMFStrategy,character}(x, y,
    exact = FALSE)
}
\arguments{
  \item{x}{an object that describes an algorithm}

  \item{y}{an object that describes a model}

  \item{...}{extra arguments to allow extension}

  \item{exact}{for logical that indicates if an algorithm is considered able to fit
only the models that it explicitly declares (\code{TRUE}), or if it should be
considered able to also fit models that extend models that it explicitly fits.}
}
\description{
  \code{canFit} is an S4 generic that tests if an algorithm can
fit a particular model.
}
\section{Methods}{
\describe{

\item{canFit}{\code{signature(x = "NMFStrategy", y = "character")}: Tells if an NMF algorithm can fit a given class of NMF models
}

\item{canFit}{\code{signature(x = "NMFStrategy", y = "NMF")}: Tells if an NMF algorithm can fit the same class of models as \code{y}
}

\item{canFit}{\code{signature(x = "character", y = "ANY")}: Tells if a registered NMF algorithm can fit a given NMF model
}

}

}
\seealso{
  Other regalgo: \code{\link{nmfAlgorithm}}
}
\keyword{methods}

