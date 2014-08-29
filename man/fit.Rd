\docType{methods}
\name{fit}
\alias{fit}
\alias{fit<-}
\alias{fit<--methods}
\alias{fit-methods}
\alias{fit,NMFfit-method}
\alias{fit<-,NMFfit,NMF-method}
\alias{fit,NMFfitX1-method}
\alias{fit,NMFfitX-method}
\alias{fit,NMFfitXn-method}
\alias{minfit}
\alias{minfit-methods}
\alias{minfit,NMFfit-method}
\alias{minfit,NMFfitX1-method}
\alias{minfit,NMFfitX-method}
\alias{minfit,NMFfitXn-method}
\title{Extracting Fitted Models}
\usage{
  fit(object, ...)

  fit(object)<-value

  minfit(object, ...)
}
\arguments{
  \item{object}{an object fitted by some algorithm, e.g. as returned by the
function \code{\link{nmf}}.}

  \item{value}{replacement value}

  \item{...}{extra arguments to allow extension}
}
\description{
  The functions \code{fit} and \code{minfit} are S4
  genetics that extract the best model object and the best
  fit object respectively, from a collection of models or
  from a wrapper object.

  \code{fit<-} sets the fitted model in a fit object. It is
  meant to be called only when developing new NMF
  algorithms, e.g. to update the value of the model stored
  in the starting point.
}
\details{
  A fit object differs from a model object in that it
  contains data about the fit, such as the initial RNG
  settings, the CPU time used, etc\ldots, while a model
  object only contains the actual modelling data such as
  regression coefficients, loadings, etc\ldots

  That best model is generally defined as the one that
  achieves the maximum/minimum some quantitative measure,
  amongst all models in a collection.

  In the case of NMF models, the best model is the one that
  achieves the best approximation error, according to the
  objective function associated with the algorithm that
  performed the fit(s).
}
\section{Methods}{
\describe{

\item{fit}{\code{signature(object = "NMFfit")}: Returns the NMF model object stored in slot \code{'fit'}.
}

\item{fit}{\code{signature(object = "NMFfitX")}: Returns the model object that achieves the lowest residual approximation
error across all the runs.

It is a pure virtual method defined to ensure \code{fit} is defined
for sub-classes of \code{NMFfitX}, which throws an error if called.
}

\item{fit}{\code{signature(object = "NMFfitX1")}: Returns the model object associated with the best fit, amongst all the
runs performed when fitting \code{object}.

Since \code{NMFfitX1} objects only hold the best fit, this method simply
returns the NMF model fitted by \code{object} -- that is stored in slot
\sQuote{fit}.
}

\item{fit}{\code{signature(object = "NMFfitXn")}: Returns the best NMF fit object amongst all the fits stored in \code{object},
i.e. the fit that achieves the lowest estimation residuals.
}


\item{fit<-}{\code{signature(object = "NMFfit", value = "NMF")}: Updates the NMF model object stored in slot \code{'fit'} with a new value.
}


\item{minfit}{\code{signature(object = "NMFfit")}: Returns the object its self, since there it is the result of a single NMF run.
}

\item{minfit}{\code{signature(object = "NMFfitX")}: Returns the fit object that achieves the lowest residual approximation
error across all the runs.

It is a pure virtual method defined to ensure \code{minfit} is defined
for sub-classes of \code{NMFfitX}, which throws an error if called.
}

\item{minfit}{\code{signature(object = "NMFfitX1")}: Returns the fit object associated with the best fit, amongst all the
runs performed when fitting \code{object}.

Since \code{NMFfitX1} objects only hold the best fit, this method simply
returns \code{object} coerced into an \code{NMFfit} object.
}

\item{minfit}{\code{signature(object = "NMFfitXn")}: Returns the best NMF model in the list, i.e. the run that achieved the lower
estimation residuals.

The model is selected based on its \code{deviance} value.
}

}

}
\keyword{methods}

