\docType{methods}
\name{nmf.equal}
\alias{nmf.equal}
\alias{nmf.equal,list,list-method}
\alias{nmf.equal,list,missing-method}
\alias{nmf.equal-methods}
\alias{nmf.equal,NMFfit,NMFfit-method}
\alias{nmf.equal,NMFfit,NMF-method}
\alias{nmf.equal,NMFfitX1,NMFfitX1-method}
\alias{nmf.equal,NMFfitX,NMF-method}
\alias{nmf.equal,NMF,NMFfit-method}
\alias{nmf.equal,NMF,NMFfitX-method}
\alias{nmf.equal,NMF,NMF-method}
\title{Testing Equality of NMF Models}
\usage{
  nmf.equal(x, y, ...)

  \S4method{nmf.equal}{NMF,NMF}(x, y, identical = TRUE,
    ...)

  \S4method{nmf.equal}{list,list}(x, y, ..., all = FALSE,
    vector = FALSE)
}
\arguments{
  \item{x}{an NMF model or an object that is associated
  with an NMF model, e.g. the result from a fit with
  \code{\link{nmf}}.}

  \item{y}{an NMF model or an object that is associated
  with an NMF model, e.g. the result from a fit with
  \code{\link{nmf}}.}

  \item{identical}{a logical that indicates if the
  comparison should be made using the function
  \code{\link{identical}} (\code{TRUE}) or
  \code{\link{all.equal}} (\code{FALSE}). See description
  for method \code{nmf.equal,NMF,NMF}.}

  \item{...}{extra arguments to allow extension, and passed
  to subsequent calls}

  \item{all}{a logical that indicates if all fits should be
  compared separately or only the best fits}

  \item{vector}{a logical, only used when \code{all=TRUE},
  that indicates if all fits must be equal for \code{x} and
  \code{y} to be declared equal, or if one wants to return
  the result of each comparison in a vector.}
}
\description{
  The function \code{nmf.equal} tests if two NMF models are
  the same, i.e. they contain -- almost -- identical data:
  same basis and coefficient matrices, as well as same
  extra parameters.
}
\details{
  \code{nmf.equal} compares two NMF models, and return
  \code{TRUE} iff they are identical acording to the
  function \code{\link{identical}} when
  \code{identical=TRUE}, or equal up to some tolerance
  acording to the function \code{\link{all.equal}}. This
  means that all data contained in the objects are
  compared, which includes at least the basis and
  coefficient matrices, as well as the extra parameters
  stored in slot \sQuote{misc}.

  If extra arguments are specified in \code{...}, then the
  comparison is performed using \code{\link{all.equal}},
  irrespective of the value of argument \code{identical}.
}
\section{Methods}{
  \describe{

  \item{nmf.equal}{\code{signature(x = "NMF", y = "NMF")}:
  Compares two NMF models.

  Arguments in \code{...} are used only when
  \code{identical=FALSE} and are passed to
  \code{all.equal}. }

  \item{nmf.equal}{\code{signature(x = "NMFfit", y =
  "NMF")}: Compares two NMF models when at least one comes
  from a NMFfit object, i.e. an object returned by a single
  run of \code{\link{nmf}}. }

  \item{nmf.equal}{\code{signature(x = "NMF", y =
  "NMFfit")}: Compares two NMF models when at least one
  comes from a NMFfit object, i.e. an object returned by a
  single run of \code{\link{nmf}}. }

  \item{nmf.equal}{\code{signature(x = "NMFfit", y =
  "NMFfit")}: Compares two fitted NMF models, i.e. objects
  returned by single runs of \code{\link{nmf}}. }

  \item{nmf.equal}{\code{signature(x = "NMFfitX", y =
  "NMF")}: Compares two NMF models when at least one comes
  from multiple NMF runs. }

  \item{nmf.equal}{\code{signature(x = "NMF", y =
  "NMFfitX")}: Compares two NMF models when at least one
  comes from multiple NMF runs. }

  \item{nmf.equal}{\code{signature(x = "NMFfitX1", y =
  "NMFfitX1")}: Compares the NMF models fitted by multiple
  runs, that only kept the best fits. }

  \item{nmf.equal}{\code{signature(x = "list", y =
  "list")}: Compares the results of multiple NMF runs.

  This method either compare the two best fit, or all fits
  separately. All extra arguments in \code{...} are passed
  to each internal call to \code{nmf.equal}. }

  \item{nmf.equal}{\code{signature(x = "list", y =
  "missing")}: Compare all elements in \code{x} to
  \code{x[[1]]}. }

  }
}
\keyword{methods}

