\docType{methods}
\name{ibterms}
\alias{bterms}
\alias{cterms}
\alias{ibasis}
\alias{ibterms}
\alias{ibterms-methods}
\alias{ibterms,NMFfit-method}
\alias{ibterms,NMFfitX-method}
\alias{ibterms,NMF-method}
\alias{ibterms,NMFstd-method}
\alias{icoef}
\alias{icterms}
\alias{icterms-methods}
\alias{icterms,NMFfit-method}
\alias{icterms,NMF-method}
\alias{icterms,NMFstd-method}
\alias{iterms}
\alias{nbterms}
\alias{ncterms}
\alias{nterms}
\title{Fixed Terms in NMF Models}
\usage{
  ibterms(object, ...)

  icterms(object, ...)

  iterms(object, ...)

  nterms(object)

  nbterms(object)

  ncterms(object)

  bterms(object)

  cterms(object)

  ibasis(object, ...)

  icoef(object, ...)
}
\arguments{
  \item{object}{NMF object}

  \item{...}{extra parameters to allow extension (currently
  not used)}
}
\description{
  Formula-based NMF models may contain fixed basis and/or
  coefficient terms. The functions documented here provide
  access to these data, which are read-only and defined
  when the model object is instantiated (e.g., see
  \code{\link[=nmfModel,formula,ANY-method]{nmfModel,formula-method}}).

  \code{ibterms}, \code{icterms} and \code{iterms}
  respectively return the indexes of the fixed basis terms,
  the fixed coefficient terms and all fixed terms, within
  the basis and/or coefficient matrix of an NMF model.

  \code{nterms}, \code{nbterms}, and \code{ncterms} return,
  respectively, the number of all fixed terms, fixed basis
  terms and fixed coefficient terms in an NMF model. In
  particular: i.e. \code{nterms(object) = nbterms(object) +
  ncterms(object)}.

  \code{bterms} and \code{cterms} return, respectively, the
  primary data for fixed basis and coefficient terms in an
  NMF model -- as stored in slots \code{bterms} and
  \code{cterms} . These are factors or numeric vectors
  which define fixed basis components, e.g., used for
  defining separate offsets for different \emph{a priori}
  groups of samples, or to incorporate/correct for some
  known covariate.

  \code{ibasis} and \code{icoef} return, respectively, the
  indexes of all latent basis vectors and estimated
  coefficients within the basis or coefficient matrix of an
  NMF model.
}
\section{Methods}{
  \describe{

  \item{ibterms}{\code{signature(object = "NMF")}: Default
  pure virtual method that ensure a method is defined for
  concrete NMF model classes. }

  \item{ibterms}{\code{signature(object = "NMFstd")}:
  Method for standard NMF models, which returns the integer
  vector that is stored in slot \code{ibterms} when a
  formula-based NMF model is instantiated. }

  \item{ibterms}{\code{signature(object = "NMFfit")}:
  Method for single NMF fit objects, which returns the
  indexes of fixed basis terms from the fitted model. }

  \item{ibterms}{\code{signature(object = "NMFfitX")}:
  Method for multiple NMF fit objects, which returns the
  indexes of fixed basis terms from the best fitted model.
  }

  \item{icterms}{\code{signature(object = "NMF")}: Default
  pure virtual method that ensure a method is defined for
  concrete NMF model classes. }

  \item{icterms}{\code{signature(object = "NMFstd")}:
  Method for standard NMF models, which returns the integer
  vector that is stored in slot \code{icterms} when a
  formula-based NMF model is instantiated. }

  \item{icterms}{\code{signature(object = "NMFfit")}:
  Method for single NMF fit objects, which returns the
  indexes of fixed coefficient terms from the fitted model.
  }

  }
}
\keyword{methods}

