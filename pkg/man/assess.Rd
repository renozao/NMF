\docType{methods}
\name{summary}
\alias{summary}
\alias{summary-methods}
\alias{summary-NMF}
\alias{summary,NMFfit-method}
\alias{summary,NMFfitX-method}
\alias{summary,NMF-method}
\title{Assessing and Comparing NMF Models}
\usage{
  summary(object, ...)

  \S4method{summary}{NMF}(object, class, target)
}
\arguments{
  \item{object}{an NMF object. See available methods in
  section \emph{Methods}.}

  \item{...}{extra arguments passed to the next
  \code{summary} method.}

  \item{class}{known classes/cluster of samples specified
  in one of the formats that is supported by the functions
  \code{\link{entropy}} and \code{\link{purity}}.}

  \item{target}{target matrix specified in one of the
  formats supported by the functions \code{\link{rss}} and
  \code{\link{evar}}}
}
\description{
  The NMF package defines \code{summary} methods for
  different classes of objects, which helps assessing and
  comparing the quality of NMF models by computing a set of
  quantitative measures, e.g. with respect to their ability
  to recover known classes and/or the original target
  matrix.

  The most useful methods are for classes
  \code{\linkS4class{NMF}}, \code{\linkS4class{NMFfit}},
  \code{\linkS4class{NMFfitX}} and
  \code{\linkS4class{NMFList}}, which compute summary
  measures for, respectively, a single NMF model, a single
  fit, a multiple-run fit and a list of heterogenous fits
  performed with the function \code{\link{nmf}}.
}
\details{
  Due to the somehow hierarchical structure of the classes
  mentionned in \emph{Description}, their respective
  \code{summary} methods call each other in chain, each
  super-class adding some extra measures, only relevant for
  objects of a specific class.
}
\section{Methods}{
  \describe{

  \item{summary}{\code{signature(object = "NMF")}: Computes
  summary measures for a single NMF model.

  The following measures are computed:

  \describe{ \item{sparseness}{Sparseness of the
  factorization computed by the function
  \code{\link{sparseness}}.} \item{entropy}{Purity of the
  clustering, with respect to known classes, computed by
  the function \code{\link{purity}}.}
  \item{entropy}{Entropy of the clustering, with respect to
  known classes, computed by the function
  \code{\link{entropy}}.} \item{RSS}{Residual Sum of
  Squares computed by the function \code{\link{rss}}.}
  \item{evar}{Explained variance computed by the function
  \code{\link{evar}}.} } }

  \item{summary}{\code{signature(object = "NMFfit")}:
  Computes summary measures for a single fit from
  \code{\link{nmf}}.

  This method adds the following measures to the measures
  computed by the method \code{summary,NMF}:

  \describe{ \item{residuals}{Residual error as measured by
  the objective function associated to the algorithm used
  to fit the model.} \item{niter}{Number of iterations
  performed to achieve convergence of the algorithm.}
  \item{cpu}{Total CPU time required for the fit.}
  \item{cpu.all}{Total CPU time required for the fit. For
  \code{NMFfit} objects, this element is always equal to
  the value in \dQuote{cpu}, but will be different for
  multiple-run fits.} \item{nrun}{Number of runs performed
  to fit the model. This is always equal to 1 for
  \code{NMFfit} objects, but will vary for multiple-run
  fits.} } }

  \item{summary}{\code{signature(object = "NMFfitX")}:
  Computes a set of measures to help evaluate the quality
  of the \emph{best fit} of the set. The result is similar
  to the result from the \code{summary} method of
  \code{NMFfit} objects. See \code{\linkS4class{NMF}} for
  details on the computed measures. In addition, the
  cophenetic correlation (\code{\link{cophcor}}) and
  \code{\link{dispersion}} coefficients of the consensus
  matrix are returned, as well as the total CPU time
  (\code{\link{runtime.all}}). }

  }
}
\examples{
# random NMF model
x <- rnmf(3, 20, 12)
summary(x)
summary(x, gl(3, 4))
summary(x, target=rmatrix(x))
summary(x, gl(3,4), target=rmatrix(x))
# generate a synthetic dataset with known classes: 50 features, 18 samples (5+5+8)
n <- 50; counts <- c(5, 5, 8);
V <- syntheticNMF(n, counts)
cl <- unlist(mapply(rep, 1:3, counts))

# perform default NMF with rank=2
x2 <- nmf(V, 2)
summary(x2, cl, V)
# perform default NMF with rank=2
x3 <- nmf(V, 3)
summary(x2, cl, V)
}
\keyword{methods}

