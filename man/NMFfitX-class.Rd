\docType{class}
\name{NMFfitX-class}
\alias{NMFfitX-class}
\title{Virtual Class to Handle Results from Multiple Runs of NMF Algorithms}
\description{
  This class defines a common interface to handle the
  results from multiple runs of a single NMF algorithm,
  performed with the \code{\link{nmf}} method.
}
\details{
  Currently, this interface is implemented by two classes,
  \code{\linkS4class{NMFfitX1}} and
  \code{\linkS4class{NMFfitXn}}, which respectively handle
  the case where only the best fit is kept, and the case
  where the list of all the fits is returned.

  See \code{\link{nmf}} for more details on the method
  arguments.
}
\section{Slots}{
  \describe{

  \item{runtime.all}{Object of class
  \code{\link[=proc.time]{proc_time}} that contains CPU
  times required to perform all the runs.}

  }
}
\section{Methods}{
  \describe{

  \item{basismap}{\code{signature(object = "NMFfitX")}:
  Plots a heatmap of the basis matrix of the best fit in
  \code{object}. }

  \item{coefmap}{\code{signature(object = "NMFfitX")}:
  Plots a heatmap of the coefficient matrix of the best fit
  in \code{object}.

  This method adds: \itemize{ \item an extra special column
  annotation track for multi-run NMF fits,
  \code{'consensus:'}, that shows the consensus cluster
  associated to each sample. \item a column sorting schema
  \code{'consensus'} that can be passed to argument
  \code{Colv} and orders the columns using the hierarchical
  clustering of the consensus matrix with average linkage,
  as returned by \code{\link{consensushc}(object)}. This is
  also the ordering that is used by default for the heatmap
  of the consensus matrix as ploted by
  \code{\link{consensusmap}}. } }

  \item{consensus}{\code{signature(object = "NMFfitX")}:
  Pure virtual method defined to ensure \code{consensus} is
  defined for sub-classes of \code{NMFfitX}. It throws an
  error if called. }

  \item{consensushc}{\code{signature(object = "NMFfitX")}:
  Compute the hierarchical clustering on the consensus
  matrix of \code{object}, or on the connectivity matrix of
  the best fit in \code{object}. }

  \item{consensusmap}{\code{signature(object = "NMFfitX")}:
  Plots a heatmap of the consensus matrix obtained when
  fitting an NMF model with multiple runs. }

  \item{cophcor}{\code{signature(object = "NMFfitX")}:
  Computes the cophenetic correlation coefficient on the
  consensus matrix of \code{object}. All arguments in
  \code{...} are passed to the method
  \code{cophcor,matrix}. }

  \item{deviance}{\code{signature(object = "NMFfitX")}:
  Returns the deviance achieved by the best fit object,
  i.e. the lowest deviance achieved across all NMF runs. }

  \item{dispersion}{\code{signature(object = "NMFfitX")}:
  Computes the dispersion on the consensus matrix obtained
  from multiple NMF runs. }

  \item{fit}{\code{signature(object = "NMFfitX")}: Returns
  the model object that achieves the lowest residual
  approximation error across all the runs.

  It is a pure virtual method defined to ensure \code{fit}
  is defined for sub-classes of \code{NMFfitX}, which
  throws an error if called. }

  \item{getRNG1}{\code{signature(object = "NMFfitX")}:
  Returns the RNG settings used for the first NMF run of
  multiple NMF runs. }

  \item{ibterms}{\code{signature(object = "NMFfitX")}:
  Method for multiple NMF fit objects, which returns the
  indexes of fixed basis terms from the best fitted model.
  }

  \item{metaHeatmap}{\code{signature(object = "NMFfitX")}:
  Deprecated method subsituted by
  \code{\link{consensusmap}}. }

  \item{minfit}{\code{signature(object = "NMFfitX")}:
  Returns the fit object that achieves the lowest residual
  approximation error across all the runs.

  It is a pure virtual method defined to ensure
  \code{minfit} is defined for sub-classes of
  \code{NMFfitX}, which throws an error if called. }

  \item{nmf.equal}{\code{signature(x = "NMFfitX", y =
  "NMF")}: Compares two NMF models when at least one comes
  from multiple NMF runs. }

  \item{NMFfitX}{\code{signature(object = "NMFfitX")}:
  Provides a way to aggregate \code{NMFfitXn} objects into
  an \code{NMFfitX1} object. }

  \item{nrun}{\code{signature(object = "NMFfitX")}: Returns
  the number of NMF runs performed to create \code{object}.

  It is a pure virtual method defined to ensure \code{nrun}
  is defined for sub-classes of \code{NMFfitX}, which
  throws an error if called.

  See \code{\link{nrun,NMFfitX-method}} for more details. }

  \item{predict}{\code{signature(object = "NMFfitX")}:
  Returns the cluster membership index from an NMF model
  fitted with multiple runs.

  Besides the type of clustering available for any NMF
  models (\code{'columns', 'rows', 'samples', 'features'}),
  this method can return the cluster membership index based
  on the consensus matrix, computed from the multiple NMF
  runs.

  See \code{\link{predict,NMFfitX-method}} for more
  details. }

  \item{residuals}{\code{signature(object = "NMFfitX")}:
  Returns the residuals achieved by the best fit object,
  i.e. the lowest residual approximation error achieved
  across all NMF runs. }

  \item{runtime.all}{\code{signature(object = "NMFfitX")}:
  Returns the CPU time required to compute all the NMF
  runs. It returns \code{NULL} if no CPU data is available.
  }

  \item{show}{\code{signature(object = "NMFfitX")}: Show
  method for objects of class \code{NMFfitX} }

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
# generate a synthetic dataset with known classes
n <- 20; counts <- c(5, 2, 3);
V <- syntheticNMF(n, counts)

# perform multiple runs of one algorithm (default is to keep only best fit)
res <- nmf(V, 3, nrun=3)
res

# plot a heatmap of the consensus matrix
\dontrun{ consensusmap(res) }

# perform multiple runs of one algorithm (keep all the fits)
res <- nmf(V, 3, nrun=3, .options='k')
res
}
\seealso{
  Other multipleNMF: \code{\link{NMFfitX1-class}},
  \code{\link{NMFfitXn-class}}
}

