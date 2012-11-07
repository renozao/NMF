\docType{class}
\name{NMFfitXn-class}
\alias{NMFfitXn-class}
\title{Structure for Storing All Fits from Multiple NMF Runs}
\description{
  This class is used to return the result from a multiple
  run of a single NMF algorithm performed with function
  \code{nmf} with option \code{keep.all=TRUE} (cf.
  \code{\link{nmf}}).
}
\details{
  It extends both classes \code{\linkS4class{NMFfitX}} and
  \code{list}, and stores the result of each run (i.e. a
  \code{NMFfit} object) in its \code{list} structure.

  IMPORTANT NOTE: This class is designed to be
  \strong{read-only}, even though all the
  \code{list}-methods can be used on its instances. Adding
  or removing elements would most probably lead to
  incorrect results in subsequent calls. Capability for
  concatenating and merging NMF results is for the moment
  only used internally, and should be included and
  supported in the next release of the package.
}
\section{Slots}{
  \describe{

  \item{.Data}{standard slot that contains the S3
  \code{list} object data. See R documentation on S3/S4
  classes for more details (e.g.,
  \code{\link{setOldClass}}).}

  }
}
\section{Methods}{
  \describe{

  \item{algorithm}{\code{signature(object = "NMFfitXn")}:
  Returns the name of the common NMF algorithm used to
  compute all fits stored in \code{object}

  Since all fits are computed with the same algorithm, this
  method returns the name of algorithm that computed the
  first fit. It returns \code{NULL} if the object is empty.
  }

  \item{basis}{\code{signature(object = "NMFfitXn")}:
  Returns the basis matrix of the best fit amongst all the
  fits stored in \code{object}. It is a shortcut for
  \code{basis(fit(object))}. }

  \item{coef}{\code{signature(object = "NMFfitXn")}:
  Returns the coefficient matrix of the best fit amongst
  all the fits stored in \code{object}. It is a shortcut
  for \code{coef(fit(object))}. }

  \item{compare}{\code{signature(object = "NMFfitXn")}:
  Compares the fits obtained by separate runs of NMF, in a
  single call to \code{\link{nmf}}. }

  \item{consensus}{\code{signature(object = "NMFfitXn")}:
  This method returns \code{NULL} on an empty object. The
  result is a matrix with several attributes attached, that
  are used by plotting functions such as
  \code{\link{consensusmap}} to annotate the plots. }

  \item{dim}{\code{signature(x = "NMFfitXn")}: Returns the
  dimension common to all fits.

  Since all fits have the same dimensions, it returns the
  dimension of the first fit. This method returns
  \code{NULL} if the object is empty. }

  \item{entropy}{\code{signature(x = "NMFfitXn", y =
  "ANY")}: Computes the best or mean entropy across all NMF
  fits stored in \code{x}. }

  \item{fit}{\code{signature(object = "NMFfitXn")}: Returns
  the best NMF fit object amongst all the fits stored in
  \code{object}, i.e. the fit that achieves the lowest
  estimation residuals. }

  \item{.getRNG}{\code{signature(object = "NMFfitXn")}:
  Returns the RNG settings used for the best fit.

  This method throws an error if the object is empty. }

  \item{getRNG1}{\code{signature(object = "NMFfitXn")}:
  Returns the RNG settings used for the first run.

  This method throws an error if the object is empty. }

  \item{minfit}{\code{signature(object = "NMFfitXn")}:
  Returns the best NMF model in the list, i.e. the run that
  achieved the lower estimation residuals.

  The model is selected based on its \code{deviance} value.
  }

  \item{modelname}{\code{signature(object = "NMFfitXn")}:
  Returns the common type NMF model of all fits stored in
  \code{object}

  Since all fits are from the same NMF model, this method
  returns the model type of the first fit. It returns
  \code{NULL} if the object is empty. }

  \item{nbasis}{\code{signature(x = "NMFfitXn")}: Returns
  the number of basis components common to all fits.

  Since all fits have been computed using the same rank, it
  returns the factorization rank of the first fit. This
  method returns \code{NULL} if the object is empty. }

  \item{nrun}{\code{signature(object = "NMFfitXn")}:
  Returns the number of runs performed to compute the fits
  stored in the list (i.e. the length of the list itself).
  }

  \item{purity}{\code{signature(x = "NMFfitXn", y =
  "ANY")}: Computes the best or mean purity across all NMF
  fits stored in \code{x}. }

  \item{runtime.all}{\code{signature(object = "NMFfitXn")}:
  If no time data is available from in slot
  \sQuote{runtime.all} and argument \code{null=TRUE}, then
  the sequential time as computed by \code{\link{seqtime}}
  is returned, and a warning is thrown unless
  \code{warning=FALSE}. }

  \item{seeding}{\code{signature(object = "NMFfitXn")}:
  Returns the name of the common seeding method used the
  computation of all fits stored in \code{object}

  Since all fits are seeded using the same method, this
  method returns the name of the seeding method used for
  the first fit. It returns \code{NULL} if the object is
  empty. }

  \item{seqtime}{\code{signature(object = "NMFfitXn")}:
  Returns the CPU time that would be required to
  sequentially compute all NMF fits stored in
  \code{object}.

  This method calls the function \code{runtime} on each fit
  and sum up the results. It returns \code{NULL} on an
  empty object. }

  \item{show}{\code{signature(object = "NMFfitXn")}: Show
  method for objects of class \code{NMFfitXn} }

  }
}
\examples{
# generate a synthetic dataset with known classes
n <- 50; counts <- c(5, 5, 8);
V <- syntheticNMF(n, counts)

# build the class factor
groups <- as.factor(do.call('c', lapply(seq(3), function(x) rep(x, counts[x]))))

# perform multiple runs of one algorithm, keeping all the fits
res <- nmf(V, 3, nrun=5, .options='k') # .options=list(keep.all=TRUE) also works
res

summary(res)
# get more info
summary(res, target=V, class=groups)

# compute/show computational times
runtime.all(res)
seqtime(res)

# plot the consensus matrix, computed on the fly
\dontrun{ consensusmap(res, annCol=groups) }
}
\seealso{
  Other multipleNMF: \code{\link{NMFfitX1-class}},
  \code{\link{NMFfitX-class}}
}

