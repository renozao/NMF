\docType{class}
\name{NMFfitX1-class}
\alias{NMFfitX1-class}
\title{Structure for Storing the Best Fit Amongst Multiple NMF Runs}
\description{
  This class is used to return the result from a multiple
  run of a single NMF algorithm performed with function
  \code{nmf} with the -- default -- option
  \code{keep.all=FALSE} (cf. \code{\link{nmf}}).
}
\details{
  It extends both classes \code{\linkS4class{NMFfitX}} and
  \code{\linkS4class{NMFfit}}, and stores a the result of
  the best fit in its \code{NMFfit} structure.

  Beside the best fit, this class allows to hold data about
  the computation of the multiple runs, such as the number
  of runs, the CPU time used to perform all the runs, as
  well as the consensus matrix.

  Due to the inheritance from class \code{NMFfit}, objects
  of class \code{NMFfitX1} can be handled exactly as the
  results of single NMF run -- as if only the best run had
  been performed.
}
\section{Slots}{
\describe{

\item{consensus}{object of class \code{matrix} used to store the
consensus matrix based on all the runs.}

\item{nrun}{an \code{integer} that contains the number of runs
performed to compute the object.}

\item{rng1}{an object that contains RNG settings used for the first
run. See \code{\link{getRNG1}}.}

}

}
\section{Methods}{
\describe{

\item{consensus}{\code{signature(object = "NMFfitX1")}: The result is the matrix stored in slot \sQuote{consensus}.
This method returns \code{NULL} if the consensus matrix is empty.
}


\item{fit}{\code{signature(object = "NMFfitX1")}: Returns the model object associated with the best fit, amongst all the
runs performed when fitting \code{object}.

Since \code{NMFfitX1} objects only hold the best fit, this method simply
returns the NMF model fitted by \code{object} -- that is stored in slot
\sQuote{fit}.
}


\item{getRNG1}{\code{signature(object = "NMFfitX1")}: Returns the RNG settings used to compute the first of all NMF runs, amongst
which \code{object} was selected as the best fit.
}


\item{minfit}{\code{signature(object = "NMFfitX1")}: Returns the fit object associated with the best fit, amongst all the
runs performed when fitting \code{object}.

Since \code{NMFfitX1} objects only hold the best fit, this method simply
returns \code{object} coerced into an \code{NMFfit} object.
}


\item{nmf.equal}{\code{signature(x = "NMFfitX1", y = "NMFfitX1")}: Compares the NMF models fitted by multiple runs, that only kept the best fits.
}


\item{nrun}{\code{signature(object = "NMFfitX1")}: Returns the number of NMF runs performed, amongst which \code{object} was
selected as the best fit.
}


\item{show}{\code{signature(object = "NMFfitX1")}: Show method for objects of class \code{NMFfitX1}
}

}

}
\examples{
\dontshow{# roxygen generated flag
options(R_CHECK_RUNNING_EXAMPLES_=TRUE)
}

# generate a synthetic dataset with known classes
n <- 20; counts <- c(5, 2, 3);
V <- syntheticNMF(n, counts)

# get the class factor
groups <- V$pData$Group

# perform multiple runs of one algorithm, keeping only the best fit (default)
#i.e.: the implicit nmf options are .options=list(keep.all=FALSE) or .options='-k'
res <- nmf(V, 3, nrun=3)
res

# compute summary measures
summary(res)
# get more info
summary(res, target=V, class=groups)

# show computational time
runtime.all(res)

# plot the consensus matrix, as stored (pre-computed) in the object
\dontrun{ consensusmap(res, annCol=groups) }
}
\seealso{
  Other multipleNMF: \code{\link{NMFfitX-class}}, \code{\link{NMFfitXn-class}}
}

