\name{NMFSet-class}

\docType{class}
\alias{NMFSet-class}

\title{Deprecated Class to store results from multiple runs of NMF algorithms}
\description{
This class is deprecated and replaced by class \code{\linkS4class{NMFfitX}} and its 
extensions.
It remains only for backward compatibility and will be defunct in the next release.
  
It extends the base class \code{list} to store the result from a multiple run 
of NMF algorithms.

The elements are of class \code{NMF}.
}

\section{Slots}{
	 \describe{
    \item{\code{consensus}:}{Object of class \code{"matrix"} used to store the 
    consensus matrix when multiple runs have been performed with option 
    \code{keep.all=FALSE}. In this case, only the best factorization is returned, 
    so the object is of length 1. However the consensus matrix across all runs 
    is still computed and stored in this slot.}
    
    \item{\code{nrun}:}{an \code{integer} that contains the number of runs when 
    NMF is performed with option \code{keep.all=FALSE}.
    
    See \code{\link{nmf}}.
    }
    
    \item{\code{runtime}:}{Object of class \code{"proc_time"} that contains 
    various measures of the time spent to perform all the runs.}
    
    \item{\code{.Data}:}{standard slot that contains the S3 \code{list} object data.
    See R documentation on S4 classes for more details.}
    
    }
}

\section{Methods}{
   All the methods for this class have been removed from the package and are 
   substituted by methods for \code{NMFfitX} objects. 
}

\author{ Renaud Gaujoux  }

\seealso{
	\code{\linkS4class{NMFfitX}}, \code{\linkS4class{NMF}} , \code{\link{nmf-methods}}, \code{\link{nmf-multiple}}  
}

\keyword{classes}
