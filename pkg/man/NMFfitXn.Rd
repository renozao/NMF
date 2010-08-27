\name{NMFfitXn-class}

\docType{class}
\alias{NMFfitXn-class}

\alias{algorithm,NMFfitXn-method}
\alias{consensus,NMFfitXn-method}
\alias{dim,NMFfitXn-method}
\alias{entropy,NMFfitXn,ANY-method}
\alias{fit,NMFfitXn-method}
\alias{nrun,NMFfitXn-method}
\alias{purity,NMFfitXn,ANY-method}
\alias{residuals,NMFfitXn-method}

\alias{seqtime,NMFfitXn-method}

\alias{show,NMFfitXn-method}

\title{Class to Store the Result from Multiple Runs of a NMF Algorithm when Keeping All the Fits}
\description{
This class is used to return the result from a multiple run of a single NMF algorithm performed 
with function \code{nmf} with option \code{keep.all=TRUE} (cf. \code{\link{nmf}}).

It extends both classes \code{\linkS4class{NMFfitX}} and \code{list}, and stores the result of each run 
(i.e. a \code{NMFfit} object) in its \code{list} structure. 

IMPORTANT NOTE: 
This class is designed to be \strong{read-only}, even though all the \code{list}-methods 
can be used on its instances. Adding or removing elements would most probably lead to incorrect results in 
subsequent calls.
Capability for concatenating and merging NMF results is for the moment only used internally, and should be 
included and supported in the next release of the package.

}

\section{Slots}{
	 \describe{
	     
    	\item{\code{.Data}:}{standard slot that contains the S3 \code{list} object data.
    	See R documentation on S4 classes for more details.}
    
    }
}

\section{Extends}{
Class \code{"\linkS4class{NMFfitX}"}, directly.
Class \code{"\linkS4class{list}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "list", distance 2.
Class \code{"\link[Biobase]{AssayData-class}"}, by class "list", distance 2.
}

\section{Validity}{

    \code{NMFfitXn} objects are designed to store results of single runs of the same NMF algorithm, 
    that have the same dimensions.
    The following checks are performed in the class validity method:
    
    \itemize{
    \item All elements must be the result of a single NMF run. That is they must be of class \code{NMFfit}, 
    and objects of class \code{NMFfitX} are not allowed.
    
    \item All elements must be the result of the same NMF algorithm.
    
    \item The dimension of the fitted problem must be the same for all elements: same dimension of the 
    target matrix, and same factorisation rank. 
    
    }
}

\section{Methods}{
  \describe{
    
    \item{algorithm}{: 
    returns the name of the common algorithm used to compute all the runs.    
    
    Since all elements in the list are results from the same algorithm, the returned name is taken from 
    the first element. The method returns \code{NULL} if the list is empty.
    }
     
    \item{consensus}{:    
    Computes the consensus matrix associated with the list of runs, i.e. the mean connectivity matrix of all 
    the fits in the list.
	It's been proposed by \emph{Brunet et al. (2004)} to help visualising and measuring the stability 
	of the clusters obtained by NMF approaches. 
	See \code{\link{consensus}}.
    }
     
    \item{dim}{: 
    returns the common dimension of the NMF problem fitted by all the runs -- based on the dimension of 
    the first fit in the list.
    See \code{\link{dim,NMF-method}} for more details on the returned value.
    }
       
    \item{entropy}{: 
    computes the mean (resp. the best) entropy of the list of NMF fit. Which value is computed 
    depends on argument \code{method}.
    See \code{\link{entropy}} for more details.
    }
        
    \item{fit}{: returns the element of the list that achieves the lowest residual approximation error.
    See \code{\link{residuals}}.
    }
    
    \item{nrun}{: 
    returns the number of NMF runs performed to compute the object, i.e. its length -- as a list.
    }
    
    \item{purity}{:  
    computes the mean (resp. the best) purity of the list of NMF fit. Which value is computed 
    depends on argument \code{method}.
    See \code{\link{purity}} for details.
    }
        
    \item{residuals}{:
    computes the mean (resp. the best) residual error of the list of NMF fit. 
    Which value is computed depends on argument \code{method}.
    See \code{\link{residuals}}  for details.
    }
    
	\item{predict}{: 
    returns the predicted cluster index based on the \emph{best} NMF factorization in the list.    
    See \code{\link{predict}} for more details.
    } 
	    
    \item{runtime.all}{\code{signature(object = "NMFfitXn", null=FALSE, warning=TRUE)}: 
    returns the CPU time used to compute all the runs in the list, as stored in slot 
    \code{runtime.all} (inherited from class \code{NMFfitX}).
    
    When the computation is performed in parallel , the result may be very different 
    from the sequential computation time returned by the \code{seqtime} method (see below).
    
    When no time data is available in slot \code{runtime.all}, the \code{runtime.all} method for class 
    \code{NMFfitXn} differs from the one defined in its parent class \code{NMFfitX}.    
    Indeed, in the case no time data is stored in slot \code{runtime.all}, setting the extra argument 
    \code{null} to \code{FALSE} (default) forces the method to return the sequential computation time 
    instead and a warning is thrown unless argument \code{warning} is \code{FALSE}. Otherwise, in such a case, 
    a call with \code{null=TRUE} would return \code{NULL}.
    }
    
    \item{seqtime}{: 
    returns the sequential CPU time spent of all the runs in the list.
    It is the sum of the CPU time used to compute each run. It returns \code{NULL} if the list is empty.   
    }
        
    \item{show}{: \code{show} method for class 
    \code{NMFfitXn}. }
		
	} % end describe
	
	Besides these above methods, class \code{NMFfitXn} inherits all the methods from class 
	\code{\linkS4class{NMFfitX}} like: \code{summary}, \code{metaHeatmap}.
}

\author{ Renaud Gaujoux \email{renaud@cbio.uct.ac.za} }

\seealso{
	\code{\linkS4class{NMFfitX}} , \code{\link{nmf-methods}}, \code{\link{nmf-multiple}}  
}
\examples{

# generate a synthetic dataset with known classes
n <- 50; counts <- c(5, 5, 8);
V <- syntheticNMF(n, counts, noise=TRUE)

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

}
\keyword{classes}

