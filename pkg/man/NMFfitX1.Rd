\name{NMFfitX1-class}

\docType{class}
\alias{NMFfitX1-class}

\alias{consensus,NMFfitX1-method}
\alias{fit,NMFfitX1-method}
\alias{minfit,NMFfitX1-method}
\alias{nrun,NMFfitX1-method}
\alias{runtime.all,NMFfitX1-method}

\alias{show,NMFfitX1-method}

\title{Class to Store the Result from Multiple Runs of a NMF Algorithm when Only the Best Fit is Kept}
\description{
This class is used to return the result from a multiple run of a single NMF algorithm performed 
with function \code{nmf} with the -- default -- option \code{keep.all=FALSE} (cf. \code{\link{nmf}}).

It extends both classes \code{\linkS4class{NMFfitX}} and \code{\linkS4class{NMFfit}}, and stores a 
the result of the best fit in its \code{NMFfit} structure.

Beside the best fit, this class allows to hold data about the computation of the multiple runs, such 
as the number of runs, the CPU time used to perform all the runs, as well as the consensus matrix.

Due to the inheritance from class \code{NMFfit}, objects of class \code{NMFfitX1} can be handled exactly 
as the results of single NMF run -- as if only the best run had been performed.
}

\section{Slots}{
  \describe{
  
    \item{\code{consensus}:}{ object of class \code{"matrix"} used to store the 
    consensus matrix based on all the runs.}
    
    \item{\code{nrun}:}{an \code{integer} that contains the number of runs performed
    to compute the object.
    }
    
    \item{\code{rng1}:}{an \code{\link{rstream}} object that contains RNG 
    settings used for the first run. See \code{\link{getRNG1}}.
    }
    
    \item{\code{runtime.all}:}{object of class \code{"proc_time"} that contains 
    various measures of the time spent to perform all the runs (inherited from \code{NMFfitX})}
  }
}

\section{Extends}{
Class \code{"\linkS4class{NMFfitX}"}, directly.
Class \code{"\linkS4class{NMFfit}"}, directly.
Class \code{"\linkS4class{NMF}"}, by class "NMFfit", distance 2.
}

\section{Validity}{
	There is currently no validity check for this class.    
}

\section{Methods}{
  \describe{
    
    \item{consensus}{:    
    returns the pre-computed consensus matrix associated with the runs. It is calculated during the NMF 
    computations itself, and is stored in slot \code{consensus}. If this slot is of length zero, then it 
    returns \code{NULL}. 
    
    It's been proposed by \emph{Brunet et al. (2004)} to help visualising and measuring the stability 
	of the clusters obtained by NMF approaches. 
	See \code{\link{consensus}}.
    }
    
    \item{fit}{: returns the NMF model of the best fit (i.e an \code{NMF} object), that is the one 
    achieving the best approximation error. Since \code{NMFfitX1} objects only store the best fit 
    of all the runs, it simply returns slot \code{fit} as inherited from \code{\linkS4class{NMFfit}}. 
    See \code{\link{fit}}.}
    
    \item{minfit}{: returns the result of the best fit as an \code{NMFfit} object. Since \code{NMFfitX1} objects only 
    store the best fit of all the runs, it simply coerce the object into a \code{\linkS4class{NMFfit}} object.
    See \code{\link{minfit}}.}
    
    \item{nrun}{: 
    returns the number of NMF runs performed to compute the object, as stored in slot \code{nrun}.
    }
       
    \item{runtime.all}{ 
    returns the CPU time used to compute all the runs in the list, as stored in slot 
    \code{runtime.all} (inherited from class \code{NMFfitX}). 
    }
            
    \item{show}{: \code{show} method for class \code{NMFfitX1}. }
		
	} % end describe
	
	Besides these above methods, class \code{NMFfitX1} inherits from all the methods from class 
	\code{\linkS4class{NMFfit}} and as such can be handled exactly as the result of a single NMF run.
}

\author{ Renaud Gaujoux  }

\seealso{
	\code{\linkS4class{NMFfitX}} , \code{\link{nmf-methods}}, \code{\link{nmf-multiple}}  
}
\examples{

# generate a synthetic dataset with known classes
n <- 50; counts <- c(5, 5, 8);
V <- syntheticNMF(n, counts, noise=TRUE)

# build the class factor
groups <- as.factor(do.call('c', lapply(seq(3), function(x) rep(x, counts[x]))))

# perform multiple runs of one algorithm, keeping only the best fit (default)
res <- nmf(V, 3, nrun=5) 
res
#NOTE: the implicit nmf options are .options=list(keep.all=FALSE) or .options='-k'

# compute summary measures
summary(res)
# get more info
summary(res, target=V, class=groups)

# show computational time
runtime.all(res)

}
\keyword{classes}

