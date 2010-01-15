\name{NMFSet-class}
\Rdversion{1.1}
\docType{class}
\alias{NMFSet-class}
\alias{algorithm,NMFSet-method}
\alias{cophcor,NMFSet-method}
\alias{compare}
\alias{compare,NMFSet-method}
\alias{compare,list-method}
\alias{connectivity,NMFSet-method}
\alias{dispersion,NMFSet-method}
\alias{entropy,NMFSet,ANY-method}
\alias{errorPlot,NMFSet-method}

\alias{featureNames,NMFSet-method}

\alias{fit,NMFSet-method}
\alias{join}
\alias{join,list-method}
\alias{metaHeatmap,NMFSet-method}
\alias{nrun}
\alias{nrun,NMFSet-method}
\alias{purity,NMFSet,ANY-method}
\alias{predict,NMFSet-method}
\alias{residuals,NMFSet-method}
\alias{rss,NMFSet-method}
\alias{runtime,NMFSet-method}

\alias{sampleNames,NMFSet-method}

\alias{show,NMFSet-method}
\alias{summary,NMFSet-method}

\title{Class to store results from multiple runs of NMF algorithms}
\description{
Class that extends base class \code{list} to store the result from a multiple run 
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
  \describe{
    \item{cophcor}{\code{signature(object = "NMFSet")}: 
    
    Computes the cophenetic correlation coefficient of the consensus matrix 
    associated to the multiple NMF runs described by \code{object}.
	It's been proposed by \emph{Brunet et al. (2004)} to measure the stability 
	of the clusters obtained by NMF approaches. 
	See \code{\link{cophcor}} for more details.
    }
    
    \item{compare}{\code{signature(x = "NMFSet")}: compute summary measures for
    each element in the list and return them in rows in a \code{data.frame}. }
    
    \item{connectivity}{\code{signature(x = "NMFSet")}: compute the consensus 
    matrix of the set of results. That is the mean of the connectivity matrices 
    of each element. See \code{\link{connectivity}}.
    }
    
    \item{dispersion}{\code{signature(object = "NMFSet")}: 
    
    Computes the dispersion coefficient of the consensus matrix associated 
    to the multiple NMF runs described by \code{object}.
	It's been proposed by \emph{Kim and Park (2007)} to measure the reproducibility 
	of the clusters. 
	See \code{\link{dispersion}} for more details.
    }
    
    \item{entropy}{\code{signature(x = "NMFSet", class = "ANY")}: computes the 
    mean entropy of the set of NMF results. See \code{\link{entropy}}.}
    
    \item{featureNames}{\code{signature(object = "NMFSet")}: 
    returns the row names of the basis matrix from the best fit of the set of results.
    If BioConductor is installed this method is defined for the generic function 
    \code{\link[Biobase]{featureNames}} from the \code{Biobase} package.
    }
    
    \item{fit}{\code{signature(object = "NMFSet")}: returns the element that 
    achieves the lowest residual approximation error. }
    
    \item{residuals}{\code{signature(object = "NMFSet")}: returns the average 
    residuals of the set of results.}
    
    \item{errorPlot}{\code{signature(x = "NMFSet")}: plot on a single graphic 
    the residuals tracks for each of the run. The tracks are normalized against their 
    maximum value, to generate tracks that lies in range \eqn{[0,1]}.      
    }
    
    \item{join}{\code{signature(x = "list")}: 
    method to create a \code{NMFSet} object from a list of 
    \code{NMFfit} or \code{NMFSet} objects. }
    
    \item{nrun}{\code{signature(object = "NMFSet")}: 
    returns the number of runs performed to create \code{object}.
    This is not necessarily the length of \code{object} that is the number of 
    NMF fits contained in \code{object}, nor always the value of slot \code{nrun}.
    
    When performing multiple NMF runs with option \code{keep.all=FALSE}, only one 
    fit is kept so \code{length(object)=1}, but the number of runs is stored 
    in slot \code{nrun}. When \code{keep.all=TRUE}, slot \code{nrun} is not set 
    at all.
    
    See \code{\link{nmf}}.
    }
    
    \item{metaHeatmap}{ Produces a heatmap of the consensus matrix using 
	function \code{\link[gplots]{heatmap.2}}. See \code{\link{metaHeatmap}}.}
	
	\item{predict}{\code{signature(object = "NMFSet")}: 
    returns a \code{factor} that gives the predicted cluster index for each sample 
    (resp. for each feature) based on the \emph{best} NMF factorization in the set 
    of results \code{object}.
    The index correspond to the basis vector that most contribute to the sample 
    (resp. to which the feature contribute the most).
    See \code{\link{predict}} details on extra arguments.
    } 
	
    \item{purity}{\code{signature(x = "NMFSet", class = "ANY")}:  computes the 
    mean purity of the set of NMF results. See \code{\link{purity}}.}
    
    \item{rss}{\code{signature(object = "NMFSet")}: 
    computes the Residual Sum of Squares (RSS) of the best factorization in the 
    set.
    See \code{\link{rss}}.
    }
    
    \item{runtime}{\code{signature(object = "NMFSet")}: 
    returns the time spent to compute all the runs.
    }
    
    \item{sampleNames}{\code{signature(object = "NMFSet")}: 
    returns the column names of the mixture coefficient matrix from the best fit 
    of the set of results.
    If BioConductor is installed this method is defined for the generic function 
    \code{\link[Biobase]{featureNames}} from the \code{Biobase} package.
    }
    
    \item{show}{\code{signature(object = "NMFSet")}: \code{show} method for class 
    \code{NMFSet}. }
	
	\item{summary}{\code{signature(x = "NMFSet")}: 
    standard generic \code{summary} method for objects of class \code{NMFSet}. 
    It computes a set of measures to evaluate the quality of the best factorization 
    of the set. The result is similar to the result from the \code{summary} method 
    of \code{NMFfit} objects. See \code{\linkS4class{NMFfit}} for details on the computed 
    measures. 
    In addition, the cophenetic correlation coefficient and the dispersion coefficient
    of the consensus matrix are returned. See methods \code{cophcor} and 
    \code{dispersion} above.
	}
	
	} % end describe
}

\references{ 

 	\emph{Metagenes and molecular pattern discovery using matrix factorization}
	Brunet, J.~P., Tamayo, P., Golub, T.~R., and Mesirov, J.~P. (2004)
	Proc Natl Acad Sci U S A
	101(12), 4164--4169.

	\emph{Sparse non-negative matrix factorizations via alternating non-negativity-constrained least squares for microarray data analysis}
	Kim, H. & Park, H. (2007)
	Bioinformatics. 
	\url{http://dx.doi.org/10.1093/bioinformatics/btm134}.

}
\author{ Renaud Gaujoux \email{renaud@cbio.uct.ac.za} }

\seealso{
	\code{\linkS4class{NMF}} , \code{\link{nmf-methods}}, \code{\link{NMF-utils}}  
}
\examples{

# generate a synthetic dataset with known classes
n <- 50; counts <- c(5, 5, 8);
V <- syntheticNMF(n, counts, noise=TRUE)

# build the class factor
groups <- as.factor(do.call('c', lapply(seq(3), function(x) rep(x, counts[x]))))

# perform multiple runs of one algorithm
res <- nmf(V, 3, nrun=20)
res

# plot a heatmap of the consensus matrix
\dontrun{metaHeatmap(res)}

# run multiple methods
res <- nmf(V, 3, list('brunet', 'lee', 'nsNMF'))
# compare the results passing prior knowledge on classes
compare(res, class=groups)

}
\keyword{classes}
