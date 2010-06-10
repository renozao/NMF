\name{NMFfitX-class}

\docType{class}
\alias{NMFfitX-class}

\alias{consensus,NMFfitX-method}

\alias{featureNames,NMFfitX-method}
\alias{sampleNames,NMFfitX-method}


\alias{nrun,NMFfitX-method}
\alias{show,NMFfitX-method}

\title{Virtual Class to Handle Results from Multiple Runs of a NMF Algorithms}
\description{
This class defines a common interface to handle the results from multiple runs of a single 
 NMF algorithm, performed with the \code{\link{nmf}} method.
 
Currently, this interface is implemented by two classes, \code{\linkS4class{NMFfitX1}} and 
\code{\linkS4class{NMFfitXn}}, which respectively handle the case where only the best fit 
is kept, and the case where the list of all the fits is returned.

See \code{\link{nmf-multiple}} for more details on the method arguments.
}

\section{Slots}{
	 \describe{
   
    \item{\code{runtime.all}:}{Object of class \code{"proc_time"} that contains 
    various measures of the time spent to perform all the runs.}
        
    }
}

\section{Methods}{
  \describe{
  
  	\item{consensus}{: 
    
    Computes the consensus matrix associated to the multiple NMF runs described by \code{object}.
	It's been proposed by \emph{Brunet et al. (2004)} to help visualising and measuring the stability 
	of the clusters obtained by NMF approaches. 
	See \code{\link{consensus}}.
	
	\emph{Technical note: this method is defined as a pure virtual method is the sense that an error is thrown 
	if it is not overloaded by the classes that implement the interface (i.e. that extends class \code{NMFfitX}).}  
    }
    
    \item{cophcor}{: 
    
    Computes the cophenetic correlation coefficient of the consensus matrix 
    associated to the multiple NMF runs described by the object.
	It's been proposed by \emph{Brunet et al. (2004)} to measure the stability 
	of the clusters obtained by NMF approaches. 
	See \code{\link{cophcor}} for more details.
    }
        
    \item{dispersion}{: 
    
    Computes the dispersion coefficient of the consensus matrix associated 
    to the multiple NMF runs described by the object.
	It's been proposed by \emph{Kim and Park (2007)} to measure the reproducibility 
	of the clusters. 
	See \code{\link{dispersion}} for more details.
    }
    
    \item{featureNames}{: 
    returns the row names of the basis matrix from the best fit of the set of results.
    
    If BioConductor is installed this method is defined for the generic function 
    \code{\link[Biobase]{featureNames}} from the \code{Biobase} package.
    }
       
    \item{fit}{: returns the element that achieves the lowest residual approximation error across
    all the runs. 
    See \code{\link{fit}} for more details.
    
    \emph{Technical note: this method is defined as a pure virtual method is the sense that an error is thrown 
	if it is not overloaded by the classes that implement the interface (i.e. that extends class \code{NMFfitX}).}  
    }
          
    %\item{join}{: method to create a \code{NMFfitX1} object from a list of \code{NMFfit} or \code{NMFSet} objects. }
    
    \item{nrun}{: 
    returns the number of runs performed to create the object.
    
    Note that because the \code{nmf} method allows to run the NMF computation keeping only the best fit, 
    \code{nrun} may return a value greater than one, while only the result of the best run is stored in 
    the object (cf. option \code{'k'} in method \code{\link{nmf}}). 
      
    See \code{\link{nmf}} and \code{\linkS4class{NMFfitX1}}.
    
    \emph{Technical note: this method is defined as a pure virtual method is the sense that an error is thrown 
	if it is not overloaded by the classes that implement the interface (i.e. that extends class \code{NMFfitX}).}
    }
    
    \item{metaHeatmap}{ Produces a heatmap of the consensus matrix using 
	function \code{\link[gplots]{heatmap.2}}. See \code{\link{metaHeatmap}}.}
		
    \item{runtime.all}{: 
    returns the total time spent to compute all the runs.
    See \code{\link{runtime.all}} for more details.
    }
    
    \item{sampleNames}{: 
    returns the column names of the mixture coefficient matrix from the best fit 
    of the set of results.
    If BioConductor is installed this method is defined for the generic function 
    \code{\link[Biobase]{featureNames}} from the \code{Biobase} package.
    }
    
    \item{show}{: \code{show} method for objects of class \code{NMFfitX}. }
	
	\item{summary}{: 
    standard generic \code{summary} method for objects of class \code{NMFfitX}. 
    It computes a set of measures to evaluate the quality of the \emph{best factorization} 
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
	\code{\link{nmf-methods}}, \code{\link{nmf-multiple}}, \code{\linkS4class{NMFfitX1}}, \code{\linkS4class{NMFfitXn}}  
}
\examples{

# generate a synthetic dataset with known classes
n <- 50; counts <- c(5, 5, 8);
V <- syntheticNMF(n, counts, noise=TRUE)

# perform multiple runs of one algorithm (default is to keep only best fit)
res <- nmf(V, 3, nrun=5)
str(res)

# perform multiple runs of one algorithm (keep all the fits)
res <- nmf(V, 3, nrun=5, .options='k')
str(res)

}
\keyword{classes}
