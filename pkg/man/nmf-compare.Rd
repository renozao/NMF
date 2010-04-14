\name{Comparing the results of different NMF runs}

\docType{methods}
\alias{nmf-compare}

\alias{as.NMFList}
\alias{as.NMFList-methods}
\alias{as.NMFList,ANY-method}

\alias{compare}
\alias{compare-methods}
\alias{compare,list-method}

\alias{plot,NMFList,missing-method}
\alias{plot,NMFList-method}

\alias{summary,NMFList-method}

\title{Comparing Results from Different NMF Runs}

\description{

This functions allow to compare the results of different NMF runs.
The results do not need to be from the same algorithm, nor even of 
the same dimension. 

}

\usage{

\S4method{as.NMFList}{ANY}(..., unlist=FALSE)

\S4method{compare}{list}(object, ..., unlist=FALSE)

\S4method{plot}{NMFList}(x, ...)

\S4method{summary}{NMFList}(object, sort.by=NULL, select=NULL, ...)

}

\arguments{
		
	\item{object}{ A \code{list} or an object of class \code{\linkS4class{NMFList}}. }	
	
	\item{select}{ the columns to be output in the result \code{data.frame}. 
	The column are given by their names (partially matched).
	The column names are the names of the summary measures returned by the 
	\code{summary} methods of the corresponding NMF results.
	}
	
	\item{sort.by}{ the sorting criteria, i.e. a partial match of a column name, 
	by which the result \code{data.frame} is sorted. 
	The sorting direction (increasing or decreasing) is computed internally depending on 
	the chosen criteria (e.g. decreasing for the cophenetic coefficient, increasing for the 
	residuals).	
	}
	
	\item{unlist}{boolean to specify if the arguments should be unlisted before wrapping them 
	into a \code{NMFList} object or comparing them.}
	
	\item{x}{ An object of class \code{\linkS4class{NMFList}}. }
		
	\item{...}{ Used to pass extra arguments to subsequent calls:
		\itemize{
			\item in \code{as.NMFList} the list of NMF results to wrap into a \code{NMFList} object.
			\item in \code{plot}: graphical parameters passed to the \code{plot} function.
			\item in \code{compare} and \code{summary}: extra arguments passed to the 
			\code{summary} method of each result object (cf. \code{\link{summary,NMF-method}}).			
		}
	}
}


\details{

\describe{
	
	\item{as.NMFList}{: wrap the arguments into a \code{\linkS4class{NMFList}} object.}
	
	\item{compare}{: shortcut for \code{summary(as.NMFList(object)} (cf. \code{summary} method below)}. 
     
     \item{plot}{: plot on a single graph the residuals tracks for each element in \code{x}. 
    See function \code{\link{nmf}} for details on how to enable the tracking of residuals.}
       
	\item{runtime}{: returns the computational time used to compute all the results in the list, 
	as stored in slot \code{runtime} of \code{object}.	
	 	
	The time is computed using the function \code{\link{system.time}} which returns object of class 
    \code{\link[=proc.time]{proc_time}}.
        
    Note that argument \code{...} is not used.
    }

	\item{summary}{: \code{summary} method for objects of class \code{NMFList}.
	
	It compute summary measures for each NMF result in the list and return them 
	in rows in a \code{data.frame}. 
	By default all the measures are included in the result, and \code{NA} values 
	are used where no data is available or the measure does not apply to the result object 
	(e.g. the dispersion for single NMF runs is not meaningful).
    This method is very useful to compare and evaluate the performance of 
    different algorithms.        
	}


}

}

\seealso{	
	\linkS4class{NMFfitX1}, \linkS4class{NMFfitXn}, \link[=summary,NMF-method]{summary}  
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

\examples{


# generate a synthetic dataset with known classes: 50 features, 18 samples (5+5+8)
n <- 50; counts <- c(5, 5, 8);
V <- syntheticNMF(n, counts, noise=TRUE)
\dontrun{metaHeatmap(V)}

# build the class factor
groups <- as.factor(do.call('c', lapply(seq(3), function(x) rep(x, counts[x]))))

# perform multiple runs of NMF (keep best only)
res <- nmf(V, 3, nrun=5)
res

# compute summary measures
summary(res)

# compute more summary measures
summary(res, target=V, class=groups)

# plot a heatmap of the consensus matrix with extra annotations
\dontrun{metaHeatmap(res, class=groups)}

# retrieve the predicted clusters of samples
predict(res)

# perform multiple runs of NMF and keep all the runs
res <- nmf(V, 3, nrun=5, .options='k')
res

# extract best fit
fit(res)

# compute/show computational times
runtime.all(res)
seqtime(res)


}


