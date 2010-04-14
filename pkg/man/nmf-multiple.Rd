\name{Handling the results of multiple NMF runs}

\docType{methods}
\alias{nmf-multiple}

\alias{consensus}
\alias{consensus-methods}

\alias{cophcor}
\alias{cophcor-methods}
\alias{cophcor,matrix-method}
\alias{cophcor,NMFfitX-method}

\alias{dispersion}
\alias{dispersion-methods}
\alias{dispersion,matrix-method}
\alias{dispersion,NMFfitX-method}

\alias{fit,NMFfitX-method}

\alias{metaHeatmap,NMFfitX-method}

\alias{nrun}
\alias{nrun-methods}

\alias{predict,NMFfitXn-method}

\alias{runtime.all}
\alias{runtime.all-methods}
\alias{runtime.all,NMFfitX-method}
\alias{runtime.all,NMFfitXn-method}

\alias{seqtime}
\alias{seqtime-methods}    
    
\alias{summary,NMFfitX-method}

\title{Handling Results from Multiple NMF Runs}

\description{
The NMF package provides an easy way to perform multiple runs of a given NMF algorithm on 
a target matrix.   

The result from the \code{\link{nmf}} method is a \code{\linkS4class{NMFfitX}} object that 
holds either all or only the best run, depending on the running options:

\code{# keep only the best run}
\code{object <- nmf(X, r, nrun=20)} 
\code{# keep all the runs}
\code{object <- nmf(X, r, nrun=20, .options='k')}

The methods documented here are used to handle such results. They are usually independent of 
the type of result and can be used without change in either situation (all runs kept or 
only the best one).

Note that when only the best result is kept, the result object conveniently inherits from 
all the methods available for single runs. Therefore it can be handled as if it had been 
computed by a single NMF run and all the methods defined for such results can be used 
(cf. \code{\linkS4class{NMFfit}} and \code{\link{NMF-utils}}).

See \code{\linkS4class{NMFfitXn}} and \code{\linkS4class{NMFfitX1}} for details on the classes 
that implement respectively the case where all the runs are kept and only the best run is kept.
}

\usage{

consensus(object, ...)

cophcor(object, ...)

dispersion(object, ...)

\S4method{fit}{NMFfitX}(object)

nrun(object, ...)

\S4method{metaHeatmap}{NMFfitX}(object, ...)

\S4method{predict}{NMFfitXn}(object, ...)

\S4method{runtime.all}{NMFfitX}(object)

\S4method{runtime.all}{NMFfitXn}(object, null=FALSE, warning=TRUE)

seqtime(object, ...)

\S4method{summary}{NMFfitX}(object, ...)

}

\arguments{
	
	\item{null}{ used in method \code{runtime.all} for \code{NMFfitXn} objects to specify if 
	the result should be \code{NULL} when the object has no time data is stored the total 
	computation time. In this case, if \code{null=FALSE} (default), the method returns the 
	sequential time (cf. \code{seqtime} below) instead of \code{NULL}. It also emits a warning 
	which can be toggle with argument \code{warning}.
	}
	
	\item{object}{ A \code{matrix} or an object that inherits from class 
	\code{\linkS4class{NMFfitX}} or \code{\linkS4class{NMFfit}} -- depending on the method. }
	
	\item{warning}{ used in method \code{runtime.all} for \code{NMFfitXn} objects to specify 
	if a warning should be emitted when the object has no time data the total computation time 
	and the sequential time is returned instead of \code{NULL} (cf. argument \code{null}).  
	}	
	
	\item{...}{ Used to pass extra arguments to subsequent calls:
		\itemize{
			\item in \code{metaHeatmap}: graphical parameters passed to function \code{\link[gplots]{heatmap.2}}
			\item in \code{predict}: extra arguments  passed to function \code{\link{predict,NMF-method}}
			\item in \code{summary}: extra arguments like \code{target} or \code{class} passed to 
			the method \code{\link{summary,NMFfit-method}}.
		}
	}
}


\details{

\describe{
	
	\item{consensus}{: 
    
    Computes the consensus matrix associated to the multiple NMF runs described by \code{object}.
    It is computed as the mean connectivity matrix of all the runs.
    
	It's been proposed by \emph{Brunet et al. (2004)} to help visualising and measuring the stability 
	of the clusters obtained by NMF approaches. 
	
	For objects of class \code{NMF} (e.g. results of a single NMF run, or NMF models), the consensus 
	matrix reduces to the connectivity matrix.
	
	Note that argument \code{...} is not used.	
	}
    
	\item{cophcor}{
	Computes the cophenetic correlation coefficient of consensus matrix \code{object}, 
	generally obtained from multiple NMF runs.
	
	The cophenetic correlation coeffificient is based on the consensus matrix 
	(i.e. the average of connectivity matrices) and was proposed by 
	\emph{Brunet et al. (2004)} to measure the stability of the clusters obtained 
	from NMF.
	
	It is defined as the Pearson correlation between the samples' distances 
	induced by the consensus matrix (seen as a similarity matrix) and their 
	cophenetic distances from a hierachical clustering based on these very 
	distances (by default an average linkage is used).
	See \emph{Brunet et al. (2004)}.
	
	Note that argument \code{...} is not used.
	}
	
	\item{dispersion}{
	Computes the dispersion coefficient of consensus matrix \code{object}, 
	generally obtained from multiple NMF runs.
	
	The dispersion coefficient is based on the consensus matrix (i.e. the average 
	of connectivity matrices) and was proposed by \emph{Kim and Park (2007)} to 
	measure the reproducibility of the clusters obtained from NMF
	.
	It is defined as:
	\deqn{\rho = \sum_{i,j=1}^n 4 (C_{ij} - \frac{1}{2})^2 .}
	, where \eqn{n} is the total number of samples.

 	We have \eqn{0 \leq \rho \leq 1} and \eqn{\rho = 1} only for a perfect 
 	consensus matrix, where all entries 0 or 1. 
 	A perfect consensus matrix is obtained only when all the connectivity 
 	matrices are the same, meaning that the algorithm gave the same 
 	clusters at each run.
 	See \emph{Kim and Park (2007)}
 	
 	Note that argument \code{...} is not used.
	}
	
	\item{fit}{: returns the element that achieves the lowest residual approximation error 
	across the runs.
	
	For \code{NMFfitX1} objects it coerces \code{object} into a \code{NMFfit} object.
	For \code{NMFfitXn} objects it builds and searches the vector of residuals of all the 
	fits and returns the one with the minimum value.
	
	Note that argument \code{...} is not used.
    }
    
	\item{metaHeatmap}{ Produces a heatmap of the consensus matrix using a 
	\code{heatmap}-like custom function, with parameters tuned for displaying 
	such result.
	
	The function used to draw the heatmap is a mixture of the function \code{heatmap.2} 
	from the \code{gplots} package, and the function \code{heatmap.plus}
	 from the \code{heatmap.plus} package. It allows to add extra annotation rows 
	 using the \code{ColSideColor} argument.
	 See \code{\link[gplots]{heatmap.2}} and \code{\link[heatmap.plus]{heatmap.plus}}.
	}

	\item{nrun}{ returns the number of NMF runs performed to compute \code{object}.
	
    In the case of a \code{NMFfitXn} object it returns its length -- as it is also a list.
    In the case of a \code{NMFfitX1} object it returns the value of its slot \code{nrun}.
    In the case of a \code{NMFfit} object it always returns 1 (this method exists to create 
    a uniform access interface to NMF results).
    
    Note that argument \code{...} is not used.
    }
    
	\item{predict}{ returns a \code{factor} that gives the predicted cluster index for 
	each sample (resp. for each feature) based on the \emph{best NMF factorization} stored in \code{object}.
	
    The index correspond to the basis vector that most contributes to the sample (resp. to which
     the feature contributes the most).
    See \code{\link{predict}} for more details.
	}
	
	\item{runtime.all}{: returns the computational time used to compute all the runs and create \code{object}. 	
	The time is computed using base function \code{\link{system.time}} which returns object of class 
    \code{\link[=proc.time]{proc_time}}.
    
    For \code{NMFfitXn} objects, there is also another time measure returned by the \code{seqtime} 
    method, which computes the sequential computational time, that is the sum of the computational 
    time used by each run. 
    
    Note that argument \code{...} is not used.
    }

	\item{seqtime}{: 
    returns the sequential CPU time spent of all the runs in the \code{object} -- which must be an 
    instance of class \code{\linkS4class{NMFfitXn}}.
    It is the sum of the CPU time used to compute each run. It returns \code{NULL} if the \code{object} 
    is empty.
    
    Note that argument \code{...} is not used.   
    }
    
	\item{summary}{: \code{summary} method for objects of class \code{NMFfitX}.
	 
    It computes a set of measures to help evaluate the quality of the \emph{best factorization} 
    of the set. The result is similar to the result from the \code{summary} method 
    of \code{NMFfit} objects. See \code{\linkS4class{NMF}} for details on the computed 
    measures. 
    In addition, the cophenetic correlation coefficient and the dispersion coefficient
    of the consensus matrix are returned, as well as the total computational time.
    See the related methods above.
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


