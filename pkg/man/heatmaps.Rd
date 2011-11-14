\name{Plots: heatmaps}
\docType{methods}

\alias{heatmap-NMF}
\alias{basismap}
\alias{basismap-methods}
\alias{basismap,NMF-method}
\alias{basismap,NMFfitX-method}

\alias{coefmap}
\alias{coefmap-methods}
\alias{coefmap,NMF-method}
\alias{coefmap,NMFfitX-method}

\alias{consensusmap}
\alias{consensusmap-methods}
\alias{consensusmap,NMF-method}
\alias{consensusmap,NMFfitX-method}
\alias{consensusmap,matrix-method}
\alias{consensusmap,NMF.rank-method}
\alias{consensusmap,list-method}

\title{Plots a Heatmap of the Basis Vector and Mixture Coefficient Matrices from an NMF model}

\description{
Produces a heatmap of the basis vectors or mixture coefficients using a 
\code{\link{aheatmap}}, with parameters tuned for displaying NMF results.

}

\usage{

\S4method{basismap}{NMF}(object, subsetRow=FALSE, tracks = 'basis', info=FALSE, ...)

\S4method{coefmap}{NMF}(object, tracks = 'basis', info=FALSE, ...)

\S4method{coefmap}{NMFfitX}(object, tracks = c('basis', 'consensus'), ...)

\S4method{consensusmap}{NMFfitX}(object, tracks = c('basis', 'consensus'), info = FALSE, ...)

}

\arguments{

	\item{object}{ An object that inherits from class \code{\linkS4class{NMF}} 
	(which includes class \code{\linkS4class{NMFfit}}) or class \code{\linkS4class{NMFfitX}}.	
	}
	
	\item{tracks}{ annotation track to add to the plot: 
		\itemize{
		\item{'basis': annotates the dominant basis component for each sample.}
		\item{'consensus': annotates the consensus clusters obtained by hierarchical 
		clustering using the consensus matrix as a similarity matrix and average 
		linkage.}
		}
	}
	
	\item{info}{ if \code{TRUE} then information about the NMF fit is added to the plot.
	If not \code{FALSE}, the this argument is passed to the heatmap plotting function 
	\code{\link{aheatmap}}. 
	}
	
	\item{subsetRow}{ 
	Argument that specifies how to filter the features that will appear in the heatmap.	
	When \code{FALSE} (default), all the features are used.
	Besides the values supported by argument \code{subsetRow} of 
	\code{\link{aheatmap}}, other possible values are:
	
	\itemize{
	\item \code{TRUE}: only the features that are basis-specific are used. 
	The default selection method is from Kim and Park (2007). Other selection methods can be specified 
	as a single \code{character} string or a custom \code{function} (cf. argument \code{method} 
	for function \code{\link{extractFeatures}}).
	
	\item a single \code{character} string that specifies the filtering method to be used to 
	select the basis-specific features that should appear in the heatmap (cf. argument \code{method} 
	for function \code{\link{extractFeatures}}).
	
	\item a single \code{numeric} value: in this case the features are sorted in decreasing order  
	according to their contribution to each basis, and only the top \code{filter} of each basis are kept.
	Multiple occurrences are removed. 
	}

	}
	
	\item{...}{ Used to pass extra parameters to the subsequent call to the internal drawing function 
	\code{heatmap.2.plus}.}
}

\details{	
Both \code{basismap} and \code{coefmap} are S4 generic functions.
}

\section{Methods}{
\describe{

\item{basismap}{\code{signature(object = "NMF")}: plots a heatmap of the basis vector matrix.
The features used in the heatmap can be filtered in various different way using the argument \code{filter}.}

\item{coefmap}{\code{signature(object = "NMF")}: plots a heatmap of the mixture coefficient matrix.}

\item{basismap, coefmap}{\code{signature(object = "NMFfitX")}:
These methods are used when \code{object} is the result of a multirun-NMF fit, as returned by the function \code{\link{nmf}} 
 when called with argument \code{nrun>1}.
It simply calls the corresponding method on the best NMF fit of 
\code{object}: \code{basismap(object, ...)} is roughly equivalent to 
\code{basismap(fit(object), ...)} and \code{coefmap(object, ...)} 
to \code{coefmap(fit(object), ...)}.
}
}

}

\seealso{ \code{\link{basis}}, \code{\link[=coef]{coef,NMF-method}}, \code{\link{consensusmap}} }

\author{ Renaud Gaujoux  }

\examples{

# generate a random target matrix
V <- rmatrix(100, 20)
res <- nmf(V, 3)
# Generates heatmaps
\dontrun{
basismap(res)
coefmap(res)
}

# filtered in different ways
\dontrun{
basismap(res, filter=TRUE)
basismap(res, filter=5)
basismap(res, filter='max')
basismap(res, filter=c(1,2,3, 10, 15, 20))
}
# the same not re-ordering the features
\dontrun{ basismap(res, filter=c(1,2,3, 10, 15, 20), Rowv=FALSE)}

## multirun NMF
res <- nmf(V, 3, nrun=3)
\dontrun{consensusmap(res)}
# compare the heatmap when there is a true signal
Vs <- syntheticNMF(100, 3, 30, noise=TRUE)
res <- nmf(Vs, 3, nrun=10)
\dontrun{consensusmap(res)} 


}

\keyword{methods}
\keyword{plot}

