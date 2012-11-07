\name{Methods from stats: predict, deviance, residuals}
\docType{methods}
\alias{stats-NMF}

\alias{deviance,NMFfit-method}
\alias{deviance,NMFfitXn-method}

%\alias{predict.nmf}
\alias{predict}
\alias{predict-methods}
\alias{predict,matrix-method}
\alias{predict,NMF-method}
\alias{predict,NMFfitX-method}

\alias{residuals}
\alias{residuals,NMFfit-method}

\title{Methods for the Interface Defined in Package stats}

\description{
The package NMF defines methods \code{deviance}, \code{predict} and \code{residuals} 
for common NMF objects returned by the function \code{\link{nmf}}.
}

\usage{

\S4method{deviance}{NMFfit}(object)

\S4method{predict}{NMF}(object, what = c("columns", "rows", "samples", "features"), prob=FALSE)

\S4method{predict}{NMFfitX}(object, what = c("columns", "rows", 'samples', 'features', 'consensus', 'cmap'), ...)

\S4method{residuals}{NMFfit}(object, track=FALSE, ...)

}

\arguments{


	\item{object}{ A \code{matrix} or an object that inherits from class 
	\code{\linkS4class{NMF}} or \code{\linkS4class{NMFfit}} -- depending on the method. }
	
	\item{prob}{ Should the probability associated with each cluster prediction be 
	computed and returned.
	}
	
	\item{track}{ if \code{TRUE}, the whole residuals track is returned. 
	Otherwise only the last residuals computed is returned.}
	
	\item{what}{ Specifies on which matrix the computation should be based:
	\itemize{
	\item \code{'columns', 'samples'}: mixture coefficients
	\item \code{'rows', 'features'}: basis components
	\item \code{'consensus', 'cmap'}: consensus matrix
	}
	See section \emph{Details} for more information.
	}
	
	\item{...}{argument passed to \code{predict,NMF} if \code{what} is not 
	\code{'consensus'}.}
}

\details{

\describe{
	
	\item{deviance}{ returns the final approximation error between an NMF model 
 	and its target matrix, as the value of the minimized objective function.
 	In the case of multiple NMF runs, it returns the approximation error of the 
 	best fit across the runs.
 	
 	If not computed by the NMF algorithm itself, the value is automatically 
 	computed at the end of the fitting process by the function \code{\link{nmf}}, 
 	using the objective function associated with the NMF algorithm.}
 
	\item{predict}{ Computes the dominant basis component for each sample (resp. feature) 
	based on its associated entries in the mixture coefficient matrix (i.e in \eqn{H}) 
	(resp. basis component matrix (i.e in \eqn{W})).
	
	When \code{what='samples'} the computation is performed on the mixture 
	coefficient matrix, or on the transposed basis matrix when \code{what='features'}.
	
	For each column, the dominant basis component is computed as the row index 
	 for which the entry is the maximum within the column.
	 
	If argument \code{prob=FALSE} (default), the result is a \code{factor}. 
	Otherwise it returns a list with two elements: element \code{predict} contains 
	the computed indexes ( as a \code{factor}) and element \code{prob} contains 
	the vector of the associated probabilities, that is the relative contribution 
	of the maximum entry within each column.
	
	}

	\item{residuals}{ returns the -- final -- residuals between the target matrix and the 
	NMF result \code{object}. They are computed using the objective function 
	associated to the NMF algorithm that returned \code{object}.
	When called with \code{track=TRUE}, the whole residuals track is returned, 
	if available. Note that method \code{\link{nmf}} does not compute the 
	residuals track, unless explicitly required.
	
	It is a S4 methods defined for the associated generic functions from package 
	\code{stats} (See \link[stats]{residuals}).
	
	Note that stricly speaking, this method does not fulfill its contract as 
	defined by the package \code{stats}, but rather acts as function \code{deviance}.  
	The will be changed in a later release to make it consistent with
	its original purpose.
	}
}

}

\seealso{	
	\code{\linkS4class{NMF}, \link[=summary,NMF-method]{summary}}  
}

\examples{


# generate random data matrix
X <- rmatrix(50, 20)

# run NMF
res <- nmf(X, 2)
res

# get the predicted column clusters
predict(res)
# get the predicted row clusters
predict(res, 'rows')

# get the final approximation error
deviance(res)
# this is the same as the residuals and no track is computed by default
residuals(res)
residuals(res, track=TRUE)
# compute the track
res <- nmf(X, 2, .options='t')
residuals(res, track=TRUE)

# perform multiple NMF runs
res <- nmf(X, 2, nrun=3)

# get the predicted consensus column clusters
predict(res, 'consensus')

}

