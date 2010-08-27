\name{NMF-class}
\Rdversion{1.1}
\docType{class}
\alias{NMF-class}

\alias{consensus,NMF-method}

\alias{dim,NMF-method}
\alias{distance}
\alias{distance,matrix,NMF-method}
\alias{distance,missing,missing-method}

\alias{featureNames}
\alias{featureNames,NMF-method}
\alias{featureNames<-}
\alias{featureNames<-,NMF-method}

\alias{fitted,NMF-method}
\alias{is.empty.nmf}
\alias{is.empty.nmf,NMF-method}

\alias{modelname}
\alias{modelname,NMF-method}

\alias{nbasis}
\alias{nbasis,NMF-method}
\alias{rnmf}
\alias{rnmf-methods}
\alias{rnmf,NMF,numeric-method}
\alias{rnmf,NMF,matrix-method}
\alias{rnmf,NMF,missing-method}

\alias{sampleNames}
\alias{sampleNames,NMF-method}
\alias{sampleNames<-}
\alias{sampleNames<-,NMF-method}
\alias{sampleNames<-,NMF,ANY-method}

\alias{show,NMF-method}
\alias{summary,NMF-method}

\title{Interface Class for Nonnegative Matrix Factorisation Models}
\description{This is a \emph{virtual class} that defines a common interface to handle 
Nonnegative Matrix Factorisation models (NMF models) in a generic way.

It provides the definition for a minimum set of generic methods that are used in 
common computations and tasks in the context of Nonnegative Matrix Factorisations.

Class \code{NMF} makes it easy to develop new models that integrates well into 
the general framework implemented by the \emph{NMF} package.

Following a few simple guidelines, new models benefit from all the functionalities 
available to built-in NMF models -- that derive themselves from class \code{NMF}. 
See section \emph{Defining new NMF models} below. 

See section \code{\linkS4class{NMFstd}}, references and links therein for details on 
the standard NMF model and its -- built-in -- extensions.

}
 
\section{Slots}{
	This class contains a single slot, that is used internally during the computations.
	
	\describe{
	 \item{\code{misc}:}{A list that is used internally to temporarily store 
	 algorithm parameters during the computation.}
    }
	
	The purpose of this class is to define a common interface for NMF models as 
	a collection of generic methods. 
	Classes that inherits from class \code{NMF} are responsible for the management of 
	data storage and the implementation of the interface's pure virtual methods.
	
}

\section{Defining new NMF models}{
	The minimum requirement to define a new NMF model that integrates into the framework 
	of the \emph{NMF} package are the followings:
	
	\itemize{
	
	\item Define a class that inherits from class \code{NMF} and implements the new model. 
	Say class \code{myNMF}.
	
	\item Implement the following S4 methods for the new class \code{myNMF}:
		\describe{
		\item{fitted}{\code{signature(object = "myNMF", value = "matrix")}: 
		Must return the estimation of the target matrix as fitted by the NMF model 
		\code{object}.
		}
		\item{basis}{\code{signature(object = "myNMF")}: 
		Must return the matrix of basis vectors (e.g. the first matrix factor in 
		the standard NMF model).
		}
		\item{basis<-}{\code{signature(object = "myNMF", value = "matrix")}: 
		Must return \code{object} with the matrix of basis vectors set to
		\code{value}.
		} 
		\item{coef}{\code{signature(object = "myNMF")}: 
		Must return the matrix of mixture coefficients (e.g. the second matrix 
		factor in the standard NMF model).
		}
		\item{coef<-}{\code{signature(object = "myNMF", value = "matrix")}:
		Must return \code{object} with the matrix of mixture coefficients set to 
		\code{value}.
		} 
		}
		The \emph{NMF} package ensures these methods are defined for classes that 
		inherits from class \code{NMF}, as the methods defined for signatures 
		\code{(object='NMF', ...)} and \code{(object='NMF', value='matrix')} 
		throw an error when called.
		
	\item Optionally, implement method \code{rnmf}(signature(object="myNMF", target="numeric")). 
	This method should fill model \code{object} (of class \code{myNMF}) with 
	random values to fit a target matrix, whose dimension is given 
	by the 2-length numeric vector \code{target}.
	
	}
	
	For concrete examples of NMF models implementations, see class 
	\code{\linkS4class{NMFstd}} and its extensions (e.g. classes 
	\code{\linkS4class{NMFOffset}} or \code{\linkS4class{NMFns}}).
	
}
 
\section{Objects from the Class}{

Strictly speaking, because class \code{NMF} is virtual, no object of class \code{NMF} can be instantiated, 
only objects from its sub-classes. However, those objects are sometimes shortly referred 
in the documentation as "\code{NMF} objects" instead of "objects that inherits from class \code{NMF}".

For built-in models or for models that inherit from the standard model class \code{NMFstd}, 
the factory method \code{nmfModel} enables to easily create valid \code{NMF} 
objects in a variety of common situations. 
See \code{\link{nmfModel}} for more details.
 
}

\section{Methods}{
  \describe{

   \item{[}{\code{signature(x = "NMF", i, j, ..., drop=FALSE)}: 
    sub-setting method for object of class \code{NMF}.
    Row subsets are applied to the basis matrix rows, while column subsets are 
    applied to the mixture coefficient matrix. A third index can be provided 
    to subset on the basis components (e.g \code{x[,,1:3]}).
    
    See \code{NMF?subset} or \code{method?'[,NMF'}.
    }

   \item{basis}{\code{signature(object = "NMF")}: 
    returns the matrix of basis vectors according to the model defined in 
    \code{object}. 
    This is a \emph{pure virtual} method that needs to be defined 
    for the sub-classes of class \code{NMF} that implements concrete models.
    See also \code{\link{basis}}.
    }

    \item{basis<-}{\code{signature(object = "NMF", value = "matrix")}: 
    sets the matrix of basis vectors in \code{object}, and returns the updated 
    object. 
    This is a \emph{pure virtual} method that needs to be defined 
    for the sub-classes of class \code{NMF} that implements concrete models.
    See also \code{\link{basis}}.
    }
    
    \item{basisnames}{\code{signature(object = "NMF")}: 
    returns the names of the columns of the matrix of basis vectors, that are also 
    always the names of the rows of the mixture coefficient matrix.
    See \code{\link{basisnames}}.
    }
    
    \item{basisnames<-}{\code{signature(object = "NMF", value='ANY')}: 
    simultaneously sets the names of the columns of the matrix of basis vectors AND 
    the rows of the mixture coefficient matrix.
    See \code{\link{basisnames<-}}.
    }

    \item{coef}{\code{signature(object = "NMF")}:
    returns the matrix of mixture coefficients according to the model defined 
    in \code{object}. 
    This is a \emph{pure virtual} method that needs to be defined 
    for the sub-classes of class \code{NMF} that implements concrete models.
    See also \code{\link{coef}}.
    }

    \item{coef<-}{\code{signature(object = "NMF", value = "matrix")}: 
    sets the matrix of mixture coefficients in \code{object}, and returns the 
    updated object. 
    This is a \emph{pure virtual} method that needs to be defined 
    for the sub-classes of class \code{NMF} that implements concrete models.
    See also \code{\link{coef<-}}.
    }

    \item{coefficients}{\code{signature(object = "NMF")}:
    This is a simple alias to method \code{coef}.
    See also \code{\link{coefficients}}.
    }

    \item{coefficients<-}{\code{signature(object = "NMF", value = "matrix")}:
    This is a simple alias to method \code{coef<-}.
    See also \code{\link{coef}}. 
    }

    \item{connectivity}{\code{signature(x = "NMF")}: 
    returns the connectivity matrix associated to the clusters based on NMF 
    factorization \code{x}.
    The connectivity matrix \eqn{C} of a clustering is the symmetric matrix that
     shows the shared membership of the samples: entry \eqn{C_{ij}} is 1 if 
     samples \eqn{i} and \eqn{j} belong to the same cluster, 0 otherwise.
    }
    
    \item{consensus}{\code{signature(object = "NMF")}: 
    returns the consensus matrix associated with multiple runs of NMF, as the  
    the mean connectivity matrix across the runs.
    For objects of class \code{NMF}, it reduces to the connectivity matrix, and 
    the method is defined to create a uniform access interface to NMF results.
    }

    \item{dim}{\code{signature(x = "NMF")}: returns a 3-length vector containing 
    the dimension of the target matrix together with the NMF factorization rank. 
    For example \code{c(2000,30,3)} for a \code{NMF} object that fits a 2000x30 
    target matrix using 3 basis vectors.}
    
    \item{dimnames}{\code{signature(x = "NMF")}: returns the dimension names of the 
    NMF model \code{x}.
    See \code{\link{dimnames-NMF}}. 
    }
    
    \item{dimnames<-}{\code{signature(x = "NMF", value)}: sets the dimension names of the 
    NMF model \code{x}.
    See \code{\link{dimnames-NMF}}.   
    }

    \item{distance}{\code{signature(target = "matrix", x = "NMF", method, ...)}: 
    returns the distance between a target matrix and a NMF model, according to 
    a given measure. 
    If both argument \code{target} and \code{x} are missing, this function returns
    the \code{function} defined by argument \code{method}. The later can either be 
    a \code{function} or a \code{character} string that correspond to a registered 
    distance metric. For the moment only the metric 'KL' and 'euclidean' are defined.
    
    %See function \code{\link{distance}} for more details.
    }

    \item{entropy}{\code{signature(x = "NMF", class = "factor")}: 
    computes the entropy of NMF model \code{x} given a priori known groups of samples. 
    See generic function \code{\link{entropy}} for more details.}

    \item{evar}{\code{signature(object = "NMF", target)}: 
    computes the explained variance of NMF model \code{object} that approximates 
    \code{target}. 
    See generic function \code{\link{evar}} for more details.}
    
    \item{featureNames}{\code{signature(object = "NMF")}: 
    returns the row names of the basis matrix.
    If BioConductor is installed this method is defined for the generic function 
    \code{\link[Biobase]{featureNames}} from the \code{Biobase} package.
    }

	\item{featureNames<-}{\code{signature(object = "NMF", value = "ANY")}: 
    sets the row names of the basis matrix. Argument \code{value} must be in a 
    format accepted by the \code{\link{rownames}} method defined for matrices.
    If BioConductor is installed this method is defined for the generic function 
    \code{\link[Biobase]{featureNames<-}} from the \code{Biobase} package.
    }
    
    \item{fitted}{\code{signature(object = "NMF")}: 
    computes the target matrix estimated by NMF model \code{object}.
    This is a \emph{pure virtual} method that needs to be implemented 
    by the sub-classes of class \code{NMF} that implements concrete models.
    }

    \item{featureScore}{\code{signature(x = "NMF")}: Computes a score for each feature 
    that reflects its specificity to one of the basis vector. 
    The definition of the score follows \emph{Kim and Park (2007)}.
    See references for more details. }

    \item{extractFeatures}{\code{signature(x = "NMF")}: extract the features 
    that are the most specific to each basis vector.
    It follows \emph{Kim and Park (2007)}'s methodology.
    See references for more details. }

    \item{is.empty.nmf}{\code{signature(object = "NMF")}: 
    Tells if \code{object} is an empty the \code{NMF} model, that is it contains 
    no data. 
    It returns \code{TRUE} if the matrices of basis vectors and mixture coefficients
    have respectively zero rows and zero columns. It returns \code{FALSE} otherwise.
    This means that an empty model can still have a non-zero number of basis vectors. 
    For example, this happens in the case of NMF models created using factory method 
    \code{\link{nmfModel}} with no initialisation for any factor matrices.
    }
    
    \item{metaHeatmap}{ Produces a heatmap of the basis or mixture 
    matrix using function \code{\link[gplots]{heatmap.2}} with parameters 
    tuned for displaying NMF results. 
    See \code{\link{metaHeatmap}} for more details.
	}
	
	\item{\code{modelname}:}{ returns the name of the model fitted by the object.
	It corresponds to the name of the S4 class of the object.}    
	
    \item{nbasis}{\code{signature(x = "NMF")}: 
    returns the number of basis vectors used in NMF model \code{x}. 
    It is the number of columns of the matrix of basis vectors.}

    \item{predict}{\code{signature(object = "NMF")}: 
    returns a \code{factor} that gives the predicted cluster index for each sample 
    (resp. for each feature) based on NMF factorization \code{object}.
    The index correspond to the basis vector that most contribute to the sample 
    (resp. to which the feature contribute the most).
    See \code{\link{predict}} details on extra arguments.
    } 

    \item{purity}{\code{signature(x = "NMF", class = "factor")}: 
    computes the purity of NMF model \code{x} given a priori known groups of samples. 
    The purity definition can be found in \emph{Kim and Park (2007)}.
    See references for more details.}

    \item{rnmf}{\code{signature(x = "NMF", target)}: 
    seeds NMF model \code{x} with random values drawn from a uniform distribution. 
    The result is a NMF model of the same class as \code{x} with basis and 
    mixture matrices filled with random values. 
    
    Argument \code{target} can be either:
    
	    \describe{
		    \item{\code{numeric}}{ it must be of length 2 (resp. 1), and give 
		    the dimension of the target matrix (resp. symmetric matrix) to fit. The result is 
		    a random NMF model.
		    } 
		    
		    \item{\code{missing}}{ it returns \code{rnmf(x, c(nrow(x), ncol(x)))}, 
		    that is a random NMF model with the same dimensions as defined in model \code{x}.
		    }
		    
		    \item{\code{matrix}}{ it returns \code{rnmf(x, dim(target))}, 
		    that is a random NMF model that fits a matrix of the same dimension as \code{target}.
		    The values are drawn within the interval \code{[0, max(target)]}.
		    }		    
	    }
	    
	    This method's version with signature \code{(object='NMF', target='numeric')} might need 
	    to be overloaded if the initialisation of the specific NMF model requires 
	    setting values for data other than the basis and mixture matrices.
	    The overloading methods must call the generic version using function 
	    \code{\link{callNextMethod}}.
    }

    \item{rss}{\code{signature(object = "NMF")}: 
    returns the Residual Sum of Squares (RSS) between the target matrix and its 
    estimation by the NMF model \code{object}. \emph{Hutchins et al. (2008)} used 
    the variation of the RSS in combination with \emph{Lee and Seung}'s algorithm 
    to estimate the correct number of basis vectors.
    See \code{\link{rss}} for details on its usage.
    
    }

	\item{sampleNames}{\code{signature(object = "NMF")}: 
    returns the column names of the mixture coefficient matrix.
    If BioConductor is installed this method is defined for the generic function 
    \code{\link[Biobase]{sampleNames}} from the \code{Biobase} package.
    }

	\item{sampleNames<-}{\code{signature(object = "NMF", value = "ANY")}: 
    sets the columns names of the basis matrix. Argument \code{value} must be in a 
    format accepted by the \code{\link{colnames}} method defined for matrices.
    If BioConductor is installed this method is defined for the generic function 
    \code{\link[Biobase]{sampleNames<-}} from the \code{Biobase} package.
    }
    
    \item{show}{\code{signature(object = "NMF")}: 
    standard generic \code{show} method for objects of class \code{NMF}. 
    It displays the model class (i.e. the name of the sub-class that implements 
    the concrete model), the dimension of the target matrix, and the number of 
    basis vectors.
    }

    \item{sparseness}{\code{signature(x = "NMF")}: 
    compute the average sparseness of the basis vectors and mixture coefficients. 
    See \emph{Hoyer (2004)} for more details. }

    \item{summary}{\code{signature(x = "NMF")}: 
    standard generic \code{summary} method for objects of class \code{NMF}. 
    It computes a set of measures to evaluate the quality of the factorization.}
	 }
}
\references{ 

	Definition of Nonnegative Matrix Factorization in its modern formulation:
	
	Lee D.D. and Seung H.S. (1999).
	Learning the parts of objects by non-negative matrix factorization. 
	\emph{Nature}, \bold{401}, 788--791.

	Historical first definition and algorithms: 
		
	Paatero, P., Tapper, U. (1994).
	Positive matrix factorization: A non-negative factor model with optimal utilization of error estimates of data values.
	\emph{Environmetrics}, \bold{2}, 111--126
	, doi:10.1002/env.3170050203.
	
	
	Reference for some utility functions:
	
	Kim, H. and Park, H. (2007).
	Sparse non-negative matrix factorizations via alternating non-negativity-constrained least squares for microarray data analysis.
	\emph{Bioinformatics}.
	
	Hoyer (2004).
	Non-negative matrix factorization with sparseness constraints.
	\emph{Journal of Machine Learning Research}, \bold{5}, 1457-1469.
	
	
}
\author{ Renaud Gaujoux \email{renaud@cbio.uct.ac.za} }

\seealso{
	Main interface to perform NMF in \code{\link{nmf-methods}}. 
	
	Built-in NMF models and factory method in \code{\link{nmfModel}}.
	
	Method \code{\link{seed}} to set NMF objects with values suitable to start 
	algorithms with. 
}
\examples{

# show all the NMF models available (i.e. the classes that inherit from class NMF)
nmfModel()
# show all the built-in NMF models available
nmfModel(builtin.only=TRUE)

# class NMF is a virtual class so cannot be instantiated: 
# the following generates an error
\dontrun{new('NMF')}

# To instantiate a NMF model, use factory method nmfModel. see ?nmfModel
nmfModel(3)
nmfModel(3, model='NMFns')


}
\keyword{classes}
