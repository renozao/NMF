\name{NMFfit-class}
\Rdversion{1.1}
\docType{class}
\alias{NMFfit-class}
\alias{algorithm}
\alias{algorithm<-}
\alias{algorithm,NMFfit-method}
\alias{algorithm<-,NMFfit,ANY-method}
\alias{basis<-,NMFfit,matrix-method}
\alias{basis,NMFfit-method}
\alias{coef,NMFfit-method}
\alias{coef<-,NMFfit,matrix-method}
\alias{distance,matrix,NMFfit-method}

\alias{$}
\alias{$<-}
\alias{$,NMFfit-method}
\alias{$<-,NMFfit-method}

\alias{fit-methods}
\alias{fit}
\alias{fit,NMFfit-method}
\alias{fit<-}
\alias{fit<-,NMFfit,NMF-method}
\alias{fitted,NMFfit-method}

\alias{modelname,NMFfit-method}

\alias{nrun,NMFfit-method}

\alias{objective}
\alias{objective<-}
\alias{objective,NMFfit-method}
\alias{objective<-,NMFfit,character-method}
\alias{objective<-,NMFfit,function-method}

\alias{residuals<-}
\alias{residuals<--methods}
\alias{residuals<-,NMFfit-method}
%\alias{run.options,NMF-method}
%\alias{run.options<-,NMF-method}
\alias{runtime}
\alias{runtime,NMFfit-method}
\alias{runtime.all,NMFfit-method}
%\alias{seed,matrix,NMF,function-method}
\alias{seeding}
\alias{seeding-methods}
\alias{seeding,NMFfit-method}
\alias{seeding<-}
\alias{seeding<-,NMFfit-method}
\alias{show,NMFfit-method}
\alias{summary,NMFfit-method}
%\alias{verbose,NMFfit-method}

\title{Base Class for to store Nonnegative Matrix Factorisation results}
\description{Base class to handle the results of general \strong{Non-negative Matrix Factorisation} 
algorithms (NMF).

It provides a general structure and generic functions to manage the results of NMF algorithms. 
It contains a slot with the fitted NMF model (see slot \code{fit}) as well as data about the methods 
and parameters used to compute the factorization.
}
 
\section{Slots}{	
	 \describe{
	 \item{\code{fit}:}{An object that inherits from class \code{"NMF"}. 
	 It contains the fitted NMF model.
	 
	 Note that class \code{"NMF"} is a virtual class. 
	 The default class for this slot is \code{NMFstd}, that implements the 
	 standard NMF model.
    }
    
    \item{\code{residuals}:}{A \code{"numeric"} vector that contains the final residuals 
    or the residuals track between the target matrix and its NMF estimate(s).
    Default value is \code{numeric()}.
    
    See method \code{\link{residuals}} for details on accessor methods 
    and main interface \code{\link{nmf}} for details on how to compute NMF with 
    residuals tracking.}
    
    \item{\code{method}:}{A single \code{"character"} string that contains the name of the algorithm used to fit 
    the model. Default value is \code{''}.}
        
    \item{\code{seed}:}{A single \code{"character"} string that contains the name of the seeding method used to seed the 
    algorithm that computed the NMF. Default value is \code{''}. 
    See \code{\link{nmf-methods}} for more details. }
    
    \item{\code{distance}:}{Either a single \code{"character"} string that 
    contains the name of the built-in objective function, or a \code{function} 
    that measures the residuals between the target matrix and its NMF estimate. }
	
    \item{\code{parameters}:}{A \code{"list"} that contains the extra parameters specific to the algorithm used 
    to fit the model. }
    
    \item{\code{runtime}:}{Object of class \code{"proc_time"} that contains 
    various measures of the time spent to fit the model.
    See \code{\link[base]{system.time}}}
    
    \item{\code{options}:}{A \code{"list"} that contains the options used to compute the object.}
    
    \item{\code{extra}:}{A \code{"list"} that contains extra miscellaneous data 
    set by the algorithm. For example it can be used to store extra parameters, 
    without the need to extend the standard NMF model.}
  }
}

\details{
The purpose of this class is to handle in a generic way the results of NMF algorithms.
Its slot \code{fit} contains the fitted NMF model as an object of class \code{\linkS4class{NMF}}.

Other slots contains data about how the factorization has been computed, such as the algorithm and seeding 
method, the computation time, the final residuals, etc\dots

Class \code{NMFfit} acts as a wrapper class for its slot \code{fit}. 
It inherits from interface class \code{\linkS4class{NMF}} defined for generic NMF models. 
Therefore, all the methods defined by this interface can be called directly on objects 
of class \code{NMFfit}. The calls are simply dispatched on slot \code{fit}, i.e. 
the results are the same as if calling the methods directly on slot \code{fit}.
 }
  
 \section{Validity checks}{ 
 
 The validity method for class \code{NMFfit} checks
 
 \itemize{
 
 \item slot \code{fit} calling the suitable validity function on this object 
 of class \code{NMF} (see \code{\linkS4class{NMF}} for more details).
 
 \item the validity of slot \code{objective} that must be either a function 
 definition or a \emph{non-empty} character string.
 }
 
 }
 
\section{Objects from the Class}{

Object of class \code{NMFfit} using the standard operator \code{\link{new}}.

However, there is usually no need to directly create such an object, as interface methods 
\code{\link{nmf}} and \code{\link{seed}} take care of it.

}

\section{Methods}{

	Class-specific methods:
  \describe{
    \item{algorithm}{\code{signature(object = "NMFfit")}: Access slot \code{method}.}    
    
    \item{algorithm<-}{\code{signature(object = "NMFfit", value="ANY")}: Set slot \code{method}.}
    
    \item{basis}{\code{signature(object = "NMFfit")}:
	 Extract the matrix of basis vectors from the fitted NMF model. 
	 It returns \code{basis(fit(object))}.
	 
	Note that this method is part of the minimum interface for NMF models, as 
    defined by class \code{NMF}. See \code{\linkS4class{NMF}}.
	}
	
	\item{basis<-}{\code{signature(object = "NMFfit", value = "matrix"}:
	Sets the matrix of basis vectors of the fitted NMF model. 
	 It calls \code{basis(fit(object), value)}.
	 
	Note that this method is part of the minimum interface for NMF models, as 
    defined by class \code{NMF}. See \code{\linkS4class{NMF}}.
	}
    
	\item{coef}{\code{signature(object = "NMFfit")}:
	 Extract the matrix of mixture coefficients from the fitted NMF model. 
	 It returns \code{coef(fit(object))}.
	 
	Note that this method is part of the minimum interface for NMF models, as 
    defined by class \code{NMF}. See \code{\linkS4class{NMF}}.
	}
	
	\item{coef<-}{\code{signature(object = "NMFfit", value = "matrix"}:
	Sets the matrix of mixture coefficients of the fitted NMF model. 
	 It calls \code{coef(fit(object), value)}.
	 
	Note that this method is part of the minimum interface for NMF models, as 
    defined by class \code{NMF}. See \code{\linkS4class{NMF}}.
	}
	
    \item{distance}{\code{signature(target = "matrix", x = "NMFfit")}: return the 
    value of the loss function given a target matrix and a NMF fit. It calls method \code{distance} 
    on slot \code{fit}. If a distance method is NOT supplied in argument \code{method}, it uses 
    the objective function from slot \code{objective}.
    %See function \code{\link{distance}} for more details.
    }
    
    \item{plot}{\code{signature(x = "NMFfit")}: plot the residuals track of the run 
    that computed object \code{x}. 
    See function \code{\link{nmf}} for details on how to enable the tracking of residuals.}

    \item{fitted}{\code{signature(object = "NMFfit")}: 
    compute the estimated target matrix according to the NMF model stored in slot
    \code{fit}. It actually dispatches the call to slot \code{fit}, returning 
    \code{fitted(fit(object))}.
    
    Note that this method is part of the minimum interface for NMF models, as 
    defined by class \code{NMF}. See \code{\linkS4class{NMF}}.
    %See \code{\link{estimate}}.
    }
    
	\item{fit}{\code{signature(object = "NMFfit")}: return the \code{NMF} object 
	stored in slot \code{fit}.
	}
	\item{fit<-}{\code{signature(object = "NMFfit", value = "NMF")}: 
	set the value of the \code{NMF} object stored in slot \code{fit}.
	}
	
	\item{\code{modelname}:}{ returns the name of the model fitted by the object.
	It corresponds to the name of the S4 class of slot \code{fit}.}
    
	\item{nrun}{\code{signature(object = "NMFfit")}: return the number of NMF runs 
	performed to compute \code{object}. In the case of a \code{NMFfit} object, this method 
	always returns 1, as it is the result of a single NMF run and the method exists 
	to create a uniform access interface to NMF results. 
	}
	
	\item{objective}{\code{signature(object = "NMFfit")}: return slot \code{distance} 
	or compute the objective value if a target is passed in argument \code{x}.}
	\item{objective<-}{\code{signature(object = "NMFfit", value = "character")}: 
	sets slot \code{distance} to a built-in distance metric identified by a 
	\code{character} string. 
	}
	
	\item{objective<-}{\code{signature(object = "NMFfit", value = "function")}: 
	sets slot \code{distance} to a custom \code{function}. The function should 
	return a single positive \code{numeric} and must take the target 
	(a \code{matrix}) as its first parameter, and an object that inherits 
	from class \code{NMF} as its second parameter. Extra parameters are passed.
	}

    \item{residuals}{\code{signature(object = "NMFfit")}: Access slot \code{residuals}. 
    See function \code{\link{residuals}} for details 
    on extra parameters.
    }

    \item{residuals<-}{\code{signature(object = "NMFfit", value="numeric")}: 
    Set slot \code{residuals} to value \code{value}. 
    }

%    \item{run.options}{\code{signature(object = "NMF")}: return the runtime 
%    options used to fit NMF model \code{object}. See \code{\link{nmf.options}}. }

%    \item{run.options<-}{\code{signature(object = "NMF")}: set the runtime 
%    options to be used to fit NMF model \code{object}. See \code{\link{nmf.options}}. }

    \item{runtime}{\code{signature(object = "NMFfit")}: return the time spent to fit 
    NMF model \code{object}. The time is computed using base function 
    \code{\link{system.time}} which returns object of class 
    \code{\link[=proc.time]{proc_time}}. }
    
    \item{runtime.all}{\code{signature(object = "NMFfit")}: return the CPU time used 
     to perform all the runs to compute the object. 
     In the case of a \code{NMFfit} object, this method is an alias to method 
     \code{runtime} (see above), as the object is the result of a single NMF run.
     The method exists to create a uniform access interface to NMF results.
    }    

    \item{seeding}{\code{signature(object = "NMFfit")}: returns the seeding method 
    used to seed the algorithm that fitted NMF model \code{object}. 
    See section \code{\link{nmf-methods}}. 
    }

    \item{seeding<-}{\code{signature(object = "NMFfit")}: sets the seeding method 
    used to seed the algorithm that fitted NMF model \code{object}. }

    \item{show}{\code{signature(object = "NMFfit")}: standard generic \code{show} 
    method for objects of class \code{NMF}. }

    \item{summary}{\code{signature(x = "NMFfit")}: standard generic \code{summary} 
    method for objects of class \code{NMF}. It returns a \code{numeric} vector 
    that contains the summary of the fitted NMF model (slot \code{fit}), plus 
    the computation time and the final residuals.}
    
    \item{$}{\code{signature(x = "NMFfit", name)}: Access element \code{name} in 
    slot \code{extra} (which is a \code{list}) partially matching argument \code{name}.
    It is equivalent to \code{x@extra[[name, exact=FALSE]]}. 
    }
    
    \item{$<-}{\code{signature(x = "NMFfit", name, value)}: Set the value in of 
    element \code{name} in slot \code{extra} (which is a \code{list}).
     Note that it does not partially match argument \code{name}. 
    }
  
    
	}
	
	Class \code{NMFfit} inherits from all the methods defined on class \code{NMF} 
	to manipulate and interpret NMF models. For those methods, class \code{NMFfit} 
	act as a wrapper class, dispatching the calls to slot \code{fit}.
	Some useful methods are: \code{dim}, \code{nbasis}, \code{predict}, \code{sparseness}. 
	See \code{\linkS4class{NMF}} for more details.
}


\author{ Renaud Gaujoux \email{renaud@cbio.uct.ac.za} }

\seealso{
	Main interface to perform NMF in \code{\link{nmf-methods}}. 
			
	Method \code{\link{seed}} to set NMF objects with values suitable to start 
	algorithms with. 
}
\examples{

# run default NMF algorithm on a random matrix
n <- 50; r <- 3; p <- 20
V <- matrix(runif(n*p), n, p)  
res <- nmf(V, r)							

# result class is NMFfit
class(res)

# show result
res

# compute summary measures
summary(res)

}
\keyword{classes}
