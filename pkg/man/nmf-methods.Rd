\name{nmf-methods}
\docType{methods}
\alias{nmf}
\alias{seed}
\alias{NMF Algorithms}
\alias{nmf-methods}
\alias{nmf,ANY,ANY,list-method}
\alias{nmf,ANY,ANY,missing-method}
\alias{nmf,data.frame,ANY,ANY-method}
\alias{nmf,ExpressionSet,ANY,ANY-method}
\alias{nmf,matrix,numeric,character-method}
\alias{nmf,matrix,numeric,function-method}
\alias{nmf,matrix,numeric,NMFStrategy-method}

\title{ Main Interface to run NMF algorithms}
\description{
 This method implements the main interface to launch NMF algorithms within the 
 framework defined in package \code{NMF}. It allows to combine NMF algorithms 
 with seeding methods. The returned object can be directly passed to 
 visualisation or comparison methods.
 
 For a tutorial on how to use the interface, please see the package's vignette: 
 \code{vignette('NMF')}  
}

\usage{

\S4method{nmf}{matrix,numeric,function}(x, rank, method, name, objective='euclidean', model='NMFstd', mixed=FALSE, ...)
\S4method{nmf}{matrix,numeric,character}(x, rank, method=nmf.getOption("default.algorithm"), ...)
\S4method{nmf}{matrix,numeric,NMFStrategy}(x, rank, method, seed=nmf.getOption('default.seed'), nrun=1, model=list(), .options=list(), ...)

}

\arguments{
  
  \item{method}{The algorithm to use to perform NMF on \code{x}.
  Different formats are allowed: \code{character}, \code{function}. 
  See section \emph{Methods} for more details on how each format is used.}
  
  \item{mixed}{Boolean that states if the algorithm requires a nonnegative input
   matrix (\code{mixed=FALSE} which is the default value) or accepts mixed sign 
   input matrices (\code{mixed=TRUE}). An error is thrown if the sign required 
   is not fulfilled. 
   This parameter is useful to plug-in algorithms such as semi-NMF, that 
   typically does not impose nonnegativity constraints on both the input and 
   the basis component matrices.
  }
   
  \item{model}{  When \code{method} is a \code{function}, argument \code{model} must be a single 
  	\code{character} string	or a \code{list} that specifies the NMF model to use. A NMF model is 
  	defined by a S4 class that extends class \code{\linkS4class{NMF}}. 
  	As a single \code{character} string, argument \code{model} must be the name of the class that 
  	defines the NMF model. 
  	As a \code{list} it contains at least the name of the class that defines the NMF model in element 
  	\code{model$class} or in the first element (if not named). The other element in the list are used 
  	to initialize the model's slots.
  	
  	When \code{method} is a single \code{character} string, argument \code{model} must be a \code{list}. 
  	It is used to initialize the slots in the NMF model associated with the NMF strategy of name 
  	\code{method}. Note that the model associated with NMF strategy cannot be changed via argument 
  	\code{model}.
  	
  	Note also that values to initialize the NMF model's slots can also be passed in \code{\dots}. 
  	However, those passed via argument \code{model} have priority over the later. 
  	This is designed to handle the situation where one wants to pass a parameter to the NMF algorithm,
  	that has the same name as a slot in the NMF model, or vice versa.
  	If a variable appears in both argument \code{model} and \code{\dots}, the former will be used to 
  	initialize the NMF model, the latter will be passed to the NMF algorithm. 
  	See code examples for an illustration of this situation. 
  	
  	}
    
  \item{name}{ A \code{character} string to be used as a name for the custom NMF algorithm.}
  
  \item{nrun}{ Used to perform multiple runs of the algorithm. It specifies 
  the number of runs to perform . 
  This argument is useful to achieve stability when using a random seeding method.}
  
  \item{objective}{ Used when \code{method} is a \code{function}. 
  It must be A \code{character} string giving the name of a built-in 
  distance method or a \code{function} to be used as the objective function. 
  It is used to compute the residuals between the target matrix and its NMF estimate.}
  
  \item{.options}{ 
  this argument is used to set some runtime options. It can be \code{list} 
  containing the named options and their values, or, in the case only boolean 
  options need to be set, a character string that specifies which options are 
  turned on or off. The string must be composed of characters that correspond to 
  a given option. Characters '+' and '-' are used to explicitly specify on and 
  off respectively. E.g. \code{.options='tv'} will toogle on options \code{track} 
  and \code{verbose}, while \code{.options='t-v'} will toogle on option \code{track} 
  and off option \code{verbose}. Note that '+' and '-' apply to all option 
  character found after them. Default behaviour assumes that \code{.options} 
  starts with a '+'.
    
  The following options are available (note the characters that correspond to 
  each option, to be used when \code{.options} is passed as a string):
  \describe{
  
  \item{debug - d}{ Toogle debug mode. Like option \code{verbose} but with more 
  information displayed.}
  
  \item{keep.all - k}{
  used when performing multiple runs (\code{nrun}>1): if toogled on, all 
  factorizations are saved and returned, otherwise only the factorization 
  achieving the minimum residuals is returned.
  }
  
%  \item{optimized - o}{
%  this option is used on multicore *nix machine only, when performing multiple runs 
%  (\code{nrun} > 1). If toggled on, the runs are performed using multiple cores 
%  thanks to package \code{link[package:doMC]{doMC}} -- which therefore needs to be 
%  installed.
%  }
  
  \item{track - t}{
  enables (resp. disables) error tracking. When \code{TRUE}, 
  the returned object's slot \code{residuals} contains the trajectory of the objective values. 
  This tracking functionality is available for all built-in algorithms.
  }
  
  \item{verbose - v}{ Toogle verbosity. If on, messages about the configuration and the 
  state of the current run(s) are displayed.}
  
  }% end describe .options
  
  }% end .options
  
  \item{rank}{The factorization rank to achieve [i.e a single positive \code{numeric}]}
   
  \item{seed}{The seeding method to use to compute the starting point passed 
  to the algorithm. 
  See section \emph{Seeding methods} %and \code{function \link{seed}} 
  for more details on the possible classes and types for argument \code{seed}.}
  
  \item{x}{The target object to estimate. It can be a \code{matrix}, 
  a \code{data.frame}, an \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}} object 
  (this requires the \code{Biobase} package to be installed). 
  See section \emph{Methods} for more details.}
  
  \item{\dots}{Extra parameters passed to the NMF algorithm's \code{run} method.  
  When there is no conflict with the slot names in the NMF model class, values for model 
  slots can also be passed in \code{...}. See argument \code{model}.
  }
}

\section{Methods}{
\describe{

\item{x = "ANY", rank = "ANY", method = "list"}{ Performs NMF on object \code{x}
 for each algorithm defined in \code{method}.}

\item{x = "ANY", rank = "ANY", method = "missing"}{ Performs default NMF 
algorithm on object \code{x}.}

\item{x = "data.frame", rank = "ANY", method = "ANY"}{ Performs NMF on a 
\code{data.frame}: the target matrix is the converted data.frame 
\code{as.matrix(x)} }

\item{x = "ExpressionSet", rank = "ANY", method = "ANY"}{ Performs NMF on an 
\code{ExpressionSet}: the target matrix is the expression matrix \code{exprs(x)}.
  
 This method requires the \code{Biobase} package to be installed.
 Special methods for bioinformatics are provided in an optional layer, 
 which is automatically loaded when the \code{Biobase} is installed.
 See \code{\link{NMF-bioc}}.
}

\item{x = "matrix", rank = "numeric", method = "character"}{ Performs NMF on a 
\code{matrix} using an algorithm whose name is given by parameter \code{method}.
 The name provided must partially match the name of a registered algorithm.
 See section \emph{Algorithms} below or the package's vignette.} 

\item{x = "matrix", rank = "numeric", method = "function"}{ Performs NMF using 
a custom algorithm defined by a \code{function}. It must have signature 
\code{(x=matrix, start=NMF, ...)} and return an object that inherits 
from class \code{NMF}. It should use its argument \code{start} as a starting 
point.}

%\item{x = "matrix", rank = "numeric", method = "NMFStrategy"}{ Performs NMF using an algorithm defined as a 
%\code{NMFStrategy} object. This version of method \code{nmf} is the one that is eventually called by all the 
%above versions.}
}

}

\value{
The returned value depends on the run mode:

  \item{Single run:}{An object that inherits from class \code{\linkS4class{NMF}}.}
  \item{Multiple runs:}{When \code{nrun > 1} or when \code{method} is a \code{list}, 
  this method returns an object of class \code{\linkS4class{NMFSet}}}
}

\section{NMF Algorithms}{

The following algorithms are available:

\describe{

\item{\code{lee}}{ Standard NMF. Based on euclidean distance, it uses simple 
multiplicative updates. See \emph{Lee and Seung (2000)}. }

\item{\code{brunet}}{ Standard NMF. Based on Kullbach-Leibler divergence, it uses simple 
multiplicative updates from \emph{Lee and Seung (2000)}, enhanced differently 
to avoid numerical underflow.
 See \emph{Brunet et al. (2004)}. }

\item{\code{lnmf}}{ Local Nonnegative Matrix Factorization. Based on a 
regularized Kullbach-Leibler divergence, it uses a modified version of 
Lee and Seung's multiplicative updates. 
See \emph{Li et al. (2001)}. }

\item{\code{nsNMF}}{ Nonsmooth NMF. Use a modified version of Lee and Seung's 
multiplicative updates for Kullbach-Leibler divergence to fit a extension 
of the standard NMF model. 
See \emph{Pascual-Montatno et al. (2006)}. }

\item{\code{offset}}{ Use a modified version of Lee and Seung's multiplicative 
updates for euclidean distance, to fit a NMF model that includes an intercept. 
See \emph{Badea (2008)}. }

\item{\code{snmf/r}, \code{snmf/l}}{ Alternating Least Square (ALS) approach 
from \emph{Kim and Park (2007)}. }

}

}

\section{Seeding methods}{
	The purpose of seeding methods is to compute initial values for the factor 
	matrices in a given NMF model. 
	This initial guess will be used as a starting point by the chosen NMF algorithm.
	
	The seeding method to use in combination with the algorithm can be passed 
	to interface \code{nmf} through argument \code{seed}.
	The following formats are supported:
	
	\describe{
	\item{a \code{character} string:}{ giving the name of a \emph{registered} 
	seeding method. The corresponding method will be called to compute 
	the starting point. }
	\item{a \code{list}:}{ giving the name of a \emph{registered} 
	seeding method and, optionally, extra parameters to pass to it.}
	\item{an object that inherits from \code{\linkS4class{NMF}}:}{ it will be 
	directly passed to the algorithm's method -- via its argument \code{start}.}
	\item{a \code{function}:}{ that computes the starting point. It must have 
	signature \code{(object=NMF, target=matrix, ...)} and return an object that 
	inherits from class \code{NMF}. Argument \code{object} should be used as a 
	template for the returned object.}
	}   
}

\author{Renaud Gaujoux \email{renaud@cbio.uct.ac.za}}

\seealso{class \code{\linkS4class{NMF}}, \code{\link{NMF-utils}}, package's vignette } 

\references{

Lee, D.~D. and Seung, H.~S. (2000).
Algorithms for non-negative matrix factorization. 
In \emph{NIPS}, 556--562.

Brunet, J.~P., Tamayo, P., Golub, T.~R., and Mesirov, J.~P. (2004).
Metagenes and molecular pattern discovery using matrix factorization. 
\emph{Proc Natl Acad Sci U S A}, \bold{101}(12), 4164--4169.

Pascual-Montano, A., Carazo, J.~M., Kochi, K., Lehmann, D., and Pascual-Marqui, R.~D. (2006).
Nonsmooth nonnegative matrix factorization (nsnmf).
\emph{IEEE transactions on pattern analysis and machine intelligence}, \bold{8}(3), 403--415.

Kim, H. and Park, H. (2007).
Sparse non-negative matrix factorizations via alternating non-negativity-constrained least squares for microarray data analysis.
\emph{Bioinformatics}.

Liviu Badea (2008).
Extracting Gene Expression Profiles Common To Colon And Pancreatic Adenocaricinoma Using Simultaneous Nonnegative Matrix Factorization.
In \emph{Pacific Symposium on Biocomputing}, \bold{13}, 279--290


S. Li, X. Hou, and H. Zhang (2001).
Learning spatially localized, parts-based representation. 
\emph{In Proc. CVPR}, 2001. 

}

\examples{

# generate a synthetic dataset with known classes: 100 features, 23 samples (10+5+8)
n <- 100; counts <- c(10, 5, 8);
V <- syntheticNMF(n, counts, noise=TRUE)

# build the class factor
groups <- as.factor(do.call('c', lapply(seq(3), function(x) rep(x, counts[x]))))

# run default algorithm
res <- nmf(V, 3)
res
summary(res, class=groups)

# run default algorithm multiple times (only keep the best fit)
res <- nmf(V, 3, nrun=20)
res
summary(res, class=groups)

# run nonsmooth NMF algorithm
res <- nmf(V, 3, 'nsNMF')
res
summary(res, class=groups)

# compare some NMF algorithms
res <- nmf(V, 3, list('brunet', 'lee', 'nsNMF'))
res
compare(res, class=groups)

# run on an ExpressionSet (requires package Biobase)
\dontrun{
data(esGolub)
nmf(esGolub, 3)
}


}

\keyword{methods}
\keyword{cluster}
\keyword{math}
\keyword{optimize}
