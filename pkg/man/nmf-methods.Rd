\name{nmf-methods: Running NMF algorithms}
\docType{methods}
\alias{nmf}
\alias{seed}
\alias{NMF Algorithms}
\alias{nmf-methods}
\alias{nmf,data.frame,ANY,ANY-method}
\alias{nmf,ExpressionSet,ANY,ANY-method}
\alias{nmf,matrix,ANY,ANY-method}
\alias{nmf,matrix,numeric,list-method}
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

\S4method{nmf}{matrix,numeric,function}(x, rank, method, name, objective='euclidean', model='NMFstd' 
, mixed=FALSE, ...)

\S4method{nmf}{matrix,numeric,character}(x, rank, method, ...)

\S4method{nmf}{matrix,numeric,NMFStrategy}(x, rank, method, seed=nmf.getOption('default.seed')
, nrun=1, model=NULL, .options=list()
, .pbackend = nmf.getOption("parallel.backend")
, .callback = NULL
, ...)

}

\arguments{
  
  \item{method}{The algorithm to use to perform NMF on \code{x}.
  Different formats are allowed: \code{character}, \code{function}.
  If missing, the method to use is retrieved from the NMF package's specific 
  options by \code{nmf.getOption("default.algorithm")} (the default built-in 
  option is \code{'brunet'}).  
  See section \emph{Methods} for more details on how each format is used.}
  
  \item{mixed}{Boolean that states if the algorithm requires a nonnegative input
   matrix (\code{mixed=FALSE} which is the default value) or accepts mixed sign 
   input matrices (\code{mixed=TRUE}). An error is thrown if the sign required 
   is not fulfilled. 
   This parameter is useful to plug-in algorithms such as semi-NMF, that 
   typically does not impose nonnegativity constraints on both the input and 
   the basis component matrices.If \code{NULL} then the NMF model is used
  }
   
  \item{model}{  
  
  	When \code{method} is a \code{function}, argument \code{model} 
  	must be either a single \code{character} string (default to 'NMFstd') or 
  	a \code{list} that specifies values for slots in the NMF model. The NMF 
  	model to be instantiated can optionally be given by its class name in the 
  	first element of the list [note: A NMF model is defined by a S4 class that 
  	extends class \code{\linkS4class{NMF}}]. If no class name is specified then 
  	the default	model is used, see \code{\linkS4class{NMFstd}}. 
  	
  	When \code{method} is a single \code{character} string or a \code{NMFStrategy}
  	object, argument \code{model} must be \code{NULL} (default) or a \code{list}.
  	Note that in this case the NMF model is defined by the NMF strategy itself 
  	and cannot be changed.
  	
  	\itemize{
  	
  	\item If a single \code{character} string, argument \code{model} must be the name
  	 of the class that defines the NMF model to be instantiated. 
  	 Arguments in \code{...} are handled in the same way as when \code{model} 
  	 is \code{NULL}, see below.
  	 
  	\item If a \code{list} \code{all and only} its elements are used to 
  	initialise the NMF model's slots.
  	  	
  	\item if \code{NULL}, the arguments in \code{...} that have the same name 
  	as slots in the NMF model associated with the NMF strategy of name 
  	\code{method} are used to initialise these slots.
  	
  	}
  	
  	\strong{Important:}
  	Values to initialise the NMF model's slots can be passed in \code{\dots}. 
  	However, if argument \code{model} is a \code{list} -- even empty -- then 
  	all and only its elements are used to initialise the model, those in \code{...} 
  	are directly passed to the algorithm.
  	
  	So to pass a parameter to the NMF algorithm, that has the same name as a 
  	slot in the NMF model, argument \code{model} MUST be a list -- possibly empty --
  	and contains all the values one wants to use for the NMF model slots.
  	
  	If a variable appears in both argument \code{model} and \code{\dots}, 
  	the former will be used to initialise the NMF model, the latter will be 
  	passed to the NMF algorithm. 
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
  
  \item{.callback}{ Used when option \code{keep.all=FALSE} (default).
  It allows to pass a callback function that is called after each run when performing 
  multiple runs (i.e. with \code{nrun>1}).
  This is useful for example if one is also interested in saving summary measures or 
  process the result of each NMF fit before it gets discarded.
  After each run, the callback function is called with only one argument, the \code{\linkS4class{NMFfit}} object that as just been fitted:
  \code{.callback(res)} 
  Therefore all other arguments should have default values.
  
  The results of the different calls to the callback function are stored in a miscellaneous 
  slot accessible by \code{result$.callback} (assuming one ran \code{result <- nmf(...)}).
  If no error occurs \code{result$.callback} contains the list "simplified" by applying the \code{\link{sapply}} function, which 
  will try to convert a list with similar components into a \code{vector}, a \code{matrix} or a \code{data.frame}.
  If any error occurs in one of the callback calls, then global computation is NOT stopped, but the error
  is still stored in \code{result$.callback}, which is then a \code{list}.  
  
  See the examples for a sample code.
  }
  
  \item{.options}{ 
  this argument is used to set some runtime options. It can be \code{list} 
  containing the named options and their values, or, in the case only boolean 
  options need to be set, a character string that specifies which options are 
  turned on or off. The string must be composed of characters that correspond to 
  a given option. Characters '+' and '-' are used to explicitly specify on and 
  off respectively. E.g. \code{.options='tv'} will toggle on options \code{track} 
  and \code{verbose}, while \code{.options='t-v'} will toggle on option \code{track} 
  and off option \code{verbose}. Note that '+' and '-' apply to all option 
  character found after them. The default behaviour is to assume that \code{.options} 
  starts with a '+'.
    
  The following options are available (note the characters that correspond to 
  each option, to be used when \code{.options} is passed as a string):
  \describe{
  
  \item{debug - d}{ Toggle debug mode. Like option \code{verbose} but with more 
  information displayed.}
  
  \item{keep.all - k}{
  used when performing multiple runs (\code{nrun}>1): if toggled on, all 
  factorizations are saved and returned, otherwise only the factorization 
  achieving the minimum residuals is returned.
  }
  
  \item{parallel - p}{
  this option is useful on multicore *nix or Mac machine only, when performing 
  multiple runs (\code{nrun} > 1). 
  If toggled on, the runs are performed using the parallel backend defined in 
  argument \code{.pbackend}. If this is set to \code{'mc'} then one tried to 
  perform the runs using multiple cores with package \code{link[package:doMC]{doMC}} 
  -- which therefore needs to be installed.
  
  Unlike option 'P' (capital 'P'), if the computation cannot be performed in parallel, 
  then it will still be carried on sequentially.
  
  
  \strong{IMPORTANT NOTE FOR MAC OS X USERS:}
  The parallel computation is based on the \code{doMC} and \code{multicore} packages, 
  so the same care should be taken as stated in the vignette of \code{doMC}:
  \emph{\dQuote{it is not safe to use doMC from R.app on Mac OS X. Instead, you should 
  use doMC from a terminal session, starting R from the command line.}}
  }
  
  \item{parallel.required - P}{ Same as \code{p}, but an error is thrown if the 
  computation cannot be performed in parallel.
  }
  
  \item{restore.seed - r}{
  used when seeding the NMF computation with a numeric seed. When \code{TRUE} (default) 
  the random seed (\code{.Random.seed}) is restored to its value as before the call to 
  the \code{nmf} function.  
  }
  
  \item{track - t}{
  enables (resp. disables) error tracking. When \code{TRUE}, 
  the returned object's slot \code{residuals} contains the trajectory of the objective values. 
  This tracking functionality is available for all built-in algorithms.
  }
  
  \item{verbose - v}{ Toggle verbosity. If on, messages about the configuration and the 
  state of the current run(s) are displayed.}
  
  }% end describe .options
  
  }% end .options
  
  \item{.pbackend}{ define the parallel backend (from the 
  \code{\link[foreach]{foreach}} package) to use when running in parallel mode.
  See options \code{p} and \code{P} in argument \code{.options}.
  Currently it accepts the following values: \code{'mc'} or a number that specifies 
  the number of cores to use, \code{'seq'} or \code{NULL} to use sequential backend.
  }
  
  \item{rank}{The factorization rank to achieve [i.e a single positive \code{numeric}]}
   
  \item{seed}{The seeding method to use to compute the starting point passed 
  to the algorithm. 
  See section \emph{Seeding methods} %and \code{function \link{seed}} 
  for more details on the possible classes and types for argument \code{seed}.}
  
  \item{x}{The target object to estimate. It can be a \code{matrix}, 
  a \code{data.frame}, an \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}} object 
  (this requires the \code{Biobase} package to be installed). 
  See section \emph{Methods} for more details.}
  
  \item{\dots}{Extra parameters passed to the NMF algorithm's \code{run} method 
  or used to initialise the NMF model slots.  
  If argument \code{model} is not supplied as a \code{list}, ANY of the arguments 
  in \code{...} that have the same name as slots in the NMF model to be instantiated 
  will be used to initialise these slots. 
  See also the \emph{Important} paragraph in argument \code{model}.
  }
}

\section{Methods}{
\describe{

\item{x = "matrix", rank = "numeric", method = "list"}{ Performs NMF on matrix 
\code{x} for each algorithm defined in the list \code{method}.}

\item{x = "data.frame", rank = "ANY", method = "ANY"}{ Performs NMF on a 
\code{data.frame}: the target matrix is the converted data.frame 
\code{as.matrix(x)} }

\item{x = "ExpressionSet", rank = "ANY", method = "ANY"}{ Performs NMF on an 
\code{ExpressionSet}: the target matrix is the expression matrix \code{exprs(x)}.
  
 This method requires the \code{Biobase} package to be installed.
 Special methods for bioinformatics are provided in an optional layer, 
 which is automatically loaded when the \code{Biobase} is installed.
 See \code{\link{bioc-NMF}}.
}

\item{x = "matrix", rank = "numeric", method = "character"}{ Performs NMF on a 
\code{matrix} using an algorithm whose name is given by parameter \code{method}.
 The name provided must partially match the name of a registered algorithm.
 See section \emph{Algorithms} below or the package's vignette for a list of the 
 implemented algorithms and their respective names.} 

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
  \item{Multiple runs, single method:}{When \code{nrun > 1} and \code{method} is 
  NOT a \code{list}, this method returns an object of class \code{\linkS4class{NMFfitX}}.}
  \item{Multiple runs, multiple methods:}{When \code{nrun > 1} and \code{method} is a \code{list}, 
  this method returns an object of class \code{\linkS4class{NMFList}}.}
}

\section{NMF Algorithms}{

All algorithms are accessible by their respective names listed below.
The following algorithms are available:

\describe{

\item{\code{brunet}}{ Standard NMF. Based on Kullbach-Leibler divergence, it uses simple 
multiplicative updates from \emph{Lee and Seung (2000)}, enhanced to avoid numerical underflow.
See \emph{Brunet et al. (2004)}.
}

\item{\code{lee}}{ Standard NMF. Based on euclidean distance, it uses simple 
multiplicative updates. See \emph{Lee and Seung (2000)}.
}

%\item{\code{lnmf}}{ Local Nonnegative Matrix Factorization. Based on a 
%regularized Kullbach-Leibler divergence, it uses a modified version of 
%Lee and Seung's multiplicative updates.
%See \emph{Li et al. (2001)}.}

\item{\code{nsNMF}}{ Nonsmooth NMF. Uses a modified version of Lee and Seung's 
multiplicative updates for Kullbach-Leibler divergence to fit a extension 
of the standard NMF model. It is meant to give sparser results.
See \emph{Pascual-Montatno et al. (2006)}.}

\item{\code{offset}}{ Uses a modified version of Lee and Seung's multiplicative 
updates for euclidean distance, to fit a NMF model that includes an intercept. 
See \emph{Badea (2008)}. }

\item{\code{pe-nmf}}{ Pattern-Expression NMF. Uses multiplicative updates to 
minimize an objective function based on the Euclidean distance and regularized 
for effective expression of patterns with basis vectors.
See \emph{Zhang et al. (2008)}. }

\item{\code{snmf/r}, \code{snmf/l}}{ Alternating Least Square (ALS) approach 
from \emph{Kim and Park (2007)}. }

}

}

\section{Optimized C++ vs. plain R}{
Lee and Seung's multiplicative updates are used by several NMF algorithms. To improve 
speed and memory usage, a C++ implementation of the specific matrix products is used 
whenever possible. It directly computes the updates for each entry in the updated matrix, 
instead of using multiple standard matrix multiplication.

The algorithms that benefit from this optimization are: 'brunet', 'lee', 'nsNMF' and 'offset'. % and 'lnmf'
However there still exists plain R versions for these methods, which implement the updates 
as standard matrix products. These are accessible by adding the prefix '.R#' to their name: 
'.R#brunet', '.R#lee', '.R#nsNMF' and '.R#offset'.
}

\section{Seeding methods}{
	The purpose of seeding methods is to compute initial values for the factor 
	matrices in a given NMF model. 
	This initial guess will be used as a starting point by the chosen NMF algorithm.
	
	The seeding method to use in combination with the algorithm can be passed 
	to interface \code{nmf} through argument \code{seed}.
	Detailed examples of how to specify the seeding method and its parameters can 
	be found in the \emph{Examples} section of this man page and in the package's 
	vignette. 
	
	Argument \code{seed} accepts the following formats:
	
	\describe{
	
	\item{a \code{character} string:}{ giving the name of a \emph{registered} 
	seeding method. The corresponding method will be called to compute 
	the starting point. }
	
	\item{a \code{list}:}{ giving the name of a \emph{registered} 
	seeding method and, optionally, extra parameters to pass to it.}
	
	\item{a single \code{numeric}:}{ that is used to seed the random number 
	generator. The value will be used in a call to the \code{\link{set.seed}} 
	function before computing a starting point with the 'random' seeding method.
	}
	
	\item{an object that inherits from \code{\linkS4class{NMF}}:}{ it should 
	contain the	data of an initialised NMF model, that is it must contain valid 
	basis and mixture coefficient matrices.
	It will be directly passed to the algorithm's method -- via its argument \code{seed}.}
	
	\item{a \code{function}:}{ that computes the starting point. It must have 
	signature \code{(object=NMF, target=matrix, ...)} and return an object that 
	inherits from class \code{NMF}. Argument \code{object} should be used as a 
	template for the returned object.}
	}   
}

\author{Renaud Gaujoux }

\seealso{class \code{\linkS4class{NMF}}, \code{\link{utils-NMF}}, package's vignette } 

\references{

Lee, D.D. and Seung, H.S. (2000).
\emph{Algorithms for non-negative matrix factorization}. 
In \emph{NIPS}, 556--562.

Brunet, J.P. et al. (2004).
\emph{Metagenes and molecular pattern discovery using matrix factorization}. 
Proc Natl Acad Sci U S A, \bold{101}(12), 4164--4169.\cr
Original MATLAB code available from:\cr
\url{http://www.broadinstitute.org/cancer/pub/nmf}

Pascual-Montano, A. et al. (2006).
\emph{Nonsmooth nonnegative matrix factorization (nsnmf)}.
IEEE transactions on pattern analysis and machine intelligence, \bold{8}(3), 403--415.

Kim, H. and Park, H. (2007).
\emph{Sparse non-negative matrix factorizations via alternating non-negativity-constrained least squares for microarray data analysis}.
Bioinformatics. 2007; 23(12):1495-502.\cr
Original MATLAB code available from:\cr
\url{http://www.cc.gatech.edu/~hpark/software/nmfsh_comb.m}\cr
\url{http://www.cc.gatech.edu/~hpark/software/fcnnls.m}

Liviu Badea (2008).
\emph{Extracting Gene Expression Profiles Common To Colon And Pancreatic Adenocaricinoma Using Simultaneous Nonnegative Matrix Factorization}.
In Pacific Symposium on Biocomputing, \bold{13}, 279--290

%S. Li, X. Hou, and H. Zhang (2001).
%Learning spatially localized, parts-based representation. 
%\emph{In Proc. CVPR}, 2001.

Zhang J, et al. (2008). 
\emph{Pattern expression nonnegative matrix factorization: algorithm and applications to blind source separation}.
Computational intelligence and neuroscience

C. Boutsidis and E. Gallopoulos (2007) 
\emph{SVD-based initialization: A head start for nonnegative matrix factorization}. 
Pattern Recognition.
doi:10.1016/j.patcog.2007.09.010\cr
Original MATLAB code available from:\cr
\url{http://www.cs.rpi.edu/~boutsc/papers/paper1/nndsvd.m}

}

\examples{

## DATA
# generate a synthetic dataset with known classes: 100 features, 23 samples (10+5+8)
n <- 100; counts <- c(10, 5, 8); p <- sum(counts) 
V <- syntheticNMF(n, counts, noise=TRUE)
dim(V)

# build the class factor
groups <- as.factor(do.call('c', lapply(seq(3), function(x) rep(x, counts[x]))))

## RUN NMF ALGORITHMS

# run default algorithm
res <- nmf(V, 3)
res
summary(res, class=groups)

# run default algorithm multiple times (only keep the best fit)
res <- nmf(V, 3, nrun=10)
res
summary(res, class=groups)

# run default algorithm multiple times keeping all the fits
res <- nmf(V, 3, nrun=10, .options='k')
res
summary(res, class=groups)

\dontrun{
## Note: one could have equivalently done
res <- nmf(V, 3, nrun=10, .options=list(keep.all=TRUE))
}

# run nonsmooth NMF algorithm
res <- nmf(V, 3, 'nsNMF')
res
summary(res, class=groups)

\dontrun{
## Note: partial match also works
nmf(V, 3, 'ns')
}

\dontrun{
# Non default values for the algorithm's parameters can be specified in '...'
res <- nmf(V, 3, 'nsNMF', theta=0.8)
}

# compare some NMF algorithms (tracking the residual error)
res <- nmf(V, 3, list('brunet', 'lee', 'nsNMF'), seed=123456, .opt='t')
res
summary(res, class=groups)
# plot the track of the residual errors
\dontrun{plot(res)}

# run on an ExpressionSet (requires package Biobase)
\dontrun{
data(esGolub)
nmf(esGolub, 3)
}

## USING SEEDING METHODS

# run default algorithm with the Non-negative Double SVD seeding method ('nndsvd')
nmf(V, 3, seed='nndsvd')

\dontrun{
## Note: partial match also works
nmf(V, 3, seed='nn')
}

# run nsNMF algorithm, fixing the seed of the random number generator 
nmf(V, 3, 'nsNMF', seed=123456)

# run default algorithm specifying the starting point following the NMF standard model
start.std <- nmfModel(W=matrix(0.5, n, 3), H=matrix(0.2, 3, p))   
nmf(V, seed=start.std)

# to run nsNMF algorithm with an explicit starting point, this one
# needs to follow the 'NMFns' model:
start.ns <- nmfModel(model='NMFns', W=matrix(0.5, n, 3), H=matrix(0.2, 3, p))   
nmf(V, seed=start.ns)
# Note: the method name does not need to be specified as it is infered from the 
# when there is only one algorithm defined for the model.

# if the model is not appropriate (as defined by the algorihtm) an error is thrown 
# [cf. the standard model doesn't include a smoothing parameter used in nsNMF] 
\dontrun{nmf(V, method='ns', seed=start.std)}

## Callback functions
# Pass a callback function to only save summary measure of each run
res <- nmf(V, 3, nrun=3, .callback=summary)
# the callback results are simplified into a matrix
res$.callback

# Pass a custom callback function
cb <- function(obj){ sparseness(obj) >= 0.5 }
res <- nmf(V, 3, nrun=3, .callback=cb)
res$.callback

# Passs a callback function which throws an error
cb <- function(){ i<-0; function(object){ i <<- i+1; if( i == 1 ) stop('SOME BIG ERROR'); summary(object) }}
res <- nmf(V, 3, nrun=3, .callback=cb())

}

\keyword{methods}
\keyword{cluster}
\keyword{math}
\keyword{optimize}
