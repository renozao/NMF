\name{nmfModel - NMF Model Factory}

\docType{methods}
\alias{nmfModel}
\alias{nmfModel-methods}
\alias{nmfModel,missing,ANY-method}
\alias{nmfModel,missing,missing-method}
\alias{nmfModel,NULL,ANY-method}
\alias{nmfModel,numeric,matrix-method}
\alias{nmfModel,numeric,missing-method}
\alias{nmfModel,numeric,numeric-method}

\title{Factory Method for NMF Models}

\description{
 \code{nmfModel} is a generic function which provides a convenient way to build 
 NMF models.
  
 It provides a unique interface to create \code{NMF} objects that can follow 
different NMF models, and is designed to resolve potential inconsistencies 
in the matrices dimensions.

}

\details{

NMF models are defined by S4 classes that inherit from class \code{\linkS4class{NMF}}.

\code{nmfModel} methods act as factory methods to help in the creation of NMF model 
objects in common situations: creating an empty model, a model with given dimensions, 
a model with dimensions compatible with a given target matrix, ...
 
All methods return an object that inherits from class \code{NMF}, except for the call 
with no argument, which lists the NMF models defined in the session (built-in and 
user-defined).

The returned \code{NMF} objects are suitable for seeding NMF algorithms via argument 
\code{seed} of the \code{\link{nmf}} method. 
In this case the factorisation rank is implicitly set by the number of columns of the 
basis vector matrix.

}

\section{Methods}{
\describe{

\item{\code{signature(rank = "missing", target = "missing", builtin.only=FALSE, ...)}}{
	When called with no argument or argument \code{builtin.only} only, the method returns 
	a character vector that contains the name of all the NMF models currently defined. 
	These are the classes that inherits from class \code{\linkS4class{NMF}}, but not 
	from \code{\linkS4class{NMFfit}}.
	If argument \code{builtin.only=TRUE}, only the models provided by the package itself 
	are returned, discarding the user-defined models.
	
	When called with extra arguments in \var{...}, then the argument \code{builtin.only} is 
	discarded and the method is equivalent to the call \code{nmfModel(0, 0, ...)}.
	See the description of the appropriate method below.
}

\item{\code{signature(rank = "numeric", target = "numeric"
			, model='NMFstd', W, H, ...)}}{

This call creates an object of class \code{model}, using the extra parameters \code{...} 
to initialise slots that are specific to the given model.
All NMF models implement get/set methods to access the matrix factors, which can be initialised 
via arguments \code{W} and \code{H}.
For example, all the built-in models derive from class \code{\linkS4class{NMFstd}}, which has two slots, 
\var{W} and \var{H}, to hold the two factors.

If only argument \code{rank} is provided, the method creates a NMF model of dimension 0x\code{rank}x0.
That is that the basis and mixture coefficient matrices, \var{W} and \var{H}, have dimension 0x\code{rank} 
and \code{rank}x0 respectively.

If target dimensions are also provided in argument \code{target}, then the method creates a 
\code{NMF} object compatible to fit a target matrix of dimension \code{target[1]}x\code{target[2]}
That is that the basis and mixture coefficient matrices, \var{W} and \var{H}, have dimension 
\code{target[1]}x\code{rank} and \code{rank}x\code{target[2]} respectively.
If no other argument is provided, these matrices are filled with NAs.

If arguments \code{W} and/or \code{H} are provided, the method creates a NMF model where 
the basis and mixture coefficient matrices, \var{W} and \var{H}, are initialised using the values
of \code{W} and/or \code{H}.

The dimensions given by \code{target}, \code{W} and \code{H}, must be compatible.
However, whenever possible, the method will reduce the dimensions to the achieve 
dimension compatibility.    

When \code{W} and \code{H} are both provided, the \code{NMF} object created is 
suitable to seed a NMF algorithm in a call to the \code{\link{nmf}} method.
Note that in this case the factorisation rank is implicitly set by the number 
of basis vectors.

}

\item{\code{signature(rank = "numeric", target = "matrix")}}{ 

This call is equivalent to \code{nmfModel(rank, dim(target), ...)}.
That is that the returned NMF object is fits a target matrix of the same 
dimension as \code{target}.}

\item{\code{signature(rank = "missing", target = "ANY", ...)}}{ 

This call is equivalent to \code{nmfModel(0, target, ...)}.
}

\item{\code{signature(rank = "NULL", target = "ANY")}}{ 

This call is equivalent to \code{nmfModel(0, target, ...)}.}

\item{\code{signature(rank = "numeric", target = "missing")}}{ 

This call is equivalent to \code{nmfModel(rank, 0, ...)}. 
}

}}
\examples{

# List all NMF models 
nmfModel()
# or list only the built-in models
nmfModel(builtin.only=TRUE)

# create a NMF object based on one random matrix: the missing matrix is deduced
# Note this only works when using factory method NMF 
n <- 50; r <- 3; 
w <- matrix(runif(n*r), n, r) 
nmfModel(W=w)

# create a NMF object based on random (compatible) matrices
p <- 20
h <- matrix(runif(r*p), r, p)
nmfModel(W=w, H=h)

# create a NMF object based on incompatible matrices: generate an error
h <- matrix(runif((r+1)*p), r+1, p)
\dontrun{new('NMFstd', W=w, H=h)}

# same thing using the factory method: dimensions are corrected and a warning 
# is thrown saying that the dimensions used are reduced 
nmfModel(W=w, H=h)

# apply default NMF algorithm to a random target matrix
V <- matrix(runif(n*p), n, p)
\dontrun{nmf(V, r)}


}
\keyword{methods}