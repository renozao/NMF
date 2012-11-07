\name{nmfModel - Factory method for NMF models}

\docType{methods}
\alias{nmfModel}
\alias{nmfModel-methods}
\alias{nmfModel,missing,ANY-method}
\alias{nmfModel,missing,missing-method}
\alias{nmfModel,NULL,ANY-method}
\alias{nmfModel,numeric,matrix-method}
\alias{nmfModel,numeric,missing-method}
\alias{nmfModel,numeric,numeric-method}
\alias{nmfModel,matrix,ANY-method}
\alias{nmfModel,matrix,matrix-method}
\alias{nmfModel,ExpressionSet,ANY-method}
\alias{nmfModel,ANY,ExpressionSet-method}

\alias{nmfModels}

\title{Factory Method for NMF Models}

\description{

\code{nmfModels} returns the types of NMF model currently defined that can be 
used to create NMF objects.
 
\code{nmfModel} is a generic function which provides a convenient way to build 
 NMF models.
  
 It provides a unique interface to create \code{NMF} objects that can follow 
different NMF models, and is designed to resolve potential inconsistencies 
in the matrices dimensions.

}

\details{

NMF model types are defined as S4 classes that inherit from class \code{\linkS4class{NMF}}.

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

\usage{

nmfModels(builtin.only=FALSE)
\S4method{nmfModel}{numeric,numeric}(rank, target, ncol=NULL, model='NMFstd', W, H, ..., force.dim=TRUE, order.basis=TRUE)

}

\arguments{

\item{builtin.only}{
If \code{TRUE}, only the models provided by the package NMF itself are returned, 
discarding the user-defined models.
}

\item{rank}{ the target factorization rank}

\item{target}{ the target matrix dimension}

\item{ncol}{ the number of columns of the target matrix. This argument is optional. 
If not missing, it is used if \code{target} is not of length 2. It takes precedence 
over the number of columns of \code{H} -- if this latter is provided --, to define 
the target number of columns.}

\item{model}{ the type of NMF model to instantiate, as the name of the respective 
S4 class that derives from class \code{\linkS4class{NMF}}.}

\item{W}{ the template basis matrix. Its dimensions must be compatible with the 
target dimensions defined by \code{target} and/or \code{ncol} and/or \code{H}.
See argument \code{force.dim} for details on how its dimensions are possibly 
reduced.}

\item{H}{ the template mixture coefficient matrix. 
Its dimensions must be compatible with the target dimensions defined by 
\code{target} and/or \code{ncol} and/or \code{W}.
See argument \code{force.dim} for details on how its dimensions are possibly 
reduced.}

\item{...}{ extra \emph{named} arguments used to initialize other slots of the 
NMF, specific to the class given in argument \code{model}.}

\item{force.dim}{ if \code{TRUE} (default), the set of dimension of the 
instantiated model \code{M} are chosen as follows:

\code{
nrow(M) = min(n[1], nrow(W))
ncol(M) = min(n[2], ncol, ncol(H))
nbasis(M) = min(rank, ncol(W), mrow(H))
}

The basis and mixture coefficient matrices are reduced accordingly, by dropping 
extra rows and/or columns.
}

\item{order.basis}{ if \code{TRUE} (default) the rows of the mixture coefficient matrix are 
ordered in the same order as the column of the basis matrix.
This is only relevant when \code{W} and \code{H} have the same column and row 
names respectively, eventually in different orders.}

}

\section{Main Factory method}{

The main factory engine of NMF models is implemented by the method with 
signature \code{numeric, numeric}.
Other factory methods provide convenient ways of creating NMF models from e.g. a 
given target matrix or known basis/coef matrices (see section \emph{Other Factory Methods}). 

This factory method creates an object of class \code{model}, using the extra parameters \code{...} 
to initialise slots that are specific to the given model.

All NMF models implement get/set methods to access the matrix factors 
(see \code{\link{basis}}), which can be initialised via arguments \code{W} and \code{H}.
For example, all the built-in models derive from class \code{\linkS4class{NMFstd}}, 
which has two slots, \var{W} and \var{H}, to hold the two factors.

If only argument \code{rank} is provided, the method creates a NMF model of dimension 0x\code{rank}x0.
That is that the basis and mixture coefficient matrices, \var{W} and \var{H}, have dimension 0x\code{rank} 
and \code{rank}x0 respectively.

If target dimensions are also provided in argument \code{target} as a 2-length vector, then the method 
creates a \code{NMF} object compatible to fit a target matrix of dimension \code{target[1]}x\code{target[2]}.
That is that the basis and mixture coefficient matrices, \var{W} and \var{H}, have dimension 
\code{target[1]}x\code{rank} and \code{rank}x\code{target[2]} respectively.
The target dimensions can also be specified using both arguments \code{target} and \code{ncol} to define the 
number of rows and the number of columns of the target matrix respectively.
If no other argument is provided, these matrices are filled with NAs.

If arguments \code{W} and/or \code{H} are provided, the method creates a NMF model where 
the basis and mixture coefficient matrices, \var{W} and \var{H}, are initialised using the values
of \code{W} and/or \code{H}.

The dimensions given by \code{target}, \code{W} and \code{H}, must be compatible.
However if \code{force.dim=TRUE}, the method will reduce the dimensions to the achieve 
dimension compatibility whenever possible.

When \code{W} and \code{H} are both provided, the \code{NMF} object created is 
suitable to seed a NMF algorithm in a call to the \code{\link{nmf}} method.
Note that in this case the factorisation rank is implicitly set by the number 
of basis vectors.
}

\section{Other Factory Methods}{
\describe{

\item{\code{signature(rank = "missing", target = "missing", ...)}}{	
	This method is equivalent to the call \code{nmfModel(0, 0, ..., force.dim=FALSE)}.	
	This means that the dimensions of the NMF model will be taken from the 
	basis and mixture coefficient slots.
	This method is used to generate an NMF model with already existing factor 
	matrices: 
	
	\code{nmfModel(W=w)} or \code{nmfModel(H=h)} or \code{nmfModel(W=w, H=h)} 
	
	Note that an error is thrown if their dimensions are not compatible.
	
	See the description of the appropriate method below.
}

\item{\code{signature(rank = "matrix", target = "matrix", ...)}}{	
	This method is equivalent to the call \code{nmfModel(0, 0, W=rank, H=target..., force.dim=FALSE)}.
	It allows for quick calls such as e.g. \code{nmfModel(w, h)}.
	This means that the dimensions of the NMF model will be taken from the 
	basis and mixture coefficient slots \code{W} and \code{H} respectively.
	Note that an error is thrown if their dimensions are not compatible.
	
	See the description of the appropriate method below.
}

\item{\code{signature(rank = "numeric", target = "matrix")}}{ 

This call is equivalent to \code{nmfModel(rank, dim(target), ...)}.
That is that the returned NMF object fits a target matrix of the same 
dimension as \code{target}.}

\item{\code{signature(rank = "missing", target = "ANY", ...)}}{ 

This call is equivalent to \code{nmfModel(0, target, ...)}.
}

\item{\code{signature(rank = "NULL", target = "ANY")}}{ 

This call is equivalent to \code{nmfModel(0, target, ...)}.}

\item{\code{signature(rank = "numeric", target = "missing")}}{ 

This call is equivalent to \code{nmfModel(rank, 0, ...)}. 
}

\item{\code{signature(rank = "matrix", target = "ANY")}}{ 

This call is equivalent to \code{nmfModel(0, target=rank, ...)}.
This allows to call the \code{nmfModel} function with arguments \code{rank} 
and \code{target} swapped. It exists for convenience: 
 
* allows typing \code{nmfModel(V)} instead of \code{nmfModel(target=V)} to create 
a model compatible with a given matrix \code{V} (i.e. of dimension \code{nrow(V), 0, ncol(V)})

* one can pass the arguments in any order (the one that comes to the user's mind first) 
and it still works as expected...
}

\item{\code{signature(rank = "matrix", target = "matrix")}}{ 

This call is equivalent to \code{nmfModel(0, 0, W=rank, H=target, ..., force.dim=FALSE)}.
This allows for a natural short cut when creating an NMF model using existing matrices 
for both factors:

\code{nmfModel(w, h)}

}

\item{\code{signature(rank = "ANY", target = "ExpressionSet")}}{
	Uses the expression matrix of \code{target}, i.e. it is equivalent to 
	\code{nmfModel(rank, exprs(target), ...)}.
}
\item{\code{signature(rank = "ExpressionSet", target = "ANY")}}{ 
	Uses the expression matrix of \code{rank}, i.e. it is equivalent to
	\code{nmfModel(exprs(rank), target, ...)}.
}

}}

\value{
For \code{nmfModels}: a character vector containing the names of the S4 classes 
that define specific NMF models.

For \code{nmfModel}: an object that inherits from class \code{\linkS4class{NMF}}.
}

\examples{

# Build an empty NMF model 
nmfModel()

# create a NMF object based on one random matrix: the missing matrix is deduced
# Note this only works when using factory method NMF 
n <- 50; r <- 3; 
w <- rmatrix(n, r) 
nmfModel(W=w)

# create a NMF object based on random (compatible) matrices
p <- 20
h <- rmatrix(r, p)
nmfModel(W=w, H=h)
# or more conveniently without the argument names
nmfModel(w, h)

# Specify the type of NMF model (e.g. 'NMFns' for non-smooth NMF)
mod <- nmfModel(w, h, model='NMFns')
mod

# One can use an NMF model as a seed (initialization) when fitting a target matrix:
# NB: when called only with the seed the rank and NMF algorithm are selected based 
# on the input NMF model. e.g. here rank is r and the algorithm nsNMF is used
# See ?nmf. 
V <- rmatrix(n, p)
nmf(V, seed=mod)

# create an empty NMF model compatible with a given target matrix
nmfModel(V)

# create a r-ranked NMF model with a given target matrix
nmfModel(r, V)

# create a r-ranked NMF model with a given target dimensions n x p as a 2-length vector
nmfModel(r, c(n,p)) # directly
nmfModel(r, dim(V)) # or from an existing matrix <=> nmfModel(r, V)
# or alternatively passing each dimension separately
nmfModel(r, n, p)

# trying to create a NMF object based on incompatible matrices generates an error
h <- rmatrix(r+1, p)
\dontrun{
new('NMFstd', W=w, H=h)
nmfModel(w, h)
nmfModel(r+1, W=w, H=h)
}

# The factory method can be force the model to match some target dimensions
nmfModel(r, W=w, H=h)
nmfModel(r, n-1, W=w, H=h)


}
\keyword{methods}