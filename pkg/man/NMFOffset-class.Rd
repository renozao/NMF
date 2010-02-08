\name{NMFOffset-class}
\Rdversion{1.1}
\docType{class}
\alias{NMFOffset-class}
\alias{fitted,NMFOffset-method}
\alias{initialize,NMFOffset-method}
\alias{offset,NMFOffset-method}
\alias{rnmf,NMFOffset,numeric-method}
\alias{show,NMFOffset-method}

\title{Nonnegative Matrix Factorization  with Offset}
\description{
Class that implements the \emph{Nonnegative Matrix Factorization with Offset} 
model, required by the NMF with Offset algorithm.

The NMF with Offset algorithm is defined by Badea (2008) as a modification of 
Lee & Seung's euclidean based NMF algorithm (see section Details and references 
below). 
It aims at obtaining 'cleaner' factor matrices, by the introduction of an offset 
matrix, explicitly modelling a feature specific baseline -- constant across samples.  

}

\section{Objects from the Class}{

Object of class \code{NMFOffset} can be created using the standard way with 
operator \code{\link{new}}

However, as for all the classes that extend class \code{\linkS4class{NMFstd}}, 
objects of class \code{NMFOffset} should be created using factory method 
\code{\link{newNMF}} :

\code{new('NMFOffset')}

\code{newNMF(model='NMFOffset')}

\code{newNMF(model='NMFOffset', W=w, offset=rep(1, nrow(w)))}

See \code{\link{newNMF}} for more details on how to use the factory method.
}

\section{Slots}{
	Class \code{NMFOffset} extends \code{\linkS4class{NMF}} adding a single slot:

	 \describe{
    \item{\code{offset}:}{A \code{"numeric"} vector to handle the common baseline 
    for each feature. Its length will always be equal to the number of features, 
    i.e. the number of rows in slot \code{W}.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{NMF}"}, directly.
}
\section{Methods}{
  \describe{
    \item{fitted}{\code{signature(object = "NMFOffset")}:
    returns the estimated target matrix according to the NMF with Offset model 
    \code{object}:
    \deqn{\hat{V} = W H + offset}
    
    Note that this method is part of the minimum interface for NMF model, as 
    defined by class \code{\linkS4class{NMF}}.
    }
    
    \item{initialize}{ \code{initialize} method for class \code{NMFOffset}. It 
    ensures consistency between slots \code{W} and slot \code{offset}. }
    
    \item{offset}{\code{signature(object = "NMFOffset")}: return the value of 
    slot \code{offset}.}
    
    \item{rnmf}{ returns the object with slots \code{W}, \code{H} and 
    \code{offset} filled with random values drawn from a uniform distribution.
    This method first calls method \code{rnmf} for \code{NMF} object to set the 
    entries of slots \code{W} and \code{H}, then sets the entries in slot 
    \code{offset} within the interval \code{[0, max(max(W),max(H))]}. }
        
    \item{show}{\code{signature(object = "NMFOffset")}: standard generic \code{show} 
    method for objects of class \code{NMFOffset}. It calls the parent class \code{show} method 
    (i.e. for class \code{NMF}) and add the value of vector \code{offset} to the display (only the first 
    five elements are displayed).}
    
    }
}

\references{ 
	Badea (2008). 
	Extracting Gene Expression Profiles Common To Colon And Pancreatic Adenocaricinoma Using Simultaneous Nonnegative Matrix Factorization. 
	In \emph{Pacific Symposium on Biocomputing}
	, \bold{13}, 279-290
}
\author{ Renaud Gaujoux \email{renaud@cbio.uct.ac.za} }

\seealso{
	\code{\linkS4class{NMF}}, \code{\link{nmf-methods}}
}
\examples{

# create a completely empty NMF object
new('NMFOffset')

# create a NMF object based on random (compatible) matrices
n <- 50; r <- 3; p <- 20
w <- matrix(runif(n*r), n, r) 
h <- matrix(runif(r*p), r, p)
newNMF(model='NMFOffset', W=w, H=h, offset=rep(0.5, nrow(w)))

# apply Nonsmooth NMF algorithm to a random target matrix
V <- matrix(runif(n*p), n, p)
\dontrun{nmf(V, r, 'offset')}


}
\keyword{classes}
