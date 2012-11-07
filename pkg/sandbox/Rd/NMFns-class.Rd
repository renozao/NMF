\name{NMFns-class: Model class for Non-smooth NMF}
\Rdversion{1.1}
\docType{class}
\alias{NMFns-class}
\alias{fitted,NMFns-method}
\alias{show,NMFns-method}
\alias{smoothing}
\alias{smoothing-methods}
\alias{smoothing,NMFns-method}

\title{Nonsmooth Nonnegative Matrix Factorization}
\description{
Class that implements the \emph{Nonsmooth Nonnegative Matrix Factorization} (nsNMF) 
model, required by the Nonsmooth NMF algorithm.

The Nonsmooth NMF algorithm is defined by Pascual-Montano et al. (2006) as a 
modification of the standard divergence based NMF algorithm (see section Details 
and references below). 
It aims at obtaining sparser factor matrices, by the introduction of a smoothing 
matrix. 

}

\usage{
\S4method{fitted}{NMFns}(object, W, H, S, ...)
\S4method{smoothing}{NMFns}(x, theta)
}

\arguments{
  \item{x}{ an object of class \code{NMFns}}
  \item{object}{ an object of class \code{NMFns}}
  \item{W}{ the \code{matrix} of basis vectors, i.e. the first matrix factor in 
  the non-smooth NMF model.}
  \item{H}{ the \code{matrix} of mixture coefficients, i.e. the third matrix factor 
  the non-smooth NMF model.}
  \item{S}{ the smoothing \code{matrix}, i.e. the middle matrix factor in 
  the non-smooth NMF model.}
  \item{theta}{ a single \code{numeric}} to be used as smoothing parameter 
  (see section Details).
  \item{...}{ extra parameters passed to method \code{smoothing}. So typically 
  used to pass a value for \code{theta}.}
}
\details{

The Nonsmooth NMF algorithm is a modification of the standard divergence based 
NMF algorithm (see \code{\linkS4class{NMF}}). 
Given a non-negative \eqn{n \times p}{n x p} matrix \eqn{V} and a factorization rank 
\eqn{r}, it fits the following model:
\deqn{V \equiv W S(\theta) H,}{V &#126; W S(theta) H,}
where:
\itemize{
\item \eqn{W} and \eqn{H} are such as in the standard model, that is 
non-negative matrices of dimension \eqn{n \times r}{n x r} and \eqn{r \times p}{r x p} 
respectively;
\item \eqn{S} is a \eqn{r \times r} square matrix whose entries depends on 
an extra parameter \eqn{0\leq \theta \leq 1} in the following way: 
\deqn{S = (1-\theta)I + \frac{\theta}{r} 11^T , }
where \eqn{I} is the identity matrix and \eqn{1} is a vector of ones.
}
 
The interpretation of S as a smoothing matrix can be explained as follows: Let \eqn{X} be a positive, nonzero, vector.
Consider the transformed vector \eqn{Y = S X}. If \eqn{\theta = 0}, then \eqn{Y = X} and no smoothing on \eqn{X} has occurred. 
However, as \eqn{\theta \to 1}{theta tends to 1}, the vector \eqn{Y} tends to the constant vector with all elements almost equal 
to the average of the elements of \eqn{X}. This is the smoothest possible vector in the sense of non-sparseness 
because all entries are equal to the same nonzero value, instead of having some values close to zero and 
others clearly nonzero.
 
}

\section{Algorithm}{

The Nonsmooth NMF algorithm uses a modified version of the multiplicative update 
equations in Lee & Seung's method for Kullback-Leibler divergence minimization. 
The update equations are modified to take into account the -- constant -- 
smoothing matrix. 
The modification reduces to using matrix \eqn{W S} instead of matrix \eqn{W} in 
the update of matrix \eqn{H}, and similarly using matrix \eqn{S H} instead of 
matrix \eqn{H} in the update of matrix \eqn{W}. 
%See \code{\link{NMF Algorithms}} for more details on the built-in NMF algorithms.
    
After matrix \eqn{W} have been updated, each of its columns is scaled so that 
it sums up to 1.
}

\section{Objects from the Class}{

Object of class \code{NMFns} can be created using the standard way with 
operator \code{\link{new}}

However, as for all the classes that extend class \code{\linkS4class{NMFstd}}, 
objects of class \code{NMFns} should be created using factory method 
\code{\link{nmfModel}} :

\code{new('NMFns', theta=0.8)}

\code{nmfModel(model='NMFns')}

\code{nmfModel(model='NMFns', theta=0.8)}

See \code{\link{nmfModel}} for more details on how to use the factory method. 
}
\section{Slots}{

	Class \code{NMFns} extends \code{\linkS4class{NMF}} adding a single slot:
 
	\describe{
    \item{\code{theta}:}{ Single \code{"numeric"} that contains the smoothing
    parameter. Default prototype value is \code{0.5} .}
  }
}
\section{Extends}{
Class \code{"\linkS4class{NMF}"}, directly.
}
\section{Methods}{
  \describe{
    \item{fitted}{\code{signature(object = "NMFns")}:
     returns the estimated target matrix according to the Nonsmooth-NMF model 
     \code{object}:
    \deqn{\hat{V} = \hat{V}(\theta) = W S(\theta) H}
    
    Note that this method is part of the minimum interface for NMF model, as 
    defined by class \code{\linkS4class{NMF}}.
    }
    
    \item{smoothing}{ returns the smoothing matrix \eqn{S(\theta)}.
    See section \emph{Details}.
    }	 
	
	\item{show}{\code{signature(object = "NMFns")}: standard generic \code{show} 
    method for objects of class \code{NMFns}. It calls the parent class \code{show} method 
    (i.e. for class \code{NMF}) and add the value of parameter \code{theta} to the display.}
    
    }
}
\references{
 
Alberto Pascual-Montano et al. (2006).
Nonsmooth Nonnegative Matrix Factorization (nsNMF).
\emph{IEEE Transactions On Pattern Analysis And Machine Intelligence}
, Vol. 28, No. 3, March 2006 403

}
\author{ Renaud Gaujoux  }

\seealso{	
	\code{\linkS4class{NMF}} , \code{\link{nmf-methods}}  
}
\examples{

# create a completely empty NMF object
new('NMFns')

# create a NMF object based on random (compatible) matrices
n <- 50; r <- 3; p <- 20
w <- rmatrix(n, r) 
h <- rmatrix(r, p)
nmfModel(model='NMFns', W=w, H=h)

# apply Nonsmooth NMF algorithm to a random target matrix
V <- rmatrix(n, p)
\dontrun{nmf(V, r, 'ns')}


}
\keyword{classes}
