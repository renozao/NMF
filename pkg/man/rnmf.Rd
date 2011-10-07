\name{rnmf-methods: Generating random NMF models}
\docType{methods}
\alias{rnmf}
\alias{rnmf-methods}
\alias{rnmf,ANY,matrix-method}
\alias{rnmf,ANY,ExpressionSet-method}
\alias{rnmf,NMF,missing-method}
\alias{rnmf,NMF,numeric-method}
\alias{rnmf,numeric,numeric-method}

\title{
Generates Random NMF Models
}

\description{
Generates NMF models with random values drawn from a uniform distribution. 
It returns an NMF model with basis and mixture coefficient matrices 
filled with random values.
The main purpose of the function \code{rnmf} is to provide a common interface to 
generate random seeds used by the \code{\link{nmf}} function.

If necessary, extensions of the standard NMF model or custom models must define a method 
"rnmf,<NMF.MODEL.CLASS>,numeric" for initializing their specific slots other 
than the basis and mixture coefficient matrices.
In order to benefit from the complete built-in interface, the overloading methods 
should call the generic version using function \code{\link{callNextMethod}}.
See for example the method \code{\link[=rnmf,NMFOffset,numeric-method]{rnmf,NMFOffset,numeric}}: 
\code{showMethods(rnmf, class='NMFOffset', include=TRUE))}. 
}

\section{Methods}{
\describe{

\item{\code{signature(x = "ANY", target = "matrix")}}{
Generates an NMF model that is compatible with the dimensions of the matrix 
\code{target}.
The entries are uniformly drawn between \code{0} and \code{max(target)}.
If \code{x} is an NMF model, the returned value is an object of the same class 
as \code{x}.
If \code{x} is a single numeric value then the returned value is a standard NMF 
model with rank \code{x} (see \code{\linkS4class{NMFstd}}).
}

\item{\code{signature(x = "ANY", target = "ExpressionSet")}}{
Generates an NMF model that is compatible with the dimensions of the target 
expression matrix \code{exprs{target}} (see \code{\link[Biobase:ExpressionSet-class]{ExpressionSet}}).

This is a shortcut for \code{rnmf(x, exprs(target), ...)}.

}

\item{\code{signature(x = "NMF", target = "missing", ...)}}{
Generates an NMF model of the same class and dimension as \code{x}.
The entries are uniformly drawn between \code{0} and \code{max} 
(optionally specified in \code{...}) that defaults to 1.
}

\item{\code{signature(x = "NMF", target = "numeric", ncol=NULL, keep.names=TRUE, ...)}}{
Generates an NMF model of the same class as \code{x}, compatible with the 
dimensions \code{target}, that can be a single or 2-length numeric vector, to 
specify a square or rectangular target matrix respectively.
The second dimension can also be passed via argument \code{ncol}.
By default the dimnames of \code{x} are set on the returned NMF model. 
This behaviour is disabled with argument \code{keep.names=FALSE}. 
See \code{\link{nmfModel}}.
The entries are uniformly drawn between \code{0} and \code{max} 
(optionally specified in \code{...}) that defaults to 1.
}

\item{\code{signature(x = "numeric", target = "numeric", ncol=NULL, ...)}}{
This is a shortcut for \code{rnmf(nmfModel(x, target, ncol)), ...)}.

It generates a standard NMF model compatible with the dimensions \code{target}, 
that can be a single or 2-length numeric vector, to specify a square or 
rectangular target matrix respectively. See \code{\link{nmfModel}}.
The entries are uniformly drawn between \code{0} and \code{max} 
(optionally specified in \code{...}) that defaults to 1.
}

}%end_describe

}%end_Methods

\value{
An NMF model, i.e. an object that inherits from class \code{\linkS4class{NMF}}.
}

\author{Renaud Gaujoux}


\seealso{
\code{\linkS4class{NMF}}, \code{\linkS4class{NMFOffset}}, \code{\link{rmatrix}}
}
\examples{

# generate a random NMF model with rank 3 that fits a 100x20 matrix  
rnmf(3, 100, 20)
\dontshow{ stopifnot( identical(dim(rnmf(3, 100, 20)), c(100L,20L,3L)) ) }
# generate a random NMF model with rank 3 that fits a 100x100 matrix
rnmf(3, 100)
\dontshow{ stopifnot( identical(dim(rnmf(3, 100)), c(100L,100L,3L)) ) }

# generate a random NMF model based on an existing NMF model
a <- nmfModel(3, 100, 20)
b <- rnmf(a)
\dontshow{ stopifnot( !nmf.equal(a,b) ) }

# generate a random NMF model of the same class as an existing NMF model but
# of different dimension
a <- nmfModel(3, 100, 20, model='NMFns')
c <- rnmf(a, 50, 10)
\dontshow{ stopifnot( identical(dim(c), c(50L,10L,3L)) ) }
\dontshow{ stopifnot( is(c, 'NMFns') ) }


}
\keyword{methods}