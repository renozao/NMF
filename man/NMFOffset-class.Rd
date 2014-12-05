\docType{class}
\name{NMFOffset-class}
\alias{initialize,NMFOffset-method}
\alias{NMFOffset-class}
\title{NMF Model - Nonnegative Matrix Factorization with Offset}
\usage{
  \S4method{initialize}{NMFOffset}(.Object, ..., offset)
}
\arguments{
  \item{offset}{optional numeric vector used to initialise slot \sQuote{offset}.}

  \item{.Object}{ An object:  see the Details section.}

  \item{...}{data to include in the new object.  Named arguments
    correspond to slots in the class definition. Unnamed arguments must
    be objects from classes that this class extends.}
}
\description{
  This class implements the \emph{Nonnegative Matrix Factorization with
Offset} model, required by the NMF with Offset algorithm.
}
\details{
  The NMF with Offset algorithm is defined by \cite{Badea (2008)} as a modification
of the euclidean based NMF algorithm from \code{Lee2001} (see section Details and
references below).
It aims at obtaining 'cleaner' factor matrices, by the introduction of an
offset matrix, explicitly modelling a feature specific baseline
-- constant across samples.
}
\section{Methods}{
\describe{

\item{fitted}{\code{signature(object = "NMFOffset")}: Computes the target matrix estimate for an NMFOffset object.

The estimate is computed as:
\deqn{ W H + offset }
}


\item{offset}{\code{signature(object = "NMFOffset")}: The function \code{offset} returns the offset vector from an NMF model
that has an offset, e.g. an \code{NMFOffset} model.
}


\item{rnmf}{\code{signature(x = "NMFOffset", target = "numeric")}: Generates a random NMF model with offset, from class \code{NMFOffset}.

The offset values are drawn from a uniform distribution between 0 and
the maximum entry of the basis and coefficient matrices, which are drawn
by the next suitable \code{\link{rnmf}} method, which is the workhorse
method \code{rnmf,NMF,numeric}.
}


\item{show}{\code{signature(object = "NMFOffset")}: Show method for objects of class \code{NMFOffset}
}

}

}
\section{Creating objects from the Class}{


Object of class \code{NMFOffset} can be created using the standard way with
operator \code{\link{new}}

However, as for all NMF model classes -- that extend class
\code{\linkS4class{NMF}}, objects of class \code{NMFOffset} should be
created using factory method \code{\link{nmfModel}} :

\code{new('NMFOffset')}

\code{nmfModel(model='NMFOffset')}

\code{nmfModel(model='NMFOffset', W=w, offset=rep(1, nrow(w)))}

See \code{\link{nmfModel}} for more details on how to use the factory
method.
}

\section{Initialize method}{

The initialize method for \code{NMFOffset} objects tries to correct the initial
value passed for slot \code{offset}, so that it is consistent with the dimensions
of the \code{NMF} model:
it will pad the offset vector with NA values to get the length equal to the
number of rows in the basis matrix.
}
\examples{
\dontshow{# roxygen generated flag
options(R_CHECK_RUNNING_EXAMPLES_=TRUE)
}

# create a completely empty NMF object
new('NMFOffset')

# create a NMF object based on random (compatible) matrices
n <- 50; r <- 3; p <- 20
w <- rmatrix(n, r)
h <- rmatrix(r, p)
nmfModel(model='NMFOffset', W=w, H=h, offset=rep(0.5, nrow(w)))

# apply Nonsmooth NMF algorithm to a random target matrix
V <- rmatrix(n, p)
\dontrun{nmf(V, r, 'offset')}

# random NMF model with offset
rnmf(3, 10, 5, model='NMFOffset')
}
\references{
  Badea L (2008). "Extracting gene expression profiles common to colon and pancreatic adenocarcinoma using simultaneous
nonnegative matrix factorization." _Pacific Symposium on Biocomputing. Pacific Symposium on Biocomputing_, *290*, pp. 267-78.
ISSN 1793-5091, <URL: http://www.ncbi.nlm.nih.gov/pubmed/18229692>.
}
\seealso{
  Other NMF-model: \code{\link{NMFns-class}}, \code{\link{NMFstd-class}}
}
\keyword{methods}

