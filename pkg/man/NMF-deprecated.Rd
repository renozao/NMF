\name{NMF-deprecated}
\Rdversion{1.1}
\docType{methods}

\alias{NMF-deprecated}
\alias{clusters}
\alias{clusters-methods}
\alias{clusters,NMF-method}
\alias{clusters,matrix-method}

\alias{errorPlot}
\alias{errorPlot,NMFfit-method}
\alias{errorPlot,ANY-method}

\alias{newNMF}

\title{Deprecated Functions in NMF Package}

\description{
These functions are provided for compatibility with older versions
of NMF package only, and may be defunct as soon as the next release.
}

\usage{

\S4method{clusters}{NMF}(x, what = c('samples', 'features'))
\S4method{errorPlot}{NMFfit}(x, ...)
newNMF(...)

}

\arguments{

	\item{what}{ Specifies on which matrix (basis components or mixture coefficients) 
	the computation should be performed.}

	\item{x}{ An object that inherits from class \code{\linkS4class{NMF}} or \code{\linkS4class{NMFfit}}. }
	
	\item{...}{ extra graphical arguments passed to the \code{plot} function.}
	
}

\details{

	\code{clusters} is a deprecated synonym for \code{\link{predict}}.
	\code{errorPlot} is a deprecated synonym for \code{\link{plot,NMFfit-method}}.
	\code{newNMF} is a deprecated synonym for \code{\link{nmfModel}}.
	
}

\keyword{internal}