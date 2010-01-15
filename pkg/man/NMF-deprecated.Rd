\name{NMF-deprecated}
\Rdversion{1.1}
\docType{methods}

\alias{NMF-deprecated}
\alias{clusters}
\alias{clusters-methods}
\alias{clusters,NMF-method}
\alias{clusters,matrix-method}

\title{Deprecated functions in NMF package}

\description{
These functions are provided for compatibility with older versions
of NMF package only, and may be defunct as soon as the next release.
}

\usage{

\S4method{clusters}{NMF}(x, what = c('samples', 'features'))

}

\arguments{

	\item{what}{ Specifies on which matrix (basis components or mixture coefficients) 
	the computation should be performed.}

	\item{x}{ An object that inherits from class \code{\linkS4class{NMF}}. }
	
}

\details{

\describe{
	
	\item{clusters}{ is a deprecated synonym for \code{\link{predict}}}
	
}

}