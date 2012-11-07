\name{Deprecated functions}
\Rdversion{1.1}
\docType{methods}

\alias{deprecated-NMF}

\alias{errorPlot}
\alias{errorPlot,NMFfit-method}
\alias{errorPlot,ANY-method}

\alias{newNMF}

\alias{metaHeatmap}
\alias{metaHeatmap-methods}
\alias{metaHeatmap,NMF-method}
\alias{metaHeatmap,NMF-method}
\alias{metaHeatmap,NMFfitX-method}

\title{Deprecated Functions in NMF Package}

\description{
These functions are provided for compatibility with older versions
of NMF package only, and may be defunct as soon as the next release.
}

\usage{

\S4method{errorPlot}{NMFfit}(x, ...)
newNMF(...)

}

\arguments{

	\item{x}{ An object that inherits from class \code{\linkS4class{NMF}} or \code{\linkS4class{NMFfit}}. }
	
	\item{...}{ extra graphical arguments passed to the \code{plot} function.}
	
}

\details{
	
	\code{errorPlot} is a deprecated synonym for \code{\link{plot,NMFfit-method}}.
	\code{newNMF} is a deprecated synonym for \code{\link{nmfModel}}.
	\code{metaHeatmap} is a deprecated synonym for \code{\link{basismap}}, \code{\link{coefmap}} or \code{\link{consensusmap}}
	, to plot a heatmap of the basis matrix, mixture coefficient matrix or consensus matrix respectively.
	
}

\keyword{internal}