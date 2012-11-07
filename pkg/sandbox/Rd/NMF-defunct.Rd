\name{Defunct functions}
\alias{defunct-NMF}
\title{Defunct Functions in NMF Package}

\alias{extra}
\alias{extra-methods}
\alias{extra,NMFfit-method}

\alias{clusters}
\alias{clusters-methods}
\alias{clusters,NMF-method}
\alias{clusters,matrix-method}

\alias{metaHeatmap,matrix-method}
\alias{has.track}

\description{

The functions or variables listed here are no longer part of package NMF as 
they are not needed (any more).

}

\details{

	\code{extra} was used to access slot \code{extra} in class \code{NMFfit}.
	It has been replaced in version 0.3 by methods \code{$} and \code{$<-}. 
	See \code{\link[=$,NMFfit-method]{$}}.
	
	\code{clusters} is a defunct synonym for \code{\link{predict}} (Defunct in 0.6).

	\code{metaHeatmap,matrix} is a defunct synonym for \code{\link{aheatmap}} that 
	provides an enhanced heatmap drawing function.
}

\keyword{internal}