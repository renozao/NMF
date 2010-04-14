\name{NMF-defunct}
\alias{NMF-defunct}
\title{Defunct Functions in NMF Package}

\alias{extra}
\alias{extra-methods}
\alias{extra,NMFfit-method}


\description{

The functions or variables listed here are no longer part of package NMF as 
they are not needed (any more).

}

\usage{

% Removed in 0.3
extra(object, name)
% </Removed in 0.3>

}

\details{

	\code{extra} was used to access slot \code{extra} in class \code{NMFfit}.
	It has been replaced in version 0.3 by methods \code{$} and \code{$<-}. 
	See \code{\link[=$,NMFfit-method]{$}}.

}

\keyword{internal}