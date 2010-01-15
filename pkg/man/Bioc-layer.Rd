\name{NMF-bioc}
\docType{methods}
\alias{NMF-bioc}

\alias{metagenes}
\alias{metagenes-methods}
\alias{metagenes,NMF-method}
\alias{metagenes<-}
\alias{metagenes<-,NMF,matrix-method}

\alias{metaprofiles}
\alias{metaprofiles-methods}
\alias{metaprofiles,NMF-method}
\alias{metaprofiles<-}
\alias{metaprofiles<-,NMF,matrix-method}

\alias{nmeta}
\alias{nmeta-methods}
\alias{nmeta,NMF-method}


\title{ Layer to use the NMF package within Bioconductor}
\description{
 The package NMF provides an optional layer for working with common objects defined in the
 Bioconductor platform.
 
 It provides:
 
 \itemize{
 
 \item alias functions whose names are more intuitive when NMF is
 applied to bioinformatics data.
 
 \item specialized vizualization methods that adapt the titles and legend using 
 bioinformatics terminology.
 
 \item functions to link the results with annotations, etc...
 
 }
}
\section{Methods}{
\describe{

	\item{metagenes}{\code{signature(object = "NMF")}: 
    returns the metagenes matrix according to the model defined in \code{object}. 
    It is a simple alias to method \code{\link{basis}}.
    }

    \item{metagenes<-}{\code{signature(object = "NMF", value = "matrix")}: 
    sets the metagenes matrix in \code{object}, and returns the updated object. 
    It is a simple alias to method \code{\link{basis<-}}.
    }

    \item{metaprofiles}{\code{signature(object = "NMF")}:
    returns the metaprofiles matrix according to the model defined in \code{object}. 
    It is a simple alias to method \code{\link{coef}}.
    }

    \item{metaprofiles<-}{\code{signature(object = "NMF", value = "matrix")}: 
    sets the metaprofiles matrix in \code{object}, and returns the updated object. 
    It is a simple alias to method \code{\link{coef<-}}.
    }
    
    \item{nmeta}{\code{signature(object = "NMF")}: 
    returns the number of metagenes use in NMF model \code{object}.
    It's an alias to \code{\link{nbasis}}.
    }
}}
\keyword{methods}
\seealso{ NMF, NMF-utils}


