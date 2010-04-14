\name{NMF - integration with Bioconductor}
\docType{methods}
\alias{NMF-bioc}

\alias{distance,ExpressionSet,NMF-method}

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
 
 \item computation functions that support \code{ExpressionSet} objects as inputs.
 
 \item alias functions whose names are more intuitive when NMF is
 applied to bioinformatics data.
 
 \item specialized vizualization methods that adapt the titles and legend using 
 bioinformatics terminology.
 
 \item functions to link the results with annotations, etc...
 
 }
}
\section{Methods}{
\describe{

    \item{distance}{\code{signature(target = "ExpressionSet", x = "NMF", method, ...)}: 
    returns the distance between the expression matrix and a NMF model, according to 
    a given measure. 
    If both argument \code{target} and \code{x} are missing, this function returns
    the \code{function} defined by argument \code{method}. The later can either be 
    a \code{function} or a \code{character} string that correspond to a registered 
    distance metric. For the moment only the metrics 'KL' and 'euclidean' are defined.
    
    %See function \code{\link{distance}} for more details.
    }

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


